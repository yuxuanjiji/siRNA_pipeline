# -*- coding: utf-8 -*-
"""毒性检测阶段运行器（Stage 6 · toxicity）。

输入：Stage2 passed 子表（rules_pass=1，默认）或任意统一记录。
处理：调用 legacy ToxicityDetector.detect_batch（importlib 隔离加载），
      子模块 A（seed g2–g7 查 4096 cell_viability 表）默认缓存广播；
      子模块 B（UGUGU/GUCCUUCAA/poly‑U/U‑rich 免疫 motif 扫描）。
      结果映射为统一 tox_*/imm_* 列（硬过滤默认关闭，软惩罚语义交由排序阶段）。
资源：cell_viability.txt 位于本包 resources/，通过环境变量 TOX_VIABILITY_TABLE
      在加载 legacy 模块前注入（legacy 模块常量在 import 时读取该环境变量）。
"""
from __future__ import annotations

import argparse
import importlib.util
import numbers
import os
from pathlib import Path

from ...common import records
from ...common.stage_io import write_manifest

LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
LEGACY_TOX = LEGACY_DIR / "toxicity_detector.py"
RESOURCES = Path(__file__).resolve().parent / "resources" / "cell_viability.txt"

# legacy 返回键 -> 统一列（bool 自动转 0/1，数值转字符串，None->''）
_TOX_KEYS = [
    "tox_seed", "tox_viability_score", "tox_viability_flag", "tox_viability_miss",
    "tox_broadcast",
    "imm_ugugu", "imm_guccuucaa",
    "imm_polyu_max_run_guide", "imm_polyu_max_run_passenger",
    "imm_polyu_high", "imm_polyu_warn",
    "imm_u_content_guide", "imm_u_content_passenger", "imm_u_content_warn",
    "imm_high_flag", "imm_warn_flag", "imm_flag", "imm_hits_detail",
]


def _load_legacy():
    # 必须先于模块导入设置：DEFAULT_TABLE_CANDIDATES 在 import 时构建并读取环境变量
    if "TOX_VIABILITY_TABLE" not in os.environ:
        os.environ["TOX_VIABILITY_TABLE"] = str(RESOURCES)
    spec = importlib.util.spec_from_file_location("_toxicity_legacy", LEGACY_TOX)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 legacy toxicity_detector：{LEGACY_TOX}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _cell(v):
    """数值/布尔 -> 字符串（None->''）。"""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, numbers.Real):
        return "%.10g" % float(v)
    return str(v)


def _map_row(row: dict, res: dict) -> dict:
    out = dict(row)
    for k in _TOX_KEYS:
        out[k] = _cell(res.get(k))
    return out


def run_toxicity(
    input_csv: str | Path,
    out_dir: str | Path,
    only_passed: bool = True,
    viability_threshold: float = 50.0,
    log=None,
) -> Path:
    """执行毒性检测；返回追加 tox_*/imm_* 列的统一记录 CSV。"""
    import logging
    log = log or logging.getLogger("sirna.toxicity")
    input_csv, out_dir = Path(input_csv), Path(out_dir)
    if not input_csv.exists():
        raise FileNotFoundError(f"输入记录不存在：{input_csv}")
    os.makedirs(out_dir, exist_ok=True)

    legacy = _load_legacy()
    rows = records.read_records(input_csv)
    if only_passed and rows and "rules_pass" in rows[0]:
        keep = [r for r in rows if r["rules_pass"] == "1"]
        log.info("toxicity: 输入 %d 行 -> rules_pass=1 的 %d 行", len(rows), len(keep))
        rows = keep
    if not rows:
        raise ValueError("无候选行进入毒性检测")
    records.validate_base(rows, "toxicity input")

    det = legacy.ToxicityDetector(viability_threshold=float(viability_threshold))
    recs = [{"guide": r["guide_checked"], "window_id": r["window_id"],
             "variant_id": r["variant_id"]} for r in rows]
    results = det.detect_batch(recs, use_cache=True)
    if len(results) != len(rows):
        raise RuntimeError("detect_batch 返回行数与输入不一致")

    out_rows = [_map_row(r, res) for r, res in zip(rows, results)]
    full_csv = out_dir / "candidates_toxicity.csv"
    records.write_records(full_csv, out_rows)

    n_tox = sum(1 for r in out_rows if r["tox_viability_flag"] == "1")
    n_imm_high = sum(1 for r in out_rows if r["imm_high_flag"] == "1")
    n_bc = sum(1 for r in out_rows if r["tox_broadcast"] == "1")
    write_manifest(out_dir, {
        "stage": "06_toxicity",
        "input_csv": str(input_csv),
        "viability_threshold": float(viability_threshold),
        "rows": len(out_rows),
        "tox_viability_flag": n_tox, "imm_high_flag": n_imm_high,
        "tox_broadcast": n_bc,
        "viability_table": str(RESOURCES),
        "cache_hits": getattr(det, "cache_hits", None),
        "cache_misses": getattr(det, "cache_misses", None),
    })
    return full_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage6 毒性检测：seed 活力查表 + 免疫 motif 扫描 -> tox_*/imm_* 列")
    ap.add_argument("--input", type=Path, required=True,
                    help="Stage2 passed 子表（candidates_passed.csv）")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--threshold", type=float, default=50.0,
                    help="viability 阈值（OligoFormer 默认 50）")
    ap.add_argument("--all-rows", action="store_true")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.toxicity")
    try:
        out = run_toxicity(args.input, args.out_dir,
                           only_passed=not args.all_rows,
                           viability_threshold=args.threshold, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("toxicity 失败：%s", e)
        return 3
    log.info("toxicity 完成：%s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
