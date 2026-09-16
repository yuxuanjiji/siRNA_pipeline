# -*- coding: utf-8 -*-
"""结构检测阶段运行器（Stage 3 · structure）。

输入：Stage2 通过子表（rules_pass=1）的统一记录（若含 rules_pass 列则自动只取通过行）。
处理：
  1. 读 CDS FASTA，按记录 cds_start 经 duplex.mrna57_from_cds 重建 57nt 上下文
     （侧翼越界窗口 -> mrna=None，走单链模式，cofold/plfold 记为空）；
  2. 调用 legacy StructureDetector（importlib 隔离加载；算法/阈值零改动），
     detect_batch(guide 列表, mrna 列表)；
  3. 把每行结果 dict 映射为统一 structure_* 列，追加 structure_pass/reject_reason/warnings。
无 ViennaRNA 时类内自动降级（python 绑定->CLI->纯 Python NN 近似），
structure_reliable 如实标注，流程不中断（与原模块降级语义一致）。
"""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

from ...common import duplex, records, seqio
from ...common.stage_io import file_sha256, write_manifest

LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
LEGACY_SD = LEGACY_DIR / "structure_detector.py"

# 结果 dict -> 统一列映射（structure 模块输出键名 -> canonical 列名）
_KEY_MAP = {
    "mfe": "mfe", "mfe_structure": "mfe_structure",
    "max_stem_length": "max_stem_length", "max_stem_gc": "max_stem_gc",
    "stem_loop_count": "stem_loop_count",
    "delta_G_5end": "delta_G_5end", "delta_G_3end": "delta_G_3end",
    "delta_deltaG_ends": "delta_deltaG_ends",
    "cofold_mfe": "cofold_mfe", "cofold_structure": "cofold_structure",
    "internal_loops": "internal_loops", "bulges": "bulges",
    "binding_energy": "binding_energy", "cofold_available": "cofold_available",
    "target_accessibility": "target_accessibility",
    "target_5end_access": "target_5end_access",
    "target_3end_access": "target_3end_access",
    "plfold_available": "plfold_available",
    "structure_reliable": "structure_reliable",  # 类输出字段(后端可靠标记)
    "pass_mfe": "pass_mfe", "pass_end_diff": "pass_end_diff",
    "pass_stem": "pass_stem", "pass_cofold": "pass_cofold",
    "pass_accessibility": "pass_accessibility",
    "reject_reason": "reject_reason", "warnings": "warnings",
}

# canonical 列中本阶段只新增这些（structure_pass 由 pass_all 派生）
CANON_OUT_STRUCT = [c for c in records.STRUCTURE_COLS
                    if c not in ("structure_pass", "pass_mfe", "pass_end_diff",
                                 "pass_stem", "pass_cofold", "pass_accessibility",
                                 "cofold_available", "plfold_available")]
# 供过滤判定保留的内部字段（不作为 canonical 追加列）
_EXTRA = ("pass_all", "structure_pass")


def _load_legacy():
    spec = importlib.util.spec_from_file_location("_structure_legacy", LEGACY_SD)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 legacy structure_detector：{LEGACY_SD}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _fmt_num(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return "%.10g" % v
    if isinstance(v, bool):
        return "1" if v else "0"
    return str(v)


def _map_result(row: dict, res: dict) -> dict:
    out = dict(row)
    for k_src, k_dst in _KEY_MAP.items():
        v = res.get(k_src)
        out[k_dst] = _fmt_num(v) if not isinstance(v, str) else v
    # 派生：pass_all -> structure_pass
    out["structure_pass"] = records.to_int01(res.get("pass_all"))
    return out


def run_structure(
    input_csv: str | Path,
    cds_fasta: str | Path,
    out_dir: str | Path,
    detector_kwargs: dict | None = None,
    only_passed: bool = True,
    log=None,
) -> Path:
    """执行结构检测；返回追加 structure_* 列的统一记录 CSV。"""
    import logging
    log = log or logging.getLogger("sirna.structure")
    input_csv, cds_fasta, out_dir = Path(input_csv), Path(cds_fasta), Path(out_dir)
    if not input_csv.exists():
        raise FileNotFoundError(f"输入记录不存在：{input_csv}")
    if not cds_fasta.exists():
        raise FileNotFoundError(f"CDS FASTA 不存在：{cds_fasta}")
    os.makedirs(out_dir, exist_ok=True)

    legacy = _load_legacy()
    rows = records.read_records(input_csv)
    if only_passed and "rules_pass" in (rows[0] if rows else {}):
        keep = [r for r in rows if r["rules_pass"] == "1"]
        log.info("structure: 输入 %d 行 -> 仅通过 rules_pass 的 %d 行", len(rows), len(keep))
        rows = keep
    if not rows:
        raise ValueError("无候选行进入结构检测")
    records.validate_base(rows, "structure input")

    header, cds = seqio.load_fasta_rna(cds_fasta)
    log.info("structure: CDS %s 长度=%d", header.split()[0][:60] if header else "-", len(cds))

    guides, mrnas = [], []
    n_no_ctx = 0
    for r in rows:
        g = r["guide_checked"]
        m57 = None
        try:
            m57 = duplex.mrna57_from_cds(cds, int(r["cds_start"]))
        except (TypeError, ValueError):
            m57 = None
        if m57 is None:
            n_no_ctx += 1
        guides.append(g)
        mrnas.append(m57)
    log.info("structure: 无 57nt 上下文(单链模式) 行数=%d/%d", n_no_ctx, len(rows))

    det = legacy.StructureDetector(**(detector_kwargs or {}))
    status = det.get_tool_status()
    log.info("structure: ViennaRNA status level=%s binding=%s rnafold=%s cofold=%s plfold=%s",
             status.get("level"), status.get("python_binding"),
             status.get("rnafold_cli"), status.get("rnacofold_cli"),
             status.get("rnaplfold_cli"))

    results = det.detect_batch(guides, mrnas)
    if len(results) != len(rows):
        raise RuntimeError("detect_batch 返回行数与输入不一致")
    out_rows = [_map_result(r, res) for r, res in zip(rows, results)]

    full_csv = out_dir / "candidates_structure.csv"
    records.write_records(full_csv, out_rows)

    n_pass = sum(1 for r in out_rows if r["structure_pass"] == "1")
    n_rel = sum(1 for r in out_rows if r.get("structure_reliable") == "1")
    write_manifest(out_dir, {
        "stage": "03_structure",
        "input_csv": str(input_csv),
        "cds_fasta": str(cds_fasta), "cds_sha256": file_sha256(cds_fasta),
        "vienna_level": status.get("level"),
        "python_binding": status.get("python_binding"),
        "rnafold_cli": status.get("rnafold_cli"),
        "rnacofold_cli": status.get("rnacofold_cli"),
        "rnaplfold_cli": status.get("rnaplfold_cli"),
        "rows": len(out_rows), "structure_pass": n_pass,
        "structure_reliable_rows": n_rel,
        "no_57nt_context_rows": n_no_ctx,
        "detector_kwargs": detector_kwargs or {},
    })
    return full_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage3 结构检测：RNAfold/cofold/plfold + 硬阈值 -> structure_* 列")
    ap.add_argument("--input", type=Path, required=True,
                    help="Stage2 passed 子表（candidates_passed.csv）")
    ap.add_argument("--cds-fasta", type=Path, required=True,
                    help="SFRP1 CDS FASTA（重建 57nt 上下文）")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--all-rows", action="store_true",
                    help="不过滤 rules_pass，处理输入全部行")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.structure")
    try:
        out = run_structure(args.input, args.cds_fasta, args.out_dir,
                            only_passed=not args.all_rows, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("structure 失败：%s", e)
        return 3
    log.info("structure 完成：%s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
