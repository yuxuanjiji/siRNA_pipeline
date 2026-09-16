# -*- coding: utf-8 -*-
"""热力学参数计算阶段运行器（Stage 4 · thermo）。

输入：Stage2 passed 子表（rules_pass=1，默认）或任意统一记录。
处理：按记录 cds_start 重建 57nt mRNA 上下文，对每条 (guide, mrna57)
      调用 legacy siRNAThermoCalculator().calculate()（PyTorch 最近邻模型，
      importlib 隔离加载，算法/默认系数零改动），把返回 dict 映射为统一 dG_* 列。

无 57nt 上下文（CDS 边界窗口）的行：热力学特征留空，并在 manifest 记数。
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import numbers
import os
from pathlib import Path

from ...common import duplex, records, seqio
from ...common.stage_io import file_sha256, write_manifest

LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
LEGACY_TC = LEGACY_DIR / "thermo_calculator.py"

# legacy 返回键 -> 统一列
_KEY_MAP = {
    "dG_total": "dG_total", "dH_total": "dH_total", "dS_total": "dS_total",
    "Tm": "Tm", "dG_core_canonical": "dG_core_canonical",
    "dG_mismatch_total": "dG_mismatch_total",
    "dG_flank_total": "dG_flank_total",
    "dG_mismatch_base": "dG_mismatch_base",
    "dG_mismatch_context": "dG_mismatch_context",
    "dG_mismatch_positional_extra": "dG_mismatch_positional_extra",
    "n_mismatch": "n_mismatch",
}
# 需要“列表/JSON 串”特殊处理
_JSON_KEYS = ("dG_mismatch_by_position",)


def _load_legacy():
    spec = importlib.util.spec_from_file_location("_thermo_legacy", LEGACY_TC)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 legacy thermo_calculator：{LEGACY_TC}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _fmt_val(v):
    """Tensor/数值 -> 字符串；None -> ''。"""
    if v is None:
        return ""
    if hasattr(v, "item"):
        try:
            v = v.item()
        except Exception:  # noqa: BLE001
            pass
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, numbers.Real):
        return "%.10g" % float(v)
    return str(v)


def _calc_row(calc, guide: str, mrna57: str | None) -> dict:
    """计算单行并映射为追加列 dict。mrna57=None -> 全空。"""
    if mrna57 is None:
        return {}
    res = calc.calculate(guide, mrna57)
    out = {}
    for k_src, k_dst in _KEY_MAP.items():
        if k_src == "n_mismatch":
            out[k_dst] = _fmt_val(res.get(k_src))
            continue
        out[k_dst] = _fmt_val(res.get(k_src))
    # 逐位惩罚：len-19 数值列表 -> JSON 串（排序/消融用；'' 表示缺失）
    vp = res.get("dG_mismatch_by_position")
    if vp is not None:
        if hasattr(vp, "item") and hasattr(vp, "tolist"):
            vp = vp.tolist()
        if isinstance(vp, (list, tuple)) and len(vp) == 19:
            out["dG_mismatch_by_position"] = json.dumps(
                [_fmt_val(x) for x in vp], ensure_ascii=False)
    return out


def run_thermo(
    input_csv: str | Path,
    cds_fasta: str | Path,
    out_dir: str | Path,
    only_passed: bool = True,
    log=None,
) -> Path:
    """执行热力学计算；返回追加 dG_* 列的统一记录 CSV。"""
    import logging
    log = log or logging.getLogger("sirna.thermo")
    input_csv, cds_fasta, out_dir = Path(input_csv), Path(cds_fasta), Path(out_dir)
    if not input_csv.exists():
        raise FileNotFoundError(f"输入记录不存在：{input_csv}")
    if not cds_fasta.exists():
        raise FileNotFoundError(f"CDS FASTA 不存在：{cds_fasta}")
    os.makedirs(out_dir, exist_ok=True)

    legacy = _load_legacy()
    rows = records.read_records(input_csv)
    if only_passed and rows and "rules_pass" in rows[0]:
        keep = [r for r in rows if r["rules_pass"] == "1"]
        log.info("thermo: 输入 %d 行 -> rules_pass=1 的 %d 行", len(rows), len(keep))
        rows = keep
    if not rows:
        raise ValueError("无候选行进入热力学计算")
    records.validate_base(rows, "thermo input")

    header, cds = seqio.load_fasta_rna(cds_fasta)
    calc = legacy.siRNAThermoCalculator()
    n_no_ctx = 0
    out_rows = []
    for r in rows:
        g = r["guide_checked"]
        m57 = None
        try:
            m57 = duplex.mrna57_from_cds(cds, int(r["cds_start"]))
        except (TypeError, ValueError):
            m57 = None
        if m57 is None:
            n_no_ctx += 1
        out_rows.append({**r, **_calc_row(calc, g, m57)})

    full_csv = out_dir / "candidates_thermo.csv"
    records.write_records(full_csv, out_rows)
    write_manifest(out_dir, {
        "stage": "04_thermo",
        "input_csv": str(input_csv),
        "cds_fasta": str(cds_fasta), "cds_sha256": file_sha256(cds_fasta),
        "rows": len(out_rows), "no_57nt_context_rows": n_no_ctx,
        "algorithm": "legacy thermo_calculator (NN 模型, 默认修正系数)",
    })
    return full_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage4 热力学：NN 杂交 ΔG/ΔH/ΔS/Tm + 错配逐位 -> dG_* 列")
    ap.add_argument("--input", type=Path, required=True,
                    help="Stage2 passed 子表（candidates_passed.csv）")
    ap.add_argument("--cds-fasta", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--all-rows", action="store_true")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.thermo")
    try:
        out = run_thermo(args.input, args.cds_fasta, args.out_dir,
                         only_passed=not args.all_rows, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("thermo 失败：%s", e)
        return 3
    log.info("thermo 完成：%s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
