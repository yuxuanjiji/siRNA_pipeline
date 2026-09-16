# -*- coding: utf-8 -*-
"""规则筛选阶段运行器（Stage 2 · rule_filter）。

输入：Stage1 统一记录 CSV（candidates.csv，含 wt+mut 全部行）。
处理：对每行的 guide_checked 调用 legacy task7 的 evaluate()（四项串联规则，算法原样），
      把判定结果追加为统一列 gc_pct/fail_*/hit_code/hit_rule/rules_pass。
输出：candidates_rules.csv（全量注释） + candidates_passed.csv（rules_pass=1 子表）。

legacy 通过 importlib 按文件加载（模块名隔离，不污染 sys.path / 不循环 import）。
"""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

from ...common import records
from ...common.stage_io import read_manifest, write_manifest

LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
LEGACY_T7 = LEGACY_DIR / "task7_sequence_rules_filter.py"

RULES_ADD_COLS = records.RULES_COLS  # 本阶段追加列


def _load_legacy():
    spec = importlib.util.spec_from_file_location("_rulefilter_legacy_t7", LEGACY_T7)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 legacy task7：{LEGACY_T7}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _ev_row(row: dict, ev: dict) -> dict:
    """把 legacy evaluate() 结果并入统一行（数值列统一字符串化）。"""
    out = dict(row)
    out.update({
        "gc_pct": ("%.2f" % float(ev["gc_pct"])),
        "fail_gc_range": records.to_int01(ev["fail_gc_range"]),
        "fail_run6_gc": records.to_int01(ev["fail_run6_gc"]),
        "fail_run5_same": records.to_int01(ev["fail_run5_same"]),
        "fail_palindrome": records.to_int01(ev["fail_palindrome"]),
        "hit_code": str(int(ev["hit_code"])),
        "hit_rule": str(ev["hit_rule"]),
        "rules_pass": records.to_int01(ev["seq_pass"]),
    })
    return out


def run_rules(candidates_csv: str | Path, out_dir: str | Path, log=None) -> tuple[Path, Path]:
    """执行规则筛选；返回 (全量注释 CSV, 通过子表 CSV)。"""
    import logging
    log = log or logging.getLogger("sirna.rule_filter")
    candidates_csv = Path(candidates_csv)
    if not candidates_csv.exists():
        raise FileNotFoundError(f"候选记录不存在：{candidates_csv}")
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    t7 = _load_legacy()
    rows = records.read_records(candidates_csv)
    records.validate_base(rows, "rule_filter input")

    total = len(rows)
    evaluated = [_ev_row(r, t7.evaluate(r["guide_checked"])) for r in rows]

    # 通过/淘汰统计
    passed = [r for r in evaluated if r["rules_pass"] == "1"]
    rejected = [r for r in evaluated if r["rules_pass"] == "0"]
    from collections import Counter
    hit_dist = Counter(r["hit_rule"] for r in rejected)
    log.info("rule_filter: total=%d passed=%d rejected=%d hit_dist=%s",
             total, len(passed), len(rejected), dict(hit_dist))

    full_csv = out_dir / "candidates_rules.csv"
    passed_csv = out_dir / "candidates_passed.csv"
    records.write_records(full_csv, evaluated)
    records.write_records(passed_csv, passed)

    write_manifest(out_dir, {
        "stage": "02_rules",
        "input_csv": str(candidates_csv),
        "total": total, "passed": len(passed), "rejected": len(rejected),
        "hit_distribution": dict(hit_dist),
        "algorithm": "legacy task7 func_filter（与 OligoFormer 对拍口径）",
    })
    return full_csv, passed_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage2 规则筛选：四项规则 -> 统一记录 rules_* 列")
    ap.add_argument("--input", type=Path, required=True, help="Stage1 candidates.csv")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.rule_filter")
    try:
        f, p = run_rules(args.input, args.out_dir, log)
    except Exception as e:  # noqa: BLE001
        log.error("rule_filter 失败：%s", e)
        return 3
    log.info("rule_filter 完成：%s / %s", f, p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
