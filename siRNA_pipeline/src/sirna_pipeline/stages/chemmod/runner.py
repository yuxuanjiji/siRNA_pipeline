# -*- coding: utf-8 -*-
"""化学修饰建议阶段运行器（Stage 8 · chemmod，第六步）。

输入：Stage7 rank 输出（含 final_rank 的排序表，通常取 Top-N）。
处理：按 `rules.build_chem_plan`（ESC-19 骨架 + 免疫掩蔽 + seed OMe 可选旋钮）对每条
      候选的 (guide_checked, rc(guide)) 生成逐位修饰建议，追加 chem_* 列写出。

设计声明：
  * 只做“建议/可解释输出”，不淘汰候选、不改序列；修饰不改碱基 → 序列级风险
    （seed 脱靶等）由第五步处理，不在本阶段补救（与任务19 论证卡 §1.4 一致）；
  * 旧目录 21 nt 占位（项目搭建/化学修饰/*）原地不动，本阶段以 19 nt 主链口径工作，
    21 nt 化作为备注交给后续成药化迭代；
  * 产出叶阶段 CSV + manifest，供报告与任务19 论证卡正式替换。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from ...common import records
from ...common.stage_io import write_manifest
from .rules import RULE_VERSION, RULES_REF, build_chem_plan

_CHEM_JSON_COLS = ("chem_guide_f_positions", "chem_guide_ome_positions",
                   "chem_guide_ps_bonds", "chem_sense_f_positions",
                   "chem_sense_ome_positions", "chem_sense_ps_bonds")
_CHEM_TEXT_COLS = ("chem_guide_mod", "chem_sense_mod")


def _fmt(v):
    if isinstance(v, bool):
        return "1" if v else "0"
    return str(v)


def _ranked_rows(rows: list[dict], top_n: int | None) -> list[dict]:
    """取 final_rank 非空行并按排名升序截断 top_n。"""
    ranked = [r for r in rows
              if str(r.get("final_rank", "")).strip() not in ("", "nan")]
    ranked.sort(key=lambda r: int(float(r["final_rank"])))
    if top_n:
        ranked = ranked[: int(top_n)]
    return ranked


def run_chemmod(
    input_csv: str | Path,
    out_dir: str | Path,
    top_n: int | None = None,
    seed_ome: bool = False,
    log=None,
) -> Path:
    """执行化学修饰建议；返回追加 chem_* 列的统一记录 CSV（叶阶段产物）。"""
    import logging
    log = log or logging.getLogger("sirna.chemmod")
    input_csv, out_dir = Path(input_csv), Path(out_dir)
    if not input_csv.exists():
        raise FileNotFoundError(f"输入记录不存在：{input_csv}")
    os.makedirs(out_dir, exist_ok=True)

    rows = records.read_records(input_csv)
    picked = _ranked_rows(rows, top_n)
    if not picked:
        raise ValueError("无带 final_rank 的行进入化学修饰（请先跑 rank 阶段）")
    log.info("chemmod: 输入 %d 行 -> 取排名 Top-N(%s) %d 行",
             len(rows), top_n or "全部", len(picked))

    out_rows = []
    for r in picked:
        guide = str(r.get("guide_checked", ""))
        plan = build_chem_plan(guide, passenger=None, seed_ome=seed_ome)
        nr = dict(r)
        for col in _CHEM_TEXT_COLS:
            nr[col] = plan["guide_mod"] if col == "chem_guide_mod" else plan["sense_mod"]
        for col in _CHEM_JSON_COLS:
            key = {"chem_guide_f_positions": "guide_f_positions",
                   "chem_guide_ome_positions": "guide_ome_positions",
                   "chem_guide_ps_bonds": "guide_ps_bonds",
                   "chem_sense_f_positions": "sense_f_positions",
                   "chem_sense_ome_positions": "sense_ome_positions",
                   "chem_sense_ps_bonds": "sense_ps_bonds"}[col]
            nr[col] = json.dumps(plan[key], ensure_ascii=False)
        nr["chem_5p_phosphate"] = _fmt(plan["guide_5p_phosphate"])
        nr["chem_notes"] = json.dumps(plan["notes"], ensure_ascii=False)
        out_rows.append(nr)

    out_csv = out_dir / "candidates_chemmod.csv"
    records.write_records(out_csv, out_rows)

    n_notes = sum(1 for r in out_rows if "免疫 motif" in (r.get("chem_notes") or ""))
    write_manifest(out_dir, {
        "stage": "08_chemmod",
        "input_csv": str(input_csv),
        "rule_version": RULE_VERSION,
        "rules_ref": RULES_REF,
        "rows_input": len(rows),
        "rows": len(out_rows),
        "top_n": top_n,
        "seed_ome": bool(seed_ome),
        "rows_with_immune_note": n_notes,
    })
    log.info("chemmod 完成：%s（%d 行）", out_csv, len(out_rows))
    return out_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage8 化学修饰：ESC-19 逐位修饰建议 -> chem_* 列")
    ap.add_argument("--input", type=Path, required=True, help="rank 输出(含 final_rank)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=None, help="只取前 N 名（默认全部）")
    ap.add_argument("--seed-ome", action="store_true",
                    help="开启 seed(g2–g8) 全 2′-OMe 旋钮")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.chemmod")
    try:
        out = run_chemmod(args.input, args.out_dir, top_n=args.top_n,
                          seed_ome=args.seed_ome, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("chemmod 失败：%s", e)
        return 3
    log.info("chemmod 完成：%s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
