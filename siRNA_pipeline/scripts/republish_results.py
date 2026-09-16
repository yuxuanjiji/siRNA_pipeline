# -*- coding: utf-8 -*-
"""重发 `outputs/results/`（大赛交付目录）：复用 orchestrator 的发布函数，口径完全一致。

背景（**重要**）：e2e 测试原本会把 rank/chemmod 产物发布到 `outputs/results/`，
导致每次跑全量单测都用**合成序列**覆盖真实交付物（已在 tests/test_pipeline.py 隔离修复）。
本脚本在修复后用真实 run 的 03/04/06 阶段产物重跑 07(rank)+08(chemmod) 并发布：

  1. 旧口径产物（reweighted，α0.4/β0.6、未分层）→ 备份到 `outputs/results/_base/`
     （由 `outputs/runs/rank_reweighted/candidates_ranked.csv` 确定性重建）；
  2. 新口径产物（含变体层）→ 发布到 `outputs/results/`；
  3. 打印可核对的摘要（行数、Top-5、是否含 SFRP1 真实窗口）。

用法：
    python scripts/republish_results.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sirna_pipeline.common import records                      # noqa: E402
from sirna_pipeline.common.results_export import export_standard_results  # noqa: E402
from sirna_pipeline.pipeline.orchestrator import (             # noqa: E402
    _export_submission, _publish_chemmod, _publish_results, _write_summary,
    build_config)
from sirna_pipeline.stages.chemmod.runner import run_chemmod    # noqa: E402
from sirna_pipeline.stages.rank.ranker import run_rank          # noqa: E402


def _backup_base(results_dir: Path, rank_reweighted: Path, top_n: int,
                 log) -> list[str]:
    """旧口径（未分层）产物备份到 results/_base/，由 reweighted run 确定性重建。"""
    base = results_dir / "_base"
    base.mkdir(parents=True, exist_ok=True)
    wrote = []
    rank_csv = rank_reweighted / "candidates_ranked.csv"
    if not rank_csv.exists():
        log.warning("旧口径 run 缺失，跳过备份：%s", rank_csv)
        return wrote
    shutil.copyfile(rank_csv, base / "rank_final.csv")
    wrote.append("rank_final.csv")
    rows = sorted((r for r in records.read_records(rank_csv) if r.get("final_rank")),
                  key=lambda r: int(r["final_rank"]))
    records.write_records(base / "rank_top.csv", rows[:top_n])
    wrote.append("rank_top.csv")
    export_standard_results(
        rank_csv, base / "results.csv", top_n=top_n,
        track="siRNA（反义寡核苷酸）靶向序列设计 — SFRP1(NM_003012)",
        rank_meta={"note": "旧口径（未分层）Top-N；由 outputs/runs/rank_reweighted 重建"})
    wrote.append("results.csv")
    cm = rank_reweighted / "candidates_chemmod.csv"
    if cm.exists():
        shutil.copyfile(cm, base / "chemmod_top.csv")
        wrote.append("chemmod_top.csv")
    return wrote


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path,
                    default=ROOT / "outputs" / "runs" / "sfrp1_oligo_on")
    ap.add_argument("--out-run", default="rank_variant_layer")
    ap.add_argument("--cfg-dir", type=Path, default=ROOT / "configs")
    args = ap.parse_args()

    import logging
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger("republish")

    cfg = build_config(args.cfg_dir, run_name=args.out_run)
    root = Path(cfg["paths"]["outputs_root"])
    results_dir = root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- 0) 备份旧口径 ----
    top_n = int((cfg.get("stages", {}).get("ranking", {}) or {}).get("top_n")
                or (cfg.get("summary", {}) or {}).get("top_n") or 50)
    wrote = _backup_base(results_dir, root / "runs" / "rank_reweighted", top_n, log)
    print("[0] 旧口径备份 -> outputs/results/_base/ %s" % wrote)

    # ---- 1) 重跑 07 rank（新口径：含变体层）----
    base = args.run_dir / "03_structure" / "candidates_passed.csv"
    extras = [args.run_dir / "04_thermo" / "candidates_thermo.csv",
              args.run_dir / "06_toxicity" / "candidates_toxicity.csv",
              args.run_dir / "06_oligoformer" / "candidates_oligoformer.csv"]
    scfg = cfg["stages"]["ranking"]
    rank_cfg = {k: v for k, v in scfg.items() if k != "enabled"}
    out_dir = root / "runs" / args.out_run / "07_ranking"
    rank_csv, meta = run_rank(base, extras, out_dir, cfg=rank_cfg, log=log)
    print("[1] rank: eligible=%d → %s" % (meta["eligible"], rank_csv))
    print("    variant_layer=%s" % json.dumps(meta["variant_layer"], ensure_ascii=False)[:200])
    _publish_results(cfg, rank_csv, meta, args.out_run, log)

    # ---- 2) 重跑 08 chemmod ----
    cm_top = scfg.get("top_n") if "chemmod" in cfg["stages"] else None
    cm_cfg = cfg["stages"].get("chemmod", {}) or {}
    cm_csv = run_chemmod(rank_csv, root / "runs" / args.out_run / "08_chemmod",
                         top_n=cm_cfg.get("top_n", cm_top or top_n),
                         seed_ome=bool(cm_cfg.get("seed_ome", False)), log=log)
    print("[2] chemmod → %s" % cm_csv)
    _publish_chemmod(cfg, cm_csv, args.out_run, log)

    # ---- 3) 导出标准结果文件（含化学修饰列）----
    _export_submission(cfg, rank_csv, meta, args.out_run, log, chem_csv=cm_csv)

    summary = {"run_name": args.out_run,
               "rank_meta": meta,
               "chemmod_csv": str(cm_csv),
               "variant_layer": meta.get("variant_layer"),
               "base_backup": wrote,
               "note": "变体层口径重发；旧口径见 outputs/results/_base/"}
    _write_summary(cfg, summary, args.out_run)

    # ---- 4) 核对 ----
    res = records.read_records(results_dir / "results.csv")
    print("\n[3] 发布核对：results.csv %d 行" % len(res))
    for r in sorted(res, key=lambda r: int(r["final_rank"]))[:5]:
        print("    #%s %-10s %-18s final=%s c_match=%s"
              % (r["final_rank"], r.get("window_id"), r.get("candidate_id"),
                 r.get("final_score"), r.get("c_match", "")))
    wids = {r.get("window_id") for r in res}
    print("    窗口样例:", sorted(wids)[:6], "...")
    # 标准结果文件表头固定（附件5）；变体层列在 rank_final/rank_top 中
    rk = records.read_records(results_dir / "rank_final.csv")
    print("    rank_final.csv %d 行；含变体层列: %s"
          % (len(rk), all(k in rk[0] for k in ("c_match", "q_window_anchor",
                                               "final_score_base", "variant_layer"))))
    print("written:", results_dir)


if __name__ == "__main__":
    main()
