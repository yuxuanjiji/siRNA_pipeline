# -*- coding: utf-8 -*-
"""变体层验收：在真实 run 上验证"WT 逐位一致 + 只改变体排序"，并量化排序变化。

验收条件（改动1 的核心承诺）：
  1. `variant_layer.enabled=True` 时，**WT 行的 score_thermo/score_oligo/penalty_total/
     final_score 与旧口径逐位一致**（不破坏冻结结果）；
  2. `enabled=False` 时应完全复现旧口径（= 冻结 run 的分数）；
  3. 变体行 final = 锚点 × C_match − penalty（逐行可核）。

用法：
    python scripts/verify_variant_layer.py
输出：
    outputs/analysis/variant_layer_verification.json
    outputs/_tmp/variant_layer/{on,off}/candidates_ranked.csv   （对照产物）
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sirna_pipeline.common import records           # noqa: E402
from sirna_pipeline.stages.rank.ranker import run_rank  # noqa: E402

_SCORE_COLS = ("score_thermo", "score_oligo", "penalty_total", "final_score")


def load_rank_cfg() -> dict:
    y = yaml.safe_load((ROOT / "configs" / "stages.yaml").read_text(encoding="utf-8"))
    return {k: v for k, v in (y.get("ranking") or {}).items() if k != "enabled"}


def is_wt(row: dict) -> bool:
    return (str(row.get("kind", "")) == "wt"
            or str(row.get("variant_id", "")).endswith("_wt"))


def key(row: dict) -> tuple:
    return (row.get("window_id"), row.get("variant_id"))


def topn(rows: list[dict], n: int) -> list[tuple]:
    ranked = sorted((r for r in rows if r.get("final_rank")),
                    key=lambda r: int(r["final_rank"]))
    return [key(r) for r in ranked[:n]]


def jaccard(a: list, b: list) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if (sa | sb) else 1.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path,
                    default=ROOT / "outputs" / "runs" / "sfrp1_oligo_on")
    ap.add_argument("--frozen", type=Path,
                    default=ROOT / "outputs" / "runs" / "rank_reweighted" /
                            "candidates_ranked.csv")
    ap.add_argument("--work", type=Path, default=ROOT / "outputs" / "_tmp" / "variant_layer")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "outputs" / "analysis" / "variant_layer_verification.json")
    args = ap.parse_args()

    base = args.run_dir / "03_structure" / "candidates_passed.csv"
    extras = [args.run_dir / "04_thermo" / "candidates_thermo.csv",
              args.run_dir / "06_toxicity" / "candidates_toxicity.csv",
              args.run_dir / "06_oligoformer" / "candidates_oligoformer.csv"]
    cfg = load_rank_cfg()

    runs = {}
    for tag, enabled in (("on", True), ("off", False)):
        c = {**cfg, "variant_layer": {**cfg.get("variant_layer", {}), "enabled": enabled}}
        csv, meta = run_rank(base, extras, args.work / tag, cfg=c)
        rows = records.read_records(csv)
        runs[tag] = {"csv": csv, "meta": meta, "rows": rows}
        print("[%s] eligible=%d  variant_layer=%s" %
              (tag, meta["eligible"], json.dumps(meta["variant_layer"], ensure_ascii=False)[:120]))

    # 仅比较**参与排序**的行（final_rank 非空）；被前序硬过滤/特征缺失的行不参与
    on = {key(r): r for r in runs["on"]["rows"] if r.get("final_rank")}
    off = {key(r): r for r in runs["off"]["rows"] if r.get("final_rank")}
    report: dict = {"run_dir": str(args.run_dir), "frozen": str(args.frozen),
                    "n_rows_total": len(runs["on"]["rows"]),
                    "n_rows_scored": len(on)}

    # ---- 1) 开关等价性：off 必须复现旧口径（自身一致性）----
    report["off_equals_base"] = sum(
        1 for k, r in off.items() if r.get("final_score") == r.get("final_score_base"))
    report["n_rows_compared"] = len(off)

    # ---- 2) WT 逐位一致：on vs off ----
    wt_keys = [k for k, r in on.items() if is_wt(r)]
    diffs = []
    for k in wt_keys:
        for col in _SCORE_COLS:
            if on[k].get(col, "") != off[k].get(col, ""):
                diffs.append({"key": k, "col": col,
                              "on": on[k].get(col, ""), "off": off[k].get(col, "")})
    report["n_wt_rows"] = len(wt_keys)
    report["n_variant_rows"] = len(on) - len(wt_keys)
    report["wt_score_diffs"] = diffs[:20]
    report["wt_bitwise_identical"] = not diffs
    print("\nWT 行 %d 条，分数列差异 %d 处 → %s"
          % (len(wt_keys), len(diffs), "逐位一致 ✓" if not diffs else "存在差异 ✗"))

    # ---- 3) 变体行逐行核对 final = q_eff × C_match − penalty（q_eff 含改善通道）----
    bad = []
    for k, r in on.items():
        if is_wt(r):
            continue
        q_eff = records.num(r.get("q_variant_effective"))
        cm = records.num(r.get("c_match"))
        pen = records.num(r.get("penalty_total")) or 0.0
        fin = records.num(r.get("final_score"))
        if None in (q_eff, cm, fin) or abs(q_eff * cm - pen - fin) > 1e-9:
            bad.append({"key": k, "q_eff": q_eff, "c_match": cm,
                        "penalty": pen, "final": fin})
    report["variant_formula_violations"] = bad[:20]
    report["variant_formula_ok"] = not bad
    print("变体行公式核对：%s（违规 %d）"
          % ("全部通过 ✓" if not bad else "有违规 ✗", len(bad)))

    # ---- 3b) 改善通道：有多少变体超过了自己窗口的 WT ----
    wt_score = {r["window_id"]: records.num(r["final_score"])
                for r in on.values() if is_wt(r)}
    above = [k for k, r in on.items()
             if not is_wt(r) and wt_score.get(r["window_id"]) is not None
             and (records.num(r["final_score"]) or 0) > wt_score[r["window_id"]]]
    report["n_variant_above_own_wt"] = len(above)
    report["variant_above_own_wt_examples"] = [
        {"variant_id": on[k]["variant_id"], "window_id": on[k]["window_id"],
         "final": on[k]["final_score"], "wt_final": wt_score[on[k]["window_id"]],
         "thermo_delta_vs_wt": on[k].get("thermo_delta_vs_wt"),
         "c_match": on[k].get("c_match"),
         "band_class": "%s|%s" % (on[k].get("c_match_band"), on[k].get("c_match_class"))}
        for k in above[:10]]
    print("改善通道：%d 个变体超过本窗口 WT（meta 自报 %d）"
          % (len(above), (runs["on"]["meta"].get("variant_layer") or {})
             .get("n_variant_above_own_wt", -1)))

    # ---- 4) 与冻结 run 的对照 ----
    if args.frozen.exists():
        fz = {key(r): r for r in records.read_records(args.frozen) if r.get("final_rank")}
        fz_wt = [k for k, r in fz.items() if is_wt(r)]
        fz_diff = [{"key": k, "col": c,
                    "on": on[k].get(c, ""), "frozen": fz[k].get(c, "")}
                   for k in fz_wt if k in on for c in _SCORE_COLS
                   if on[k].get(c, "") != fz[k].get(c, "")]
        report["frozen_wt_compared"] = len([k for k in fz_wt if k in on])
        report["frozen_wt_score_diffs"] = fz_diff[:20]
        report["frozen_wt_bitwise_identical"] = not fz_diff
        report["jaccard_top20_vs_frozen"] = jaccard(topn(runs["on"]["rows"], 20),
                                                    topn(list(fz.values()), 20))
        report["jaccard_top50_vs_frozen"] = jaccard(topn(runs["on"]["rows"], 50),
                                                    topn(list(fz.values()), 50))
        new_top = sorted((r for r in runs["on"]["rows"] if r.get("final_rank")),
                         key=lambda r: int(r["final_rank"]))[:10]
        report["new_top10"] = [{"rank": int(r["final_rank"]), "window_id": r["window_id"],
                                "variant_id": r["variant_id"],
                                "final_score": r["final_score"],
                                "final_score_base": r.get("final_score_base", ""),
                                "c_match": r.get("c_match", ""),
                                "c_match_band": r.get("c_match_band", ""),
                                "c_match_class": r.get("c_match_class", "")}
                               for r in new_top]
        print("与冻结 run：WT 差异 %d 处 | Top20 Jaccard %.3f | Top50 Jaccard %.3f"
              % (len(fz_diff), report["jaccard_top20_vs_frozen"],
                 report["jaccard_top50_vs_frozen"]))

    # ---- 5) C_match 分布（审计：谁被压制）----
    var = [r for r in on.values() if not is_wt(r)]
    cm_vals = sorted(records.num(r.get("c_match")) or 0.0 for r in var)
    report["c_match_stats"] = {
        "n": len(var), "min": cm_vals[0] if cm_vals else None,
        "median": cm_vals[len(cm_vals) // 2] if cm_vals else None,
        "max": cm_vals[-1] if cm_vals else None,
        "by_band_class": dict(Counter("%s|%s" % (r.get("c_match_band", "?"),
                                                r.get("c_match_class", "?"))
                                      for r in var)),
    }
    print("C_match：n=%d min=%.4f median=%.4f max=%.4f"
          % (len(var), cm_vals[0], cm_vals[len(cm_vals) // 2], cm_vals[-1]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nwritten:", args.out)


if __name__ == "__main__":
    main()
