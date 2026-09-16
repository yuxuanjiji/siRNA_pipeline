# -*- coding: utf-8 -*-
"""权重/α-β 稳健性扫描（任务2）：在真实 SFRP1 候选集上扰动权重与 α/β，看交付榜稳定性。

基准 = outputs/runs/rank_reweighted/candidates_ranked.csv（新默认 α0.4/β0.6，已验证）。
扰动：
  A) 单权重 ±20%（feat_mfe / feat_ddg_ends / feat_dG_duplex 逐个）
  B) 同时随机扰动：每个权重 U(0.8,1.2)×，200 次抽样
  C) α/β 扫描：β∈{0.4,0.5,0.6,0.7,0.8}（α=1-β）
  D) 参考：把 feat_dG_seed 权重恢复到 0.2 / 0.8
指标：top-20/top-50 Jaccard、rank1 是否变化、全表 Spearman(与基准排名)。
输出：outputs/analysis/weight_robustness.json
"""
from __future__ import annotations

import csv
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "outputs" / "runs" / "rank_reweighted" / "candidates_ranked.csv"
OUT = ROOT / "outputs" / "analysis" / "weight_robustness.json"
BASE_W = {"feat_mfe": 1.0, "feat_ddg_ends": 1.0, "feat_dG_duplex": 1.2, "feat_dG_seed": 0.0}
DIR = {"feat_mfe": 1, "feat_ddg_ends": 1, "feat_dG_duplex": -1, "feat_dG_seed": -1}
COL = {"feat_mfe": "feat_mfe", "feat_ddg_ends": "feat_ddg_ends",
       "feat_dG_duplex": "feat_dG_duplex", "feat_dG_seed": "feat_dG_seed"}
ALPHA, BETA = 0.4, 0.6


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


rows = []
with open(SRC, encoding="utf-8-sig", newline="") as fh:
    for r in csv.DictReader(fh):
        if not str(r.get("final_rank", "")).strip():
            continue
        feats = {}
        ok = True
        for k, c in COL.items():
            v = num(r.get(c))
            if v is None:
                ok = False
                break
            feats[k] = v
        ol = num(r.get("score_oligo"))
        pen = num(r.get("penalty_total")) or 0.0
        if not ok or ol is None:
            continue
        rows.append({"id": r["variant_id"], "win": r["window_id"], "kind": r["kind"],
                     "feats": feats, "oligo": ol, "pen": pen})

print("eligible rows =", len(rows))

# 归一化（与 ranker 一致：方向统一后 min-max，在 eligible 全集上）
norm = {}
for k in BASE_W:
    vals = [(r["feats"][k] if DIR[k] > 0 else -r["feats"][k]) for r in rows]
    lo, hi = min(vals), max(vals)
    span = hi - lo
    norm[k] = [((v - lo) / span if span else 0.5) for v in vals]


def score(weights, beta):
    a = 1.0 - beta
    wsum = sum(weights.values()) or 1.0
    out = []
    for i, r in enumerate(rows):
        th = sum(weights[k] * norm[k][i] for k in weights) / wsum
        out.append(a * th + beta * r["oligo"] - r["pen"])
    return out


def ranking(scores):
    order = sorted(range(len(rows)), key=lambda i: -scores[i])
    return [rows[i]["id"] for i in order]


def rank_corr(a_ids, b_ids):
    ra = {v: i for i, v in enumerate(a_ids)}
    rb = {v: i for i, v in enumerate(b_ids)}
    n = len(a_ids)
    x = [ra[v] for v in a_ids]
    y = [rb[v] for v in a_ids]
    mx, my = sum(x) / n, sum(y) / n
    num_ = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    den = (sum((xi - mx) ** 2 for xi in x) * sum((yi - my) ** 2 for yi in y)) ** 0.5
    return num_ / den if den else 0.0


def jac(a, b, k):
    A, B = set(a[:k]), set(b[:k])
    return len(A & B) / len(A | B)


base_scores = score(BASE_W, BETA)
base_ids = ranking(base_scores)
result = {"n_rows": len(rows), "base_weights": BASE_W, "base_alpha_beta": [ALPHA, BETA],
          "base_top1": base_ids[0], "scenarios": []}


def record(name, weights, beta):
    ids = ranking(score(weights, beta))
    entry = {"scenario": name, "weights": weights, "beta": beta, "alpha": round(1 - beta, 3),
             "jaccard_top20_vs_base": round(jac(ids, base_ids, 20), 4),
             "jaccard_top50_vs_base": round(jac(ids, base_ids, 50), 4),
             "spearman_vs_base": round(rank_corr(base_ids, ids), 4),
             "top1": ids[0], "rank1_changed": ids[0] != base_ids[0]}
    result["scenarios"].append(entry)
    print("%-34s β=%.1f J20=%.3f J50=%.3f ρ=%.3f top1=%s%s"
          % (name, beta, entry["jaccard_top20_vs_base"], entry["jaccard_top50_vs_base"],
             entry["spearman_vs_base"], entry["top1"], "  ←rank1变了" if entry["rank1_changed"] else ""))


for k in BASE_W:
    if BASE_W[k] == 0:
        continue
    for f in (0.8, 1.2):
        w = dict(BASE_W); w[k] = round(BASE_W[k] * f, 4)
        record("%s×%.1f" % (k, f), w, BETA)

rng = random.Random(42)
j20, j50, sp, flip = [], [], [], 0
for _ in range(200):
    w = {k: (round(BASE_W[k] * rng.uniform(0.8, 1.2), 4) if BASE_W[k] > 0 else 0.0) for k in BASE_W}
    ids = ranking(score(w, BETA))
    j20.append(jac(ids, base_ids, 20)); j50.append(jac(ids, base_ids, 50))
    sp.append(rank_corr(base_ids, ids)); flip += (ids[0] != base_ids[0])
result["random_perturbation_200"] = {
    "jaccard_top20_mean": round(sum(j20) / len(j20), 4),
    "jaccard_top20_min": round(min(j20), 4),
    "jaccard_top50_mean": round(sum(j50) / len(j50), 4),
    "jaccard_top50_min": round(min(j50), 4),
    "spearman_mean": round(sum(sp) / len(sp), 4),
    "spearman_min": round(min(sp), 4),
    "rank1_changed_ratio": round(flip / len(j20), 3),
}
print("随机扰动(200):", result["random_perturbation_200"])

for beta in (0.4, 0.5, 0.6, 0.7, 0.8):
    if abs(beta - BETA) > 1e-9:
        record("beta=%.1f" % beta, BASE_W, beta)
record("seed权重恢复0.2", {**BASE_W, "feat_dG_seed": 0.2}, BETA)
record("seed权重恢复0.8（旧默认）", {**BASE_W, "feat_dG_seed": 0.8}, BETA)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print("written:", OUT)
