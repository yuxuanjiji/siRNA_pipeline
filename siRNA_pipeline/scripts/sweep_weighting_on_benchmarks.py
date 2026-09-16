# -*- coding: utf-8 -*-
"""在四个基准集上做"特征权重 × α/β"网格实验，为 ranker 默认值提供证据。

口径与 task14 完全一致：
  * 热力分 = 四特征 min-max（按 direction 统一越大越好）后加权平均（权重可配）；
  * DL 分 = OligoFormer predictions（pred）在同数据集内 min-max；
  * 终分 = (1-β)·热力分 + β·DL 分（α=1-β，与 ranker 的 alpha/beta 同义）；
  * 指标 = 与实测 label 的 Spearman(per-dataset)、跨集 pooled Spearman、
    以及最终候选配置的 bootstrap 95% CI 与置换 p。

用法：python scripts/sweep_weighting_on_benchmarks.py
输出：outputs/analysis/weighting_sweep.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import run_experiments as R  # noqa: E402

FEATS = ("Hu", "Taka", "Mix", "Simone")
BETAS = (0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0)
# (预设名, {feat_mfe, feat_ddg_ends, feat_dG_duplex, feat_dG_seed})
PRESETS = {
    "current_1.0_1.0_1.2_0.8": {"feat_mfe": 1.0, "feat_ddg_ends": 1.0, "feat_dG_duplex": 1.2, "feat_dG_seed": 0.8},
    "equal_1_1_1_1": {"feat_mfe": 1.0, "feat_ddg_ends": 1.0, "feat_dG_duplex": 1.0, "feat_dG_seed": 1.0},
    "mfe_dominant_1.8_0.4_0.4_0.2": {"feat_mfe": 1.8, "feat_ddg_ends": 0.4, "feat_dG_duplex": 0.4, "feat_dG_seed": 0.2},
    "mfe_heavy_2.0_0.5_1.0_0.0": {"feat_mfe": 2.0, "feat_ddg_ends": 0.5, "feat_dG_duplex": 1.0, "feat_dG_seed": 0.0},
    "no_seed_1.0_1.0_1.2_0.0": {"feat_mfe": 1.0, "feat_ddg_ends": 1.0, "feat_dG_duplex": 1.2, "feat_dG_seed": 0.0},
    "duplex_heavy_0.5_0.5_2.0_0.0": {"feat_mfe": 0.5, "feat_ddg_ends": 0.5, "feat_dG_duplex": 2.0, "feat_dG_seed": 0.0},
    "mfe_only": {"feat_mfe": 1.0, "feat_ddg_ends": 0.0, "feat_dG_duplex": 0.0, "feat_dG_seed": 0.0},
}
SOURCE = {"feat_mfe": ("MFE_guide",), "feat_ddg_ends": ("delta_deltaG_ends",),
          "feat_dG_duplex": ("dG_total",), "feat_dG_seed": ("dG_seed", "dG_seed_only_canonical")}
DIRECTION = {"feat_mfe": 1, "feat_ddg_ends": 1, "feat_dG_duplex": -1, "feat_dG_seed": -1}


def rank_corr(a, b) -> float | None:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 3 or len(np.unique(a)) < 2 or len(np.unique(b)) < 2:
        return None
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den else None


def load_dataset(name: str):
    feats = {r_["siRNA"] + "\u0001" + r_["mRNA"]: r_ for r_ in R._read(ROOT / "outputs" / "analysis" / "features" / f"{name}.csv")}
    preds = {p["siRNA"] + "\u0001" + p["mRNA"]: p for p in R._read(ROOT / "OligoFormer部分" / f"{name}_predictions.csv")}
    keys = [k for k in feats if k in preds]
    rows = []
    for k in keys:
        f, p = feats[k], preds[k]
        rows.append({"label": float(f["label"]), "dl": float(p["pred"]),
                     **{name_: f.get(name_ if name_ != "feat_dG_seed" else "dG_seed", "") for name_ in FEATS}})
    # 直接取四特征列（dG_seed 若空则退 dG_seed_only_canonical）
    for r, k in zip(rows, keys):
        f = feats[k]
        for feat, cands in SOURCE.items():
            val = ""
            for c in cands:
                v = f.get(c, "")
                if str(v).strip() not in ("", "nan"):
                    val = float(v); break
            r[feat] = val
    return rows


def thermo_scores(rows, preset):
    norms = {}
    for feat, w in preset.items():
        if w == 0:
            norms[feat] = np.zeros(len(rows))
            continue
        vals = np.array([r[feat] for r in rows], dtype=float)
        if DIRECTION[feat] < 0:
            vals = -vals
        lo, hi = vals.min(), vals.max()
        norms[feat] = (vals - lo) / (hi - lo) if hi > lo else np.full(len(vals), 0.5)
    wsum = sum(preset.values()) or 1.0
    return sum(norms[f] * w / wsum for f, w in preset.items())


def main() -> None:
    data = {n: load_dataset(n) for n in FEATS}
    for n, rows in data.items():
        print("%-7s rows=%d" % (n, len(rows)))
    result = {"datasets": {n: len(r) for n, r in data.items()}, "grid": []}
    for pname, preset in PRESETS.items():
        th = {n: thermo_scores(rows, preset) for n, rows in data.items()}
        dl = {}
        for n, rows in data.items():
            v = np.array([r["dl"] for r in rows])
            lo, hi = v.min(), v.max()
            dl[n] = (v - lo) / (hi - lo) if hi > lo else np.full(len(v), 0.5)
        for beta in BETAS:
            per, labels_all, scores_all = {}, [], []
            for n, rows in data.items():
                y = np.array([r["label"] for r in rows])
                s = (1 - beta) * th[n] + beta * dl[n]
                per[n] = rank_corr(y, s)
                labels_all.append((y - y.mean()) / (y.std() + 1e-12))
                scores_all.append((s - s.mean()) / (s.std() + 1e-12))
            pooled = rank_corr(np.concatenate(labels_all), np.concatenate(scores_all))
            vals = [v for v in per.values() if v is not None]
            entry = {"preset": pname, "beta": beta, "alpha": round(1 - beta, 3),
                     "per_dataset": {k: (round(v, 4) if v is not None else None) for k, v in per.items()},
                     "mean_spearman": round(float(np.mean(vals)), 4),
                     "worst_spearman": round(float(np.min(vals)), 4),
                     "pooled_spearman": round(pooled, 4) if pooled is not None else None}
            result["grid"].append(entry)

    grid = result["grid"]
    by_mean = sorted(grid, key=lambda e: -e["mean_spearman"])[:6]
    by_worst = sorted(grid, key=lambda e: -e["worst_spearman"])[:4]
    print("\n== 按 mean Spearman 前 6 ==")
    for e in by_mean:
        print("  %-28s β=%.1f mean=%.3f worst=%.3f pooled=%.3f | %s"
              % (e["preset"], e["beta"], e["mean_spearman"], e["worst_spearman"],
                 e["pooled_spearman"] or float("nan"),
                 " ".join("%s=%.3f" % (k, v) for k, v in e["per_dataset"].items() if v is not None)))
    print("\n== 按 worst-dataset 前 4（稳健性）==")
    for e in by_worst:
        print("  %-28s β=%.1f mean=%.3f worst=%.3f" % (e["preset"], e["beta"], e["mean_spearman"], e["worst_spearman"]))

    # 纯 DL 与当前默认对照
    cur = [e for e in grid if e["preset"] == "current_1.0_1.0_1.2_0.8" and abs(e["beta"] - 0.2) < 1e-9][0]
    pure = [e for e in grid if e["preset"] == "current_1.0_1.0_1.2_0.8" and abs(e["beta"] - 1.0) < 1e-9][0]
    result["reference"] = {"current_default_alpha0.8_beta0.2": cur, "pure_dl_beta1.0": pure}
    (ROOT / "outputs" / "analysis" / "weighting_sweep.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nwritten: outputs/analysis/weighting_sweep.json")


if __name__ == "__main__":
    main()
