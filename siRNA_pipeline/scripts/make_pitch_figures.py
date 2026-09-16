# -*- coding: utf-8 -*-
"""答辩极简版图：一张图 = 一个信息点，图上几乎不留字（数值/注解/轴标题全部交给 PPT 文本框）。

与 make_ppt_figures.py（带标注版）的区别：
  * 去掉所有 n=、数值标签、ρ 标注、数据来源行、轴标题、注释性文字；
  * 只保留：坐标轴刻度、图例（≤3 项）、以及 1–2 个绝对必要的数字；
  * 字号放大到答辩可读（刻度 14、图例 14、大数字 20+）。

产出（outputs/analysis/figures/pitch_*.png，同步到 PPT用图/）：
  P05 四套数据集效能   P06 正对照回收   P07 β 敏感性
  P08 分层前后         P09 识别范围热图 P09b 错配数量效应

用法：python scripts/make_pitch_figures.py
"""
from __future__ import annotations

import ast
import json
import math
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np              # noqa: E402
import pandas as pd             # noqa: E402

import sys
ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
FIG = AN / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))
from sirna_pipeline.stages.rank import ranker  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 200,
    "axes.grid": True, "grid.alpha": 0.22,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 14, "axes.labelsize": 14, "xtick.labelsize": 14, "ytick.labelsize": 14,
    "legend.fontsize": 14,
})
C_TH, C_DL, C_CO, C_RED = "#E07A5F", "#4C78A8", "#2E7D32", "#C62828"
J = lambda p: json.loads(Path(p).read_text(encoding="utf-8"))  # noqa: E731


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


def p05() -> None:
    t14 = J(AN / "task14_benchmark.json")
    sets = ["Hu", "Taka", "Mix", "Simone"]
    keys = [("pure_thermo", "热力学", C_TH), ("pure_dl", "深度学习", C_DL),
            ("thermo_plus_aux", "本方案", C_CO)]
    fig, ax = plt.subplots(figsize=(9.6, 5.2), constrained_layout=True)
    x, w = np.arange(4), 0.26
    for i, (k, name, c) in enumerate(keys):
        v = [t14["datasets"][ds]["metrics"][k]["spearman"] for ds in sets]
        ax.bar(x + (i - 1) * w, v, w, label=name, color=c)
    ax.axvspan(2.5, 3.5, color="#FFF1DC", zorder=0)
    ax.set_xticks(x, sets)
    ax.axhline(0, color="#90A4AE", lw=0.8)
    ax.set_ylim(-0.10, 0.80)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.set_ylabel("ρ")
    save(fig, "pitch_p05")


def p06() -> None:
    d = pd.read_csv(AN / "task17_positive_controls.csv", dtype=str, encoding="utf-8-sig")
    m = d[d["matched"] == "1"].copy()
    m["pct"] = m["rank_percentile_current"].astype(float) * 100
    m = m.sort_values("pct")
    fig, ax = plt.subplots(figsize=(9.6, 4.2), constrained_layout=True)
    ys = np.arange(len(m))[::-1]
    ax.barh(ys, m["pct"].tolist(), color=C_CO, height=0.5)
    ax.axvline(50, color=C_RED, ls="--", lw=2)
    ax.set_yticks(ys, m["id"].tolist())
    ax.set_xlim(0, 55)
    ax.set_xlabel("排名百分位（%）")
    save(fig, "pitch_p06")


def p07() -> None:
    rob = J(AN / "weight_robustness.json")
    betas, j20 = [], []
    for s in rob["scenarios"]:
        if str(s.get("scenario", "")).startswith("beta="):
            betas.append(float(s["beta"])); j20.append(s["jaccard_top20_vs_base"])
    betas.append(rob["base_alpha_beta"][1]); j20.append(1.0)
    o = np.argsort(betas)
    fig, ax = plt.subplots(figsize=(8.4, 5.0), constrained_layout=True)
    ax.plot([betas[i] for i in o], [j20[i] for i in o], "o-", color=C_CO, lw=3, ms=10)
    ax.scatter([0.6], [1.0], s=260, facecolor="none", edgecolor=C_RED, lw=2.6, zorder=5)
    ax.set_xlabel("β")
    ax.set_ylabel("Top-20 一致性")
    ax.set_ylim(0.25, 1.06)
    save(fig, "pitch_p07")


def p08() -> None:
    d = pd.read_csv(ROOT / "outputs" / "results" / "rank_final.csv", dtype=str)
    d = d[d["final_rank"].notna() & (d["final_rank"] != "")].copy()
    d["wt"] = d["variant_id"].str.endswith("_wt")
    ref = {}
    for wid, g in d.groupby("window_id"):
        w = g[g["wt"]]
        if len(w):
            ref[wid] = float(w["final_score"].iloc[0])
    v = d[(~d["wt"]) & d["window_id"].isin(ref)].copy()
    xs = [float(a) for a in v["final_score_base"]]
    ys = [float(b) for b in v["final_score"]]
    up = [b > ref[w] for b, w in zip(ys, v["window_id"])]
    fig, ax = plt.subplots(figsize=(7.6, 6.4), constrained_layout=True)
    ax.scatter([a for a, u in zip(xs, up) if not u], [b for b, u in zip(ys, up) if not u],
               s=16, alpha=0.35, color="#90A4AE", edgecolor="none")
    ax.scatter([a for a, u in zip(xs, up) if u], [b for b, u in zip(ys, up) if u],
               s=24, alpha=0.9, color=C_CO, edgecolor="none")
    ax.plot([0.15, 0.9], [0.15, 0.9], "--", color="#546E7A", lw=1.2)
    ax.set_xlim(0.15, 0.9); ax.set_ylim(0.15, 0.9)
    ax.set_xlabel("改造前终分"); ax.set_ylabel("改造后终分")
    save(fig, "pitch_p08")


def p09() -> None:
    lam = dict(ranker.DEFAULT_RANK_CFG["variant_layer"]["lambda_bands"])
    cost = np.zeros((3, 19))
    for j, p in enumerate(range(1, 20)):
        l = float(lam.get(ranker.variant_band(p), 0.0))
        for i, (cls, fixed) in enumerate([("GU", 0.5), ("transition", None),
                                          ("transversion", None)]):
            tf = fixed if fixed is not None else (
                (1.0 if cls == "transition" else 0.7) if p >= 12
                else (0.85 if cls == "transition" else 1.0))
            cost[i, j] = l * tf
    fig, ax = plt.subplots(figsize=(11.0, 4.2), constrained_layout=True)
    im = ax.imshow(cost, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=0.35)
    ax.set_yticks(range(3), ["GU", "转换", "颠换"])
    ax.set_xticks(range(19), range(1, 20))
    ax.set_xlabel("引导链位置")
    ax.grid(False)
    fig.colorbar(im, ax=ax, pad=0.012)
    save(fig, "pitch_p09")


def _pl(x):
    if x in ("[]", "", None) or (isinstance(x, float) and math.isnan(x)):
        return []
    try:
        return ast.literal_eval(x)
    except (ValueError, SyntaxError):
        return []


def p09b() -> None:
    d = pd.read_csv(AN / "mismatch" / "mismatch_curated_v2.csv", dtype=str,
                    encoding="utf-8-sig")
    d = d[d["model_eligible_v2"] == "True"].copy()
    d["mc"] = pd.to_numeric(d["mismatch_count"], errors="coerce")
    old = pd.to_numeric(d["declared_mismatch_count"], errors="coerce")
    wt = d[d["mc"] == 0].set_index("group_id")["label"].astype(float).to_dict()
    m = d[(d["mc"] > 0) & d["group_id"].map(wt).notna()].copy()
    m["ret"] = m["label"].astype(float) / m["group_id"].map(wt).astype(float)
    m["k"] = m["verified_positions"].apply(lambda x: len(_pl(x))).clip(1, 2)
    g = m.groupby("k")["ret"]
    fig, ax = plt.subplots(figsize=(6.6, 5.4), constrained_layout=True)
    mus = [float(g.mean()[k]) for k in (1, 2)]
    ax.bar(["1 个错配", "2 个错配"], mus, color=[C_DL, C_TH], width=0.5)
    for i, mu in enumerate(mus):
        ax.text(i, mu + 0.03, "%.2f" % mu, ha="center", fontsize=22, fontweight="bold")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("效力保留率")
    save(fig, "pitch_p09b")


def main() -> int:
    for fn in (p05, p06, p07, p08, p09, p09b):
        fn()
    dst = ROOT.parent / "PPT用图"
    if dst.exists():
        for stem, name in {"pitch_p05": "P05_四基准集效能.png",
                           "pitch_p06": "P06_正对照回收.png",
                           "pitch_p07": "P07_β敏感性.png",
                           "pitch_p08": "P08_分层前后.png",
                           "pitch_p09": "P09_识别范围.png",
                           "pitch_p09b": "P09b_错配数量效应.png"}.items():
            shutil.copyfile(FIG / ("%s.png" % stem), dst / name)
        print("synced ->", dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
