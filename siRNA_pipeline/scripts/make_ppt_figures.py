# -*- coding: utf-8 -*-
"""生成"PPT 直用"版图：图上只留数据与坐标轴，**不写标题/不写脚注**，文字留给 PPT 文本框。

设计原则（回应"文字太多 / 字叠在一起 / 不好修改"）：
  1. 无 suptitle、无长注释、无脚注 → 全部可编辑文字交给 PPT；
  2. 每张图 ≤3 个元素；标签字号 ≥11；用 constrained_layout 防重叠；
  3. 逻辑修正：
     - 错配图不再写"模型能预判"（原图 ρ=−0.17 与标题矛盾）→ 改为**实测**的
       "错配数量效应"（1 vs 2 错配保留率）与"逐位保留率"，只陈述数据；
     - 架构图不再混用两个数据集的口径（原图 82.5% 来自错配集、散点来自生产 run）
       → 全部只用生产 run，并给出同一数据上的符号一致率；
     - 识别范围图去掉与标题碰撞的"增/容/脆"文字行，改为图例色块。

产出（outputs/analysis/figures/ppt_*.png，同步到 展示图片/PPT用图）：
  ppt_p05_efficacy      四套公开数据集效能
  ppt_p06_controls      文献正对照回收
  ppt_p07_robustness    权重稳健性（β 扫描 + 200 次扰动）
  ppt_p08_architecture  分层架构证据（两信号方向不一致 + 分层前后）
  ppt_p09_recognition   逐位错配代价 + 位置×类型代价
  ppt_p09b_mismatch     实测：错配数量效应 + 逐位保留率（替代原自相矛盾那张）

用法：python scripts/make_ppt_figures.py
"""
from __future__ import annotations

import ast
import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np              # noqa: E402
import pandas as pd             # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
FIG = AN / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))
from sirna_pipeline.stages.rank import ranker  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 200,
    "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 12, "axes.labelsize": 12.5, "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5, "legend.fontsize": 11.5,
})
C_TH, C_DL, C_CO, C_RED = "#E07A5F", "#4C78A8", "#2E7D32", "#C62828"
J = lambda p: json.loads(Path(p).read_text(encoding="utf-8"))  # noqa: E731


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


# --------------------------------------------------------------------------- P05
def p05() -> None:
    t14 = J(AN / "task14_benchmark.json")
    sets = ["Hu", "Taka", "Mix", "Simone"]
    labels = ["Hu\nn=2361", "Taka\nn=702", "Mix\nn=464", "Simone\nn=322"]
    keys = [("pure_thermo", "纯热力学", C_TH), ("pure_dl", "纯深度学习", C_DL),
            ("thermo_plus_aux", "组合（本方案）", C_CO)]
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    x, w = np.arange(4), 0.26
    for i, (k, name, c) in enumerate(keys):
        v, lo, hi = [], [], []
        for ds in sets:
            m = t14["datasets"][ds]["metrics"][k]
            v.append(m["spearman"])
            lo.append(max(0.0, m["spearman"] - m["ci95_low"]))
            hi.append(max(0.0, m["ci95_high"] - m["spearman"]))
        pos = x + (i - 1) * w
        ax.bar(pos, v, w, label=name, color=c, yerr=[lo, hi], capsize=4,
               error_kw={"elinewidth": 1.2, "ecolor": "#455A64"})
        for p, val in zip(pos, v):
            ax.text(p, val + 0.025, "%.2f" % val, ha="center", fontsize=10.5)
    ax.axvspan(2.5, 3.5, color="#FFF3E0", zorder=0)
    ax.text(3, -0.155, "外部集", ha="center", fontsize=11.5, color="#E65100",
            fontweight="bold")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Spearman ρ")
    ax.set_ylim(-0.20, 0.85)
    ax.axhline(0, color="#90A4AE", lw=0.8)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    save(fig, "ppt_p05_efficacy")


# --------------------------------------------------------------------------- P06
def p06() -> None:
    d = pd.read_csv(AN / "task17_positive_controls.csv", dtype=str, encoding="utf-8-sig")
    m = d[d["matched"] == "1"].copy()
    m["pct"] = m["rank_percentile_current"].astype(float) * 100
    m = m.sort_values("pct")
    names, vals = [], []
    for _, r in m.iterrows():
        nm = str(r["source"]).split(",")[0].split(" ")[0]
        names.append("%s  %s" % (r["id"], nm))
        vals.append(r["pct"])
    fig, ax = plt.subplots(figsize=(9.6, 4.4), constrained_layout=True)
    ys = np.arange(len(vals))[::-1]
    ax.barh(ys, vals, color=[C_CO if v <= 5 else (C_DL if v <= 25 else "#90A4AE")
                             for v in vals], height=0.55)
    for y, v, (_, r) in zip(ys, vals, m.iterrows()):
        ax.text(v + 0.8, y, "第 %s 名" % r["final_rank_current"], va="center", fontsize=11)
    ax.axvline(50, color=C_RED, ls="--", lw=1.6)
    ax.text(51, float(np.mean(ys)), "随机期望 50%", color=C_RED, fontsize=11.5, va="center")
    ax.set_yticks(ys, names)
    ax.set_xlim(0, 60)
    ax.set_xlabel("在 4218 条候选中的排名百分位（%）")
    save(fig, "ppt_p06_controls")


# --------------------------------------------------------------------------- P07
def p07() -> None:
    rob = J(AN / "weight_robustness.json")
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.6), constrained_layout=True)
    ax = axes[0]
    betas, j20 = [], []
    for s in rob["scenarios"]:
        if str(s.get("scenario", "")).startswith("beta="):
            betas.append(float(s["beta"])); j20.append(s["jaccard_top20_vs_base"])
    base = rob["base_alpha_beta"][1]
    betas.append(base); j20.append(1.0)
    o = np.argsort(betas)
    ax.plot([betas[i] for i in o], [j20[i] for i in o], "o-", color=C_CO, lw=2.6, ms=8)
    ax.scatter([base], [1.0], s=220, facecolor="none", edgecolor=C_RED, lw=2.2, zorder=5)
    ax.set_xlabel("β（深度学习分权重）")
    ax.set_ylabel("Top-20 与默认口径一致性")
    ax.set_ylim(0.25, 1.06)

    ax = axes[1]
    rp = rob["random_perturbation_200"]
    names = ["Top-20", "Top-50", "全体 ρ"]
    means = [rp["jaccard_top20_mean"], rp["jaccard_top50_mean"], rp["spearman_mean"]]
    mins = [rp["jaccard_top20_min"], rp["jaccard_top50_min"], rp["spearman_min"]]
    bars = ax.bar(names, means, color=[C_CO, C_DL, "#7E57C2"], width=0.55)
    for b, mn, mean in zip(bars, mins, means):
        cx = b.get_x() + b.get_width() / 2
        ax.plot([cx, cx], [mn, mean], color="#37474F", lw=2.2)
        ax.plot([cx], [mn], marker="_", color="#37474F", ms=15, mew=2)
        ax.text(cx, mean + 0.02, "%.2f" % mean, ha="center", fontsize=11)
    ax.set_ylim(0.6, 1.10)
    ax.set_ylabel("200 次随机扰动后一致性\n（纵轴自 0.6 起）")
    save(fig, "ppt_p07_robustness")


# --------------------------------------------------------------------------- P08
def p08() -> None:
    # 直接读已发布交付表（含 final_score_base 与变体层诊断列），不依赖 outputs/_tmp 中间产物
    d = pd.read_csv(ROOT / "outputs" / "results" / "rank_final.csv", dtype=str)
    d = d[d["final_rank"].notna() & (d["final_rank"] != "")].copy()
    d["wt"] = d["variant_id"].str.endswith("_wt")
    ref = {}
    for wid, g in d.groupby("window_id"):
        w = g[g["wt"]]
        if len(w):
            ref[wid] = (float(w["score_thermo"].iloc[0]), float(w["score_oligo"].iloc[0]),
                        float(w["final_score"].iloc[0]))
    var = d[(~d["wt"]) & d["window_id"].isin(ref)].copy()
    var["dth"] = [float(t) - ref[w][0] for t, w in zip(var["score_thermo"], var["window_id"])]
    var["dol"] = [float(o) - ref[w][1] for o, w in zip(var["score_oligo"], var["window_id"])]
    above = [float(f) > ref[w][2] for f, w in zip(var["final_score"], var["window_id"])]
    agree = float(np.mean((var["dol"] > 0) == (var["dth"] > 0)))
    rho = float(var["dol"].corr(var["dth"], method="spearman"))

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), constrained_layout=True)
    ax = axes[0]
    ax.scatter(var["dol"], var["dth"], s=12, alpha=0.35, color=C_DL, edgecolor="none")
    ax.axhline(0, color="#90A4AE", lw=0.9); ax.axvline(0, color="#90A4AE", lw=0.9)
    ax.set_xlabel("Δ 深度学习分（变体 − 同窗口野生型）")
    ax.set_ylabel("Δ 热力学分（变体 − 同窗口野生型）")
    ax.text(0.03, 0.95, "ρ = %.2f\n方向一致率 %.2f" % (rho, agree), transform=ax.transAxes,
            va="top", fontsize=12)
    ax = axes[1]
    xs = [float(v) for v in var["final_score_base"]]
    ys2 = [float(v) for v in var["final_score"]]
    ax.scatter(xs, ys2, s=12, alpha=0.35, color="#90A4AE", edgecolor="none")
    ax.scatter([x for x, a in zip(xs, above) if a], [y for y, a in zip(ys2, above) if a],
               s=18, alpha=0.85, color=C_CO, edgecolor="none")
    ax.scatter([ref[w][2] for w in ref], [ref[w][2] for w in ref], marker="D", s=55,
               color="#37474F", zorder=4, edgecolor="white", linewidth=0.6)
    lim = (0.15, 0.9)
    ax.plot(lim, lim, "--", color="#546E7A", lw=1.1)
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("改造前终分")
    ax.set_ylabel("改造后终分")
    ax.text(0.03, 0.95, "绿点 = 越过本窗口野生型的变体（%d 个）\n黑菱形 = 野生型（逐位不变）"
            % sum(above), transform=ax.transAxes, va="top", fontsize=11.5)
    save(fig, "ppt_p08_architecture")


# --------------------------------------------------------------------------- P09
def p09() -> None:
    lam = dict(ranker.DEFAULT_RANK_CFG["variant_layer"]["lambda_bands"])
    pos = list(range(1, 20))
    lams = [float(lam.get(ranker.variant_band(p), 0.0)) for p in pos]
    classes = [("GU", 0.5), ("transition", None), ("transversion", None)]
    cost = np.zeros((3, 19))
    for j, p in enumerate(pos):
        for i, (cls, fixed) in enumerate(classes):
            tf = fixed if fixed is not None else (
                (1.0 if cls == "transition" else 0.7) if p >= 12
                else (0.85 if cls == "transition" else 1.0))
            cost[i, j] = lams[p - 1] * tf

    fig, axes = plt.subplots(2, 1, figsize=(12.2, 7.2), constrained_layout=True,
                             gridspec_kw={"height_ratios": [1, 1.05]})
    ax = axes[0]
    ax.bar(pos, lams, color=C_RED, width=0.6)
    ax.set_xticks(pos)
    ax.set_ylabel("λ（错配代价系数）")
    ax.set_ylim(0, 0.38)
    ax.text(0.99, 0.92, "灰柱 = 无数据、取文献先验", transform=ax.transAxes,
            ha="right", va="top", fontsize=11.5, color="#546E7A")
    ax.bar([p for p in pos if ranker.variant_band(p) in ("g1", "p3")],
           [lams[p - 1] for p in pos if ranker.variant_band(p) in ("g1", "p3")],
           color="#CFD8DC", width=0.6)

    ax = axes[1]
    im = ax.imshow(cost, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=0.35)
    ax.set_yticks(range(3), [c for c, _ in classes])
    ax.set_xticks(range(19), pos)
    ax.set_xlabel("引导链位置（1-based，5′ 端起算）")
    ax.grid(False)
    fig.colorbar(im, ax=ax, pad=0.01).set_label("代价 1 − C_match")
    save(fig, "ppt_p09_recognition")


# -------------------------------------------------------------------------- P09b
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
    d["pos"] = d["verified_positions"].apply(_pl)
    wt = d[d["mc"] == 0].set_index("group_id")["label"].astype(float).to_dict()
    m = d[(d["mc"] > 0) & d["group_id"].map(wt).notna()].copy()
    m["wtl"] = m["group_id"].map(wt)
    m["ret"] = m["label"].astype(float) / m["wtl"].astype(float)
    m["n_mm"] = m["pos"].apply(len)

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.8), constrained_layout=True)
    ax = axes[0]
    g = m.groupby("n_mm")["ret"]
    ns = [int(g.size().get(k, 0)) for k in (1, 2)]
    means = [float(g.mean().get(k, np.nan)) for k in (1, 2)]
    sem = [float(g.std().get(k, np.nan) / max(1, ns[i]) ** 0.5) for i, k in enumerate((1, 2))]
    ax.bar(["单错配", "双错配"], means, yerr=sem, capsize=6, color=[C_DL, C_TH], width=0.5,
           error_kw={"elinewidth": 1.2, "ecolor": "#45524F"})
    for i, (mu, n) in enumerate(zip(means, ns)):
        ax.text(i, mu + 0.03, "%.2f\nn=%d" % (mu, n), ha="center", fontsize=11.5)
    ax.set_ylabel("实测效力保留率（label / 同家族 WT label）")
    ax.set_ylim(0, 1.05)

    ax = axes[1]
    s1 = m[m["n_mm"] == 1].copy()
    s1["p"] = s1["pos"].apply(lambda x: int(x[0]))
    cnt = s1.groupby("p")["ret"]
    ps = sorted(cnt.groups)
    ax.bar([str(p) for p in ps], [float(cnt.mean()[p]) for p in ps], color=C_DL, width=0.62)
    for i, p in enumerate(ps):
        ax.text(i, float(cnt.mean()[p]) + 0.02, "n=%d" % int(cnt.size()[p]),
                ha="center", fontsize=10)
    ax.axhline(1.0, color="#90A4AE", lw=0.9)
    ax.set_xlabel("错配位置（引导链，1-based）")
    ax.set_ylabel("实测保留率")
    ax.set_ylim(0, 1.12)
    save(fig, "ppt_p09b_mismatch")


def main() -> int:
    for fn in (p05, p06, p07, p08, p09, p09b):
        fn()
    # 同步到 展示图片/ 与 PPT用图/
    dest_dirs = [ROOT.parent / "展示图片", ROOT.parent / "PPT用图"]
    mapping = {
        "ppt_p05_efficacy": "P05_性能证据1_四基准集效能.png",
        "ppt_p06_controls": "P06_性能证据2_正对照回收.png",
        "ppt_p07_robustness": "P07_稳健性与规范性.png",
        "ppt_p08_architecture": "P08_架构设计_分层.png",
        "ppt_p09_recognition": "P09_识别范围与安全性.png",
        "ppt_p09b_mismatch": "P09b_错配实测影响.png",
    }
    for d in dest_dirs:
        if not d.exists():
            continue
        for stem, name in mapping.items():
            shutil.copyfile(FIG / ("%s.png" % stem), d / name)
        print("synced ->", d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
