# -*- coding: utf-8 -*-
"""验证图表生成（task14–16 · 以 Spearman 为核心，含 95% CI 与显著性标注）。

数据源（由 scripts/run_experiments.py 产出，勿手工改数）：
    outputs/analysis/task14_benchmark.json     # 纯DL / 纯热力(ranker口径) / 热力+DL
    outputs/analysis/task16_ablation.json      # 错配集：逐特征 Spearman、权重预设、MFE 权重扫描

输出（outputs/analysis/figures/）：
    validation_task14_spearman.png/.pdf        # 三基准 × 三种方法（误差棒=95% CI）
    validation_task15_scores_spearman.png/.pdf # 错配集：旧口径 vs 新口径对比
    validation_feature_spearman.png/.pdf       # 逐特征 Spearman
    validation_mfe_weight_sweep.png/.pdf       # MFE 权重扫描（CV 与全样本）
    validation_dashboard.png/.pdf              # 2×2 汇总

用法：python scripts/plot_validation_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AN = Path(__file__).resolve().parents[1] / "outputs" / "analysis"
FIG = AN / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150

C_DL, C_THERMO, C_COMBO = "#4C78A8", "#E07A5F", "#4C956C"
C_OLD, C_NEW, C_WARN = "#B0BEC5", "#2E7D32", "#C62828"


def load(name: str) -> dict:
    return json.loads((AN / name).read_text(encoding="utf-8"))


def err(stats: dict):
    lo, hi, rho = stats.get("ci95_low"), stats.get("ci95_high"), stats.get("spearman")
    if lo is None or hi is None or rho is None:
        return None
    return [[max(0.0, rho - lo)], [max(0.0, hi - rho)]]


def mark(stats: dict) -> str:
    p = stats.get("p_perm")
    if p is None:
        return ""
    return "**" if p < 0.01 else ("*" if p < 0.05 else "ns")


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


# ---------------------------------------------------------------------------
# 图 1：task14 三基准（纯DL / 纯热力·ranker口径 / 热力+DL）
# ---------------------------------------------------------------------------
def fig_task14(ax=None):
    data = load("task14_benchmark.json")
    ab = None
    for ds in data.get("datasets", {}).values():
        if ds.get("combo_alpha_beta"):
            ab = ds["combo_alpha_beta"]
            break
    combo_label = ("热力+DL(α%.1f/β%.1f)" % (ab[0], ab[1])) if ab else "热力+DL(ranker 默认)"
    methods = [("pure_dl", "纯 DL", C_DL), ("pure_thermo", "纯热力(ranker口径)", C_THERMO),
               ("thermo_plus_aux", combo_label, C_COMBO)]
    names = [k for k in ("Hu", "Taka", "Mix", "Simone") if k in data["datasets"]]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
    width = 0.26
    xs = range(len(names))
    for i, (key, label, color) in enumerate(methods):
        vals, errs, labels = [], [], []
        for n in names:
            s = data["datasets"][n]["metrics"].get(key) or {}
            vals.append(s.get("spearman") if s.get("spearman") is not None else 0.0)
            errs.append(err(s) or [[0.0], [0.0]])
            labels.append(mark(s))
        pos = [x + (i - 1) * width for x in xs]
        bars = ax.bar(pos, vals, width, label=label, color=color, edgecolor="#263238", linewidth=0.6)
        for b, e, lab in zip(bars, errs, labels):
            if e != [[0.0], [0.0]]:
                ax.errorbar(b.get_x() + b.get_width() / 2, b.get_height(),
                            yerr=e, fmt="none", ecolor="#263238", capsize=3, linewidth=1.0)
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02, lab,
                    ha="center", va="bottom", fontsize=8)
    ax.axhline(0, color="#333", linewidth=0.8)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(["%s\n(n=%d)" % (n, data["datasets"][n]["rows"]) for n in names])
    ax.set_ylabel("Spearman $\\rho$（与实验效率）")
    ax.set_title("task14 三基准：三种口径的排序一致性（误差棒=95% CI）", fontweight="bold")
    ax.legend(fontsize=8, ncol=3, loc="upper left", frameon=False)
    ax.grid(axis="y", alpha=0.25)
    if own:
        save(fig, "validation_task14_spearman")
    return ax


# ---------------------------------------------------------------------------
# 图 2：task15 错配集 —— 旧口径 vs 新口径
# ---------------------------------------------------------------------------
def fig_task15(ax=None):
    abl = load("task16_ablation.json")
    fs = abl.get("feature_spearman", {})
    wp = abl.get("weight_presets", {})
    rows = [
        ("纯 OligoFormer（旧）", fs.get("dl_score"), C_OLD),
        ("纯热力·等权（旧）", fs.get("S_thermo"), C_OLD),
        ("原综合 0.7热+0.3奥（旧）", fs.get("S_combo"), C_OLD),
        ("纯热力·ranker 口径", fs.get("thermo_score"), C_THERMO),
        ("MFE 单特征", fs.get("MFE_guide"), C_NEW),
        ("MFE 主导 1.8/0.4/0.4/0.2", wp.get("mfe_dominant"), C_NEW),
        ("no_seed 1.4/0.6/0.6", wp.get("no_seed"), C_NEW),
        ("dG_seed 单特征", fs.get("dG_seed"), C_WARN),
    ]
    rows = [r for r in rows if r[1] and r[1].get("spearman") is not None]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
    labels = [r[0] for r in rows]
    vals = [r[1]["spearman"] for r in rows]
    errs = [err(r[1]) or [[0.0], [0.0]] for r in rows]
    colors = [r[2] for r in rows]
    y = range(len(rows))
    ax.barh(list(y), vals, color=colors, edgecolor="#263238", linewidth=0.6, height=0.62)
    for i, (v, e, r) in enumerate(zip(vals, errs, rows)):
        if e != [[0.0], [0.0]]:
            ax.errorbar(v, i, xerr=e, fmt="none", ecolor="#263238", capsize=3, linewidth=1.0)
        ax.text(v + 0.012, i, "%.3f %s" % (v, mark(r[1])), va="center", fontsize=8)
    ax.axvline(0, color="#333", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    n = abl.get("rows", 0)
    ax.set_xlabel("Spearman $\\rho$（预测分 vs 实验效率，n=%d）" % n)
    ax.set_title("task15 错配集：旧口径 vs 新口径（误差棒=95% CI）", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    if own:
        save(fig, "validation_task15_scores_spearman")
    return ax


# ---------------------------------------------------------------------------
# 图 3：逐特征 Spearman
# ---------------------------------------------------------------------------
def fig_features(ax=None):
    abl = load("task16_ablation.json")
    fs = abl.get("feature_spearman", {})
    zh = {"MFE_guide": "MFE 引导链自折叠", "dl_score": "OligoFormer",
          "dG_total": "双链 ΔG", "dG_seed": "seed 结合能（现有定义）",
          "delta_deltaG_ends": "末端 ΔΔG", "GC_content": "GC 含量",
          "thermo_score": "热力·ranker 口径", "S_thermo": "热力·等权（旧）",
          "S_combo": "原综合（旧）"}
    items = [(zh.get(k, k), v) for k, v in fs.items() if v.get("spearman") is not None]
    items.sort(key=lambda t: t[1]["spearman"])
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
    vals = [v["spearman"] for _, v in items]
    errs = [err(v) or [[0.0], [0.0]] for _, v in items]
    colors = [C_NEW if v.get("p_perm") is not None and v["p_perm"] < 0.05 else C_OLD
              for _, v in items]
    y = list(range(len(items)))
    ax.barh(y, vals, color=colors, edgecolor="#263238", linewidth=0.6, height=0.62)
    for i, (v, e) in enumerate(zip(vals, errs)):
        if e != [[0.0], [0.0]]:
            ax.errorbar(v, i, xerr=e, fmt="none", ecolor="#263238", capsize=3, linewidth=1.0)
        ax.text(v + (0.012 if v >= 0 else -0.012), i,
                "%.3f %s" % (v, mark(items[i][1])),
                va="center", ha="left" if v >= 0 else "right", fontsize=8)
    ax.axvline(0, color="#333", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([t[0] for t in items], fontsize=9)
    ax.set_xlabel("Spearman $\\rho$（n=%d；绿色=置换检验 p<0.05）" % abl.get("rows", 0))
    ax.set_title("task15/16 逐特征预测力（误差棒=95% CI）", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    if own:
        save(fig, "validation_feature_spearman")
    return ax


# ---------------------------------------------------------------------------
# 图 4：MFE 权重扫描
# ---------------------------------------------------------------------------
def fig_sweep(ax=None):
    abl = load("task16_ablation.json")
    sw = abl.get("mfe_weight_sweep") or {}
    cv = sw.get("cv_mean_spearman") or {}
    full = sw.get("full_sample") or {}
    if not cv:
        return None
    a = sorted(float(k) for k in cv)
    cv_v = [cv.get(str(k)) if str(k) in cv else cv.get("%g" % k) for k in a]
    full_v = [full.get(str(k), {}).get("spearman") for k in a]
    best = sw.get("recommended_a_cv")
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(a, cv_v, "-o", color=C_NEW, label="CV 平均 Spearman（K=5）")
    ax.plot(a, full_v, "--s", color=C_THERMO, label="全样本 Spearman", markersize=4)
    if best is not None:
        ax.axvline(best, color="#263238", linestyle=":", linewidth=0.9)
        ax.text(best, max(v for v in cv_v if v is not None) + 0.01,
                "CV 最优 a=%.1f" % best, fontsize=8, ha="center")
    ax.axhline(0, color="#333", linewidth=0.8)
    ax.set_xlabel("a = 非 MFE 热力权重（a=0 为纯 MFE；a=1 为去掉 MFE）")
    ax.set_ylabel("Spearman $\\rho$")
    ax.set_title("task16 MFE 权重扫描（错配集）", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    if own:
        save(fig, "validation_mfe_weight_sweep")
    return ax


def fig_pair():
    """与原汇报版式一致：左=task14 三基准，右=task15 错配集（旧 vs 新）。"""
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 4.6))
    fig_task14(axes[0])
    fig_task15(axes[1])
    fig.tight_layout()
    save(fig, "validation_task14_task15_pair")


def main() -> int:
    fig_pair()
    fig_task14()
    fig_task15()
    fig_features()
    fig_sweep()
    # 2×2 汇总
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 9.0))
    fig_task14(axes[0][0])
    fig_task15(axes[0][1])
    fig_features(axes[1][0])
    if fig_sweep(axes[1][1]) is None:
        axes[1][1].axis("off")
    fig.suptitle("siRNA 计算验证结果（口径修正版 · 2026-09-11）", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, "validation_dashboard")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
