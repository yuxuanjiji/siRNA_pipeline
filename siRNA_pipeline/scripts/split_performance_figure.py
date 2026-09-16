# -*- coding: utf-8 -*-
"""把 2x2 的「性能总图」拆成两张可直接上 PPT 的独立图。

原 `fig_performance_robustness.png` 是 2x2：A 四基准集效能 / B β 敏感性 / C 权重扰动 / D 漏斗。
P05 只用 A、P07 只用 B+C+D，于是此前把同一文件复制了两份、要在 PPT 里手工裁 —— 本脚本消除这一步：

  outputs/analysis/figures/fig_perf_panelA.png     ← 仅四基准集效能（P05 用）
  outputs/analysis/figures/fig_perf_panelBCD.png   ← 仅稳健性三项 + 漏斗（P07 用）

用法：python scripts/split_performance_figure.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np              # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
FIG = AN / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 200
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

C_TH, C_DL, C_CO = "#E07A5F", "#4C78A8", "#2E7D32"
SETS = ["Hu", "Taka", "Mix", "Simone"]
LABELS = ["Hu\n(n=2361)", "Taka\n(n=702)", "Mix\n(n=464)", "Simone\n(n=322)"]


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


def panel_a(ax) -> None:
    t14 = json.loads((AN / "task14_benchmark.json").read_text(encoding="utf-8"))
    keys = [("pure_thermo", "纯热力学分", C_TH), ("pure_dl", "纯 OligoFormer", C_DL),
            ("thermo_plus_aux", "组合分 α0.4/β0.6", C_CO)]
    x, w = np.arange(len(SETS)), 0.26
    for i, (k, name, c) in enumerate(keys):
        vals, lo, hi = [], [], []
        for ds in SETS:
            m = t14["datasets"][ds]["metrics"][k]
            vals.append(m["spearman"])
            lo.append(max(0.0, m["spearman"] - m["ci95_low"]))
            hi.append(max(0.0, m["ci95_high"] - m["spearman"]))
        pos = x + (i - 1) * w
        ax.bar(pos, vals, w, label=name, color=c, yerr=[lo, hi], capsize=4,
               error_kw={"elinewidth": 1.2, "ecolor": "#455A64"})
        for p, v in zip(pos, vals):
            ax.text(p, v + 0.022, "%.3f" % v, ha="center", va="bottom", fontsize=9.5)
    ax.axvspan(2.5, 3.5, color="#FFF3E0", zorder=0)
    ax.text(3, -0.115, "外部验证集（未参与调参）", ha="center", fontsize=10,
            color="#E65100", fontweight="bold")
    ax.set_xticks(x, LABELS)
    ax.set_ylabel("Spearman ρ（沉默效率预测）", fontsize=11)
    ax.set_ylim(-0.17, 0.80)
    ax.axhline(0, color="#90A4AE", lw=0.8)
    ax.legend(frameon=False, fontsize=10, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.16))
    ax.set_title("四套公开数据集：组合分均高于单用热力学或单用深度学习（含 95% CI）",
                 fontsize=12, loc="left", pad=34)


def panel_b(ax) -> None:
    rob = json.loads((AN / "weight_robustness.json").read_text(encoding="utf-8"))
    betas, j20, rho = [], [], []
    for s in rob["scenarios"]:
        if str(s.get("scenario", "")).startswith("beta="):
            betas.append(float(s["beta"]))
            j20.append(s["jaccard_top20_vs_base"])
            rho.append(s["spearman_vs_base"])
    base = rob["base_alpha_beta"][1]
    betas.append(base); j20.append(1.0); rho.append(1.0)
    o = np.argsort(betas)
    betas = [betas[i] for i in o]; j20 = [j20[i] for i in o]; rho = [rho[i] for i in o]
    ax.plot(betas, j20, "o-", color=C_CO, lw=2.4, ms=8, label="Top-20 Jaccard")
    ax.plot(betas, rho, "s--", color=C_DL, lw=2.2, ms=7, label="全体排序 Spearman")
    ax.scatter([base], [1.0], s=200, facecolor="none", edgecolor="#C62828", lw=2.2, zorder=5)
    ax.annotate("采用 β=0.6\n（四基准集网格最优）", xy=(base, 1.0), xytext=(0.41, 0.62),
                fontsize=10, color="#C62828",
                arrowprops={"arrowstyle": "->", "color": "#C62828", "lw": 1.3})
    ax.set_xlabel("β（深度学习辅助分权重）", fontsize=10.5)
    ax.set_ylabel("与默认口径的一致性", fontsize=10.5)
    ax.set_ylim(0.25, 1.08)
    ax.set_title("β 是最敏感参数，且选在证据最优处", fontsize=11.5, loc="left")
    ax.legend(frameon=False, fontsize=9.5, loc="lower right")


def panel_c(ax) -> None:
    rob = json.loads((AN / "weight_robustness.json").read_text(encoding="utf-8"))
    rp = rob["random_perturbation_200"]
    names = ["Top-20\nJaccard", "Top-50\nJaccard", "全体\nSpearman"]
    means = [rp["jaccard_top20_mean"], rp["jaccard_top50_mean"], rp["spearman_mean"]]
    mins = [rp["jaccard_top20_min"], rp["jaccard_top50_min"], rp["spearman_min"]]
    bars = ax.bar(names, means, color=[C_CO, C_DL, "#7E57C2"], width=0.55)
    for b, mn, mean in zip(bars, mins, means):
        cx = b.get_x() + b.get_width() / 2
        ax.plot([cx, cx], [mn, mean], color="#37474F", lw=2)
        ax.plot([cx], [mn], marker="_", color="#37474F", ms=14, mew=2)
        ax.text(cx, mean + 0.012, "均值 %.3f\n最差 %.3f" % (mean, mn),
                ha="center", va="bottom", fontsize=9.5)
    ax.set_ylim(0.6, 1.10)
    ax.set_ylabel("与默认权重下排序的一致性", fontsize=10.5)
    ax.set_title("200 次随机权重扰动：结论稳定\n（rank1 在 %.0f%% 的扰动下会变 → 只报 Top-N）"
                 % (100 * rp["rank1_changed_ratio"]), fontsize=11.5, loc="left")


def _funnel_counts():
    """漏斗各级行数：优先读 run 产物；产物已归档/清理时退回交付文档记录的数值。"""
    run = ROOT / "outputs" / "runs" / "sfrp1_oligo_on"
    try:
        def nrows(p: Path) -> int:
            return sum(1 for _ in p.open(encoding="utf-8-sig")) - 1

        return [("01 全窗口", nrows(run / "01_generation" / "candidates.csv")),
                ("02 规则硬过滤", nrows(run / "02_rules" / "candidates_passed.csv")),
                ("03 结构硬过滤", nrows(run / "03_structure" / "candidates_passed.csv")),
                ("07 可排序候选", json.loads((run / "07_ranking" / "manifest.json")
                                             .read_text(encoding="utf-8"))["eligible"]),
                ("最终交付 Top-50", 50)], "run sfrp1_oligo_on 产物"
    except (FileNotFoundError, KeyError):
        return [("01 全窗口", 14832), ("02 规则硬过滤", 6524), ("03 结构硬过滤", 4417),
                ("07 可排序候选", 4218), ("最终交付 Top-50", 50)], "已归档 run 记录"


def panel_d(ax) -> None:
    stages, src = _funnel_counts()
    ys = np.arange(len(stages))[::-1]
    vals = [v for _, v in stages]
    ax.barh(ys, vals, color=["#B0BEC5", "#90A4AE", "#78909C", C_DL, C_CO], height=0.6)
    for y, (name, v), prev in zip(ys, stages, [None] + vals[:-1]):
        tag = "%d 条" % v if prev is None else "%d 条（保留 %.0f%%）" % (v, 100 * v / prev)
        ax.text(v * 1.05, y, tag, va="center", fontsize=9.5)
    ax.set_yticks(ys, [n for n, _ in stages])
    ax.set_xscale("log")
    ax.set_xlim(30, 60000)
    ax.set_xlabel("候选数（log 轴）　数据来源：%s" % src, fontsize=10.5)
    ax.set_title("端到端过滤漏斗：每级条件可审计", fontsize=11.5, loc="left")


def main() -> int:
    fig, ax = plt.subplots(figsize=(9.4, 5.6))
    panel_a(ax)
    fig.tight_layout()
    save(fig, "fig_perf_panelA")

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))
    panel_b(axes[0]); panel_c(axes[1]); panel_d(axes[2])
    fig.tight_layout()
    save(fig, "fig_perf_panelBCD")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
