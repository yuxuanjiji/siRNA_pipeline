# -*- coding: utf-8 -*-
"""答辩用示意图 4 张（用 matplotlib 画，风格与 pitch 图一致：字少、色块、无多余装饰）。

  S02_SFRP1通路      Wnt → FZD → β-catenin → 毛囊生长期；SFRP1 ⊣ Wnt；我们的干预点
  S03_六步流程       六个环节 + 箭头 + 下方候选数条
  S04_方法公式       终分公式 + α/β 取值（一行小字）
  S11_湿实验方案     三项实验 + 各自读出指标

产出：outputs/analysis/figures/sch_*.png，并同步到 PPT用图/
用法：python scripts/make_schematic_figures.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "outputs" / "analysis" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 200

GREEN, BLUE, ORANGE, GREY, DARK = "#2E7D32", "#4C78A8", "#E07A5F", "#90A4AE", "#37474F"


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


def box(ax, x, y, w, h, text, fc, tc="white", fs=15, r=0.12):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=%s" % r,
                                linewidth=0, facecolor=fc))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)


def arrow(ax, x1, y1, x2, y2, color=DARK, lw=2.4):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw,
                            "shrinkA": 0, "shrinkB": 0})


def canvas(w, h, xlim, ylim):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.axis("off")
    return fig, ax


def sch_pathway() -> None:
    fig, ax = canvas(11.0, 5.2, (0, 11), (0, 5.2))
    y = 3.5
    box(ax, 0.3, y, 1.5, 0.95, "Wnt", BLUE)
    box(ax, 2.6, y, 1.5, 0.95, "FZD", BLUE)
    box(ax, 4.9, y, 2.0, 0.95, "β-catenin", BLUE)
    box(ax, 7.7, y, 3.0, 0.95, "毛囊生长期", GREEN)
    for x1, x2 in ((1.8, 2.6), (4.1, 4.9), (6.9, 7.7)):
        arrow(ax, x1, y + 0.48, x2, y + 0.48)
    # SFRP1 抑制 Wnt
    box(ax, 0.3, 0.7, 1.9, 0.95, "SFRP1", ORANGE)
    arrow(ax, 1.25, 1.65, 1.05, y - 0.02, color=ORANGE)
    ax.plot([0.78, 1.32], [y - 0.10, y - 0.10], color=ORANGE, lw=3)
    ax.text(2.5, 1.15, "siRNA 敲低 SFRP1", fontsize=15, color=ORANGE, va="center")
    # 上调提示
    ax.text(0.55, 3.05, "脱发者表达上调", fontsize=12, color=GREY)
    save(fig, "sch_s02_pathway")


def sch_pipeline() -> None:
    fig, ax = canvas(13.0, 5.0, (0, 13), (0, 5.0))
    steps = ["① 生成候选", "② 规则过滤", "③ 结构过滤", "④ 热力学", "⑤ 脱靶·毒性", "⑥ 综合排序"]
    w, gap, y, h = 1.86, 0.24, 3.0, 1.05
    for i, s in enumerate(steps):
        x = 0.25 + i * (w + gap)
        box(ax, x, y, w, h, s, GREEN if i == 5 else BLUE)
        if i:
            arrow(ax, x - gap, y + h / 2, x, y + h / 2)
    box(ax, 0.25 + 6 * (w + gap), y, 1.5, h, "＋化学修饰", ORANGE, fs=13)
    box(ax, 0.25, 1.1, 12.5, 0.95,
        "14,832  →  6,524  →  4,417  →  4,218  →  Top-50", GREY, fs=16)
    save(fig, "sch_s03_pipeline")


def sch_formula() -> None:
    fig, ax = canvas(11.0, 4.2, (0, 11), (0, 4.2))
    ax.text(5.5, 3.0, "终分 = α · 热力学分  ＋  β · 深度学习分  −  脱靶 / 毒性惩罚",
            ha="center", va="center", fontsize=20, color=DARK)
    box(ax, 3.6, 1.35, 3.8, 0.85, "α = 0.4      β = 0.6", BLUE, fs=17)
    ax.text(5.5, 0.6, "权重由四套公开数据集网格搜索确定", ha="center", fontsize=13, color=GREY)
    save(fig, "sch_s04_formula")


def sch_experiments() -> None:
    fig, ax = canvas(13.0, 4.8, (0, 13), (0, 4.8))
    items = [("双荧光素酶报告基因", "读出：IC50、突变/WT 区分度"),
             ("体外 AGO2 切割", "读出：单周转切割速率"),
             ("RNA-seq 脱靶谱", "读出：脱靶位点数")]
    w, gap, y, h = 4.0, 0.5, 0.8, 2.6
    for i, (head, ro) in enumerate(items):
        x = 0.25 + i * (w + gap)
        box(ax, x, y, w, h, "", "#ECEFF1")
        box(ax, x, y + h - 0.95, w, 0.95, head, BLUE, fs=15)
        ax.text(x + w / 2, y + 0.85, ro, ha="center", va="center", fontsize=14, color=DARK)
    ax.text(6.5, 4.15, "下一步实验（已设计，待做）", ha="center", fontsize=16, color=DARK)
    save(fig, "sch_s11_experiments")


def main() -> int:
    sch_pathway(); sch_pipeline(); sch_formula(); sch_experiments()
    dst = ROOT.parent / "PPT用图"
    if dst.exists():
        for stem, name in {"sch_s02_pathway": "S02_SFRP1通路.png",
                           "sch_s03_pipeline": "S03_六步流程.png",
                           "sch_s04_formula": "S04_方法公式.png",
                           "sch_s11_experiments": "S11_湿实验方案.png"}.items():
            shutil.copyfile(FIG / ("%s.png" % stem), dst / name)
        print("synced ->", dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
