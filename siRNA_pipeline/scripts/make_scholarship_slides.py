# -*- coding: utf-8 -*-
"""奖学金答辩两页版式稿（16:9，2400x1350）：色带 + 卡片 + 大数字 + 漏斗条 + 真实图嵌入。

  S_奖学金页1_问题与方法.png   顶部色带 + 左侧通路图 + 右侧三张卡片 + 底部候选漏斗
  S_奖学金页2_结果与贡献.png   顶部色带 + 左侧正对照回收图 + 右侧三张数据卡 + 底部信息条

用途：① 直接当整页背景图用；② 或照它排版、把文字换成自己的（文字都在图上，位置已定）。
用法：python scripts/make_scholarship_slides.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
OUT = AN / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

GREEN, GREEN_L, BLUE, ORANGE, GREY, DARK = "#2E7D32", "#E8F5E9", "#4C78A8", "#E07A5F", "#90A4AE", "#263238"
FUNNEL = ["#B0BEC5", "#90A4AE", "#78909C", "#4C78A8", "#2E7D32"]


def canvas(w=16, h=9):
    fig = plt.figure(figsize=(w, h), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 16); ax.set_ylim(0, 9); ax.axis("off")
    return fig, ax


def card(ax, x, y, w, h, bar, title, body, tc=GREEN, fs_t=17, fs_b=13.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                linewidth=0, facecolor="#FFFFFF", zorder=2))
    ax.add_patch(FancyBboxPatch((x, y), 0.12, h, boxstyle="round,pad=0.0,rounding_size=0.06",
                                linewidth=0, facecolor=bar, zorder=3))
    ax.text(x + 0.42, y + h - 0.42, title, fontsize=fs_t, color=tc, fontweight="bold",
            va="center", zorder=4)
    ax.text(x + 0.42, y + h * 0.30, body, fontsize=fs_b, color=DARK, va="center", zorder=4)


def band(ax, y, h, fc, text, tc, fs=25, align="left"):
    ax.add_patch(plt.Rectangle((0, y), 16, h, facecolor=fc, linewidth=0, zorder=1))
    x = 0.55 if align == "left" else 8
    ax.text(x, y + h / 2, text, fontsize=fs, color=tc, va="center", zorder=4,
            ha="left" if align == "left" else "center")


def place_image(fig, path: Path, rect):
    """把图片等比缩放进 rect=[l,b,w,h]（画布数据坐标；画布 16x9 英寸、坐标 0..16 / 0..9，
    故 1 数据单位 = 1 英寸，可按英寸算等比后换算回 add_axes 的比例坐标）。"""
    img = plt.imread(path)
    ih, iw = img.shape[0], img.shape[1]
    ar = iw / ih
    l, b, w, h = rect
    w2, h2 = (h * ar, h) if (w / h) > ar else (w, w / ar)
    l2, b2 = l + (w - w2) / 2, b + (h - h2) / 2
    ax = fig.add_axes([l2 / 16.0, b2 / 9.0, w2 / 16.0, h2 / 9.0])
    ax.imshow(img)
    ax.axis("off")


def slide1(fig, ax) -> None:
    band(ax, 8.25, 0.75, GREEN, "脱发毛囊的「刹车分子」SFRP1  →  我们为它设计 siRNA",
         "white", fs=26)
    place_image(fig, OUT / "sch_s02_pathway.png", [0.30, 3.30, 6.20, 4.60])
    ax.text(3.35, 3.05, "毛囊再生通路与干预点", fontsize=13, color=GREY, ha="center")
    card(ax, 6.85, 6.30, 8.75, 1.25, BLUE, "科学问题",
         "1.5 万条候选，选出真正有效的「那一条」", tc=BLUE)
    card(ax, 6.85, 4.90, 8.75, 1.25, GREEN, "方法",
         "热力学 4 特征 ＋ 深度学习效率模型加权融合（α0.4 / β0.6），叠脱靶·毒性·免疫过滤")
    card(ax, 6.85, 3.50, 8.75, 1.25, ORANGE, "产出",
         "一条可复现流水线 →  50 条可直接合成的候选", tc=ORANGE)
    ax.text(0.45, 2.62, "候选漏斗（每级过滤条件可审计）", fontsize=13.5, color=GREY)
    labels = ["14,832\n全部窗口", "6,524\n规则过滤", "4,417\n结构过滤", "4,218\n可排序", "50\n交付清单"]
    x, w, gap = 0.45, 2.86, 0.28
    for i, (lab, c) in enumerate(zip(labels, FUNNEL)):
        ax.add_patch(FancyBboxPatch((x, 1.05), w, 1.30,
                                    boxstyle="round,pad=0.02,rounding_size=0.10",
                                    linewidth=0, facecolor=c, zorder=2))
        ax.text(x + w / 2, 1.70, lab, ha="center", va="center", fontsize=14.5,
                color="white", zorder=4, linespacing=1.5)
        if i < 4:
            ax.annotate("", xy=(x + w + gap - 0.04, 1.70), xytext=(x + w + 0.02, 1.70),
                        arrowprops={"arrowstyle": "-|>", "color": GREY, "lw": 2.2})
        x += w + gap


def slide2(fig, ax) -> None:
    band(ax, 8.25, 0.75, GREEN_L,
         "文献已发表有效的 5 条 siRNA，全部排进候选库前 25%", GREEN, fs=26)
    place_image(fig, ROOT.parent / "PPT用图" / "P06_正对照回收.png", [0.30, 2.10, 8.90, 5.80])
    ax.text(4.75, 1.72, "红虚线 = 随机期望（50%）　　横轴 = 排名百分位", fontsize=13,
            color=GREY, ha="center")
    cards = [("5 / 5", "文献有效 siRNA 进入前 25%"),
             ("13.9%", "平均排名百分位（随机期望 50%）"),
             ("p ≈ 0.0014", "Irwin–Hall 检验，n = 5")]
    y = 6.35
    for num, lab in cards:
        ax.add_patch(FancyBboxPatch((9.30, y), 6.30, 1.80,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    linewidth=0, facecolor="#FFFFFF", zorder=2))
        ax.text(9.75, y + 1.08, num, fontsize=33, color=GREEN, fontweight="bold", zorder=4)
        ax.text(9.75, y + 0.44, lab, fontsize=13, color=DARK, zorder=4)
        y -= 2.10
    ax.add_patch(plt.Rectangle((0, 0.30), 16, 1.15, facecolor="#ECEFF1", linewidth=0, zorder=1))
    ax.text(0.55, 1.08, "我的贡献：〔数据构建 · 模型实现 · 流水线工程 · 结果复现，按实际填 2–3 项〕",
            fontsize=15, color=DARK, va="center", zorder=4)
    ax.text(0.55, 0.62, "可复现：固定随机种子 42 · 标准交付文件 · 任何人可复跑",
            fontsize=13, color=GREY, va="center", zorder=4)
    ax.text(15.45, 0.62, "竞赛：〔级别〕〔奖项〕｜本人排名：〔第 N 完成人〕｜〔年-月〕",
            fontsize=13, color=GREY, ha="right", va="center", zorder=4)


def main() -> int:
    fig, ax = canvas(); slide1(fig, ax)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("S_奖学金页1_问题与方法.%s" % ext))
    plt.close(fig)
    fig, ax = canvas(); slide2(fig, ax)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("S_奖学金页2_结果与贡献.%s" % ext))
    plt.close(fig)
    dst = ROOT.parent / "PPT用图"
    if dst.exists():
        for stem, name in {"S_奖学金页1_问题与方法": "S_奖学金页1_问题与方法.png",
                           "S_奖学金页2_结果与贡献": "S_奖学金页2_结果与贡献.png"}.items():
            shutil.copyfile(OUT / ("%s.png" % stem), dst / name)
        print("synced ->", dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
