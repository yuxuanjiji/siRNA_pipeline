# -*- coding: utf-8 -*-
"""答辩用性能图统一生成脚本（16:9, 300dpi, 中文）"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import spearmanr
from scipy.stats import pearsonr

# ---------- 全局风格 ----------
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14
plt.rcParams["axes.edgecolor"] = "#444"
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["grid.linestyle"] = "--"

C_NAVY = "#1f3a6b"
C_BLUE = "#4a7bb7"
C_LBLUE = "#9bb8d9"
C_ORANGE = "#e07b39"
C_GREEN = "#4a8c5f"
C_GRAY = "#9aa0a6"
C_RED = "#c0392b"

OUT = r"D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\展示图片\答辩版"
os.makedirs(OUT, exist_ok=True)
PRED = r"D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\siRNA_pipeline\OligoFormer部分"

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved:", name)

# ============================================================
# 图1：预测 vs 实测散点（Hu 集，门面图）
# ============================================================
df = pd.read_csv(os.path.join(PRED, "Hu_predictions.csv"))
rho, p = spearmanr(df["label"], df["pred"])
pcc, pp = pearsonr(df["label"], df["pred"])
fig, ax = plt.subplots(figsize=(10, 6.2))
ax.scatter(df["label"], df["pred"], s=14, alpha=0.35, color=C_BLUE, edgecolors="none")
lim = [0, 1.02]
ax.plot(lim, lim, "--", color=C_GRAY, lw=1.5, label="y = x（完美预测）")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("实验测得沉默效率", fontsize=15)
ax.set_ylabel("OligoFormer 预测沉默效率", fontsize=15)
ax.set_title("预测引擎选型复现：OligoFormer（Bai et al., Bioinformatics 2024）在 Hu 公开基准集上的水平",
             fontsize=15, weight="bold", pad=12)
ax.text(0.04, 0.95,
        f"Spearman ρ = {rho:.3f}    Pearson r = {pcc:.3f}    n = 2361\n"
        "本团队以 OligoFormer 为效率预测底座，\n"
        "在其上叠加组合分、变体层与化学修饰交付",
        transform=ax.transAxes, fontsize=12, va="top",
        bbox=dict(boxstyle="round,pad=0.45", fc="#fff8ec", ec=C_ORANGE))
ax.legend(loc="lower right", fontsize=12)
save(fig, "fig1_预测vs实测_散点.png")

# ============================================================
# 图2：跨数据集泛化（4 集 × 3 模型）
# ============================================================
datasets = ["Hu\n(n=2361)", "Taka\n(n=702)", "Mix\n(n=464)", "Simone 外部集\n(n=322)"]
thermo = [0.285, -0.021, 0.292, 0.105]
dl     = [0.644,  0.565, 0.689, 0.305]
comb   = [0.660,  0.590, 0.707, 0.325]
x = np.arange(len(datasets)); w = 0.26
fig, ax = plt.subplots(figsize=(11, 6.2))
b1 = ax.bar(x-w, thermo, w, label="纯热力学分", color=C_LBLUE)
b2 = ax.bar(x,   dl,     w, label="纯 OligoFormer", color=C_BLUE)
b3 = ax.bar(x+w, comb,   w, label="组合分（主0.4/辅0.6）", color=C_NAVY)
for bars in (b1,b2,b3):
    for b in bars:
        v = b.get_height()
        ax.text(b.get_x()+b.get_width()/2, v + (0.02 if v>=0 else -0.05),
                f"{v:.2f}", ha="center", fontsize=10)
ax.axhline(0, color="#333", lw=0.8)
ax.axvspan(2.5, 3.5, color="#fff4e6", alpha=0.6, zorder=0)
ax.text(3, 0.86, "未参与任何调参\n的独立外部集", ha="center", color=C_ORANGE, fontsize=11, weight="bold")
ax.set_xticks(x); ax.set_xticklabels(datasets, fontsize=12)
ax.set_ylabel("Spearman ρ（沉默效率预测）", fontsize=14)
ax.set_ylim(-0.15, 0.95)
ax.set_title("跨数据集预测性能：组合分在全部数据集上稳定优于单一支柱", fontsize=15, weight="bold", pad=10)
ax.legend(fontsize=11, loc="upper left")
save(fig, "fig2_跨数据集泛化.png")

# ============================================================
# 图3：正对照回收
# ============================================================
ctrls = [
    ("17-4b Suzuki 2008", 75, 1.8),
    ("17-2 Wang 2024",    182, 4.3),
    ("17-4a Suzuki 2008",  767, 18.2),
    ("17-4c Suzuki 2008", 857, 20.3),
    ("17-1 Broadley 2026",1057, 25.1),
]
names = [c[0] for c in ctrls][::-1]
pct   = [c[2] for c in ctrls][::-1]
ranks = [c[1] for c in ctrls][::-1]
colors = [C_GREEN if p<=10 else (C_BLUE if p<=25 else C_GRAY) for p in pct]
fig, ax = plt.subplots(figsize=(11, 5.6))
bars = ax.barh(names, pct, color=colors, height=0.6)
ax.axvline(50, color=C_RED, ls="--", lw=2)
ax.text(50.5, 4.4, "随机期望 = 50%", color=C_RED, fontsize=12)
for b, r, p_ in zip(bars, ranks, pct):
    ax.text(p_+0.8, b.get_y()+b.get_height()/2, f"rank {r}/4218（前 {p_:.1f}%）",
            va="center", fontsize=11)
ax.set_xlim(0, 70)
ax.set_xlabel("在 4218 条候选中的排名百分位（越靠左越好）", fontsize=13)
ax.set_title("5 条文献已发表有效 siRNA：4 条进入前 25%、2 条进入前 10%（平均 13.9%，p≈1.4e-3）",
             fontsize=14, weight="bold", pad=10)
ax.grid(axis="x"); ax.grid(axis="y", visible=False)
save(fig, "fig3_正对照回收.png")

# ============================================================
# 图4：错配位置效应（Holen 单错配，U 形曲线）
# ============================================================
holen = pd.read_csv(r"D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\artifacts\holen2005_gki785.csv")
# 只保留单错配（突变位置是单个整数，不含逗号）
single = holen[~holen["突变位置"].astype(str).str.contains(",")].copy()
single["pos"] = single["突变位置"].astype(int)
single = single[(single["pos"]>=1) & (single["pos"]<=19)]
fig, ax = plt.subplots(figsize=(11, 5.8))
# 散点
ax.scatter(single["pos"], single["沉默效率"], s=40, alpha=0.5, color=C_BLUE, zorder=3)
# 分箱均值
g = single.groupby("pos")["沉默效率"].agg(["mean","count"]).reset_index()
ax.plot(g["pos"], g["mean"], "-o", color=C_ORANGE, lw=2.5, ms=8, label="位置平均沉默效率", zorder=4)
# 高亮中心切割区
ax.axvspan(9, 11, color="#fdecea", alpha=0.5, zorder=0)
ax.text(10, 0.95, "中心切割区\n(9–11 nt)\n对错配最敏感", ha="center", color=C_RED, fontsize=11, weight="bold")
ax.set_xlabel("guide 链上错配位置（5′→3′）", fontsize=14)
ax.set_ylabel("错配后相对沉默效率", fontsize=14)
ax.set_xticks(range(1,20))
ax.set_ylim(0, 1.05)
ax.set_title("错配位置效应：中心切割区效率骤降，两端错配几乎不影响", fontsize=15, weight="bold", pad=10)
ax.legend(loc="lower left", fontsize=12)
save(fig, "fig4_错配位置U形.png")

# ============================================================
# 图5：端到端漏斗
# ============================================================
stages = ["01 全窗口\n19nt 候选", "02 规则硬过滤", "03 结构过滤", "07 可排序候选", "最终交付\nTop-50"]
counts = [14832, 6524, 4417, 4218, 50]
keep   = ["", "保留 44%", "保留 68%", "保留 95%", "保留 1%"]
fig, ax = plt.subplots(figsize=(11, 5.2))
ypos = np.arange(len(stages))[::-1]
colors_fun = [C_LBLUE, C_LBLUE, C_BLUE, C_BLUE, C_NAVY]
for i,(c,k,s) in enumerate(zip(counts, keep, stages)):
    ax.barh(ypos[i], np.log10(c)+0.5, color=colors_fun[i], height=0.6)
    ax.text(0.1, ypos[i], f"{c:,} 条", va="center", fontsize=13, color="white" if i>=2 else "#222", weight="bold")
    ax.text(np.log10(c)+0.7, ypos[i], k, va="center", fontsize=11, color="#555")
ax.set_yticks(ypos); ax.set_yticklabels(stages, fontsize=12)
ax.set_xticks([])
ax.set_xlim(0, 5.5)
ax.set_title("端到端过滤漏斗：14832 → 50，每一级条件可审计", fontsize=15, weight="bold", pad=10)
ax.grid(axis="x", visible=False); ax.grid(axis="y", visible=False)
save(fig, "fig5_漏斗.png")

print("\n全部完成，输出目录:", OUT)
