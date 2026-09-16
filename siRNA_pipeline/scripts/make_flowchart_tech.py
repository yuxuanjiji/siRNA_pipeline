# -*- coding: utf-8 -*-
"""答辩用：深色科技风项目流程图（16:9, 1920x1080@300dpi）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

# 配色（深色科技风）
BG_TOP    = "#0a1628"   # 深蓝底
BG_BOT    = "#13294b"
CARD      = "#1b2d4d"   # 卡片
CARD_EDGE = "#3b6ea5"
CYAN      = "#00d4ff"   # 主强调
CYAN_DIM  = "#4a90b8"
ORANGE    = "#ff9f43"
GREEN     = "#2ecc71"
TEXT      = "#eaf2fb"
TEXT_DIM  = "#9bb3d0"

fig, ax = plt.subplots(figsize=(13.33, 7.5))
ax.set_xlim(0, 133.33); ax.set_ylim(0, 75); ax.axis("off")

# 渐变背景
grad = np.linspace(0, 1, 256).reshape(-1, 1)
ax.imshow(grad, extent=[0,133.33,0,75], aspect="auto",
          cmap=matplotlib.colors.LinearSegmentedColormap.from_list("bg", [BG_TOP, BG_BOT]), zorder=-10)

def glow_card(x, y, w, h, text, fc=CARD, ec=CARD_EDGE, fs=10, tc=TEXT, weight="normal", lw=1.3):
    # 外发光（两层）
    for dx, dy, a in [(0,0,0.0)]:
        pass
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.9",
                       fc=fc, ec=ec, lw=lw, zorder=2)
    ax.add_patch(p)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center",
            fontsize=fs, color=tc, weight=weight, zorder=3, linespacing=1.35)

def arrow(x1, y1, x2, y2, color=CYAN_DIM):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>",
                                 mutation_scale=15, color=color, lw=1.8, zorder=1))

# 标题
ax.text(66.6, 72, "SFRP1 siRNA 计算设计流水线", ha="center",
        fontsize=20, weight="bold", color="white")
ax.text(66.6, 68.5, "纯计算 · 零湿实验 · 全链可审计", ha="center",
        fontsize=11, color=CYAN)

# 主链（第一行：5 步）
y1, h = 36, 13
steps = [
    (2,   "① 候选生成\nSFRP1 CDS\n全窗口 19nt\n14,832 条"),
    (27,  "② 规则过滤\nGC/长度/连续GC\n→ 6,524"),
    (52,  "③ 结构预测\nViennaRNA\nMFE/ΔΔG\n→ 4,417"),
    (77,  "④ 效率打分\n热力学×3 +\nOligoFormer 底座\n(引用 Bai 2024)", CARD, CYAN),
    (102, "⑤ 安全过滤\n脱靶/毒性\n→ 4,218"),
]
w = 22
for i, s in enumerate(steps):
    x = s[0]; txt = s[1]
    ec = s[3] if len(s)>3 else CARD_EDGE
    glow_card(x, y1, w, h, txt, ec=ec, fs=10.5, weight="bold" if len(s)>3 else "normal")
    if i < len(steps)-1:
        arrow(x+w, y1+h/2, steps[i+1][0], y1+h/2)

# 第二行（⑥⑦⑧）
y2, h2 = 15, 12
glow_card(77, y2, w, h2, "⑥ 组合排序\nα0.4 热力 + β0.6 DL\n− 惩罚项", ec=ORANGE, fs=10.5)
arrow(113, y1, 113, y2+h2)

glow_card(52, y2, w, h2, "⑦ 变体层\n错配 c_match\n(文献先验)", ec=ORANGE, fs=10.5)
arrow(77, y2+h2/2, 74, y2+h2/2)

glow_card(2, y2, w, h2, "⑧ 化学修饰交付\n2′-F/OMe/PS\n→ Top-50", ec=CYAN, fs=10.5, weight="bold")
arrow(52, y2+h2/2, 24, y2+h2/2)

# 顶部验证条
ax.text(66.6, 64.5, "独立验证证据（不进主链）", ha="center", fontsize=11, color=CYAN)
ev = [
    (2,   "公开基准 ρ\nHu 0.644 / Taka 0.565\nMix 0.707 / Simone 0.325"),
    (35,  "盲排回收\n4/5 文献有效 siRNA\n进前 25% (p≈1.4e-3)"),
    (68,  "稳健性\n200 次扰动\nTop-50 Jaccard 0.93"),
    (101, "可复现\nseed=42 / CPU\n参数全可追溯"),
]
for x, t in ev:
    glow_card(x, 55, 26, 8, t, fc="#14273f", ec=CARD_EDGE, fs=9.5)

# 底部定位
glow_card(10, 2, 113, 7,
          "定位：零湿实验交付「可直接送合成、每步可审计」的设计文档  |  最终效率需细胞实验验证（后续工作）",
          fc="#1a2a44", ec=ORANGE, fs=11.5, weight="bold")

out = r"D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\展示图片\答辩版\fig0_项目流程图_科技风.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG_TOP)
print("saved", out)
