# -*- coding: utf-8 -*-
"""答辩用：项目流程图（16:9, 300dpi）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

NAVY="#1f3a6b"; BLUE="#4a7bb7"; LBLUE="#cfe0f0"; ORANGE="#e07b39"; GREEN="#4a8c5f"; GRAY="#9aa0a6"; CREAM="#fff4e6"

fig, ax = plt.subplots(figsize=(14, 7.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc, ec=NAVY, fs=11, tc="#222", weight="normal"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.8",
                       fc=fc, ec=ec, lw=1.4)
    ax.add_patch(p)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, color=tc, weight=weight)

def arrow(x1,y1,x2,y2,color="#555"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,color=color,lw=1.6))

# 标题
ax.text(50, 96, "面向 SFRP1 的 siRNA 计算设计流水线（纯计算 / 零湿实验 / 全链可审计）",
        ha="center", fontsize=15, weight="bold", color=NAVY)

# 主链（从左到右）
y = 55; h = 11
box(2,  y, 13, h, "① 候选生成\nSFRP1 CDS\n全窗口 19nt 滑窗\n14,832 条", LBLUE)
box(17, y, 12, h, "② 规则硬过滤\nGC% / 长度 /\n连续 GC / 位置\n→ 6,524 条", LBLUE)
box(31, y, 12, h, "③ 结构预测\nViennaRNA\nMFE / ΔΔG ends\n→ 4,417 条", LBLUE)
box(45, y, 14, h, "④ 效率打分\n热力学特征 ×3\n+ OligoFormer 底座\n（Bai 2024, 引用）", "#dbe9d8", ec=GREEN)
box(61, y, 13, h, "⑤ 安全过滤\n脱靶 / 毒性 /\n免疫 motif\n→ 4,218 条", LBLUE)
box(76, y, 12, h, "⑥ 组合排序\nα0.4 热力 +\nβ0.6 DL\n− 惩罚项", "#fde8d4", ec=ORANGE)

for x1,x2 in [(15,17),(29,31),(43,45),(59,61),(74,76)]:
    arrow(x1, y+h/2, x2, y+h/2)

# 变体层（⑥下方）
box(76, 38, 12, 9, "⑦ 变体层\n错配 c_match\n（文献先验）", "#fde8d4", ec=ORANGE, fs=10)
arrow(82, 55, 82, 47)

# 输出
box(76, 20, 12, 11, "⑧ 化学修饰\n2′-F / 2′-OMe / PS\n→ Top-50 交付", NAVY, tc="white", weight="bold")
arrow(82, 38, 82, 31)

# 左侧：验证证据（旁挂）
ax.text(25, 88, "性能与可复现性证据（不进主链，独立验证）", fontsize=12, weight="bold", color=NAVY)
box(2, 74, 22, 10, "公开基准复现\nHu ρ=0.644 / Taka 0.565 / Mix 0.707\nSimone 外部集 ρ=0.325", "#eef3fa", fs=10)
box(26, 74, 22, 10, "SFRP1 盲排回收\n4/5 文献有效 siRNA\n进前 25%（p≈1.4e-3）", "#eef3fa", fs=10)
box(50, 74, 22, 10, "稳健性\n200 次权重扰动\nTop-50 Jaccard 均值 0.93", "#eef3fa", fs=10)
box(74, 74, 22, 10, "可复现性\nseed=42 / CPU /\n全部参数与中间产物可追溯", "#eef3fa", fs=10)

# 底部定位条
box(8, 5, 84, 8,
    "定位：零湿实验条件下交付「可直接送合成、每一步可审计」的候选清单  |  最终效率需细胞实验验证（后续工作）",
    CREAM, ec=ORANGE, fs=12, weight="bold")

fig.savefig(r"D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\展示图片\答辩版\fig0_项目流程图.png",
            dpi=300, bbox_inches="tight", facecolor="white")
print("saved fig0_项目流程图.png")
