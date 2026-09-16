# -*- coding: utf-8 -*-
"""双链几何与坐标约定（duplex geometry）。

统一口径（与各旧模块注释一致，勿改）：
  * 碱基 RNA 字母；序列 5'->3' 书写。
  * 每 19nt 靶窗口 target(mRNA 5'->3')：guide = rc(target)，sense = target。
  * guide 位点 g1..g19 自 5' 端起 1-based；反平行配对 guide[p](0-based) ↔ target[18-p]。
  * mRNA 57nt = 19 上游 + 19 靶 + 19 下游，靶窗落在 [19:38]。
  * seed6（OligoFormer 毒性口径）= g2..g7；seed7（热力学/任务6口径）= g2..g8。
"""
from __future__ import annotations

# 任务6 突变位点（引导链 5'->3' 1-based）
MUTATION_SITES: tuple[int, ...] = (1, 12, 17, 18, 19)
WINDOW_LEN = 19
FLANK = 19                      # 侧翼长度（57nt 上下文 = 19+19+19）
MRNA57_LEN = 57

# seed 区定义（1-based 位点集合；两口径并存，调用方按 config 选择）
SEED6_POS: tuple[int, ...] = (2, 3, 4, 5, 6, 7)      # g2..g7
SEED7_POS: tuple[int, ...] = (2, 3, 4, 5, 6, 7, 8)   # g2..g8


def paired_mrna_nt(target: str, guide_pos: int) -> str:
    """与 guide 第 guide_pos 位(1-based)反平行配对的 mRNA 靶窗碱基。

    等价 task6._paired_mrna_nt：返回 target[19-guide_pos]（0-based）。
    target 必须为 19nt。
    """
    if len(target) != WINDOW_LEN:
        raise ValueError(f"靶窗长度必须为 {WINDOW_LEN}，当前 {len(target)}")
    if not (1 <= guide_pos <= WINDOW_LEN):
        raise ValueError(f"guide 位点须在 1..{WINDOW_LEN}：{guide_pos}")
    return target[WINDOW_LEN - guide_pos]


def mrna57_from_cds(cds_rna: str, cds_start: int) -> str | None:
    """由 CDS 全长与窗口起点重建 57nt 上下文。

    cds_start 为靶窗在 CDS 上的 1-based 起点（即 t5 表 cds_start）。
    返回 19 上游+19 靶+19 下游；任一侧翼越界时返回 None（由调用方决定处理策略：
    通常记该行结构/热力特征缺失，见 docs/design/phase2 §2.3）。
    """
    n = len(cds_rna)
    w0 = cds_start - 1                      # 靶窗 0-based 起点
    if w0 < 0 or w0 + WINDOW_LEN > n:
        return None
    b = w0 - FLANK
    e = w0 + WINDOW_LEN + FLANK
    if b < 0 or e > n:
        return None
    return cds_rna[b:e]


def seed_positions(seed_kind: str = "seed7") -> tuple[int, ...]:
    """按口径返回 seed 1-based 位点集合：seed6(g2-7) / seed7(g2-8)。"""
    if seed_kind == "seed6":
        return SEED6_POS
    if seed_kind == "seed7":
        return SEED7_POS
    raise ValueError(f"未知 seed 口径：{seed_kind}（可选 seed6/seed7）")


def guide_seed(guide: str, seed_kind: str = "seed6") -> str:
    """提取引导链 seed 子串。seed6: guide[1:7]（复刻 OligoFormer -tox）；seed7: guide[1:8]。"""
    if len(guide) != WINDOW_LEN:
        raise ValueError(f"guide 须为 {WINDOW_LEN}nt：{len(guide)}")
    if seed_kind == "seed6":
        return guide[1:7]
    if seed_kind == "seed7":
        return guide[1:8]
    raise ValueError(f"未知 seed 口径：{seed_kind}")
