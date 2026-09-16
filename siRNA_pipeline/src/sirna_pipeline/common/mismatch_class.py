# -*- coding: utf-8 -*-
"""错配分类：细粒度 WC/GU/PP/YY/PY（thermo/calibrate/validate 使用的语义）。

等价来源：
  - 项目搭建/热力学参数计算/calibrate_mismatch.py::classify
  - 项目搭建/热力学参数计算/validate_mismatch_dataset.py::classify
细粒度用于热力学惩罚查表；粗粒度 WC/GU/MM 见 seqio.pair_type。
"""
from __future__ import annotations

_COMP = {"A": "U", "U": "A", "C": "G", "G": "C"}


def mismatch_class(guide_nt: str, mrna_nt: str) -> str:
    """guide 碱基与 mRNA 靶窗碱基（反平行列）的配对分类。

    返回：WC / GU / PP(嘌呤-嘌呤) / YY(嘧啶-嘧啶) / PY(嘌呤-嘧啶)。
    """
    g, m = guide_nt.upper(), mrna_nt.upper()
    if g == _COMP.get(m, ""):
        return "WC"
    if {g, m} == {"G", "U"}:
        return "GU"
    pur_g, pur_m = g in "AG", m in "AG"
    if pur_g and pur_m:
        return "PP"
    if (not pur_g) and (not pur_m):
        return "YY"
    return "PY"


def column_classes(guide: str, window: str) -> list[str]:
    """逐列（guide 0-based p ↔ window[18-p]）分类，返回 19 元素列表。"""
    if len(guide) != len(window):
        raise ValueError(f"guide 与 window 长度不一致：{len(guide)} vs {len(window)}")
    return [mismatch_class(g, window[len(window) - 1 - p])
            for p, g in enumerate(guide)]
