# -*- coding: utf-8 -*-
"""化学修饰规则引擎（第六步 · chemmod）。

按行业金标准（临床获批 siRNA 的“增强稳定化学”ESC 骨架逻辑，见
docs/design/chemmod_rules.md 证据表）对 19 nt 主链给出逐位修饰建议。

默认骨架（ESC-19，把 21 nt 临床配比缩放到 19 nt 主链）：
  * guide（引导链/反义链）：全链默认 2′-OMe；偶数位 {2,4,…,18} 换 2′-F；
    g1、g19 保持 2′-OMe；端部硫代磷酸酯(PS) 于键 {1,2,18}；5′ 端加磷酸(5′-P)。
  * passenger（有义链/乘客链）：全链 2′-OMe；{7,9,11,13,15,17} 换 2′-F；
    端部 PS 于键 {1,2,18}。
  叠加掩蔽规则（优先级高于默认配比）：
  * 免疫掩蔽：引导链/乘客链命中 U 串(≥3 连续 U)、UGUGU、GUCCUUCAA 免疫 motif
    的 U 位强制 2′-OMe（哪怕位于偶数 F 位）——2′-OMe 拮抗 TLR7/8（Judge 2006 等）；
  * 可选旋钮 seed_ome=True：seed(g2–g8) 全 2′-OMe（降低 miRNA 样脱靶，代价部分效力），
    默认 False（本管道脱靶已由第五步软惩罚处理，双保险按需开）。

工程口径：
  * 只改“修饰”，不改碱基序列 → seed 脱靶/序列属性全部不变；
  * passenger 序列按 rc(guide)（完全互补口径）生成，并在 chem_notes 标注与
    “sense=target(1-错配)”口径的差异（docs/design/phase2 §2.4 #3）；
  * 21 nt 化（加 2 nt 悬垂、临床全长骨架）不在本引擎内做，以备注给出后续动作；
  * 本模块不触碰旧目录 `项目搭建/化学修饰/`（21 nt 占位）的任何文件。
"""
from __future__ import annotations

import re

COMP = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}

# 免疫 motif（与毒性阶段任务10 口径一致）
IMMUNE_MOTIFS = ("UGUGU", "GUCCUUCAA")
POLY_U_MIN_RUN = 3            # ≥3 连续 U 才做免疫掩蔽（避免普通 UU 干扰 F/OMe 交替）

# ESC-19 默认布点（1-based 位置）
GUIDE_F_POSITIONS = tuple(range(2, 20, 2))     # 2,4,...,18
GUIDE_PS_BONDS = (1, 2, 18)
SENSE_F_POSITIONS = (7, 9, 11, 13, 15, 17)
SENSE_PS_BONDS = (1, 2, 18)

RULE_VERSION = "esc19-v1"
RULES_REF = ("ESC 骨架（2′-OMe/2′-F 交替 + 端部 PS + guide 5′-P）为临床 siRNA 金标准配比；"
             "免疫掩蔽 2′-OMe 依据 Judge 2006 / Robbins 2009；U-rich 与 UGUGU/GUCCUUCAA "
             "为序列依赖免疫刺激 motif（Judge 2005）。详见 docs/design/chemmod_rules.md。")


def rc(seq: str) -> str:
    """反向互补（RNA；T→U）。"""
    return "".join(COMP[b] for b in reversed(seq.upper().replace("T", "U")))


def find_motif_positions(seq: str, motifs=IMMUNE_MOTIFS) -> set[int]:
    """返回 motif 覆盖的 1-based 位置集合（跨 motif 可重叠计）。"""
    out: set[int] = set()
    for m in motifs:
        for mch in re.finditer(re.escape(m), seq):
            for k in range(mch.start(), mch.end()):
                out.add(k + 1)
    return out


def find_polyu_positions(seq: str, min_run: int = POLY_U_MIN_RUN) -> set[int]:
    """返回长度≥min_run 的连续 U 串内所有 1-based 位置。"""
    out: set[int] = set()
    for mch in re.finditer(r"U{%d,}" % min_run, seq):
        for k in range(mch.start(), mch.end()):
            out.add(k + 1)
    return out


def _immune_mask(seq: str) -> tuple[set[int], list[str]]:
    """U 位免疫掩蔽集 + 命中说明。掩蔽仅针对 motif/U 串中的 U 碱基。"""
    hits = find_motif_positions(seq) | find_polyu_positions(seq)
    mask = {p for p in hits if seq[p - 1] == "U"}
    notes = []
    for m in IMMUNE_MOTIFS:
        if m in seq:
            notes.append(f"命中免疫 motif {m} → 其 U 位 2′-OMe 掩蔽")
    if any(len(r) >= POLY_U_MIN_RUN for r in re.findall(r"U+", seq)):
        notes.append(f"命中连续 U 串(≥{POLY_U_MIN_RUN}) → U 位 2′-OMe 掩蔽")
    return mask, notes


def mod_notation(seq: str, f_positions: set[int], ps_bonds: set[int],
                 prefix5p: bool = False) -> str:
    """逐位修饰记法：mX=2′-OMe-X，fX=2′-F-X；token 间 's' 表示 PS 键。

    与旧目录 `_mod_notation.py` 记法一致；ps_bonds 为 1-based“键”号
    （键 k = 第 k 与 k+1 位之间）。
    """
    tokens = []
    for i, nt in enumerate(seq, 1):
        p = "f" if i in f_positions else "m"
        tokens.append(p + nt)
    parts = []
    for i in range(len(seq)):
        parts.append(tokens[i])
        if (i + 1) in ps_bonds:
            parts.append("s")
    s = " ".join(parts)
    return ("5′-P " + s) if prefix5p else s


def build_chem_plan(guide: str, passenger: str | None = None,
                    seed_ome: bool = False, n: int = 19) -> dict:
    """对 (guide, passenger) 生成修饰建议计划 dict。

    guide/passenger 长度必须为 n（默认 19）；passenger=None 时按 rc(guide) 生成。
    返回字段见模块 docstring 与 runner 映射。
    """
    if len(guide) != n:
        raise ValueError(f"guide 长度须为 {n}：{len(guide)}")
    if passenger is None:
        passenger = rc(guide)
    if len(passenger) != n:
        raise ValueError(f"passenger 长度须为 {n}：{len(passenger)}")

    g_f = set(GUIDE_F_POSITIONS) & set(range(1, n + 1))
    s_f = set(SENSE_F_POSITIONS) & set(range(1, n + 1))

    notes: list[str] = []
    # 免疫掩蔽（U 位强制 OMe）
    g_mask, g_notes = _immune_mask(guide)
    s_mask, s_notes = _immune_mask(passenger)
    g_f -= g_mask
    s_f -= s_mask
    notes += g_notes
    notes += [f"[passenger]{x}" for x in s_notes]

    # seed OMe 旋钮（仅对 guide 生效）
    if seed_ome:
        seed = set(range(2, 9)) & set(range(1, n + 1))
        g_f -= seed
        notes.append("旋钮 seed_ome=True：seed(g2–g8) 全 2′-OMe（降 miRNA 样脱靶，代价部分效力）")

    notes.insert(0, f"ESC-19 骨架：guide 偶数位 2′-F、passenger {sorted(SENSE_F_POSITIONS)} 位 2′-F，"
                    "其余 2′-OMe；端部 PS")
    if passenger == rc(guide):
        notes.append("passenger 按 rc(guide) 完全互补口径生成（sense=target 的 1-错配口径见 "
                     "docs/design/phase2 §2.4 #3）")
    notes.append("21 nt 化（2 nt 悬垂/临床全长）属后续成药迭代（任务19），本表为 19 nt 主链建议")

    return {
        "guide": guide,
        "passenger": passenger,
        "guide_f_positions": sorted(g_f),
        "guide_ome_positions": sorted(set(range(1, n + 1)) - g_f),
        "guide_ps_bonds": sorted(set(GUIDE_PS_BONDS) & set(range(1, n))),
        "sense_f_positions": sorted(s_f),
        "sense_ome_positions": sorted(set(range(1, n + 1)) - s_f),
        "sense_ps_bonds": sorted(set(SENSE_PS_BONDS) & set(range(1, n))),
        "guide_mod": mod_notation(guide, g_f, set(GUIDE_PS_BONDS), prefix5p=True),
        "sense_mod": mod_notation(passenger, s_f, set(SENSE_PS_BONDS), prefix5p=False),
        "guide_5p_phosphate": True,          # RISC 装载必需；sense 端不加 5′-P
        "notes": notes,
    }
