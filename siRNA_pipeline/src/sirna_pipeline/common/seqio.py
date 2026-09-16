# -*- coding: utf-8 -*-
"""共享基础工具（sequence IO / 几何约定）。

本文件是各旧模块内重复实现的“序列工具”收敛位，算法与以下原始文件逐条等价：
  - 项目搭建/候选序列生成/seq_utils.py
  - 项目搭建/初级规则筛选/task7_sequence_rules_filter.py（normalize/gc_pct/rc）
原始文件保持不变且仍为权威；本模块只读等价复刻 + 少量纯函数补充。
"""
from __future__ import annotations

import re
from pathlib import Path

RNA_ALPHABET = "AUGC"
_COMP_RNA = str.maketrans("AUGC", "UACG")  # A<->U, G<->C

WC_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
GU_PAIRS = {("G", "U"), ("U", "G")}


def to_rna(seq: str) -> str:
    """DNA -> RNA（T->U，大小写均处理）；已为 RNA 时原样返回（大写化）。"""
    return str(seq).replace("T", "U").replace("t", "u").upper()


def normalize(seq: str) -> str:
    """清洗：去空白、大写、T->U（不含长度/字母校验，由 assert_rna 负责）。"""
    return "".join(str(seq).split()).upper().replace("T", "U")


def assert_rna(seq: str, label: str = "seq") -> str:
    """校验为合法 A/U/G/C 序列；否则抛中文 ValueError。返回原串。"""
    bad = sorted({c for c in seq if c not in RNA_ALPHABET})
    if bad:
        raise ValueError(f"{label} 含非 A/U/G/C 字符：{bad}")
    return seq


def rc_rna(rna: str) -> str:
    """RNA 反向互补（5'->3'）。输入须为 A/U/G/C。"""
    assert_rna(rna, "rna")
    return rna.translate(_COMP_RNA)[::-1]


def gc_pct(seq: str) -> float:
    """GC 含量（%），0~100。与 task7 gc_pct 一致。"""
    return (seq.count("G") + seq.count("C")) / len(seq) * 100.0


def pair_type(guide_nt: str, mrna_nt: str) -> str:
    """guide 某位碱基与其反平行配对的 mRNA 碱基的配对类型。

    与 seq_utils.pair_type 等价：'WC'=完全互补；'GU'=G:U 摆动；'MM'=其它错配。
    """
    if (guide_nt, mrna_nt) in WC_PAIRS:
        return "WC"
    if (guide_nt, mrna_nt) in GU_PAIRS:
        return "GU"
    return "MM"


def load_fasta_rna(path: str | Path) -> tuple[str, str]:
    """读取 FASTA（单条记录，取第一条）返回 (header, RNA 序列)。

    与 seq_utils.load_fasta_rna 等价：只保留 ACGTU 字母，DNA 自动转 RNA；
    含歧义碱基抛出 ValueError，防止静默污染候选库。
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    header, blocks = "", []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if not header:
                header = line
            else:  # 多记录 FASTA 只取第一条（本项目 mRNA 为单条 CDS）
                break
            continue
        if header:
            blocks.append(line.upper())
    joined = "".join(blocks)
    clean = re.sub(r"[^ACGTU]", "", joined)
    if len(clean) != len(joined):
        raise ValueError(
            f"FASTA 含非 ACGTU 字符（歧义碱基/数字等）：共 {len(joined)} 字符，"
            f"保留 {len(clean)} 字符。请先清洗序列文件。"
        )
    return header, to_rna(clean)


def sliding_windows(rna_seq: str, window_len: int = 19) -> list[dict]:
    """1nt 滑窗。与 seq_utils.sliding_windows 等价。

    返回 [{'start','end','target','guide','sense','gc_pct'}, ...]
    start/end 为 1-based 坐标；target=mRNA 5'->3' 窗口；
    guide=rc(target)；sense=target。
    """
    n = len(rna_seq)
    if n < window_len:
        raise ValueError(f"序列长度 {n} < 窗口长度 {window_len}")
    rows = []
    for i in range(n - window_len + 1):
        target = rna_seq[i:i + window_len]
        guide = rc_rna(target)
        gc = (target.count("G") + target.count("C")) / window_len
        rows.append({
            "start": i + 1,
            "end": i + window_len,
            "target": target,
            "guide": guide,
            "sense": target,
            "gc_pct": round(gc * 100.0, 2),
        })
    return rows


def read_csv_rows(path: str | Path) -> list[dict]:
    """读 UTF-8-sig CSV 为行字典列表（本项目全管道唯一 CSV 读取规范）。"""
    import csv
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def write_csv(path: str | Path, fieldnames: list[str], rows: list[dict]) -> int:
    """UTF-8-sig（带 BOM，Excel 中文不乱码）写 CSV，返回行数。

    与 seq_utils.write_csv 等价（None -> 空串）。
    """
    import csv
    rows = [{k: ("" if v is None else v) for k, v in r.items()} for r in rows]
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
