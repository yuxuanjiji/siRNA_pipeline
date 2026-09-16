# -*- coding: utf-8 -*-
"""候选序列生成模块共用工具（任务 5 / 任务 6）。

约定（与全项目其它模块一致，务必核对后复用）：
* 碱基一律用 RNA 字母书写（A/U/G/C）。
* mRNA 序列按 5'->3' 方向存储；FASTA 为 DNA（T）时先转 RNA（U）。
* 每个 19 nt 靶窗口（target，mRNA 5'->3'）：
    - 引导链（guide / antisense，5'->3'）= target 的反向互补 rc(target)；
    - 有义链（sense / passenger，5'->3'）= 与 mRNA 窗口同向，序列 = target。
* 引导链位点编号 g1..g19 从 5' 端起（1-based）：g1 = guide 第 1 个碱基（5' 端），
  g19 = guide 第 19 个碱基（3' 端）。
* 双链反平行配对寄存器：guide 第 k 位（g_k，1-based）与 mRNA 靶窗口第 (20-k) 位
  （1-based，即 target[19-k] 0-based）碱基配对。
* SFRP1 参考：NM_003012.5 CDS（数据文件 数据集\\SFRP1-mRNA.txt，header location=315..1259），
  CDS 全长 945 nt；CDS 坐标 1-based，与转录本 NM_003012.5 坐标相差 +314。
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

RNA_ALPHABET = "AUGC"
_COMP_RNA = str.maketrans("AUGC", "UACG")   # A<->U, G<->C

WC_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
GU_PAIRS = {("G", "U"), ("U", "G")}


def to_rna(seq: str) -> str:
    """DNA -> RNA（仅替换 T->U；已为 RNA 时原样返回）。"""
    return seq.replace("T", "U").replace("t", "u")


def rc_rna(rna: str) -> str:
    """RNA 反向互补（输入须为 A/U/G/C）。"""
    return rna.translate(_COMP_RNA)[::-1]


def load_fasta_rna(path: str | Path) -> tuple[str, str]:
    """读取 FASTA，返回 (header, RNA 序列)。

    只保留 ACGTU 字母（忽略其它字符/空白）；DNA 自动转 RNA。
    若含非 ACGTU 的歧义碱基（N/R/Y...）会抛出 ValueError，防止静默污染候选库。
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
            else:  # 多记录 FASTA 时只取第一条（本项目 mRNA 文件为单条 CDS）
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
    """1 nt 滑窗。返回 [{'start','end','target','guide','sense','gc_pct'}, ...]。

    start/end 为 1-based CDS（或输入序列）坐标；target 为 mRNA 5'->3' 19 nt，
    guide 为其反向互补（5'->3'），sense 与 target 同向。
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


def pair_type(guide_nt: str, mrna_nt: str) -> str:
    """判断 guide 某位碱基与配对的 mRNA 碱基的配对类型。

    'WC'=完全互补；'GU'=G:U 摆动；'MM'=其它错配。
    """
    if (guide_nt, mrna_nt) in WC_PAIRS:
        return "WC"
    if (guide_nt, mrna_nt) in GU_PAIRS:
        return "GU"
    return "MM"


def replace_others(base: str, alphabet: str = RNA_ALPHABET) -> list[str]:
    """返回该位点可替换的其它 3 个碱基（按 AUGC 固定顺序，保证结果确定可复现）。"""
    return [b for b in alphabet if b != base]


def write_csv(path: str | Path, fieldnames: list[str], rows: list[dict]) -> int:
    """UTF-8-sig（带 BOM，Excel 打开中文不乱码）写 CSV，返回行数。"""
    rows = [{k: ("" if v is None else v) for k, v in r.items()} for r in rows]
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def default_fasta() -> Path:
    """.../生科挑战赛/数据集/SFRP1-mRNA.txt（脚本位于 项目搭建/候选序列生成 时）。"""
    here = Path(__file__).resolve()          # .../候选序列生成/seq_utils.py
    root = here.parents[2]                   # .../生科挑战赛
    return root / "数据集" / "SFRP1-mRNA.txt"


def default_outdir() -> Path:
    """候选序列生成 文件夹本身。"""
    return Path(__file__).resolve().parent
