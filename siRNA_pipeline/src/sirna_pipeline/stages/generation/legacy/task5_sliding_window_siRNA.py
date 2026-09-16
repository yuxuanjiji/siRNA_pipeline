# -*- coding: utf-8 -*-
"""任务 5：1 nt 滑窗脚本——从 SFRP1 mRNA 生成全部 19 nt 完全互补 siRNA。

对 mRNA（NM_003012.5 CDS，数据集\\SFRP1-mRNA.txt，945 nt）做 1 nt 步长的 19 nt
滑窗；每个窗口产出 3 条序列，写入 SFRP1_task5_...csv：
  * target_mRNA_19   : mRNA 5'->3' 靶窗口（RNA 字母，U）
  * guide_antisense_19 : 引导链（反义链）5'->3' = rc(target)，即与靶窗口完全互补的
                         siRNA 引导序列（后续任务中的“siRNA 序列”即指它）
  * sense_strand_19  : 有义链（过客链）5'->3'，与 mRNA 窗口同向（=target）
完全互补性验证：guide = reverse_complement(target) 逐窗口断言。

用法:
    python task5_sliding_window_siRNA.py [--fasta 数据集\\SFRP1-mRNA.txt]
                                        [--out 输出CSV路径] [--window 19]

运行环境：仅标准库（Python >= 3.9）。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from seq_utils import (  # noqa: E402
    default_fasta,
    default_outdir,
    load_fasta_rna,
    rc_rna,
    sliding_windows,
    write_csv,
)

OUT_NAME = "SFRP1_task5_sliding_window_19nt_fullcomp.csv"
FIELDNAMES = [
    "window_id",          # 窗口编号 W0001..W0927（= CDS 上 1..927 的 19 nt 窗口）
    "cds_start",          # CDS 1-based 起点（=CDS 坐标，转录本坐标 = cds_start+314）
    "cds_end",            # CDS 1-based 终点
    "nm_003012_start",    # NM_003012.5 转录本坐标起点（= cds_start+314）
    "nm_003012_end",      # NM_003012.5 转录本坐标终点（= cds_end+314）
    "target_mRNA_19",     # mRNA 5'->3' 靶窗口（RNA 字母）
    "guide_antisense_19", # 引导链(反义)5'->3'，与靶窗口完全互补 = rc(target)
    "sense_strand_19",    # 有义链(过客)5'->3'，与 mRNA 同向
    "window_gc_pct",      # 靶窗口 GC%
]


def build_rows(windows, cds_len: int) -> list[dict]:
    """滑窗结果 -> 输出行（附坐标、ID）。"""
    rows = []
    for k, w in enumerate(windows, start=1):
        s, e = w["start"], w["end"]
        rows.append({
            "window_id": f"W{k:04d}",
            "cds_start": s,
            "cds_end": e,
            "nm_003012_start": s + (315 - 1),   # CDS 1 -> NM_003012.5 坐标 315
            "nm_003012_end": e + (315 - 1),
            "target_mRNA_19": w["target"],
            "guide_antisense_19": w["guide"],
            "sense_strand_19": w["sense"],
            "window_gc_pct": w["gc_pct"],
        })
    return rows


def verify(windows) -> None:
    """完全互补性 + 唯一性自检。"""
    n_bad = sum(1 for w in windows if rc_rna(w["target"]) != w["guide"])
    if n_bad:
        raise AssertionError(f"{n_bad} 个窗口 guide != rc(target)，完全互补性不成立！")
    seqs = [w["guide"] for w in windows]
    assert len(seqs) == len(set(seqs)), "存在重复引导链序列（滑窗步长应为 1 nt）"
    bad_len = [w["guide"] for w in windows if len(w["guide"]) != 19]
    if bad_len:
        raise AssertionError(f"存在非 19 nt 引导链：{bad_len[:3]}")
    for w in windows:
        for ch in w["guide"] + w["target"]:
            if ch not in "AUGC":
                raise AssertionError(f"含非 RNA 字母：{ch}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="任务5：SFRP1 mRNA 1nt 滑窗生成 19nt 完全互补 siRNA")
    ap.add_argument("--fasta", type=Path, default=default_fasta(), help="SFRP1 mRNA FASTA（默认读取数据集\\SFRP1-mRNA.txt）")
    ap.add_argument("--out", type=Path, default=default_outdir() / OUT_NAME, help="输出 CSV 路径")
    ap.add_argument("--window", type=int, default=19, help="滑窗长度（默认 19）")
    args = ap.parse_args(argv)

    header, rna_seq = load_fasta_rna(args.fasta)
    print(f"[输入] {args.fasta}")
    print(f"[参考] {header.splitlines()[0][:120]}")
    print(f"[序列] mRNA(RNA) 长度 = {len(rna_seq)} nt；窗口长度 = {args.window} nt")

    windows = sliding_windows(rna_seq, window_len=args.window)
    verify(windows)
    rows = build_rows(windows, len(rna_seq))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = write_csv(args.out, FIELDNAMES, rows)
    print(f"[输出] {args.out}  共 {n} 行（{len(rna_seq) - args.window + 1} 个窗口 × 1 引导链）")

    print("\n[样例] 前 3 行：")
    for r in rows[:3]:
        print(f"  {r['window_id']}  CDS {r['cds_start']}-{r['cds_end']}  "
              f"target={r['target_mRNA_19']}  guide={r['guide_antisense_19']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
