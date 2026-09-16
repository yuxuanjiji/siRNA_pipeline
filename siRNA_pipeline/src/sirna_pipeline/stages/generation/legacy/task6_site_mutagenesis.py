# -*- coding: utf-8 -*-
"""任务 6：定点突变脚本——在引导链 g1/g12/g17/g18/g19 五个位点做单位点替换，
每个完全互补窗口扩展出 15 条突变体，输出完整候选库（错配变体库）。

突变规则（按任务要求 & 团队确认口径）：
  * 位点：引导链 5'->3' 1-based 编号中的 g1、g12、g17、g18、g19；
  * 每个位点分别替换为“其它 3 个碱基”（A/U/G/C 中除去原碱基）；
  * 每次仅突变一个位点（单点突变），与原始完全互补引导链仅差 1 nt；
  * 因此每个 19 nt 靶窗口生成 5 × 3 = 15 条突变体，全部窗口共
    windows × 15 行。完全互补“野生型”引导链由任务 5 输出单独保留，不入本库。

设计说明 / 口径标注：
  * 突变后与 mRNA 靶窗口在该列形成的配对类型（pair_type）:
      'GU' = G:U 摆动；'MM' = 其它错配（替换碱基排除了原完全互补碱基，故必非 WC）。
  * 文档原话“g1 限定 G:U 摆动”默认不筛（团队已确认：g1 与其余 4 位点同样替换
    其它 3 碱基，全保留）。如需按文档字面只保留摆动替换，加参数 --g1-wobble-only。
  * 五位点均位于 seed 区(g2-g8)之外 ⇒ 同一窗口 15 条变体 seed 完全相同
    （下游脱靶/毒性模块可“窗口级算一次、广播 15 条”，见任务 9/10 说明）。
  * 顺序固定、无随机数 ⇒ 结果完全可复现。

用法:
    python task6_site_mutagenesis.py
        [--input 任务5输出CSV]           默认读取本文件夹的 SFRP1_task5_...csv
        [--fasta SFRP1 mRNA]            当 --input 不存在时直接由 FASTA 重建窗口
        [--out 输出CSV]                 默认 本文件夹/SFRP1_task6_site_mutagenesis_library.csv
        [--g1-wobble-only]              可选：g1 只保留能与 mRNA 形成 G:U 摆动的替换
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from seq_utils import (  # noqa: E402
    default_fasta,
    default_outdir,
    load_fasta_rna,
    pair_type,
    replace_others,
    sliding_windows,
    write_csv,
)
from task5_sliding_window_siRNA import OUT_NAME as TASK5_OUT  # noqa: E402

MUTATION_SITES = [1, 12, 17, 18, 19]      # 引导链 5'->3' 位点（g1..g19 1-based）
OUT_NAME = "SFRP1_task6_site_mutagenesis_library.csv"
FIELDNAMES = [
    "window_id",          # 对应任务 5 的窗口编号（W0001..）
    "cds_start", "cds_end",
    "target_mRNA_19",     # mRNA 5'->3' 靶窗口（RNA）
    "guide_wt_19",        # 原始完全互补引导链（=任务5 guide_antisense_19）
    "variant_no",         # 窗口内序号 v01..v15（固定顺序，便于核对）
    "variant_id",         # 唯一变体名：W0001_g1:A>U
    "site",               # g1/g12/g17/g18/g19
    "position",           # 突变位点（引导链 1-based 数值）
    "wt_nt", "mut_nt",    # 原碱基 -> 突变碱基
    "guide_mut_19",       # 突变后引导链 5'->3'（RNA，19 nt）
    "paired_mRNA_nt",     # 与突变位点配对的 mRNA 碱基（= target[19-position]）
    "pair_type",          # 'GU'=G:U 摆动, 'MM'=其它错配（WC 必不会出现）
]


def _paired_mrna_nt(target: str, pos: int) -> str:
    """guide 第 pos 位(1-based) 配对的 mRNA 靶窗口碱基。

    反平行寄存器：guide 5'->3' 第 pos 位 与 target[19-pos] (0-based) 配对。
    """
    return target[19 - pos]


def expand_window(window_row: dict, g1_wobble_only: bool) -> list[dict]:
    """对单个窗口（完全互补 guide）生成 15 条单点突变体。"""
    target = window_row["target_mRNA_19"]
    wt = window_row["guide_antisense_19"]
    if len(wt) != 19:
        raise ValueError(f"{window_row['window_id']}: 引导链长度 {len(wt)} != 19")

    rows = []
    no = 0
    for pos in MUTATION_SITES:                     # g1 -> g12 -> g17 -> g18 -> g19
        wt_nt = wt[pos - 1]
        mrna_nt = _paired_mrna_nt(target, pos)
        for mut_nt in replace_others(wt_nt):       # 其它 3 碱基，顺序固定
            if g1_wobble_only and pos == 1 and pair_type(mut_nt, mrna_nt) != "GU":
                continue                            # g1 只保留 G:U 摆动替换（可选）
            no += 1
            mut = wt[: pos - 1] + mut_nt + wt[pos:]
            rows.append({
                "window_id": window_row["window_id"],
                "cds_start": window_row["cds_start"],
                "cds_end": window_row["cds_end"],
                "target_mRNA_19": target,
                "guide_wt_19": wt,
                "variant_no": f"v{no:02d}",
                "variant_id": f"{window_row['window_id']}_g{pos}:{wt_nt}>{mut_nt}",
                "site": f"g{pos}",
                "position": pos,
                "wt_nt": wt_nt,
                "mut_nt": mut_nt,
                "guide_mut_19": mut,
                "paired_mRNA_nt": mrna_nt,
                "pair_type": pair_type(mut_nt, mrna_nt),
            })
    return rows


def read_task5_windows(csv_path: Path) -> list[dict]:
    """读取任务 5 输出 CSV 为窗口行（转成通用字段名）。"""
    with open(csv_path, encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def windows_from_fasta(fasta: Path) -> list[dict]:
    """无任务 5 CSV 时由 FASTA 重建窗口（与任务 5 同口径）。"""
    header, rna_seq = load_fasta_rna(fasta)
    windows = sliding_windows(rna_seq, window_len=19)
    rows = []
    for k, w in enumerate(windows, start=1):
        s, e = w["start"], w["end"]
        rows.append({
            "window_id": f"W{k:04d}",
            "cds_start": s,
            "cds_end": e,
            "nm_003012_start": s + 314,
            "nm_003012_end": e + 314,
            "target_mRNA_19": w["target"],
            "guide_antisense_19": w["guide"],
            "sense_strand_19": w["sense"],
            "window_gc_pct": w["gc_pct"],
        })
    return rows


def expected_per_window(wt: str, target: str, g1_wobble_only: bool) -> dict:
    """复算某窗口在给定模式下的逐位点期望：{pos: 允许的替换碱基集合}。"""
    exp = {}
    for pos in MUTATION_SITES:
        wt_nt = wt[pos - 1]
        alts = set(replace_others(wt_nt))
        if g1_wobble_only and pos == 1:
            mrna_nt = _paired_mrna_nt(target, pos)
            alts = {b for b in alts if pair_type(b, mrna_nt) == "GU"}
        exp[pos] = alts
    return exp


def verify_library(rows: list[dict], n_windows: int, g1_wobble_only: bool) -> None:
    """库级自检：窗口数/每窗口条数/单点突变/位点覆盖/碱基替换正确。"""
    if not rows:
        raise AssertionError("候选库为空！")
    per_window: dict[str, int] = {}
    for r in rows:
        per_window[r["window_id"]] = per_window.get(r["window_id"], 0) + 1
        wt, mut = r["guide_wt_19"], r["guide_mut_19"]
        assert len(wt) == len(mut) == 19, "长度 != 19"
        diff = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
        assert len(diff) == 1, f"{r['variant_id']}: 非单点突变，差异位 {diff}"
        pos = int(r["position"])
        assert diff[0] == pos - 1, f"{r['variant_id']}: 差异位点与 position 不符"
        assert wt[pos - 1] == r["wt_nt"] and mut[pos - 1] == r["mut_nt"]
        assert r["mut_nt"] != r["wt_nt"] and r["mut_nt"] in "AUGC"
        # 突变换掉的是原完全互补碱基 -> 该列绝不可能是 WC
        assert r["pair_type"] in ("GU", "MM"), r["variant_id"]
        assert r["paired_mRNA_nt"] == _paired_mrna_nt(r["target_mRNA_19"], pos)
        assert r["guide_mut_19"] == wt[: pos - 1] + r["mut_nt"] + wt[pos:]
    assert len(per_window) == n_windows, "窗口数与输入不符"

    # 每个窗口：条数 = 期望值；位点替换碱基集 = 期望集合；15 条序列互不相同
    from collections import defaultdict
    per_win_data: dict[str, dict] = defaultdict(lambda: {"seqs": set(), "wt": None,
                                                         "target": None,
                                                         "alts": defaultdict(set)})
    for r in rows:
        d = per_win_data[r["window_id"]]
        d["seqs"].add(r["guide_mut_19"])
        d["wt"] = r["guide_wt_19"]
        d["target"] = r["target_mRNA_19"]
        d["alts"][int(r["position"])].add(r["mut_nt"])
    for wid, d in per_win_data.items():
        exp = expected_per_window(d["wt"], d["target"], g1_wobble_only)
        exp_n = sum(len(s) for s in exp.values())
        assert len(d["seqs"]) == exp_n, f"{wid}: 条数 {len(d['seqs'])} != 期望 {exp_n}"
        for pos in MUTATION_SITES:
            assert d["alts"][pos] == exp[pos], \
                f"{wid}: 位点 g{pos} 替换碱基集 {sorted(d['alts'][pos])} != 期望 {sorted(exp[pos])}"
    mode = "g1 仅 G:U 摆动" if g1_wobble_only else "g1 不过滤(3 碱基全保留)"
    print("[校验] 全部通过：窗口数=%d，总行数=%d，每窗口条数与位点替换集均等于期望"
          "（模式：%s），无重复序列。" % (n_windows, len(rows), mode))
    print("[校验] 全部通过：窗口数=%d，总行数=%d，每窗口恰 15 条单点突变、无重复序列。"
          % (n_windows, len(rows)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="任务6：g1/g12/g17/g18/g19 单点突变扩库（15 条/窗口）")
    ap.add_argument("--input", type=Path, default=default_outdir() / TASK5_OUT,
                    help="任务5输出 CSV（默认本文件夹 SFRP1_task5_...csv）")
    ap.add_argument("--fasta", type=Path, default=default_fasta(),
                    help="--input 不存在时改由该 FASTA 直接重建窗口")
    ap.add_argument("--out", type=Path, default=default_outdir() / OUT_NAME,
                    help="输出完整候选库 CSV 路径")
    ap.add_argument("--g1-wobble-only", action="store_true",
                    help="可选：g1 位只保留能与 mRNA 形成 G:U 摆动的替换（默认不筛）")
    args = ap.parse_args(argv)

    if args.input.exists():
        print(f"[输入] 任务5窗口表：{args.input}")
        windows = read_task5_windows(args.input)
    else:
        print(f"[输入] 未找到 {args.input}，改由 FASTA 直接重建窗口")
        windows = windows_from_fasta(args.fasta)
    n_windows = len(windows)
    if n_windows == 0:
        raise SystemExit("错误：输入窗口数为 0。请先运行 task5_sliding_window_siRNA.py。")

    rows: list[dict] = []
    for w in windows:
        rows.extend(expand_window(w, g1_wobble_only=args.g1_wobble_only))

    verify_library(rows, n_windows, g1_wobble_only=args.g1_wobble_only)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = write_csv(args.out, FIELDNAMES, rows)
    from collections import Counter
    pair_c = Counter(r["pair_type"] for r in rows)
    site_c = Counter(r["site"] for r in rows)
    print(f"[输出] 完整候选库：{args.out}  共 {n} 行 = {n_windows} 窗口 × 15 条/窗口")
    print(f"[统计] 配对类型分布：{dict(pair_c)}")
    print(f"[统计] 各位点变体数：{dict(site_c)}")
    print("\n[样例] 首个窗口前 3 条：")
    for r in [x for x in rows if x["window_id"] == rows[0]["window_id"]][:3]:
        print(f"  {r['variant_id']:<24} guide={r['guide_mut_19']}  "
              f"mRNA配对碱基={r['paired_mRNA_nt']}  {r['pair_type']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
