# -*- coding: utf-8 -*-
"""BLAST 脱靶检测模块（脱靶验证第三层：近全长序列同源筛选）
================================================================
来源：本项目自建层（算法 Altschul et al. 1990 J Mol Biol 215:403-410；
      实现 NCBI BLAST+ 2.14.1，Camacho et al. 2009 BMC Bioinformatics 10:421；
      参数依据 NCBI BLAST Command Line Applications User Manual）
数据库：GENCODE v46 蛋白编码全转录本（约 11 万 ENST，含 5'/3' UTR）

正式判据（tolerant 容错口径，2026-09-06 第三轮修正定稿）：
  一条候选判为"脱靶"当且仅当其任一 BLAST 命中同时满足：
    1) 命中基因不含靶基因名（默认 SFRP1）
    2) 反向互补方向命中（sstart > send，即转录本含 rc(guide)，guide 可结合）
    3) 正确配对数 (length - mismatch) >= 17
    4) 错配数 mismatch <= 1
  可选加固：--require-no-gaps（额外要求 gapopen == 0；正式归档口径未启用）

用法（Python 3.8+，仅需 pandas + 本地 BLAST+ 可执行文件）：
  python blast_offtarget.py --input candidates.csv --out-dir out --prefix sfrp1
  可选：--make-db transcripts.fa --db-prefix mydb   # 从 FASTA 重建 BLAST 库

输入：CSV（须含 siRNA 列，19nt RNA 含 U）或 FASTA（标题行为候选名）
输出（--out-dir 下）：
  <prefix>_queries.fa           实际送检的 FASTA（U->T）
  <prefix>_hits_all.tsv         全部命中（>=16nt，审计留痕，与归档口径一致）
  <prefix>_flagged.csv          判脱靶的候选
  <prefix>_clean.csv            通过的候选
  <prefix>_stats.json           参数与计数
"""
import argparse
import json
import os
import subprocess
import sys

import pandas as pd

DEFAULT_BLAST_BIN = r"E:\大创2026\blast\ncbi-blast-2.14.1+\bin"
DEFAULT_DB = r"E:\blastdb\gencode_v46_pc_simple"
DEFAULT_TARGET_GENE = "SFRP1"

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


# ---------- 输入读取 ----------

def read_candidates(path, column="siRNA"):
    """返回候选 DataFrame（一列 name + 一列 seq，DNA 字母）。"""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".fa", ".fasta", ".fna"):
        names, seqs = [], []
        name = None
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    name = line[1:].split()[0] if len(line) > 1 else "q%d" % len(names)
                    names.append(name)
                    seqs.append("")
                elif name is not None:
                    seqs[-1] += line
        df = pd.DataFrame({"name": names, "seq": [s.upper().replace("U", "T") for s in seqs]})
    else:
        tbl = pd.read_csv(path)
        if column not in tbl.columns:
            raise SystemExit("ERROR: input CSV 缺少列 %r（现有列：%s）" % (column, list(tbl.columns)))
        seqs = tbl[column].astype(str).str.upper().str.replace("U", "T")
        df = pd.DataFrame({"name": ["q%d" % i for i in range(len(tbl))], "seq": seqs})
        df = pd.concat([tbl.reset_index(drop=True), df], axis=1)
    # 19nt 校验（警告不阻断）
    bad = df[df["seq"].str.len() != 19]
    if len(bad):
        print("WARNING: %d 条候选长度 != 19nt，仍照常送检" % len(bad), file=sys.stderr)
    if df["seq"].duplicated().any():
        print("WARNING: 存在重复序列，重复者按同名处理", file=sys.stderr)
    return df


def write_fasta(df, out_fa):
    with open(out_fa, "w", encoding="utf-8") as fh:
        for i, row in df.iterrows():
            fh.write(">%s\n%s\n" % (row["name"], row["seq"]))


# ---------- BLAST ----------

def run_blastn(query_fa, out_tsv, blast_bin, db, threads=8):
    """与归档运行完全同参：blastn-short / word_size 7 / evalue 1000 / 12 列输出。"""
    blastn = os.path.join(blast_bin, "blastn.exe")
    if not os.path.exists(blastn):
        blastn = os.path.join(blast_bin, "blastn")  # Linux
    cmd = [blastn, "-db", db, "-query", query_fa, "-out", out_tsv,
           "-outfmt", "6 " + " ".join(BLAST_COLS),
           "-task", "blastn-short", "-word_size", "7", "-evalue", "1000",
           "-num_threads", str(threads)]
    print("[BLAST]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def make_db(transcripts_fa, db_prefix, blast_bin):
    makeblastdb = os.path.join(blast_bin, "makeblastdb.exe")
    if not os.path.exists(makeblastdb):
        makeblastdb = os.path.join(blast_bin, "makeblastdb")
    subprocess.run([makeblastdb, "-in", transcripts_fa, "-dbtype", "nucl",
                    "-out", db_prefix], check=True)
    print("[makeblastdb] done ->", db_prefix)


# ---------- 判据 ----------

def apply_tolerant(hits, df, target_gene, min_paired, max_mismatch, require_no_gaps):
    """返回 (flag 布尔 Series[与 df 同序], 命中明细 DataFrame)。"""
    if len(hits) == 0:
        return pd.Series(False, index=df.index), hits

    non = hits[~hits["sseqid"].str.contains(target_gene, case=False, na=False)].copy()
    non["minus"] = non["sstart"] > non["send"]
    non["paired"] = non["length"] - non["mismatch"]
    ok = (non["minus"]) & (non["paired"] >= min_paired) & (non["mismatch"] <= max_mismatch)
    if require_no_gaps:
        ok &= non["gapopen"] == 0
    offending = non[ok]

    flag = df["name"].isin(set(offending["qseqid"]))
    return flag, offending


# ---------- 主流程 ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="候选 CSV（含 siRNA 列）或 FASTA")
    ap.add_argument("--column", default="siRNA", help="CSV 中序列列名（默认 siRNA）")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="blast_offtarget")
    ap.add_argument("--blast-bin", default=DEFAULT_BLAST_BIN)
    ap.add_argument("--db", default=DEFAULT_DB, help="BLAST 库前缀")
    ap.add_argument("--target-gene", default=DEFAULT_TARGET_GENE)
    ap.add_argument("--min-paired", type=int, default=17)
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--require-no-gaps", action="store_true",
                    help="额外要求 gapopen==0（正式归档口径未启用）")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--make-db", dest="make_db", default=None,
                    help="传入转录本 FASTA，先建库到 --db 指定前缀")
    args = ap.parse_args()

    if args.make_db:
        make_db(args.make_db, args.db, args.blast_bin)

    os.makedirs(args.out_dir, exist_ok=True)
    df = read_candidates(args.input, args.column)
    print("candidates:", len(df))

    fa = os.path.join(args.out_dir, args.prefix + "_queries.fa")
    write_fasta(df, fa)

    raw_tsv = os.path.join(args.out_dir, args.prefix + "_blast_raw.tsv")
    run_blastn(fa, raw_tsv, args.blast_bin, args.db, args.threads)

    hits = pd.read_csv(raw_tsv, sep="\t", header=None, names=BLAST_COLS)
    hits = hits[hits["length"] >= 16].copy()          # 审计留痕口径（与归档一致）
    hits_tsv = os.path.join(args.out_dir, args.prefix + "_hits_all.tsv")
    hits.to_csv(hits_tsv, sep="\t", index=False)
    print("raw hits (>=16nt):", len(hits))

    flag, offending = apply_tolerant(hits, df, args.target_gene,
                                     args.min_paired, args.max_mismatch,
                                     args.require_no_gaps)

    flagged = df[flag].drop(columns=["name"])
    clean = df[~flag].drop(columns=["name"])
    f_csv = os.path.join(args.out_dir, args.prefix + "_flagged.csv")
    c_csv = os.path.join(args.out_dir, args.prefix + "_clean.csv")
    flagged.to_csv(f_csv, index=False)
    clean.to_csv(c_csv, index=False)

    off_genes = (offending["sseqid"].str.split("|").str[-1]
                 .value_counts().head(10) if len(offending) else {})
    stats = {
        "criterion": "tolerant v3 (minus & paired>=%d & mismatch<=%d%s)"
                     % (args.min_paired, args.max_mismatch,
                        " & gapopen==0" if args.require_no_gaps else ""),
        "blast_task": "blastn-short", "word_size": 7, "evalue": 1000,
        "db": args.db, "target_gene": args.target_gene,
        "n_candidates": int(len(df)),
        "n_hits_all_ge16": int(len(hits)),
        "n_flagged_offtarget": int(flag.sum()),
        "n_clean": int(~flag.sum()),
        "top_offtarget_genes": {str(k): int(v) for k, v in off_genes.items()},
    }
    j = os.path.join(args.out_dir, args.prefix + "_stats.json")
    with open(j, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=2)

    print("\n=== 判定结果（tolerant 口径）===")
    print("候选 %d | 命中(>=16nt) %d 行 | 判脱靶 %d | 通过 %d"
          % (len(df), len(hits), int(flag.sum()), int((~flag).sum())))
    if len(off_genes):
        print("脱靶基因 top:", dict(off_genes))
    print("saved:", f_csv)
    print("saved:", c_csv)
    print("saved:", hits_tsv)
    print("saved:", j)


if __name__ == "__main__":
    main()
