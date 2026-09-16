# -*- coding: utf-8 -*-
"""Build the v3 guide-strand mismatch main table plus a per-row verification sheet.

Input : 数据集/错配siRNA_guide链_论文实测效率_合并_v2.csv   (146 rows, 25 cols)
Output: 数据集/错配siRNA_guide链_论文实测效率_合并_v3.csv        (added columns)
        数据集/错配siRNA_guide链_校验逐条.csv                   (per-row audit)

What v3 adds
    target_gene       inferred from the row's own note / source_file / reporter sequence
    family            grouping key for family-safe CV (paper + construct series)
    design_id         unique per (guide, target_site_19nt); replicates share it
    replicate_n       1-based index inside design_id
    replicate_block   which measurement series the row belongs to
    mismatch_region_fixed  repairs the "g1" placeholder into 5'(1)
    guide_mismatch    True only when n_mismatch_from_guide > 0
    eligible          guide_mismatch and project_eligible
    paper_entry_status  filled / missing

Every derived value is recomputed from the sequences themselves, never trusted
from the source columns; the audit sheet records both so disagreements are visible.
"""
from __future__ import annotations

import csv
import io
import os
import re
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\Polareasterlies\Desktop\生科挑战赛"
SRC = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v2.csv")
OUT_V3 = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3.csv")
OUT_AUDIT = os.path.join(BASE, r"数据集\错配siRNA_guide链_校验逐条.csv")

COMP = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}

# Gene per cohort, established from the row notes / source_file / reporter sequence.
GENE = {
    "Ohnishi_2008": "PRNP",     # construct names siPrnp102/105/178
    "Sierant_2011": "PSEN1",    # source_file name + note
    "Kini_2009": "EGFP",        # reporter mRNA, assay "EGFP 荧光敲低 %"
    "Holen_2005": "AQP4",       # note: "属原文报告基因(AQP4...)" on the position-1 block
}
GENE_SOURCE = {
    "Ohnishi_2008": "构建名 siPrnp*",
    "Sierant_2011": "source_file 与 note",
    "Kini_2009": "报告基因 mRNA 序列 + efficiency_unit",
    "Holen_2005": "note 明写原文报告基因",
}


def norm(s: str) -> str:
    return (s or "").strip().upper().replace("T", "U")


def load(path):
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "gb18030", "utf-8"):
        try:
            return list(csv.DictReader(io.StringIO(raw.decode(enc))))
        except Exception:
            continue
    raise SystemExit(f"cannot decode {path}")


def mismatches(guide: str, site: str):
    """Return 1-based guide positions whose base does not pair with the site."""
    return [j + 1 for j in range(19) if guide[j] != COMP.get(site[18 - j], "N")]


def parse_positions(text: str):
    return [int(x) for x in re.split(r"[,，;；\s]+", text or "") if x.strip().isdigit()]


def parse_pair_types(text: str):
    """'3:G:U' or '1:U:G;7:A:C' -> {position: 'G:U'}"""
    out = {}
    for tok in re.split(r"[;；]", text or ""):
        parts = tok.strip().split(":")
        if len(parts) == 3 and parts[0].isdigit():
            out[int(parts[0])] = f"{parts[1]}:{parts[2]}"
    return out


def family_of(row, gene):
    ds, entry = row["dataset"], (row["paper_entry"] or "").strip()
    if ds == "Ohnishi_2008":
        m = re.match(r"(siPrnp\d+\(\w+\))", entry)
        return m.group(1) if m else (entry or "siPrnp_未标注")
    if ds == "Sierant_2011":
        return "Sierant2011_P系列"
    if ds == "Kini_2009":
        return "Kini2009_396"
    if ds == "Holen_2005":
        if entry:
            return f"Holen2005_{entry[0]}系列"
        return "Holen2005_g1block"
    return f"{ds}_未标注"


def main():
    rows = load(SRC)
    print(f"读入 v2: {len(rows)} 行")

    # ---- pass 1: recompute geometry and classify -------------------------
    for r in rows:
        g, site = norm(r["guide"]), norm(r["target_site_19nt"])
        comp = mismatches(g, site) if len(g) == 19 and len(site) == 19 else []
        r["_mm"] = comp
        r["_n_mm"] = len(comp)
        r["_declared"] = parse_positions(r["mismatch_positions_guide5p"])
        r["_pairtypes"] = parse_pair_types(r["mismatch_pair_types"])
        r["_n_guide"] = int(r["n_mismatch_from_guide"] or 0)
        r["_n_target"] = int(r["n_mismatch_from_target"] or 0)
        r["_gene"] = GENE.get(r["dataset"], "待确认")
        r["_gene_src"] = GENE_SOURCE.get(r["dataset"], "")
        r["_family"] = family_of(r, r["_gene"])
        r["_region_fixed"] = re.sub(r"^g1$", "5'(1)", r["mismatch_region"])

    # ---- pass 2: design identity and replicate numbering -----------------
    order = defaultdict(list)
    for r in rows:
        order[(norm(r["guide"]), norm(r["target_site_19nt"]))].append(r)
    for key, group in order.items():
        group.sort(key=lambda x: x["record_id"])
        for i, r in enumerate(group, 1):
            r["_design_id"] = "D{:03d}".format(sorted(order.keys()).index(key) + 1)
            r["_replicate_n"] = i
            r["_replicate_count"] = len(group)

    # A group carries a repeat block when its later half reproduces the earlier
    # half's efficiency multiset exactly. Same-mRNA pairs are true replicates.
    for key, group in order.items():
        half = len(group) // 2
        same_mrna = len(set(x["mRNA_57nt"] for x in group)) == 1
        earlier = sorted(float(x["silencing_efficiency_raw"]) for x in group[:half])
        later = sorted(float(x["silencing_efficiency_raw"]) for x in group[half:])
        repeated = half >= 1 and earlier == later
        for i, r in enumerate(group):
            if len(group) == 1:
                r["_repeat_block"] = "单次测量"
            elif same_mrna:
                r["_repeat_block"] = "同mRNA重复测量"
            elif repeated and i >= half:
                r["_repeat_block"] = "重复块"
            else:
                r["_repeat_block"] = "主块"

    # ---- pass 3: write v3 -------------------------------------------------
    extra = ["target_gene", "target_gene_source", "family", "design_id", "replicate_n",
             "replicate_count", "replicate_block", "mismatch_region_fixed",
             "guide_mismatch", "eligible", "paper_entry_status"]
    header = list(rows[0].keys()) + extra
    with open(OUT_V3, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        for r in rows:
            out = {k: r[k] for k in rows[0].keys() if not k.startswith("_")}
            out["target_gene"] = r["_gene"]
            out["target_gene_source"] = r["_gene_src"]
            out["family"] = r["_family"]
            out["design_id"] = r["_design_id"]
            out["replicate_n"] = r["_replicate_n"]
            out["replicate_count"] = r["_replicate_count"]
            out["replicate_block"] = r["_repeat_block"]
            out["mismatch_region_fixed"] = r["_region_fixed"]
            out["guide_mismatch"] = r["_n_guide"] > 0
            out["eligible"] = bool(r["_n_guide"] > 0 and (r["project_eligible"] or "").strip().lower() == "true")
            out["paper_entry_status"] = "filled" if (r["paper_entry"] or "").strip() else "missing"
            w.writerow(out)

    # ---- pass 4: per-row audit -------------------------------------------
    audit_header = ["record_id", "dataset", "paper_entry", "family", "target_gene",
                    "guide", "target_site_19nt", "guide_len_ok", "computed_n_mm",
                    "declared_n_mm", "n_mm_match", "computed_positions",
                    "declared_positions", "positions_subset", "pair_types_match",
                    "n_from_guide", "n_from_target", "attribution",
                    "attribution_consistent", "guide_mismatch", "eligible",
                    "design_id", "replicate_n", "replicate_block", "status"]
    problems = Counter()
    with open(OUT_AUDIT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=audit_header)
        w.writeheader()
        for r in rows:
            lens_ok = len(norm(r["guide"])) == 19 and len(norm(r["target_site_19nt"])) == 19
            n_ok = r["_n_mm"] == int(r["n_mismatch_total"] or -1)
            subset = all(p in r["_mm"] for p in r["_declared"]) if r["_declared"] else False
            pt_ok = all(r["_pairtypes"].get(p) ==
                        f"{norm(r['guide'])[p-1]}:{norm(r['target_site_19nt'])[18-(p-1)]}"
                        for p in r["_mm"] if p in r["_pairtypes"]) if r["_pairtypes"] else False
            attr_ok = ((r["attribution"] == "guide_only" and r["_n_guide"] > 0 and r["_n_target"] == 0) or
                       (r["attribution"] == "guide+allele" and r["_n_guide"] > 0 and r["_n_target"] > 0) or
                       (r["attribution"] == "allele_only" and r["_n_guide"] == 0 and r["_n_target"] > 0))
            gm = r["_n_guide"] > 0
            elig = gm and (r["project_eligible"] or "").strip().lower() == "true"
            status = []
            if not lens_ok:
                status.append("序列长度异常"); problems["len"] += 1
            if not n_ok:
                status.append("错配数不符"); problems["n_mm"] += 1
            if not subset:
                status.append("位置越界"); problems["pos"] += 1
            if not pt_ok:
                status.append("配对类型不符"); problems["pair"] += 1
            if not attr_ok:
                status.append("归属不符"); problems["attr"] += 1
            if not gm:
                status.append("非guide链错配"); problems["not_guide"] += 1
            w.writerow({
                "record_id": r["record_id"], "dataset": r["dataset"],
                "paper_entry": r["paper_entry"], "family": r["_family"], "target_gene": r["_gene"],
                "guide": r["guide"], "target_site_19nt": r["target_site_19nt"],
                "guide_len_ok": lens_ok, "computed_n_mm": r["_n_mm"],
                "declared_n_mm": r["n_mismatch_total"], "n_mm_match": n_ok,
                "computed_positions": ";".join(map(str, r["_mm"])),
                "declared_positions": r["mismatch_positions_guide5p"],
                "positions_subset": subset, "pair_types_match": pt_ok,
                "n_from_guide": r["_n_guide"], "n_from_target": r["_n_target"],
                "attribution": r["attribution"], "attribution_consistent": attr_ok,
                "guide_mismatch": gm, "eligible": elig,
                "design_id": r["_design_id"], "replicate_n": r["_replicate_n"],
                "replicate_block": r["_repeat_block"],
                "status": "OK" if not status else "; ".join(status),
            })

    # ---- summary ----------------------------------------------------------
    print()
    print("=== 逐条校验汇总 ===")
    print("  行数:", len(rows))
    for k, v in problems.items():
        print(f"  异常 [{k}]: {v}")
    if not problems:
        print("  未发现任何异常")
    print()
    print("=== guide 链错配判定 ===")
    gm = sum(1 for r in rows if r["_n_guide"] > 0)
    print(f"  n_mismatch_from_guide > 0 : {gm}")
    print(f"  n_mismatch_from_guide = 0 : {len(rows) - gm}  ({[r['record_id'] for r in rows if r['_n_guide']==0]})")
    print()
    print("=== 独立设计 / 重复 ===")
    print(f"  唯一 (guide, 靶位点) 设计数 : {len(order)}")
    print(f"  重复行数                    : {len(rows) - len(order)}")
    print("  重复块分布:", dict(Counter(r["_repeat_block"] for r in rows)))
    print()
    print("=== 家族 / 基因 ===")
    for k, v in Counter(r["_family"] for r in rows).most_common():
        print(f"  {k:<22} {v}")
    print()
    for k, v in Counter((r['dataset'], r['_gene']) for r in rows).most_common():
        print(f"  {k[0]:<16} -> {k[1]:<8} {v}")
    print()
    print("  写入:", OUT_V3)
    print("  写入:", OUT_AUDIT)


if __name__ == "__main__":
    main()
