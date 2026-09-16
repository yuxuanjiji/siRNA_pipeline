# -*- coding: utf-8 -*-
"""v3 -> v3.1: neutral provenance labels and two family boundaries.

Applied after review feedback:
  1. Do not privilege the first-seen block. Rows that belong to a repeated
     measurement group are labelled "重复-待核" on both sides; which mRNA
     version is authentic is undecidable from the table alone.
  2. Replace the ad-hoc paper_entry-prefix families with the agreed boundary
     (experimental block x genetic locus), and also emit the coarser variant
     the reviewer proposed, so both can be compared side by side.

Input : 数据集/错配siRNA_guide链_论文实测效率_合并_v3.csv
Output: 数据集/错配siRNA_guide链_论文实测效率_合并_v3_1.csv
"""
from __future__ import annotations

import csv
import io
import os
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\Polareasterlies\Desktop\生科挑战赛"
SRC = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3.csv")
OUT = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_1.csv")

# block x locus  ->  merged w+m because both sit on the identical 19nt site and
# the identical 57nt context; Ohnishi stays split because its four constructs
# share no target site at all.
LOCUS_FAMILY = {
    "Holen2005_w系列": "Holen2005_主座",
    "Holen2005_m系列": "Holen2005_主座",
    "Holen2005_g1block": "Holen2005_g1报告基因块",
}
# the coarser view: one family per paper's reporter system
COARSE = {"Ohnishi_2008": "Ohnishi2008_PRNP体系"}


def norm(s: str) -> str:
    return (s or "").strip().upper().replace("T", "U")


def main():
    raw = open(SRC, "rb").read()
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    print("读入 v3:", len(rows), "行")

    groups = defaultdict(list)
    for r in rows:
        groups[(norm(r["guide"]), norm(r["target_site_19nt"]))].append(r)

    for r in rows:
        key = (norm(r["guide"]), norm(r["target_site_19nt"]))
        group = sorted(groups[key], key=lambda x: x["record_id"])
        n = len(group)
        idx = [x["record_id"] for x in group].index(r["record_id"]) + 1
        same_mrna = len(set(x["mRNA_57nt"] for x in group)) == 1
        if n == 1:
            r["_prov"] = "单次测量"
        elif same_mrna:
            r["_prov"] = "同mRNA重复测量"
        else:
            r["_prov"] = "重复-待核"        # both halves, no privileging
        r["_rep_n"], r["_rep_count"] = idx, n
        r["_fam_locus"] = LOCUS_FAMILY.get(r["family"], r["family"])
        r["_fam_coarse"] = COARSE.get(r["dataset"], r["_fam_locus"])

    cols = list(rows[0].keys())
    for c in ("provenance", "family_cv", "family_coarse"):
        if c not in cols:
            cols.append(c)
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            out = {k: v for k, v in r.items() if not k.startswith("_")}
            out["provenance"] = r["_prov"]
            out["family_cv"] = r["_fam_locus"]
            out["family_coarse"] = r["_fam_coarse"]
            w.writerow(out)

    print()
    print("=== provenance（不再区分主块/重复块）===")
    for k, v in Counter(r["_prov"] for r in rows).most_common():
        print("  {:<20} {}".format(k, v))
    dup = [r for r in rows if r["_prov"] == "重复-待核"]
    print("  '重复-待核' 行数: {}  (涉及 {} 个设计)".format(
        len(dup), len(set((norm(r['guide']), norm(r['target_site_19nt'])) for r in dup))))

    print()
    print("=== family_cv：实验块 × 基因座（8 族）===")
    seen = {}
    for r in rows:
        seen.setdefault(r["_fam_locus"], set()).add((norm(r["guide"]), norm(r["target_site_19nt"])))
    for k, v in sorted(seen.items(), key=lambda x: -len(x[1])):
        print("  {:<26} 设计数={}".format(k, len(v)))
    print("  合计设计数:", sum(len(v) for v in seen.values()), "| 族数:", len(seen))

    print()
    print("=== family_coarse：审阅者版本（5 族）===")
    seen2 = {}
    for r in rows:
        seen2.setdefault(r["_fam_coarse"], set()).add((norm(r["guide"]), norm(r["target_site_19nt"])))
    for k, v in sorted(seen2.items(), key=lambda x: -len(x[1])):
        print("  {:<26} 设计数={}".format(k, len(v)))
    print("  合计设计数:", sum(len(v) for v in seen2.values()), "| 族数:", len(seen2))
    print()
    print("写入:", OUT)


if __name__ == "__main__":
    main()
