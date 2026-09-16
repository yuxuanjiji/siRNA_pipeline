# -*- coding: utf-8 -*-
"""v3.3 -> v3.4: drop the falsified family column, give the remaining family
columns an explicit three-level hierarchy, and freeze the schema.

`family_cv` (the 8-family version built from exact-string target-site
intersection) is falsified: best-offset alignment shows three of the four
Ohnishi constructs are shifted windows of ONE sequence region, so using that
column as a fold unit would leak the mRNA context across folds. It is removed.

Levels kept:
    family_construct  9  construction / experimental block label
    family_locus      6  real sequence regions after best-offset alignment
    family_coarse     5  paper-level reporter system (recommended CV unit)
"""
from __future__ import annotations

import csv
import io
import os
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\Polareasterlies\Desktop\生科挑战赛"
SRC = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_3.csv")
OUT = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_4.csv")

DROP = ["family_cv"]


def main():
    rows = list(csv.DictReader(io.StringIO(open(SRC, "rb").read().decode("utf-8-sig"))))
    print("读入 v3.3:", len(rows), "行 /", len(rows[0]), "列")

    present = [c for c in DROP if c in rows[0]]
    print("删除已证伪的列:", present)
    cols = [c for c in rows[0] if c not in DROP]
    cols = ["family_construct" if c == "family" else c for c in cols]
    for r in rows:
        r["family_construct"] = r.pop("family")
        for c in DROP:
            r.pop(c, None)

    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    back = list(csv.DictReader(io.StringIO(open(OUT, "rb").read().decode("utf-8-sig"))))
    hdr = list(back[0].keys())
    print()
    print("=== 输出核对 ===")
    print("  行数:", len(back), "| 列数:", len(hdr))
    print("  残留 _ 列 / 全空列:", [c for c in hdr if c.startswith("_")] or
          [c for c in hdr if not any((r[c] or "").strip() for r in back)] or "无")
    print("  非满填充:", {c: sum(1 for r in back if (r[c] or "").strip()) for c in hdr
                        if sum(1 for r in back if (r[c] or "").strip()) < len(back)})

    print()
    print("=== 三级族列 ===")
    for col in ("family_construct", "family_locus", "family_coarse"):
        fams = {}
        for r in back:
            fams.setdefault(r[col], set()).add((r["guide"], r["target_site_19nt"]))
        print("  {:<18} 族数={} 设计数合计={}".format(col, len(fams), sum(len(v) for v in fams.values())))
        for k, v in sorted(fams.items(), key=lambda x: -len(x[1])):
            print("        {:<24} {}".format(k, len(v)))

    print()
    print("=== 断言 ===")
    assert len(set(r["design_id"] for r in back)) == 123
    assert len(set(r["family_locus"] for r in back)) == 6
    assert len(set(r["family_coarse"] for r in back)) == 5
    assert len(set(r["family_construct"] for r in back)) == 9
    assert "family_cv" not in hdr
    assert sum(1 for r in back if r["eligible_guide_strand_mismatch"] == "True") == 142
    print("  ✔ 123 设计 / 6 locus / 5 coarse / 9 construct / 142 eligible / 无 family_cv")
    print()
    print("写入:", OUT)


if __name__ == "__main__":
    main()
