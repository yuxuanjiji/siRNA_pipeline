# -*- coding: utf-8 -*-
"""v3.1 -> v3.2: drop internal columns, make the normalisation scope explicit.

Two defects found in the v3.1 final-readiness audit:
  1. Twenty internal `_`-prefixed scratch columns leaked into the CSV because
     the header was taken from the row dicts after they had been annotated.
     They are empty in every row and must not ship.
  2. `silencing_efficiency_norm` had been min-max normalised *per cohort*
     (Ohnishi and Sierant both span exactly [0,1]), so it is not comparable
     across papers. It is renamed to make that scope explicit, and a genuinely
     global normalisation is added alongside it. `silencing_efficiency_raw`
     stays the primary modelling target.

Input : 数据集/错配siRNA_guide链_论文实测效率_合并_v3_1.csv
Output: 数据集/错配siRNA_guide链_论文实测效率_合并_v3_2.csv
"""
from __future__ import annotations

import csv
import io
import os
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\Polareasterlies\Desktop\生科挑战赛"
SRC = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_1.csv")
OUT = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_2.csv")


def main():
    rows = list(csv.DictReader(io.StringIO(open(SRC, "rb").read().decode("utf-8-sig"))))
    print("读入 v3.1:", len(rows), "行")

    junk = [c for c in rows[0] if c.startswith("_")]
    print("剔除内部空列:", len(junk))
    keep = [c for c in rows[0] if not c.startswith("_")]
    for r in rows:
        for c in junk:
            r.pop(c, None)

    # rename the per-cohort normalised column; keep its values untouched
    keep = ["silencing_efficiency_norm_cohort" if c == "silencing_efficiency_norm" else c
            for c in keep]
    for r in rows:
        r["silencing_efficiency_norm_cohort"] = r.pop("silencing_efficiency_norm")

    # genuinely global scaling over all 146 rows, for comparison
    raws = [float(r["silencing_efficiency_raw"]) for r in rows]
    lo, hi = min(raws), max(raws)
    avg = sum(raws) / len(raws)
    var = sum((x - avg) ** 2 for x in raws) / len(raws)
    sd = var ** 0.5
    for r in rows:
        v = float(r["silencing_efficiency_raw"])
        r["silencing_efficiency_norm_global"] = round((v - lo) / (hi - lo), 6)

    out_cols = keep + ["silencing_efficiency_norm_global"]
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols)
        w.writeheader()
        w.writerows(rows)

    print("全局归一化参数: min={:.3f} max={:.3f} mean={:.3f} sd={:.3f}".format(lo, hi, avg, sd))
    print()
    print("=== 输出核对 ===")
    back = list(csv.DictReader(io.StringIO(open(OUT, "rb").read().decode("utf-8-sig"))))
    hdr = list(back[0].keys())
    print("  行数:", len(back), "| 列数:", len(hdr))
    print("  残留 _ 列:", [c for c in hdr if c.startswith("_")] or "无")
    empty = [c for c in hdr if sum(1 for r in back if (r[c] or "").strip()) == 0]
    print("  全空列:", empty or "无")
    print("  非满填充列:", {c: sum(1 for r in back if (r[c] or '').strip()) for c in hdr
                          if sum(1 for r in back if (r[c] or '').strip()) < len(back)} or "无")
    print()
    print("=== 三个目标列的取值范围对照 ===")
    for c in ("silencing_efficiency_raw", "silencing_efficiency_norm_cohort",
              "silencing_efficiency_norm_global"):
        v = [float(r[c]) for r in back]
        print("  {:<38} [{:.3f}, {:.3f}]".format(c, min(v), max(v)))
    print()
    print("  按数据集看 norm_cohort 是否仍为各 cohort 满量程:")
    for ds in sorted(set(r["dataset"] for r in back)):
        v = [float(r["silencing_efficiency_norm_cohort"]) for r in back if r["dataset"] == ds]
        print("     {:<16} [{:.3f}, {:.3f}]".format(ds, min(v), max(v)))
    print()
    print("写入:", OUT)


if __name__ == "__main__":
    main()
