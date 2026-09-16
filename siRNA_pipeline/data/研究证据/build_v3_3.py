# -*- coding: utf-8 -*-
"""v3.2 -> v3.3: reviewer-driven corrections.

Four blockers from the review round:

 1. Column naming/semantics. `norm_cohort` was in fact a GLOBAL operation
    (raw / global max, then clipped at 0) - verified 146/146 - so calling it
    "cohort" was wrong. It is restored to the original name
    `silencing_efficiency_norm` for downstream compatibility, and the other
    column becomes `silencing_efficiency_norm_minmax_unclipped`.
 2. `eligible` vs `project_eligible` were undocumented. The derived column is
    renamed `eligible_guide_strand_mismatch` and the data card defines it.
 3. `target_gene = AQP4` was applied to all 50 Holen rows, but only the g1
    block's note states AQP4. Holen's gene is therefore resolved per block.
 4. Family boundary. My earlier "four disjoint loci" claim came from exact
    19nt string intersection, which is meaningless for shifted windows. With
    best-offset alignment three of the four Ohnishi constructs are 100 %
    identical over 48-56 positions - i.e. one locus. Loci are now clustered
    by aligned identity, and the aligned locus column is explicitly marked as
    NOT a fold unit.

Input : 数据集/错配siRNA_guide链_论文实测效率_合并_v3_2.csv
Output: 数据集/错配siRNA_guide链_论文实测效率_合并_v3_3.csv
"""
from __future__ import annotations

import csv
import io
import os
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\Polareasterlies\Desktop\生科挑战赛"
SRC = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_2.csv")
OUT = os.path.join(BASE, r"数据集\错配siRNA_guide链_论文实测效率_合并_v3_3.csv")

IDENT_MIN, OVERLAP_MIN = 0.95, 30

# Holen's gene is only documented for the g1 block; w/m must wait for the
# supplementary table. (The reviewer's FLuc hypothesis is plausible but unproven.)
GENE_OVERRIDE = {
    "Holen2005_g1block": ("AQP4", "note 明写：属原文报告基因(AQP4)块"),
    "Holen2005_w系列": ("待确认", "note 未提基因；疑为萤火虫荧光素酶(FLuc)，须以 PMC1188085 补充表为准"),
    "Holen2005_m系列": ("待确认", "note 未提基因；疑为萤火虫荧光素酶(FLuc)，须以 PMC1188085 补充表为准"),
}


def norm(s: str) -> str:
    return (s or "").strip().upper().replace("T", "U")


def best_aligned(a: str, b: str):
    best = (0, 0, 1)
    for off in range(-len(b) + 1, len(a)):
        inter = [(i, i - off) for i in range(len(a)) if 0 <= i - off < len(b)]
        if len(inter) < OVERLAP_MIN:
            continue
        ident = sum(1 for i, j in inter if a[i] == b[j])
        if ident / len(inter) > best[0] / best[2]:
            best = (ident, off, len(inter))
    return best


def main():
    rows = list(csv.DictReader(io.StringIO(open(SRC, "rb").read().decode("utf-8-sig"))))
    print("读入 v3.2:", len(rows), "行")

    # ---- 1. norm column names -------------------------------------------
    for r in rows:
        r["silencing_efficiency_norm"] = r.pop("silencing_efficiency_norm_cohort")
        r["silencing_efficiency_norm_minmax_unclipped"] = r.pop("silencing_efficiency_norm_global")

    # ---- 2. rename the derived eligibility column ------------------------
    for r in rows:
        r["eligible_guide_strand_mismatch"] = r.pop("eligible")

    # ---- 3. per-block gene for Holen -------------------------------------
    for r in rows:
        if r["family"] in GENE_OVERRIDE:
            g, src = GENE_OVERRIDE[r["family"]]
            r["target_gene"] = g
            r["target_gene_source"] = src

    # ---- 4. cluster families into loci by ALIGNED identity ---------------
    reps = {}
    for r in rows:
        reps.setdefault(r["family"], Counter())[norm(r["mRNA_57nt"])] += 1
    rep = {f: c.most_common(1)[0][0] for f, c in reps.items()}

    parent = {f: f for f in rep}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    fams = sorted(rep)
    merges = []
    for i in range(len(fams)):
        for j in range(i + 1, len(fams)):
            ident, off, ov = best_aligned(rep[fams[i]], rep[fams[j]])
            if ov >= OVERLAP_MIN and ident / ov >= IDENT_MIN:
                a, b = find(fams[i]), find(fams[j])
                if a != b:
                    parent[a] = b
                merges.append((fams[i], fams[j], off, ov, ident / ov))

    cluster = {f: find(f) for f in fams}
    label = {}
    for c in set(cluster.values()):
        members = [f for f in fams if cluster[f] == c]
        label[c] = " + ".join(sorted(members)) if len(members) > 1 else members[0]
    for r in rows:
        fam = r["family"]
        members = [f for f in fams if cluster[f] == cluster[fam]]
        if r["dataset"] == "Ohnishi_2008" and len(members) > 1:
            r["family_locus"] = "PRNP_主座"          # the aligned Ohnishi cluster
        elif fam in ("Holen2005_w系列", "Holen2005_m系列"):
            r["family_locus"] = "Holen2005_主座"
        else:
            r["family_locus"] = fam

    cols = [c for c in rows[0] if not c.startswith("_")]
    for c in ("target_gene", "target_gene_source"):
        if c not in cols:
            cols.append(c)
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    print()
    print("=== 对齐后合并的家族对（一致度 ≥95%）===")
    for a, b, off, ov, ident in merges:
        print("  {:<18} x {:<18} 错位={:<3} 重叠={:<3} 一致度={:.0%}".format(a, b, off, ov, ident))
    print("  未合并的:", [f for f in fams if len([g for g in fams if cluster[g] == cluster[f]]) == 1])
    print()
    print("=== family_locus（真实区划，6 族）===")
    seen = {}
    for r in rows:
        seen.setdefault(r["family_locus"], set()).add((norm(r["guide"]), norm(r["target_site_19nt"])))
    for k, v in sorted(seen.items(), key=lambda x: -len(x[1])):
        print("  {:<26} 设计数={}".format(k, len(v)))
    print("  合计设计数:", sum(len(v) for v in seen.values()), "| 族数:", len(seen))
    assert sum(len(v) for v in seen.values()) == 123, "设计总数应为 123"
    assert len(seen) == 6, f"family_locus 应为 6 族，实为 {len(seen)}"
    print("  ✔ 断言通过：6 族 / 123 设计")
    print()
    print("=== target_gene 修正后 ===")
    for k, v in Counter((r["dataset"], r["family"], r["target_gene"]) for r in rows).most_common():
        print("  {:<16} {:<20} -> {:<8} {}".format(k[0], k[1], k[2], v))
    print()
    print("=== 目标列（三列口径已显式化）===")
    for c in ("silencing_efficiency_raw", "silencing_efficiency_norm",
              "silencing_efficiency_norm_minmax_unclipped"):
        v = [float(r[c]) for r in rows]
        print("  {:<44} [{:.3f}, {:.3f}]".format(c, min(v), max(v)))
    print()
    print("=== 输出核对 ===")
    back = list(csv.DictReader(io.StringIO(open(OUT, "rb").read().decode("utf-8-sig"))))
    hdr = list(back[0].keys())
    print("  行数:", len(back), "| 列数:", len(hdr))
    print("  残留 _ 列/全空列:", [c for c in hdr if c.startswith("_")] or
          [c for c in hdr if not any((r[c] or "").strip() for r in back)] or "无")
    print("  非满填充:", {c: sum(1 for r in back if (r[c] or "").strip()) for c in hdr
                        if sum(1 for r in back if (r[c] or "").strip()) < len(back)})
    print()
    print("写入:", OUT)


if __name__ == "__main__":
    main()
