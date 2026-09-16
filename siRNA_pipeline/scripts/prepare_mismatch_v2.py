# -*- coding: utf-8 -*-
"""Prepare mismatch dataset v2: geometry repair + verified annotations + split-safe ids.

修复 prepare_mismatch_dataset.py 暴露的问题（research-only，不改生产输出）：

  1. 几何：Holen 等 21nt 双链被截成 19nt 后 guide/target 不再互补。这里改为**二维对齐搜索**：
     在 mRNA 上的窗口偏移 w_off 与 guide 起始偏移 g_off 上枚举，取"反平行不配对位数最少"的
     组合（同分时优先与论文注释一致的、偏移更小的）。
  2. 注释：以序列实算的错配位点/类型为准（verified_*）；论文注释保留原列并记录是否一致。
     Kini 的 "Guide 1" 这类无法解析的注释 → 用实算结果补齐（不再当成 0 错配）。
  3. 主键：row_id 全局唯一化（dataset|原id），避免跨数据集的键冲突。
  4. 划分：外层 Leave-One-Dataset-Out；内层 GroupKFold(by family group)；
     另附 Leave-One-Family-Out 折；并断言"任何 group 不同时出现在 train/test"。

输出：outputs/analysis/mismatch/{mismatch_curated_v2.csv, mismatch_splits_v2.json, mismatch_qc_v2.json}
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupKFold

COMP = {"A": "U", "U": "A", "C": "G", "G": "C"}
SOURCE_FILES = (
    "Ohnishi_2008_78_mutated_siRNA_normalized.csv",
    "Holen_2005_18_mutated_siRNA_final_flanked.csv",
    "Kini_2009_terminal_mismatch_siRNA.csv",
    "mismatch_wt_references.csv",
)
ANN_RE = re.compile(r"([ACGU]):([ACGU])\((\d+)\)", re.IGNORECASE)


def norm(value: str) -> str:
    return re.sub(r"\s+", "", str(value).upper()).replace("T", "U")


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("AUCG", "UAGC"))[::-1]


def mismatch_sites(guide19: str, window19: str) -> tuple[list[int], list[str]]:
    """guide[p-1] ↔ window[19-p]（反平行）；返回 1-based 位置与类型。"""
    pos, types = [], []
    for p in range(1, 20):
        gb, wb = guide19[p - 1], window19[19 - p]
        if gb != COMP.get(wb, ""):
            pos.append(p)
            types.append(f"{gb}:{wb}")
    return pos, types


def best_alignment(guide_raw: str, mrna: str, target_raw: str,
                   ann_positions: list[int]) -> dict | None:
    """枚举 (窗口偏移, guide 偏移)，取不配对数最少者；返回最优 19nt 双链与元数据。"""
    start = mrna.find(target_raw)
    if start < 0:
        start = 19 if len(mrna) >= 57 else -1
    if start < 0:
        return None
    cand = []
    for w_off in range(max(0, start - 3), min(len(mrna) - 19, start + 3) + 1):
        window = mrna[w_off:w_off + 19]
        if len(window) != 19:
            continue
        for g_off in range(0, max(1, len(guide_raw) - 19 + 1)):
            guide = guide_raw[g_off:g_off + 19]
            if len(guide) != 19:
                continue
            pos, types = mismatch_sites(guide, window)
            score = (len(pos),                      # 不配对数最少
                     0 if (ann_positions and sorted(pos) == sorted(ann_positions)) else 1,
                     abs(w_off - start), g_off)    # 偏移更小优先
            cand.append((score, w_off, g_off, guide, window, pos, types))
    if not cand:
        return None
    cand.sort(key=lambda x: x[0])
    score, w_off, g_off, guide, window, pos, types = cand[0]
    return {"w_off": w_off, "g_off": g_off, "guide": guide, "window": window,
            "positions": pos, "types": types, "n_unpaired": len(pos)}


def main() -> None:
    root = Path(__file__).parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path,
                    default=root.parent / "数据集" / "最终数据" / "错配数据集")
    ap.add_argument("--out-dir", type=Path,
                    default=root / "outputs" / "analysis" / "mismatch")
    args = ap.parse_args()

    frames = []
    digest = hashlib.sha256()
    for name in SOURCE_FILES:
        path = args.input_dir / name
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
        frame = pd.read_csv(path, dtype=str, encoding="utf-8-sig").fillna("")
        frame["source_file"] = name
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)

    rows, unresolved = [], []
    for _, row in raw.iterrows():
        ann = [{"pos": int(p), "type": f"{a.upper()}:{b.upper()}"}
               for a, b, p in ANN_RE.findall(str(row["mismatch"]))]
        ann_positions = sorted({x["pos"] for x in ann})
        guide_raw = norm(row["sirna_sequence"])
        target_raw = norm(row["target_sequence"])
        mrna_raw = norm(row["mRNA"])

        align = best_alignment(guide_raw, mrna_raw, target_raw, ann_positions)
        dataset = str(row["dataset"])
        uid = "%s|%s" % (dataset, row["id"])
        if align is None or align["w_off"] < 19 or align["w_off"] + 38 > len(mrna_raw):
            unresolved.append({"uid": uid, "sirna_id": row["sirna_id"]})
            continue
        mrna57 = mrna_raw[align["w_off"] - 19: align["w_off"] + 38]
        if len(mrna57) != 57 or mrna57[19:38] != align["window"]:
            unresolved.append({"uid": uid, "sirna_id": row["sirna_id"]})
            continue

        verified_pos = sorted(align["positions"])
        verified_types = align["types"]
        reason = ""
        if not ann_positions:
            reason = "annotation_unparsed_use_sequence"
        elif verified_pos != ann_positions:
            reason = "annotation_disagrees_use_sequence"

        family = re.sub(r"-.*$", "", str(row["sirna_id"]))
        rows.append({
            "uid": uid,
            "row_id": str(row["id"]),
            "dataset": dataset,
            "source_file": row["source_file"],
            "sirna_id": row["sirna_id"],
            "family": family,
            "group_id": "%s|%s" % (dataset, family),
            "guide": align["guide"],
            "target": align["window"],
            "mRNA": mrna57,
            "label": float(row["silencing_efficiency_norm"]),
            "mismatch_count": len(verified_pos),
            "declared_mismatch_count": int(float(row["mismatch_count"] or 0)),
            "annotation_raw": row["mismatch"],
            "annotation_positions": json.dumps(ann_positions),
            "verified_positions": json.dumps(verified_pos),
            "verified_types": json.dumps(verified_types),
            "guide_offset": align["g_off"],
            "window_offset": align["w_off"],
            "mRNA_window_offset": align["w_off"],
            "repair_note": reason,
        })

    data = pd.DataFrame(rows)
    # ---- 几何可信资格：实测错配数 ≤3（单/双/三错配）且与论文声明数一致或声明为空 ----
    data["declared_vs_verified_ok"] = (
        (data["declared_mismatch_count"] == 0) | (data["declared_mismatch_count"] == data["mismatch_count"])
    )
    data["model_eligible_v2"] = (data["mismatch_count"] <= 3) & data["declared_vs_verified_ok"]
    excluded = data.loc[~data["model_eligible_v2"], ["uid", "dataset", "sirna_id",
                                                     "mismatch_count", "declared_mismatch_count"]]
    data = data.reset_index(drop=True)
    # ---- 划分（仅对几何可信子集）：外层 LOSO；内层 GroupKFold；另附 leave-one-family-out ----
    elig = data[data["model_eligible_v2"]].reset_index(drop=True)
    splits = {"outer_loso": {}, "inner_group_kfold": {}, "leave_one_family_out": [],
              "eligible_rows": len(elig), "eligible_groups": int(elig["group_id"].nunique()),
              "eligible_by_dataset": elig.groupby("dataset").size().to_dict(),
              "source_hash": digest.hexdigest()}
    for ds in sorted(elig["dataset"].unique()):
        te = elig[elig["dataset"] == ds]
        tr = elig[elig["dataset"] != ds]
        inner = []
        n_splits = min(5, tr["group_id"].nunique())
        if n_splits >= 2:
            gkf = GroupKFold(n_splits=n_splits)
            for fold, (a, b) in enumerate(gkf.split(tr, groups=tr["group_id"])):
                inner.append({"fold": fold, "train_uids": tr.iloc[a]["uid"].tolist(),
                              "valid_uids": tr.iloc[b]["uid"].tolist()})
        splits["outer_loso"][ds] = {"train_uids": tr["uid"].tolist(),
                                    "test_uids": te["uid"].tolist(),
                                    "test_groups": int(te["group_id"].nunique())}
        splits["inner_group_kfold"][ds] = inner
        # 组完整性断言
        assert not (set(tr["group_id"]) & set(te["group_id"])), "group 泄漏：%s" % ds
    for g in sorted(elig["group_id"].unique()):
        va = elig[elig["group_id"] == g]
        splits["leave_one_family_out"].append(
            {"group": g, "train_uids": elig[elig["group_id"] != g]["uid"].tolist(),
             "valid_uids": va["uid"].tolist()})

    qc = {
        "total_rows_all": len(data),
        "model_eligible_v2_rows": int(data["model_eligible_v2"].sum()),
        "excluded_by_geometry": excluded.to_dict(orient="records"),
        "by_dataset": data.groupby("dataset").size().to_dict(),
        "eligible_by_dataset": data[data["model_eligible_v2"]].groupby("dataset").size().to_dict(),
        "groups": int(data["group_id"].nunique()),
        "eligible_groups": int(data[data["model_eligible_v2"]]["group_id"].nunique()),
        "families_by_dataset": data.groupby("dataset")["family"].nunique().to_dict(),
        "repaired_rows": int((data["repair_note"] != "").sum()),
        "repair_notes": data.loc[data["repair_note"] != "", "repair_note"].value_counts().to_dict(),
        "offset_shift_guide": data["guide_offset"].value_counts().to_dict(),
        "offset_shift_window": (data["window_offset"] - 19).value_counts().to_dict(),
        "mismatch_count_hist": data["mismatch_count"].value_counts().sort_index().to_dict(),
        "unresolved_rows": unresolved,
        "source_hash": digest.hexdigest(),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.out_dir / "mismatch_curated_v2.csv", index=False, encoding="utf-8-sig")
    (args.out_dir / "mismatch_splits_v2.json").write_text(
        json.dumps(splits, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.out_dir / "mismatch_qc_v2.json").write_text(
        json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
