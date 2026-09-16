"""Prepare the updated mismatch collection without changing production outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupKFold


MISMATCH_RE = re.compile(r"([ACGU]):([ACGU])\((\d+)\)", re.IGNORECASE)
SOURCE_FILES = (
    "Ohnishi_2008_78_mutated_siRNA_normalized.csv",
    "Holen_2005_18_mutated_siRNA_final_flanked.csv",
    "Kini_2009_terminal_mismatch_siRNA.csv",
    "mismatch_wt_references.csv",
)


def normalize_sequence(value: str) -> str:
    return re.sub(r"\s+", "", str(value).upper()).replace("T", "U")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("AUCG", "UAGC"))[::-1]


def parse_annotation(value: str) -> list[dict[str, object]]:
    return [
        {"guide_base": left.upper(), "target_base": right.upper(), "position": int(position)}
        for left, right, position in MISMATCH_RE.findall(str(value))
    ]


def compare_pair(guide: str, target: str) -> list[dict[str, object]]:
    if len(guide) != len(target):
        return []
    return [
        {"guide_base": guide[index], "target_base": target[index], "position": index + 1}
        for index, (guide_base, target_base) in enumerate(zip(guide, target))
        if guide_base != target_base
    ]


def source_hash(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def prepare(input_dir: Path) -> tuple[pd.DataFrame, dict]:
    frames = []
    paths = [input_dir / name for name in SOURCE_FILES]
    for path in paths:
        frame = pd.read_csv(path, dtype=str, encoding="utf-8-sig").fillna("")
        frame["source_file"] = path.name
        frames.append(frame)

    rows = []
    for _, row in pd.concat(frames, ignore_index=True).iterrows():
        guide_raw = normalize_sequence(row["sirna_sequence"])
        target_raw = normalize_sequence(row["target_sequence"])
        mrna_raw = normalize_sequence(row["mRNA"])
        normalization_action = "none"
        if len(guide_raw) >= 19 and len(mrna_raw) >= 57 and mrna_raw.find(target_raw) == 19:
            guide = guide_raw[:19]
            mrna = mrna_raw[:57]
            target = mrna[19:38]
            if len(guide_raw) != 19 or len(target_raw) != 19 or len(mrna_raw) != 57:
                normalization_action = "truncate_to_oligoformer_19_57_contract"
        else:
            guide, target, mrna = guide_raw, target_raw, mrna_raw
        annotation = parse_annotation(row["mismatch"])
        reasons = []
        geometry_ok = len(guide) == 19 and len(target) == 19 and len(mrna) == 57
        if not geometry_ok:
            reasons.append("geometry_not_19nt_57nt")
        if mrna.find(target) != 19:
            reasons.append("target_not_at_mrna_19_38")

        mismatch_count = int(row["mismatch_count"] or 0)
        if not annotation and mismatch_count > 0:
            reasons.append("mismatch_annotation_unparsed")

        oligoformer_eligible = geometry_ok and mrna.find(target) == 19
        model_eligible = oligoformer_eligible
        dataset = str(row["dataset"])
        family = re.sub(r"-.*$", "", str(row["sirna_id"]))
        rows.append(
            {
                "row_id": str(row["id"]),
                "dataset": dataset,
                "source_file": row["source_file"],
                "reference_source": row.get("reference_source", ""),
                "sirna_id": row["sirna_id"],
                "siRNA": guide,
                "target": target,
                "mRNA": mrna,
                "normalization_action": normalization_action,
                "label": float(row["silencing_efficiency_norm"]),
                "mismatch_count": mismatch_count,
                "mismatch_annotation": row["mismatch"],
                "mismatch_positions": json.dumps(sorted(int(x["position"]) for x in annotation)),
                "mismatch_types": json.dumps([f"{x['guide_base']}:{x['target_base']}" for x in annotation]),
                "oligoformer_eligible": oligoformer_eligible,
                "family": family,
                "group_id": f"{dataset}|{family}" if model_eligible else "",
                "model_eligible": model_eligible,
                "exclusion_reason": ";".join(reasons),
            }
        )

    result = pd.DataFrame(rows)
    eligible = result[result["model_eligible"]].copy()
    splits = {"outer_loso": {}, "inner_group_kfold": {}, "eligible_rows": len(eligible), "eligible_groups": eligible["group_id"].nunique()}
    for dataset in sorted(eligible["dataset"].unique()):
        test = eligible[eligible["dataset"] == dataset]
        train = eligible[eligible["dataset"] != dataset]
        inner = []
        n_splits = min(5, train["group_id"].nunique())
        if n_splits >= 2:
            splitter = GroupKFold(n_splits=n_splits)
            for fold, (train_idx, valid_idx) in enumerate(splitter.split(train, groups=train["group_id"])):
                inner.append({"fold": fold, "train_row_ids": train.iloc[train_idx]["row_id"].tolist(), "valid_row_ids": train.iloc[valid_idx]["row_id"].tolist()})
        splits["outer_loso"][dataset] = {"train_row_ids": train["row_id"].tolist(), "test_row_ids": test["row_id"].tolist()}
        splits["inner_group_kfold"][dataset] = inner
    return result, {"source_hash": source_hash(paths), "splits": splits}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).parents[2] / "数据集" / "最终数据" / "错配数据集")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame, metadata = prepare(args.input_dir)
    frame.to_csv(args.output_dir / "mismatch_curated.csv", index=False, encoding="utf-8-sig")
    (args.output_dir / "mismatch_splits.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    qc = {
        "total_rows": len(frame),
        "model_eligible_rows": int(frame["model_eligible"].sum()),
        "excluded_rows": int((~frame["model_eligible"]).sum()),
        "excluded_reasons": frame.loc[~frame["model_eligible"], "exclusion_reason"].value_counts().to_dict(),
        "rows_by_dataset": frame.groupby("dataset").size().to_dict(),
        "eligible_rows_by_dataset": frame[frame["model_eligible"]].groupby("dataset").size().to_dict(),
        "source_hash": metadata["source_hash"],
    }
    (args.output_dir / "mismatch_qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
