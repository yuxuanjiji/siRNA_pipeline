"""Prepare OligoFormer-compatible CSV/FASTA inputs for the research run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def write_fasta(sequences: list[str], output: Path, prefix: str) -> dict[str, str]:
    mapping = {}
    with output.open("w", encoding="utf-8") as handle:
        for index, sequence in enumerate(sorted(set(sequences))):
            name = f"{prefix}_{index}"
            mapping[sequence] = name
            handle.write(f">{name}\n{sequence}\n")
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated_td.csv")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parents[1] / "OligoFormer部分" / "data" / "mismatch_finetune")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input, dtype={"row_id": str})
    data = data[data["model_eligible"]].copy()
    group_wt = data[data["mismatch_count"] == 0].groupby("group_id")["label"].first().to_dict()
    family_wt = data[data["mismatch_count"] == 0].groupby("family")["label"].first().to_dict()
    data["wt_label"] = data["group_id"].map(group_wt).fillna(data["family"].map(family_wt))
    if data["wt_label"].isna().any():
        raise ValueError("Every row must have an auditable WT label before OligoFormer fine-tuning")
    data["delta_label"] = data["label"] - data["wt_label"]
    siRNA_names = write_fasta(data["siRNA"].tolist(), args.output_dir / "siRNA.fa", "siRNA")
    mRNA_names = write_fasta(data["mRNA"].tolist(), args.output_dir / "mRNA.fa", "mRNA")
    data["siRNA_name"] = data["siRNA"].map(siRNA_names)
    data["mRNA_name"] = data["mRNA"].map(mRNA_names)
    data["y"] = (data["label"] >= data["label"].median()).astype(int)
    columns = ["row_id", "dataset", "family", "siRNA", "mRNA", "label", "y", "td", "wt_label", "delta_label", "group_id", "siRNA_name", "mRNA_name"]
    data.to_csv(args.output_dir / "mismatch_finetune.csv", columns=columns, index=False, encoding="utf-8-sig")
    print({"rows": len(data), "unique_siRNA": len(siRNA_names), "unique_mRNA": len(mRNA_names), "output": str(args.output_dir)})


if __name__ == "__main__":
    main()