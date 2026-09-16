"""Join the verified Birmingham rows with their reconstructed 57-nt target windows."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matched", type=Path, required=True)
    parser.add_argument("--filled-c", type=Path, required=True)
    parser.add_argument("--filled-d", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    matched = pd.read_csv(args.matched, dtype=str)
    filled = pd.concat([pd.read_csv(args.filled_c, dtype=str), pd.read_csv(args.filled_d, dtype=str)], ignore_index=True)
    filled = filled.drop_duplicates(subset=["row_list_idx"], keep="first")
    joined = matched.merge(filled[["row_list_idx", "mRNA_57nt"]], on="row_list_idx", how="left", validate="one_to_one")
    joined["mRNA"] = joined["mRNA_57nt"]
    joined["沉默效率"] = joined[["log2_rep1", "log2_rep2"]].apply(lambda row: ";".join(value for value in row if pd.notna(value) and str(value).strip()), axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output, index=False, encoding="utf-8-sig")
    print({"matched_rows": len(matched), "joined_rows": len(joined), "mRNA_rows": int(joined["mRNA"].notna().sum()), "output": str(args.output)})


if __name__ == "__main__":
    main()