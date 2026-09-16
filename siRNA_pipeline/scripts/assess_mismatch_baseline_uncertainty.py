"""Bootstrap and permutation checks for the saved LOSO baseline predictions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def correlation(labels: np.ndarray, predictions: np.ndarray) -> float | None:
    if len(labels) < 3 or np.unique(labels).size < 2 or np.unique(predictions).size < 2:
        return None
    value = pd.Series(labels).corr(pd.Series(predictions), method="spearman")
    return None if pd.isna(value) else float(value)


def summarize(labels: np.ndarray, predictions: np.ndarray, rng: np.random.Generator, bootstrap: int, permutations: int) -> dict:
    observed = correlation(labels, predictions)
    if observed is None:
        return {"observed_spearman": None, "bootstrap_ci95": None, "permutation_p": None}
    boot = []
    for _ in range(bootstrap):
        indices = rng.integers(0, len(labels), len(labels))
        value = correlation(labels[indices], predictions[indices])
        if value is not None:
            boot.append(value)
    null = []
    for _ in range(permutations):
        value = correlation(rng.permutation(labels), predictions)
        if value is not None:
            null.append(value)
    p_value = (1 + sum(abs(value) >= abs(observed) for value in null)) / (len(null) + 1)
    return {
        "observed_spearman": observed,
        "bootstrap_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if boot else None,
        "permutation_p_two_sided": float(p_value),
        "bootstrap_replicates": len(boot),
        "permutation_replicates": len(null),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "baseline_loso.json")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "baseline_uncertainty.json")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--permutations", type=int, default=5000)
    args = parser.parse_args()
    report = json.loads(args.input.read_text(encoding="utf-8"))
    rng = np.random.default_rng(20260911)
    result = {"input": str(args.input), "bootstrap": args.bootstrap, "permutations": args.permutations, "models": {}}
    for model_name, model in report["models"].items():
        result["models"][model_name] = {}
        for dataset, metrics in model["outer_loso"].items():
            rows = metrics["predictions"]
            labels = np.asarray([row["label"] for row in rows], dtype=float)
            predictions = np.asarray([row["prediction"] for row in rows], dtype=float)
            result["models"][model_name][dataset] = summarize(labels, predictions, rng, args.bootstrap, args.permutations)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    for model_name, datasets in result["models"].items():
        for dataset, summary in datasets.items():
            print(model_name, dataset, summary)


if __name__ == "__main__":
    main()