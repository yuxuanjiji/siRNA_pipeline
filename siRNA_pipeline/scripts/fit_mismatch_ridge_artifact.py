"""Fit and export the selected mismatch Ridge research artifact."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)
FEATURE_NAMES = ["mismatch_count", "parsed_mismatch_count", "band_1_2", "band_3_8", "band_9_12", "band_13_19"] + [f"position_{position}" for position in range(1, 20)]


def feature_row(row: pd.Series) -> list[float]:
    positions = ast.literal_eval(row["mismatch_positions"])
    values = [float(row["mismatch_count"]), float(len(positions))]
    values.extend(float(sum(start <= position <= end for position in positions)) for start, end in ((1, 2), (3, 8), (9, 12), (13, 19)))
    values.extend(float(position in positions) for position in range(1, 20))
    return values


def select_alpha(data: pd.DataFrame) -> tuple[float, dict]:
    x = np.asarray([feature_row(row) for _, row in data.iterrows()])
    y = data["label"].to_numpy()
    groups = data["group_id"].to_numpy()
    splitter = GroupKFold(n_splits=min(5, data["group_id"].nunique()))
    scores = {}
    for alpha in ALPHAS:
        fold_scores = []
        for train_idx, valid_idx in splitter.split(x, y, groups):
            model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
            model.fit(x[train_idx], y[train_idx])
            prediction = model.predict(x[valid_idx])
            value = pd.Series(y[valid_idx]).corr(pd.Series(prediction), method="spearman")
            if not pd.isna(value):
                fold_scores.append(float(value))
        scores[str(alpha)] = fold_scores
    means = {alpha: float(np.mean(values)) if values else float("-inf") for alpha, values in scores.items()}
    selected = max(means, key=means.get)
    return float(selected), {"fold_scores": scores, "mean_spearman": means}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_ridge_artifact.json")
    args = parser.parse_args()
    data = pd.read_csv(args.input, dtype={"row_id": str})
    data = data[data["model_eligible"]].reset_index(drop=True)
    alpha, selection = select_alpha(data)
    x = np.asarray([feature_row(row) for _, row in data.iterrows()])
    y = data["label"].to_numpy()
    model = make_pipeline(StandardScaler(), Ridge(alpha=alpha)).fit(x, y)
    scaler = model.named_steps["standardscaler"]
    ridge = model.named_steps["ridge"]
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    artifact = {
        "model": "Ridge",
        "feature_set": "position_only",
        "feature_names": FEATURE_NAMES,
        "alpha": alpha,
        "scaler_mean": scaler.mean_.tolist(),
        "scaler_scale": scaler.scale_.tolist(),
        "coefficients": ridge.coef_.tolist(),
        "intercept": float(ridge.intercept_),
        "training_rows": len(data),
        "training_groups": int(data["group_id"].nunique()),
        "source_hash": digest,
        "selection": selection,
        "status": "research_only_not_connected_to_production",
    }
    args.output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: artifact[key] for key in ("model", "alpha", "training_rows", "training_groups", "source_hash", "status")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()