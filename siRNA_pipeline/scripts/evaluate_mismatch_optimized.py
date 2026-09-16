"""Evaluate WT-paired mismatch-effect Ridge models without changing production."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)
BASES = "AUCG"
TYPE_NAMES = [f"{left}:{right}" for left in BASES for right in BASES]


def parse_list(value: str) -> list:
    return ast.literal_eval(value)


def feature_row(row: pd.Series) -> list[float]:
    positions = [int(position) for position in parse_list(row["mismatch_positions"])]
    types = parse_list(row["mismatch_types"])
    values = [float(len(positions))]
    values.extend(float(sum(start <= position <= end for position in positions)) for start, end in ((1, 2), (3, 8), (9, 12), (13, 19)))
    values.extend(float(position in positions) for position in range(1, 20))
    values.extend(float(types.count(type_name)) for type_name in TYPE_NAMES)
    return values


def spearman(labels: np.ndarray, predictions: np.ndarray) -> float | None:
    if len(labels) < 2 or np.unique(labels).size < 2 or np.unique(predictions).size < 2:
        return None
    value = pd.Series(labels).corr(pd.Series(predictions), method="spearman")
    return None if pd.isna(value) else float(value)


def add_references(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    group_wt = data[data["mismatch_count"] == 0].groupby("group_id")["label"].first().to_dict()
    family_wt = data[data["mismatch_count"] == 0].groupby("family")["label"].first().to_dict()
    data["wt_label"] = data["group_id"].map(group_wt)
    data["wt_label"] = data["wt_label"].fillna(data["family"].map(family_wt))
    if data["wt_label"].isna().any():
        raise ValueError("Some eligible records have no auditable WT reference")
    data["delta_label"] = data["label"] - data["wt_label"]
    return data


def fit_model(train: pd.DataFrame, alpha: float):
    x = np.asarray([feature_row(row) for _, row in train.iterrows()])
    model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
    model.fit(x, train["delta_label"].to_numpy())
    return model


def choose_alpha(train: pd.DataFrame) -> tuple[float, dict]:
    x = np.asarray([feature_row(row) for _, row in train.iterrows()])
    y = train["delta_label"].to_numpy()
    groups = train["group_id"].to_numpy()
    splitter = GroupKFold(n_splits=min(5, train["group_id"].nunique()))
    scores = {}
    for alpha in ALPHAS:
        fold_scores = []
        for train_idx, valid_idx in splitter.split(x, y, groups):
            model = fit_model(train.iloc[train_idx], alpha)
            prediction = model.predict(x[valid_idx])
            score = spearman(y[valid_idx], prediction)
            if score is not None:
                fold_scores.append(score)
        scores[str(alpha)] = fold_scores
    means = {key: float(np.mean(value)) if value else float("-inf") for key, value in scores.items()}
    selected = max(means, key=means.get)
    return float(selected), {"fold_scores": scores, "mean_spearman": means}


def evaluate(input_path: Path, output_path: Path) -> dict:
    data = pd.read_csv(input_path, dtype={"row_id": str})
    data = add_references(data[data["model_eligible"]].reset_index(drop=True))
    report = {"input": str(input_path), "rows": len(data), "model": "wt_paired_delta_ridge", "feature_count": len(feature_row(data.iloc[0])), "outer_loso": {}}
    for dataset in sorted(data["dataset"].unique()):
        train = data[data["dataset"] != dataset].reset_index(drop=True)
        test = data[data["dataset"] == dataset].reset_index(drop=True)
        alpha, selection = choose_alpha(train)
        model = fit_model(train, alpha)
        delta_prediction = model.predict(np.asarray([feature_row(row) for _, row in test.iterrows()]))
        raw_prediction = test["wt_label"].to_numpy() + delta_prediction
        raw_label = test["label"].to_numpy()
        delta_label = test["delta_label"].to_numpy()
        report["outer_loso"][dataset] = {
            "train_rows": len(train),
            "test_rows": len(test),
            "selected_alpha": alpha,
            "inner_selection": selection,
            "delta_spearman": spearman(delta_label, delta_prediction),
            "raw_spearman": spearman(raw_label, raw_prediction),
            "raw_mae": float(mean_absolute_error(raw_label, raw_prediction)),
            "raw_rmse": float(np.sqrt(mean_squared_error(raw_label, raw_prediction))),
            "predictions": [{"row_id": row_id, "label": float(label), "wt_label": float(wt), "delta_label": float(delta), "delta_prediction": float(pred), "raw_prediction": float(raw)} for row_id, label, wt, delta, pred, raw in zip(test["row_id"], raw_label, test["wt_label"], delta_label, delta_prediction, raw_prediction)],
            "interpretation": "exploratory_small_n" if len(test) < 10 else "external_loso",
        }
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "optimized_delta_ridge_loso.json")
    args = parser.parse_args()
    report = evaluate(args.input, args.output)
    for dataset, metrics in report["outer_loso"].items():
        print(dataset, {key: metrics[key] for key in ("test_rows", "selected_alpha", "delta_spearman", "raw_spearman", "raw_mae", "raw_rmse", "interpretation")})


if __name__ == "__main__":
    main()