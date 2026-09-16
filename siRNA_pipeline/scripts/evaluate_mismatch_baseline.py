"""Evaluate a small, non-deep baseline on the curated mismatch collection."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)


def feature_row(row: pd.Series, include_thermo: bool = False) -> list[float]:
    positions = ast.literal_eval(row["mismatch_positions"])
    features = [float(row["mismatch_count"]), float(len(positions))]
    for start, end in ((1, 2), (3, 8), (9, 12), (13, 19)):
        features.append(float(sum(start <= position <= end for position in positions)))
    features.extend(float(position in positions) for position in range(1, 20))
    if include_thermo:
        features.extend(float(row[column]) for column in ("dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content"))
    return features


def spearman(y_true: np.ndarray, prediction: np.ndarray) -> float | None:
    if len(y_true) < 2 or pd.Series(y_true).nunique() < 2 or pd.Series(prediction).nunique() < 2:
        return None
    return float(pd.Series(y_true).corr(pd.Series(prediction), method="spearman"))


def fit_model(train: pd.DataFrame, alpha: float, model_name: str, include_thermo: bool):
    estimator = Ridge(alpha=alpha) if model_name == "ridge" else ElasticNet(alpha=alpha, l1_ratio=0.2, max_iter=10000)
    model = make_pipeline(StandardScaler(), estimator)
    model.fit(np.asarray([feature_row(row, include_thermo) for _, row in train.iterrows()]), train["label"].to_numpy())
    return model


def select_alpha(train: pd.DataFrame, model_name: str, include_thermo: bool) -> tuple[float, dict]:
    scores = {}
    n_splits = min(5, train["group_id"].nunique())
    if n_splits < 2:
        return 1.0, scores
    splitter = GroupKFold(n_splits=n_splits)
    x = np.asarray([feature_row(row, include_thermo) for _, row in train.iterrows()])
    y = train["label"].to_numpy()
    groups = train["group_id"].to_numpy()
    for alpha in ALPHAS:
        fold_scores = []
        for train_idx, valid_idx in splitter.split(x, y, groups):
            model = fit_model(train.iloc[train_idx], alpha, model_name, include_thermo)
            pred = model.predict(x[valid_idx])
            score = spearman(y[valid_idx], pred)
            if score is not None:
                fold_scores.append(score)
        scores[str(alpha)] = fold_scores
    means = {alpha: (float(np.mean(values)) if values else float("-inf")) for alpha, values in scores.items()}
    selected = max(means, key=means.get)
    return float(selected), {"fold_scores": scores, "mean_spearman": means}


def evaluate(input_path: Path, output_path: Path) -> dict:
    data = pd.read_csv(input_path, dtype={"row_id": str})
    data = data[data["model_eligible"]].reset_index(drop=True)
    report = {"input": str(input_path), "rows": len(data), "models": {}}
    feature_sets = {"position_only": False}
    if {"dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content"}.issubset(data.columns):
        feature_sets["position_plus_thermo"] = True
    for feature_name, include_thermo in feature_sets.items():
      for model_name in ("ridge", "elasticnet"):
        model_report = {"feature_set": feature_name, "outer_loso": {}}
        for dataset in sorted(data["dataset"].unique()):
            train = data[data["dataset"] != dataset].reset_index(drop=True)
            test = data[data["dataset"] == dataset].reset_index(drop=True)
            alpha, selection = select_alpha(train, model_name, include_thermo)
            model = fit_model(train, alpha, model_name, include_thermo)
            prediction = model.predict(np.asarray([feature_row(row, include_thermo) for _, row in test.iterrows()]))
            y_true = test["label"].to_numpy()
            model_report["outer_loso"][dataset] = {
                "train_rows": len(train),
                "test_rows": len(test),
                "test_groups": int(test["group_id"].nunique()),
                "selected_alpha": alpha,
                "inner_selection": selection,
                "mae": float(mean_absolute_error(y_true, prediction)),
                "rmse": float(np.sqrt(mean_squared_error(y_true, prediction))),
                "spearman": spearman(y_true, prediction),
                "predictions": [{"row_id": row_id, "label": float(label), "prediction": float(pred)} for row_id, label, pred in zip(test["row_id"], y_true, prediction)],
                "interpretation": "exploratory_only_n4" if len(test) < 5 else "external_loso",
            }
        report["models"][f"{feature_name}_{model_name}"] = model_report
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "baseline_loso.json")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = evaluate(args.input, args.output)
    for model_name, result in report["models"].items():
        print(model_name)
        for dataset, metrics in result["outer_loso"].items():
            print(dataset, {key: metrics[key] for key in ("test_rows", "selected_alpha", "mae", "rmse", "spearman", "interpretation")})


if __name__ == "__main__":
    main()