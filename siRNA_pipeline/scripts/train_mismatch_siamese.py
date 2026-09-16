"""Train/evaluate the research-only WT/mismatch Siamese CNN."""

from __future__ import annotations

import argparse
import ast
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import GroupKFold
from torch import nn

from mismatch_siamese import MismatchSiameseCNN


SEEDS = (11, 23, 37)
EPOCH_OPTIONS = (50, 100)
BASES = "AUCG"
TYPE_NAMES = [f"{left}:{right}" for left in BASES for right in BASES]
THERMO_NAMES = ("dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content")


def encode(sequence: str) -> np.ndarray:
    mapping = {"A": 0, "U": 1, "C": 2, "G": 3}
    output = np.zeros((4, 19), dtype=np.float32)
    for index, base in enumerate(sequence):
        output[mapping[base], index] = 1.0
    return output


def add_references(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    group_wt = data[data["mismatch_count"] == 0].groupby("group_id")["label"].first().to_dict()
    family_wt = data[data["mismatch_count"] == 0].groupby("family")["label"].first().to_dict()
    data["wt_sequence"] = data["group_id"].map(
        data[data["mismatch_count"] == 0].groupby("group_id")["siRNA"].first().to_dict()
    )
    data["wt_sequence"] = data["wt_sequence"].fillna(data["family"].map(data[data["mismatch_count"] == 0].groupby("family")["siRNA"].first().to_dict()))
    data["wt_label"] = data["group_id"].map(group_wt).fillna(data["family"].map(family_wt))
    if data[["wt_sequence", "wt_label"]].isna().any().any():
        raise ValueError("Missing auditable WT reference")
    data["delta_label"] = data["label"] - data["wt_label"]
    return data


def auxiliary_row(row: pd.Series) -> list[float]:
    positions = [int(value) for value in ast.literal_eval(row["mismatch_positions"])]
    types = ast.literal_eval(row["mismatch_types"])
    values = [float(len(positions))]
    values.extend(float(sum(start <= position <= end for position in positions)) for start, end in ((1, 2), (3, 8), (9, 12), (13, 19)))
    values.extend(float(position in positions) for position in range(1, 20))
    values.extend(float(types.count(type_name)) for type_name in TYPE_NAMES)
    values.extend(float(row[name]) for name in THERMO_NAMES)
    return values


def tensors(data: pd.DataFrame) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    wt = torch.tensor(np.stack([encode(value) for value in data["wt_sequence"]]))
    mutant = torch.tensor(np.stack([encode(value) for value in data["siRNA"]]))
    auxiliary = torch.tensor(np.asarray([auxiliary_row(row) for _, row in data.iterrows()], dtype=np.float32))
    labels = torch.tensor(data["delta_label"].to_numpy(dtype=np.float32))
    return wt, mutant, auxiliary, labels


def ranking_loss(prediction: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    difference = prediction[:, None] - prediction[None, :]
    target = torch.sign(labels[:, None] - labels[None, :])
    mask = target != 0
    if not torch.any(mask):
        return prediction.new_tensor(0.0)
    return torch.relu(0.05 - target[mask] * difference[mask]).mean()


def fit(train: pd.DataFrame, epochs: int, seed: int) -> MismatchSiameseCNN:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    wt, mutant, auxiliary, labels = tensors(train)
    model = MismatchSiameseCNN(auxiliary.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-2)
    huber = nn.SmoothL1Loss()
    model.train()
    for _ in range(epochs):
        prediction = model(wt, mutant, auxiliary)
        loss = huber(prediction, labels) + 0.25 * ranking_loss(prediction, labels)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return model.eval()


def spearman(labels: np.ndarray, prediction: np.ndarray) -> float | None:
    if len(labels) < 2 or np.unique(labels).size < 2 or np.unique(prediction).size < 2:
        return None
    value = pd.Series(labels).corr(pd.Series(prediction), method="spearman")
    return None if pd.isna(value) else float(value)


def choose_epochs(train: pd.DataFrame) -> tuple[int, dict]:
    wt, mutant, auxiliary, labels = tensors(train)
    scores = {}
    splitter = GroupKFold(n_splits=min(5, train["group_id"].nunique()))
    for epochs in EPOCH_OPTIONS:
        fold_scores = []
        for train_idx, valid_idx in splitter.split(wt, labels, train["group_id"]):
            model = fit(train.iloc[train_idx], epochs, 101)
            with torch.no_grad():
                prediction = model(wt[valid_idx], mutant[valid_idx], auxiliary[valid_idx]).numpy()
            score = spearman(labels[valid_idx].numpy(), prediction)
            if score is not None:
                fold_scores.append(score)
        scores[str(epochs)] = fold_scores
    means = {key: float(np.mean(value)) if value else float("-inf") for key, value in scores.items()}
    return int(max(means, key=means.get)), {"fold_scores": scores, "mean_spearman": means}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated_td.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_siamese_loso.json")
    args = parser.parse_args()
    data = add_references(pd.read_csv(args.input, dtype={"row_id": str}))
    report = {"model": "WT-mismatch Siamese CNN", "rows": len(data), "auxiliary_dim": len(auxiliary_row(data.iloc[0])), "loss": "Huber + 0.25 pairwise hinge", "outer_loso": {}}
    for dataset in sorted(data["dataset"].unique()):
        train = data[data["dataset"] != dataset].reset_index(drop=True)
        test = data[data["dataset"] == dataset].reset_index(drop=True)
        epochs, selection = choose_epochs(train)
        wt, mutant, auxiliary, labels = tensors(test)
        metrics = []
        for seed in SEEDS:
            model = fit(train, epochs, seed)
            with torch.no_grad():
                delta_prediction = model(wt, mutant, auxiliary).numpy()
            raw_prediction = test["wt_label"].to_numpy() + delta_prediction
            raw_label = test["label"].to_numpy()
            metrics.append({"seed": seed, "delta_spearman": spearman(labels.numpy(), delta_prediction), "raw_spearman": spearman(raw_label, raw_prediction), "raw_mae": float(np.mean(np.abs(raw_label - raw_prediction))), "raw_rmse": float(np.sqrt(np.mean((raw_label - raw_prediction) ** 2)))})
        report["outer_loso"][dataset] = {"train_rows": len(train), "test_rows": len(test), "selected_epochs": epochs, "inner_selection": selection, "metrics": metrics}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    for dataset, metrics in report["outer_loso"].items():
        print(dataset, metrics["selected_epochs"], metrics["metrics"])


if __name__ == "__main__":
    main()