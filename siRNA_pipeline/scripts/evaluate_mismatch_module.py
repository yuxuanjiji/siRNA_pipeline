"""Train the legacy mismatch CNN only where an explicit WT guide exists."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


class MismatchModule(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 16, kernel_size=3, padding=1)
        self.pool = nn.MaxPool1d(2, 2)
        self.fc1 = nn.Linear(16 * 9 * 2, 64)
        self.fc2 = nn.Linear(64, 1)

    def forward(self, base, variant):
        base = self.pool(torch.relu(self.conv1(base))).flatten(1)
        variant = self.pool(torch.relu(self.conv1(variant))).flatten(1)
        return torch.sigmoid(self.fc2(torch.relu(self.fc1(torch.cat((base, variant), dim=1))))).squeeze(1)


def encode(sequence: str) -> np.ndarray:
    mapping = {"A": 0, "U": 1, "C": 2, "G": 3}
    result = np.zeros((4, 19), dtype=np.float32)
    for index, base in enumerate(sequence):
        result[mapping[base], index] = 1.0
    return result


def attach_wt(data: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    references = {}
    for group_id, group in data.groupby("group_id"):
        wt = group[group["mismatch_count"] == 0]
        if len(wt) == 1:
            references[group_id] = wt.iloc[0]["siRNA"]
        elif len(wt) > 1 and wt["siRNA"].nunique() == 1:
            references[group_id] = wt.iloc[0]["siRNA"]
    data = data.copy()
    data["family"] = data["sirna_id"].astype(str).str.replace(r"-.*$", "", regex=True)
    family_references = {}
    for family, family_rows in data[data["mismatch_count"] == 0].groupby("family"):
        if family_rows["siRNA"].nunique() == 1:
            family_references[family] = family_rows.iloc[0]["siRNA"]
    result = data.copy()
    result["wt_siRNA"] = result["group_id"].map(references)
    result["wt_reference_origin"] = np.where(result["wt_siRNA"].notna(), "same_group", "")
    family_mask = result["wt_siRNA"].isna()
    result.loc[family_mask, "wt_siRNA"] = result.loc[family_mask, "family"].map(family_references)
    result.loc[family_mask & result["wt_siRNA"].notna(), "wt_reference_origin"] = "same_family"
    result = result[result["wt_siRNA"].notna()].copy()
    audit = {
        "total_groups": int(data["group_id"].nunique()),
        "groups_with_same_group_wt": len(references),
        "families_with_wt": len(family_references),
        "excluded_groups_without_wt_reference": sorted(set(data["group_id"]) - set(result["group_id"])),
        "rows_with_wt_reference": len(result),
        "rows_by_reference_origin": result["wt_reference_origin"].value_counts().to_dict(),
    }
    return result, audit


def tensors(frame: pd.DataFrame) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    base = torch.tensor(np.stack([encode(value) for value in frame["wt_siRNA"]]))
    variant = torch.tensor(np.stack([encode(value) for value in frame["siRNA"]]))
    label = torch.tensor(frame["label"].to_numpy(dtype=np.float32))
    return base, variant, label


def train_one(train: pd.DataFrame, test: pd.DataFrame, seed: int, epochs: int = 80) -> dict:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    model = MismatchModule()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-2)
    loss_fn = nn.SmoothL1Loss()
    base, variant, label = tensors(train)
    loader = DataLoader(TensorDataset(base, variant, label), batch_size=min(8, len(train)), shuffle=True)
    model.train()
    for _ in range(epochs):
        for batch_base, batch_variant, batch_label in loader:
            loss = loss_fn(model(batch_base, batch_variant), batch_label)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
    model.eval()
    test_base, test_variant, test_label = tensors(test)
    with torch.no_grad():
        prediction = model(test_base, test_variant).numpy()
    truth = test_label.numpy()
    spearman = pd.Series(truth).corr(pd.Series(prediction), method="spearman") if len(test) > 1 else None
    return {"mae": float(mean_absolute_error(truth, prediction)), "rmse": float(np.sqrt(mean_squared_error(truth, prediction))), "spearman": None if pd.isna(spearman) else float(spearman), "predictions": [{"row_id": row_id, "label": float(label), "prediction": float(pred)} for row_id, label, pred in zip(test["row_id"], truth, prediction)]}


def evaluate(input_path: Path, output_path: Path, seeds: tuple[int, ...]) -> dict:
    raw = pd.read_csv(input_path, dtype={"row_id": str})
    data, audit = attach_wt(raw[raw["model_eligible"]].copy())
    report = {"input": str(input_path), "audit": audit, "outer_loso": {}}
    for dataset in sorted(raw["dataset"].unique()):
        train = data[data["dataset"] != dataset]
        test = data[data["dataset"] == dataset]
        result = {"train_rows": len(train), "test_rows": len(test), "seeds": list(seeds), "metrics": []}
        if len(train) and len(test):
            result["metrics"] = [train_one(train, test, seed) for seed in seeds]
        if len(test) == 0:
            result["interpretation"] = "missing_wt_reference"
        elif len(train) < 5:
            result["interpretation"] = "insufficient_training_rows"
        elif len(test) < 5:
            result["interpretation"] = "exploratory_only_n4"
        else:
            result["interpretation"] = "external_loso"
        report["outer_loso"][dataset] = result
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_curated.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parents[1] / "outputs" / "analysis" / "mismatch" / "mismatch_module_loso.json")
    args = parser.parse_args()
    report = evaluate(args.input, args.output, (11, 23, 37))
    print(json.dumps({"audit": report["audit"], "outer_loso": {key: {k: value[k] for k in ("train_rows", "test_rows", "interpretation")} for key, value in report["outer_loso"].items()}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()