"""Train a frozen OligoFormer representation with a small delta-efficiency head."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import GroupKFold
from torch import nn


ROOT = Path(__file__).parents[1]
OLIGO_SCRIPTS = ROOT / "OligoFormer部分" / "scripts"
DATA_ROOT = ROOT / "OligoFormer部分" / "data" / "mismatch_finetune"
FM_ROOT = ROOT / "OligoFormer部分" / "data" / "RNAFM"
sys.path.insert(0, str(OLIGO_SCRIPTS))
from model import Oligo  # noqa: E402


DEVICE = torch.device("cpu")
SEEDS = (11, 23, 37)
EPOCH_OPTIONS = (20, 50, 100)


def one_hot(sequence: str) -> torch.Tensor:
    mapping = {"A": 0, "U": 1, "C": 2, "G": 3}
    result = np.zeros((1, 1, len(sequence), 5), dtype=np.float32)
    for index, base in enumerate(sequence):
        result[0, 0, index, mapping[base]] = 1.0
    return torch.tensor(result)


def read_fasta(path: Path) -> dict[str, str]:
    result = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
        elif name:
            result[name] = line.strip()
    return result


def load_model(weight: Path) -> Oligo:
    model = Oligo().to(DEVICE)
    model.load_state_dict(torch.load(weight, map_location=DEVICE, weights_only=False))
    model.eval()
    for parameter in model.parameters():
        parameter.requires_grad = False
    return model


def extract_merge(model: Oligo, row: pd.Series, siRNA_names: dict[str, str], mRNA_names: dict[str, str]) -> torch.Tensor:
    siRNA = one_hot(row["siRNA"])
    mRNA = one_hot(row["mRNA"])
    siRNA_fm = torch.tensor(np.load(FM_ROOT / "mismatch_finetune_siRNA" / "representations" / f"{siRNA_names[row['siRNA']]}.npy"), dtype=torch.float32).unsqueeze(0)
    mRNA_fm = torch.tensor(np.load(FM_ROOT / "mismatch_finetune_mRNA" / "representations" / f"{mRNA_names[row['mRNA']]}.npy"), dtype=torch.float32).unsqueeze(0)
    td = torch.tensor([[float(value) for value in str(row["td"]).split(",")]], dtype=torch.float32)
    with torch.no_grad():
        siRNA_encoded, _ = model.siRNA_encoder(siRNA)
        mRNA_encoded, _ = model.mRNA_encoder(mRNA)
        siRNA_fm = model.siRNA_avgpool(siRNA_fm).view(1, -1)
        mRNA_fm = model.mRNA_avgpool(mRNA_fm).view(1, -1)
        merge = torch.cat([model.flatten(siRNA_encoded), model.flatten(mRNA_encoded), model.flatten(siRNA_fm), model.flatten(mRNA_fm), td], dim=-1)
    return merge.squeeze(0)


def spearman(labels: np.ndarray, prediction: np.ndarray) -> float | None:
    if len(labels) < 2 or np.unique(labels).size < 2 or np.unique(prediction).size < 2:
        return None
    value = pd.Series(labels).corr(pd.Series(prediction), method="spearman")
    return None if pd.isna(value) else float(value)


def fit_head(features: torch.Tensor, labels: torch.Tensor, epochs: int, seed: int) -> nn.Linear:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    head = nn.Linear(features.shape[1], 1)
    optimizer = torch.optim.AdamW(head.parameters(), lr=1e-3, weight_decay=1e-2)
    loss_fn = nn.SmoothL1Loss()
    for _ in range(epochs):
        prediction = head(features).squeeze(1)
        loss = loss_fn(prediction, labels)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return head.eval()


def choose_epochs(features: torch.Tensor, labels: torch.Tensor, groups: np.ndarray) -> tuple[int, dict]:
    splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    scores = {}
    for epochs in EPOCH_OPTIONS:
        values = []
        for train_idx, valid_idx in splitter.split(features, labels, groups):
            head = fit_head(features[train_idx], labels[train_idx], epochs, 101)
            prediction = head(features[valid_idx]).detach().numpy().ravel()
            score = spearman(labels[valid_idx].numpy(), prediction)
            if score is not None:
                values.append(score)
        scores[str(epochs)] = values
    means = {key: (float(np.mean(value)) if value else float("-inf")) for key, value in scores.items()}
    return int(max(means, key=means.get)), {"fold_scores": scores, "mean_spearman": means}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight", type=Path, default=ROOT / "OligoFormer部分" / "model" / "best_model.pth")
    parser.add_argument("--input", type=Path, default=DATA_ROOT / "mismatch_finetune.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "analysis" / "mismatch" / "oligoformer_head_loso.json")
    args = parser.parse_args()
    data = pd.read_csv(args.input, dtype={"row_id": str})
    siRNA_names = {sequence: name for name, sequence in read_fasta(DATA_ROOT / "siRNA.fa").items()}
    mRNA_names = {sequence: name for name, sequence in read_fasta(DATA_ROOT / "mRNA.fa").items()}
    model = load_model(args.weight)
    features = torch.stack([extract_merge(model, row, siRNA_names, mRNA_names) for _, row in data.iterrows()])
    labels = torch.tensor(data["delta_label"].to_numpy(dtype=np.float32))
    report = {"rows": len(data), "feature_dim": features.shape[1], "backbone_frozen": True, "device": str(DEVICE), "outer_loso": {}}
    for dataset in sorted(data["dataset"].unique()):
        train_mask = data["dataset"].to_numpy() != dataset
        test_mask = ~train_mask
        train = data[train_mask].reset_index(drop=True)
        test = data[test_mask].reset_index(drop=True)
        train_features = features[torch.tensor(np.flatnonzero(train_mask))]
        test_features = features[torch.tensor(np.flatnonzero(test_mask))]
        train_labels = labels[torch.tensor(np.flatnonzero(train_mask))]
        test_labels = labels[torch.tensor(np.flatnonzero(test_mask))]
        epochs, selection = choose_epochs(train_features, train_labels, train["group_id"].to_numpy())
        seed_metrics = []
        for seed in SEEDS:
            head = fit_head(train_features, train_labels, epochs, seed)
            delta_prediction = head(test_features).detach().numpy().ravel()
            raw_prediction = test["wt_label"].to_numpy() + delta_prediction
            raw_label = test["label"].to_numpy()
            seed_metrics.append({"seed": seed, "delta_spearman": spearman(test_labels.numpy(), delta_prediction), "raw_spearman": spearman(raw_label, raw_prediction), "raw_mae": float(np.mean(np.abs(raw_label - raw_prediction)))})
        report["outer_loso"][dataset] = {"train_rows": len(train), "test_rows": len(test), "selected_epochs": epochs, "inner_selection": selection, "metrics": seed_metrics}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("rows", "feature_dim", "backbone_frozen", "device")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()