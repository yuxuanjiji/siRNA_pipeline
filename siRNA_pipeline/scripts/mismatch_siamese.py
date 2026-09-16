"""Small WT/mismatch Siamese CNN for research-only efficiency deltas."""

from __future__ import annotations

import torch
from torch import nn


class MismatchSiameseCNN(nn.Module):
    def __init__(self, auxiliary_dim: int):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Conv1d(4, 16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(16, 32, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.head = nn.Sequential(
            nn.Linear(32 * 4 + auxiliary_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, 1),
        )

    def encode(self, sequence: torch.Tensor) -> torch.Tensor:
        return self.encoder(sequence).squeeze(-1)

    def forward(self, wt: torch.Tensor, mutant: torch.Tensor, auxiliary: torch.Tensor) -> torch.Tensor:
        wt_features = self.encode(wt)
        mutant_features = self.encode(mutant)
        difference = mutant_features - wt_features
        absolute_difference = torch.abs(difference)
        merged = torch.cat((wt_features, mutant_features, difference, absolute_difference, auxiliary), dim=1)
        return self.head(merged).squeeze(1)