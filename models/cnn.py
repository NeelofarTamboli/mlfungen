"""Minimal 1D CNN that predicts acceptor usage from one-hot sequence.

Plain CNN for the short windows (80 / 400 bp), as in implementation_plan.md:
    * input: one-hot sequence, shape (batch, 4, L); N bases are all zeros
      (see ``src.extract_sequences.one_hot_encode``, then transpose to 4 x L);
    * a few "same"-padded conv layers, so every position keeps a feature vector;
    * one output per site, read at the centre position (the junction; the
      canonical AG sits at ``[L//2-2 : L//2]``);
    * loss: binary cross-entropy on the 0-1 usage value (SpliceAI / Pangolin),
      optionally weighted per site (e.g. by read support).

SpliceAI-style dilated residual blocks for 2 kb / 10 kb are not here yet.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class SpliceCNN(nn.Module):
    """Stack of Conv1d -> BatchNorm -> ReLU layers with a centre-position readout."""

    def __init__(self, channels: int = 32, kernel_size: int = 11, n_layers: int = 3, dropout: float = 0.1):
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd so 'same' padding stays centred")
        layers, in_ch = [], 4
        for _ in range(n_layers):
            layers += [
                nn.Conv1d(in_ch, channels, kernel_size, padding=kernel_size // 2),
                nn.BatchNorm1d(channels),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            in_ch = channels
        self.body = nn.Sequential(*layers)
        self.head = nn.Conv1d(channels, 1, kernel_size=1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Map (batch, 4, L) one-hot input to (batch,) logits at the centre position."""
        logits = self.head(self.body(x)).squeeze(1)  # (batch, L)
        return logits[:, x.shape[-1] // 2]


def usage_loss(logits: torch.Tensor, usage: torch.Tensor, weight: torch.Tensor | None = None) -> torch.Tensor:
    """Binary cross-entropy between predicted logits and 0-1 usage targets."""
    return F.binary_cross_entropy_with_logits(logits, usage.float(), weight=weight)


def predict_usage(model: nn.Module, x: torch.Tensor) -> torch.Tensor:
    """Predicted usage in [0, 1]."""
    return torch.sigmoid(model(x))
