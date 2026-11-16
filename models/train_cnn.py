"""Train a minimal short-context CNN for acceptor usage."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error
from torch.utils.data import DataLoader, TensorDataset

from models.cnn import SpliceCNN, predict_usage, usage_loss
from src.extract_sequences import load_context, one_hot_encode


def target_values(sites: pd.DataFrame, choice: str) -> np.ndarray:
    if choice == "mean":
        return sites[["usage_rep1", "usage_rep2"]].mean(axis=1).to_numpy(dtype=np.float32)
    return sites[f"usage_{choice}"].to_numpy(dtype=np.float32)


def metrics(model: SpliceCNN, x: torch.Tensor, y: torch.Tensor, device: torch.device) -> dict[str, float]:
    model.eval()
    with torch.no_grad():
        prediction = predict_usage(model, x.to(device)).cpu().numpy()
    actual = y.numpy()
    return {
        "mae": float(mean_absolute_error(actual, prediction)),
        "rmse": float(mean_squared_error(actual, prediction) ** 0.5),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", default="data/processed/person_b_acceptors_all.splits.tsv.gz")
    ap.add_argument("--sequences", default="data/processed/person_b_acceptors_all.sequences.npy")
    ap.add_argument("--context", type=int, choices=[80, 400], required=True)
    ap.add_argument("--target", choices=["rep1", "rep2", "mean"], default="mean")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--learning-rate", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out", default=None, help="model checkpoint path")
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str, "strand": str})
    codes = load_context(args.sequences, args.context)
    if len(codes) != len(sites):
        raise ValueError(f"sequence rows ({len(codes)}) do not match site rows ({len(sites)})")
    y = target_values(sites, args.target)
    valid = np.isfinite(y)
    if not valid.all():
        sites, codes, y = sites.loc[valid].reset_index(drop=True), codes[valid], y[valid]

    x = torch.from_numpy(one_hot_encode(codes, dtype=np.float32).transpose(0, 2, 1).copy())
    targets = torch.from_numpy(y)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = SpliceCNN().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)
    train_mask = sites["split"].to_numpy() == "train"
    validation_mask = sites["split"].to_numpy() == "validation"
    test_mask = sites["split"].to_numpy() == "test"
    if not all(mask.any() for mask in (train_mask, validation_mask, test_mask)):
        raise ValueError("train, validation, and test must each contain at least one site")

    train_data = TensorDataset(x[train_mask], targets[train_mask])
    loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True)
    best_validation = float("inf")
    best_state = None
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        for batch_x, batch_y in loader:
            optimizer.zero_grad()
            loss = usage_loss(model(batch_x.to(device)), batch_y.to(device))
            loss.backward()
            optimizer.step()
            losses.append(loss.item())
        validation = metrics(model, x[validation_mask], targets[validation_mask], device)
        print(f"epoch {epoch:03d} train_loss={np.mean(losses):.5f} validation_mae={validation['mae']:.5f}", flush=True)
        if validation["mae"] < best_validation:
            best_validation = validation["mae"]
            best_state = {key: value.cpu().clone() for key, value in model.state_dict().items()}

    if best_state is None:
        raise RuntimeError("training did not produce a checkpoint")
    model.load_state_dict(best_state)
    output = Path(args.out or f"results/cnn_{args.context}bp_{args.target}.pt")
    output.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"model_state": model.state_dict(), "context": args.context, "target": args.target}, output)
    summary = {
        "context": args.context,
        "target": args.target,
        "device": str(device),
        "n_train": int(train_mask.sum()),
        "n_validation": int(validation_mask.sum()),
        "n_test": int(test_mask.sum()),
        "validation": metrics(model, x[validation_mask], targets[validation_mask], device),
        "test": metrics(model, x[test_mask], targets[test_mask], device),
        "checkpoint": str(output),
    }
    summary_path = output.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
