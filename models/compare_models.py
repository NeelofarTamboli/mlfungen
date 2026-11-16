"""Compare trained CNNs with the local baselines and the pretrained SpliceAI / Pangolin scores.

Test sets:
    * ``test_all``: every test chromosome (chr1, chr3, chr5, chr7, chr9); CNNs and local baselines.
    * ``test_chr1_chr3``: the chr1 + chr3 test sites (23,474), the rows on which SpliceAI and Pangolin
      were scored (``models.pretrained_baselines``). Those per-site scores are not stored here, so
      their metrics are read from ``results/pretrained_metrics_test_chr1_chr3.tsv``. The site IDs and
      usage values of these rows are identical in the chr1-4 and genome-wide tables, so the metrics are
      computed on the same sites and the same target.

The CNNs and local baselines are all trained on the genome-wide training chromosomes (validation chr2,
chr4). Metrics come from ``analysis.evaluate.evaluate_predictions`` (all sites and usage 0.1-0.9).

    python -m models.compare_models --checkpoints results/cnn_80bp_mean.pt
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from analysis.evaluate import evaluate_predictions
from models.cnn import SpliceCNN, predict_usage
from src.config import REPO_ROOT
from src.extract_sequences import load_context, one_hot_encode

PRETRAINED_CHROMS = ["chr1", "chr3"]


def cnn_predictions(checkpoint: Path, sequences: Path, batch_size: int = 1024) -> tuple[str, np.ndarray]:
    """Predicted usage for every row of the sequence array, using a train_cnn checkpoint."""
    saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
    context = int(saved["context"])
    model = SpliceCNN()
    model.load_state_dict(saved["model_state"])
    model.eval()
    codes = load_context(sequences, context)
    out = np.empty(len(codes), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(codes), batch_size):
            batch = one_hot_encode(np.asarray(codes[start:start + batch_size]), dtype=np.float32).transpose(0, 2, 1)
            out[start:start + batch_size] = predict_usage(model, torch.from_numpy(batch.copy())).numpy()
    return f"prediction_cnn_{context}bp", out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoints", nargs="+", required=True, help="train_cnn .pt files")
    ap.add_argument("--baseline-predictions", default="data/processed/baseline_predictions_all.tsv.gz")
    ap.add_argument("--sequences", default="data/processed/person_b_acceptors_all.sequences.npy")
    ap.add_argument("--pretrained-metrics", default="results/pretrained_metrics_test_chr1_chr3.tsv")
    ap.add_argument("--predictions-out", default="data/processed/model_predictions_all.tsv.gz")
    ap.add_argument("--out", default="results/model_comparison.tsv")
    a = ap.parse_args()
    root = lambda p: Path(p) if Path(p).is_absolute() else REPO_ROOT / p

    table = pd.read_csv(root(a.baseline_predictions), sep="\t", dtype={"chrom": str, "strand": str})
    if not (table["sequence_row"].to_numpy() == np.arange(len(table))).all():
        raise ValueError("baseline prediction rows are not in sequence-array order")
    for checkpoint in a.checkpoints:
        column, values = cnn_predictions(root(checkpoint), root(a.sequences))
        if len(values) != len(table):
            raise ValueError(f"{checkpoint}: {len(values)} predictions for {len(table)} sites")
        table[column] = values
    columns = [c for c in table if c.startswith("prediction_")]

    frames = []
    for name, rows in (("test_all", table), ("test_chr1_chr3", table[table["chrom"].isin(PRETRAINED_CHROMS)])):
        m = evaluate_predictions(rows, columns, target_column="target", split="test")
        m.insert(0, "test_set", name)
        m["source"] = "this run (trained on genome-wide training chromosomes)"
        frames.append(m)

    pretrained = root(a.pretrained_metrics)
    if pretrained.exists():
        p = pd.read_csv(pretrained, sep="\t")
        p = p[~p["model"].isin(["Training mean", "Ridge regression", "Gradient boosting"])]  # chr2-trained versions
        expected = frames[1].groupby("subset")["n_sites"].first()
        mismatch = p[p["n_sites"] != p["subset"].map(expected)]
        if len(mismatch):
            raise ValueError(f"pretrained metrics are on different sites:\n{mismatch}")
        frames.append(pd.DataFrame({
            "test_set": "test_chr1_chr3", "split": "test", "site_type": "acceptor", "subset": p["subset"],
            "model": p["model"], "n_sites": p["n_sites"], "pearson": p["pearson"], "spearman": p["spearman"],
            "r2": p["r2"], "mae": p["mae"], "source": [f"{pretrained.name}: {t}" for t in p["score_type"]],
        }))

    out = pd.concat(frames, ignore_index=True)
    root(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(root(a.out), sep="\t", index=False, float_format="%.4f")
    keep = ["site_id", "chrom", "split", "target", *columns]
    table[keep].to_csv(root(a.predictions_out), sep="\t", index=False, float_format="%.6g")
    print(out.drop(columns=["split", "site_type", "source"]).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"metrics -> {a.out}\npredictions -> {a.predictions_out}")


if __name__ == "__main__":
    main()
