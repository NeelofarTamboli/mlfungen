"""Evaluate usage predictions on all and alternative acceptor sites.

Metrics are Pearson correlation, Spearman correlation, R-squared and MAE.
The alternative subset is defined from the observed target as 0.1 <= usage <=
0.9. Final reporting should use only the held-out test chromosomes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import mean_absolute_error, r2_score

from src.config import load_config, project_path


def regression_metrics(y_true: Iterable[float], y_pred: Iterable[float]) -> dict[str, float | int]:
    """Compute safe regression metrics after dropping non-finite pairs."""
    y_true, y_pred = np.asarray(y_true, dtype=float), np.asarray(y_pred, dtype=float)
    keep = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true, y_pred = y_true[keep], y_pred[keep]
    out: dict[str, float | int] = {"n_sites": int(len(y_true))}
    if not len(y_true):
        return {**out, "pearson": np.nan, "spearman": np.nan, "r2": np.nan, "mae": np.nan}
    variable_true = len(y_true) >= 2 and np.ptp(y_true) > 0
    variable_pred = len(y_pred) >= 2 and np.ptp(y_pred) > 0
    out["pearson"] = float(stats.pearsonr(y_true, y_pred).statistic) if variable_true and variable_pred else np.nan
    out["spearman"] = float(stats.spearmanr(y_true, y_pred).statistic) if variable_true and variable_pred else np.nan
    out["r2"] = float(r2_score(y_true, y_pred)) if variable_true and len(y_true) >= 2 else np.nan
    out["mae"] = float(mean_absolute_error(y_true, y_pred))
    return out


def evaluate_predictions(
    predictions: pd.DataFrame,
    prediction_columns: Iterable[str],
    target_column: str = "target",
    split: str = "test",
) -> pd.DataFrame:
    """Return tidy metrics for every model on all and alternative sites."""
    data = predictions[predictions["split"] == split] if "split" in predictions else predictions
    rows = []
    groups = [("all", data), ("alternative_0.1_to_0.9", data[data[target_column].between(0.1, 0.9, inclusive="both")])]
    for model_col in prediction_columns:
        for subset, group in groups:
            rows.append({
                "split": split,
                "site_type": "acceptor",
                "subset": subset,
                "model": model_col.removeprefix("prediction_"),
                **regression_metrics(group[target_column], group[model_col]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions", default=str(project_path(cfg, "data_processed") / "baseline_predictions_chr1-4.tsv.gz"))
    ap.add_argument("--out", default=str(project_path(cfg, "results") / "baseline_metrics_chr1-4.tsv"))
    ap.add_argument("--target-column", default="target")
    ap.add_argument("--split", default="test")
    args = ap.parse_args()
    predictions = pd.read_csv(args.predictions, sep="\t")
    prediction_columns = [column for column in predictions if column.startswith("prediction_")]
    if not prediction_columns:
        raise ValueError("prediction table has no prediction_* columns")
    metrics = evaluate_predictions(predictions, prediction_columns, args.target_column, args.split)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(metrics.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(f"metrics -> {out}")


if __name__ == "__main__":
    main()
