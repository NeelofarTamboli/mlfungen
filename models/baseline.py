"""Local sequence baselines for quantitative acceptor usage.

Features are computed from an oriented 100 bp window: the MaxEntScan 23-mer
acceptor score, GC fraction, and normalized 2-mer/3-mer counts. Models are a
training-set mean, ridge regression (alpha selected on validation MAE), and a
histogram gradient-boosting regressor.

The handoff keeps both biological-replicate labels. This script therefore
requires an explicit ``--target`` choice instead of silently deciding how to
combine them. ``--target mean`` is useful for a clearly labelled preliminary
analysis; ``rep1`` and ``rep2`` keep either measurement separate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from analysis.evaluate import evaluate_predictions
from src.config import load_config, project_path
from src.extract_sequences import decode_codes, load_context

_DNA = frozenset("ACGT")


def _kmer_names(k_values: Iterable[int]) -> list[str]:
    return ["kmer_" + "".join(chars) for k in k_values for chars in itertools.product("ACGT", repeat=k)]


def _kmer_frequencies(sequence: str, k: int) -> dict[str, float]:
    counts: dict[str, int] = {}
    valid = 0
    for start in range(len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if set(kmer) <= _DNA:
            counts[kmer] = counts.get(kmer, 0) + 1
            valid += 1
    return {kmer: count / valid for kmer, count in counts.items()} if valid else {}


def local_features(codes: np.ndarray, k_values: Iterable[int] = (2, 3)) -> tuple[np.ndarray, list[str]]:
    """Build MaxEntScan, GC and normalized k-mer features from oriented windows."""
    from maxentpy import maxent

    if codes.ndim != 2 or codes.shape[1] < 40:
        raise ValueError("codes must be a 2D array with enough context around the junction")
    k_values = tuple(k_values)
    names = ["maxentscan_acceptor_23mer", "gc_fraction", *_kmer_names(k_values)]
    name_to_col = {name: i for i, name in enumerate(names)}
    features = np.zeros((len(codes), len(names)), dtype=np.float32)
    features[:, 0] = np.nan
    matrix = maxent.load_matrix3()
    half = codes.shape[1] // 2
    for i, row in enumerate(codes):
        sequence = decode_codes(row)
        canonical = sequence[half - 20 : half + 3]
        if len(canonical) == 23 and set(canonical) <= _DNA:
            features[i, 0] = maxent.score3(canonical, matrix=matrix)
        valid = row < 4
        features[i, 1] = float(np.isin(row[valid], [1, 2]).mean()) if valid.any() else np.nan
        for k in k_values:
            for kmer, frequency in _kmer_frequencies(sequence, k).items():
                features[i, name_to_col[f"kmer_{kmer}"]] = frequency
    return features, names


def _target(sites: pd.DataFrame, choice: str) -> np.ndarray:
    if choice == "mean":
        return sites[["usage_rep1", "usage_rep2"]].mean(axis=1).to_numpy(dtype=float)
    column = {"rep1": "usage_rep1", "rep2": "usage_rep2"}[choice]
    return sites[column].to_numpy(dtype=float)


def fit_baselines(
    sites: pd.DataFrame,
    features: np.ndarray,
    target_choice: str,
    random_state: int = 17,
) -> tuple[pd.DataFrame, dict]:
    """Fit baselines on train, select ridge alpha on validation, predict every row."""
    if "split" not in sites:
        raise ValueError("site table needs a split column")
    y = _target(sites, target_choice)
    train = (sites["split"] == "train").to_numpy() & np.isfinite(y)
    validation = (sites["split"] == "validation").to_numpy() & np.isfinite(y)
    test = (sites["split"] == "test").to_numpy() & np.isfinite(y)
    if not train.any() or not validation.any() or not test.any():
        raise ValueError("train, validation and test must each contain at least one finite target")

    train_mean = float(y[train].mean())
    ridge_scores: dict[float, float] = {}
    ridge_models = {}
    for alpha in (0.01, 0.1, 1.0, 10.0, 100.0):
        model = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), Ridge(alpha=alpha, solver="lsqr"))
        model.fit(features[train], y[train])
        ridge_scores[alpha] = float(mean_absolute_error(y[validation], model.predict(features[validation])))
        ridge_models[alpha] = model
    best_alpha = min(ridge_scores, key=ridge_scores.get)
    ridge = ridge_models[best_alpha]

    gradient_boosting = make_pipeline(
        SimpleImputer(strategy="median"),
        HistGradientBoostingRegressor(
            max_iter=200,
            learning_rate=0.05,
            max_leaf_nodes=15,
            l2_regularization=1.0,
            early_stopping=False,
            random_state=random_state,
        ),
    )
    # Some macOS/Python combinations cannot report physical (rather than
    # logical) core count. Joblib falls back correctly; keep that environment
    # warning out of the scientific output.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Could not find the number of physical cores.*")
        gradient_boosting.fit(features[train], y[train])

    predictions = sites.copy()
    predictions["target"] = y
    predictions["target_definition"] = target_choice
    predictions["prediction_mean"] = train_mean
    predictions["prediction_ridge"] = np.clip(ridge.predict(features), 0.0, 1.0)
    predictions["prediction_gradient_boosting"] = np.clip(gradient_boosting.predict(features), 0.0, 1.0)
    settings = {
        "target": target_choice,
        "n_train": int(train.sum()),
        "n_validation": int(validation.sum()),
        "n_test": int(test.sum()),
        "training_target_mean": train_mean,
        "ridge_validation_mae_by_alpha": {str(key): value for key, value in ridge_scores.items()},
        "selected_ridge_alpha": best_alpha,
        "random_state": random_state,
    }
    return predictions, settings


def main() -> None:
    cfg = load_config()
    processed, results = project_path(cfg, "data_processed"), project_path(cfg, "results")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default=str(processed / "person_b_acceptors_chr1-4.splits.tsv.gz"))
    ap.add_argument("--sequences", default=str(processed / "person_b_acceptors_chr1-4.sequences.npy"))
    ap.add_argument("--target", required=True, choices=["rep1", "rep2", "mean"])
    ap.add_argument("--local-width", type=int, default=100)
    ap.add_argument("--predictions-out", default=str(processed / "baseline_predictions_chr1-4.tsv.gz"))
    ap.add_argument("--metrics-out", default=str(results / "baseline_metrics_chr1-4.tsv"))
    ap.add_argument("--settings-out", default=str(results / "baseline_settings_chr1-4.json"))
    ap.add_argument("--features-out", default=str(processed / "baseline_features_chr1-4.npz"))
    ap.add_argument("--random-state", type=int, default=17)
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str, "strand": str})
    codes = load_context(args.sequences, args.local_width)
    if len(codes) != len(sites):
        raise ValueError(f"sequence rows ({len(codes)}) do not match site rows ({len(sites)})")
    features, feature_names = local_features(codes)
    features_path = Path(args.features_out)
    features_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(features_path, X=features, feature_names=np.asarray(feature_names))
    predictions, settings = fit_baselines(sites, features, args.target, args.random_state)

    predictions_path, metrics_path = Path(args.predictions_out), Path(args.metrics_out)
    predictions_path.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(predictions_path, sep="\t", index=False, compression="gzip" if predictions_path.suffix == ".gz" else None)
    prediction_columns = [column for column in predictions if column.startswith("prediction_")]
    metrics = evaluate_predictions(predictions, prediction_columns)
    metrics.to_csv(metrics_path, sep="\t", index=False, float_format="%.6g")
    settings.update({"local_width": args.local_width, "features": feature_names})
    settings_path = Path(args.settings_out)
    settings_path.parent.mkdir(parents=True, exist_ok=True)
    settings_path.write_text(json.dumps(settings, indent=2) + "\n")
    print(metrics.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(f"features -> {features_path}")
    print(f"predictions -> {predictions_path}")
    print(f"metrics -> {metrics_path}")


if __name__ == "__main__":
    main()
