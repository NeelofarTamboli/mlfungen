"""Run the reproducible Person B sequence-and-baseline pipeline.

This is the one user-facing entry point for the modules in ``src/``, ``models/``
and ``analysis/``. Large/reproducible artifacts are kept under git-ignored
``data/``; only code and compact QC/metric summaries go under ``results/``.

Current chr1-4 handoff (preliminary target/split):
    python -m src.run_part_b

Later full-genome handoff, for example:
    python -m src.run_part_b \
      --sites results/person_b_acceptors_all.tsv \
      --tag all \
      --validation chr2 chr4 \
      --target mean

The runner is resumable: an existing stage output is skipped. Pass ``--force``
to recompute selected stages, and ``--steps`` to run only part of the pipeline.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

from src.config import REPO_ROOT, load_config

STAGES = ("download", "extract", "split", "baseline", "evaluate")


def run(command: list[str]) -> None:
    """Run a pipeline module with the active Python interpreter."""
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default="results/person_b_acceptors_chr1-4.tsv", help="portable input TSV path")
    ap.add_argument("--tag", default="chr1-4", help="short label used in generated file names")
    ap.add_argument("--validation", nargs="+", default=["chr4"], help="validation chromosomes")
    ap.add_argument("--target", choices=["rep1", "rep2", "mean"], default="mean")
    ap.add_argument("--steps", nargs="+", choices=STAGES, default=list(STAGES))
    ap.add_argument("--force", action="store_true", help="recompute outputs that already exist")
    args = ap.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.tag):
        ap.error("--tag may contain only letters, numbers, dots, underscores and hyphens")

    sites = Path(args.sites)
    if not sites.is_absolute():
        sites = REPO_ROOT / sites
    if not sites.exists():
        ap.error(f"site table does not exist: {sites}")

    tag = args.tag
    raw_fasta = REPO_ROOT / cfg["reference"]["genome_fasta"]
    prefix = Path("data/processed") / f"person_b_acceptors_{tag}"
    sequence_path = Path(f"{prefix}.sequences.npy")
    metadata_path = Path(f"{prefix}.metadata.tsv.gz")
    splits_path = Path(f"{prefix}.splits.tsv.gz")
    features_path = Path("data/processed") / f"baseline_features_{tag}.npz"
    predictions_path = Path("data/processed") / f"baseline_predictions_{tag}.tsv.gz"
    qc_path = Path("results") / f"sequence_qc_{tag}.json"
    split_counts_path = Path("results") / f"split_counts_{tag}.tsv"
    split_settings_path = Path("results") / f"split_settings_{tag}.json"
    metrics_path = Path("results") / f"baseline_metrics_{tag}.tsv"
    baseline_settings_path = Path("results") / f"baseline_settings_{tag}.json"

    python = sys.executable

    if "download" in args.steps:
        if args.force or not raw_fasta.exists():
            run([python, "-m", "src.download_data", "fasta"])
        else:
            print(f"skip download (exists: {raw_fasta.relative_to(REPO_ROOT)})")

    if "extract" in args.steps:
        extract_outputs = (sequence_path, metadata_path, qc_path)
        if args.force or any(not (REPO_ROOT / path).exists() for path in extract_outputs):
            run([
                python, "-m", "src.extract_sequences",
                "--sites", str(sites),
                "--out-prefix", str(prefix),
                "--qc-out", str(qc_path),
            ])
        else:
            print(f"skip extract (exists: {sequence_path})")

    if "split" in args.steps:
        split_outputs = (splits_path, split_counts_path, split_settings_path)
        if args.force or any(not (REPO_ROOT / path).exists() for path in split_outputs):
            run([
                python, "-m", "src.make_splits",
                "--sites", str(metadata_path),
                "--out", str(splits_path),
                "--summary-out", str(split_counts_path),
                "--settings-out", str(split_settings_path),
                "--validation", *args.validation,
            ])
        else:
            print(f"skip split (exists: {splits_path})")

    if "baseline" in args.steps:
        baseline_outputs = (features_path, predictions_path, metrics_path, baseline_settings_path)
        if args.force or any(not (REPO_ROOT / path).exists() for path in baseline_outputs):
            run([
                python, "-m", "models.baseline",
                "--sites", str(splits_path),
                "--sequences", str(sequence_path),
                "--target", args.target,
                "--features-out", str(features_path),
                "--predictions-out", str(predictions_path),
                "--metrics-out", str(metrics_path),
                "--settings-out", str(baseline_settings_path),
            ])
        else:
            print(f"skip baseline (exists: {predictions_path})")

    if "evaluate" in args.steps:
        if args.force or not (REPO_ROOT / metrics_path).exists():
            run([
                python, "-m", "analysis.evaluate",
                "--predictions", str(predictions_path),
                "--out", str(metrics_path),
            ])
        else:
            print(f"skip evaluate (exists: {metrics_path})")

    print("Person B pipeline complete.")
    print(f"Compact results: results/sequence_qc_{tag}.json, results/split_counts_{tag}.tsv, "
          f"results/baseline_metrics_{tag}.tsv")


if __name__ == "__main__":
    main()
