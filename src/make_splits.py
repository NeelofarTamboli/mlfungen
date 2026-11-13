"""Assign splice sites to chromosome-disjoint train, validation and test sets.

Test chromosomes come from ``config.yaml`` (chr1, chr3, chr5, chr7, chr9).
Validation chromosomes must be supplied in the config or with ``--validation``;
all remaining chromosomes become training chromosomes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

from src.config import load_config, project_path


def assign_splits(
    sites: pd.DataFrame,
    test_chromosomes: Iterable[str],
    validation_chromosomes: Iterable[str],
) -> pd.DataFrame:
    """Return a copy with a chromosome-derived ``split`` column."""
    if "chrom" not in sites:
        raise ValueError("site table needs a chrom column")
    test, validation = set(test_chromosomes), set(validation_chromosomes)
    overlap = test & validation
    if overlap:
        raise ValueError(f"test and validation chromosomes overlap: {sorted(overlap)}")
    out = sites.copy()
    out["split"] = "train"
    out.loc[out["chrom"].isin(validation), "split"] = "validation"
    out.loc[out["chrom"].isin(test), "split"] = "test"
    return out


def split_summary(sites: pd.DataFrame) -> pd.DataFrame:
    """Count rows by split and chromosome."""
    return (
        sites.groupby(["split", "chrom"], observed=True)
        .size()
        .rename("n_sites")
        .reset_index()
        .sort_values(["split", "chrom"])
    )


def main() -> None:
    cfg = load_config()
    default_input = project_path(cfg, "data_processed") / "person_b_acceptors_chr1-4.metadata.tsv.gz"
    default_output = project_path(cfg, "data_processed") / "person_b_acceptors_chr1-4.splits.tsv.gz"
    results = project_path(cfg, "results")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", default=str(default_input))
    ap.add_argument("--out", default=str(default_output))
    ap.add_argument("--summary-out", default=str(results / "split_counts_chr1-4.tsv"))
    ap.add_argument("--settings-out", default=str(results / "split_settings_chr1-4.json"))
    ap.add_argument("--test", nargs="+", default=cfg["splits"]["test_chromosomes"])
    ap.add_argument("--validation", nargs="+", default=cfg["splits"]["validation_chromosomes"])
    args = ap.parse_args()
    if not args.validation:
        ap.error("choose at least one validation chromosome with --validation (the config value is still null)")

    sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str, "strand": str})
    split_sites = assign_splits(sites, args.test, args.validation)
    counts = split_sites["split"].value_counts()
    missing = [name for name in ("train", "validation", "test") if counts.get(name, 0) == 0]
    if missing:
        raise ValueError(f"empty split(s) {missing}; choose chromosomes that leave data for every split")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    split_sites.to_csv(out, sep="\t", index=False, compression="gzip" if out.suffix == ".gz" else None)
    summary = split_summary(split_sites)
    summary_path = Path(args.summary_out)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    settings_path = Path(args.settings_out)
    settings_path.write_text(json.dumps({"test_chromosomes": args.test, "validation_chromosomes": args.validation}, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"split table -> {out}")


if __name__ == "__main__":
    main()
