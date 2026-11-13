"""Join the two replicates' internal-acceptor tables into one label table.

First regression target: internal ACCEPTOR usage only (config `labels.target_site_type`); donor rows
stay in the QC outputs because McSplicer donor values are near-binary (docs/mcsplicer_notes.md).

This step only lines the replicates up. It does not average them, filter on read support, or pick
chromosomes: those decisions (config `labels.combine_replicates`, `labels.min_read_support`,
`splits.*`) are left open until the QC in analysis/ has been looked at.

A site is the same site in both replicates when (gene_id, chrom, locus, strand) match. Sites seen in
only one replicate are kept, with empty values for the other, and flagged in `in_rep1` / `in_rep2`.

Output: data/processed/acceptor_labels.tsv, one row per site, columns
    site_id, gene_id, chrom, locus, strand, site_class,
    usage_rep1/2, bootstrap_sd_rep1/2, read_support_rep1/2,        <- the label columns
    n_bootstraps_rep1/2   (bootstrap_sd is empty where no bootstraps were run; the full genome was run without them)
    has_read_support_rep1/2   read_support > 0
    also_transcript_start     internal acceptor that is also a transcript start in another transcript (kept, but see in_primary_qc)
    in_primary_qc             in both replicates, nonzero support in both, not also_transcript_start (no other filter)
    spliced_in_rep1/2, contiguous_in_rep1/2, n_incompat_warnings_rep1/2,
    in_rep1, in_rep2, in_both, index_rep1/2, raw_csv_rep1/2        <- QC flags and trace back to raw output

With --handoff it also writes the compact table for the sequence-extraction step (primary-QC sites only):
    site_id, gene_id, chrom, locus, strand, usage_rep1/2, read_support_rep1/2, also_transcript_start, in_primary_qc

Usage:
    python -m src.build_dataset [--rep1 data/processed/acceptors_rep1.tsv] [--rep2 ...] [--out ...] [--chroms ...] [--handoff PATH]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.config import load_config, project_path

KEY = ["gene_id", "chrom", "locus", "strand"]
PER_REP = ["usage", "bootstrap_sd", "n_bootstraps", "read_support", "spliced_in", "contiguous_in", "n_incompat_warnings", "index", "raw_csv",
           "also_transcript_start"]


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "strand": str})
    if df.duplicated(KEY).any():
        raise ValueError(f"{path}: duplicated sites, cannot join")
    return df[KEY + PER_REP]


def build(rep1: Path, rep2: Path, chroms: list[str] | None = None) -> pd.DataFrame:
    a, b = load(rep1), load(rep2)
    if chroms:
        a, b = a[a.chrom.isin(chroms)], b[b.chrom.isin(chroms)]
    df = a.merge(b, on=KEY, how="outer", suffixes=("_rep1", "_rep2"), indicator=True)
    df["in_rep1"] = df["_merge"].isin(["both", "left_only"])
    df["in_rep2"] = df["_merge"].isin(["both", "right_only"])
    df["in_both"] = df["_merge"] == "both"
    df["site_class"] = "acceptor"
    df["also_transcript_start"] = df["also_transcript_start_rep1"].combine_first(df["also_transcript_start_rep2"]).astype(bool)
    for r in ("rep1", "rep2"):
        df[f"has_read_support_{r}"] = df[f"read_support_{r}"].fillna(0) > 0
    # primary QC set: in both replicates, reads at the acceptor in both, and not a transcript start as well
    df["in_primary_qc"] = df["in_both"] & df["has_read_support_rep1"] & df["has_read_support_rep2"] & ~df["also_transcript_start"]
    df["site_id"] = df["chrom"] + ":" + df["locus"].astype(str) + ":" + df["strand"] + ":" + df["gene_id"]
    cols = (["site_id", "gene_id", "chrom", "locus", "strand", "site_class",
             "usage_rep1", "usage_rep2", "bootstrap_sd_rep1", "bootstrap_sd_rep2", "n_bootstraps_rep1", "n_bootstraps_rep2",
             "read_support_rep1", "read_support_rep2", "has_read_support_rep1", "has_read_support_rep2", "also_transcript_start", "in_primary_qc",
             "spliced_in_rep1", "spliced_in_rep2", "contiguous_in_rep1", "contiguous_in_rep2",
             "n_incompat_warnings_rep1", "n_incompat_warnings_rep2", "in_rep1", "in_rep2", "in_both",
             "index_rep1", "index_rep2", "raw_csv_rep1", "raw_csv_rep2"])
    return df[cols].sort_values(["chrom", "locus", "strand", "gene_id"]).reset_index(drop=True)


HANDOFF_COLUMNS = ["site_id", "gene_id", "chrom", "locus", "strand", "usage_rep1", "usage_rep2",
                   "read_support_rep1", "read_support_rep2", "also_transcript_start", "in_primary_qc"]


def handoff(df: pd.DataFrame) -> pd.DataFrame:
    """Compact table for the sequence-extraction step: primary-QC sites only, labels and support, nothing raw."""
    out = df.loc[df.in_primary_qc, HANDOFF_COLUMNS].copy()
    return out.round({"usage_rep1": 6, "usage_rep2": 6, "read_support_rep1": 2, "read_support_rep2": 2})


def main() -> None:
    proc = project_path(load_config(), "data_processed")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep1", default=str(proc / "acceptors_rep1.tsv"))
    ap.add_argument("--rep2", default=str(proc / "acceptors_rep2.tsv"))
    ap.add_argument("--out", default=str(proc / "acceptor_labels.tsv"))
    ap.add_argument("--chroms", nargs="+", help="only these chromosomes (e.g. chr1 chr2 chr3 chr4)")
    ap.add_argument("--handoff", help="also write the compact primary-QC table for the sequence step to this path")
    a = ap.parse_args()
    df = build(Path(a.rep1), Path(a.rep2), a.chroms)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)
    if a.handoff:
        h = handoff(df)
        h.to_csv(a.handoff, sep="\t", index=False)
        print(f"handoff table: {len(h)} primary-QC sites -> {a.handoff}")
    print(f"{len(df)} internal acceptors: {int(df.in_both.sum())} in both replicates, "
          f"{int((df.in_rep1 & ~df.in_rep2).sum())} only rep1, {int((df.in_rep2 & ~df.in_rep1).sum())} only rep2 -> {a.out}")


if __name__ == "__main__":
    main()
