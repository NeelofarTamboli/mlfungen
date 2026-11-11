"""Combine per-replicate usage tables into the final labelled site table.

First regression target: internal ACCEPTOR usage only (config `labels.target_site_type`).
Donor rows are kept for QC/reporting but not in the regression dataset, because McSplicer
donor values are near-binary in this release (docs/mcsplicer_notes.md).

Planned steps:
    1. Load the parsed tables for replicate 1 and replicate 2.
    2. Join on (chrom, locus, strand, site type).
    3. Measure replicate agreement (Pearson and Spearman on usage, acceptors and, for the
       report, donors) at several read-support thresholds. This is roughly the ceiling for
       a model, so it is saved to `results/` and reported next to every model result.
    4. Only after looking at that agreement, decide whether and how to combine the
       replicates (config `labels.combine_replicates` is unset on purpose).
    5. Apply the read-support filter (config `labels.min_read_support`, unset on purpose)
       and report how many acceptors pass, per chromosome.
    6. Write the site table to `data/processed/`.

TODO:
    * Depends on parse_mcsplicer.py.
    * Choose the read-support threshold and the combination rule from the real data.
    * Decide the output format (TSV or Parquet) and columns.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented yet."""
    raise NotImplementedError("Requires parsed McSplicer tables (parse_mcsplicer.py).")


if __name__ == "__main__":
    main()
