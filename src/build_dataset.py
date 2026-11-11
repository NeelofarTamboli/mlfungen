"""Combine per-replicate usage tables into the final labelled site table.

Planned steps:
    1. Load the parsed McSplicer tables for replicate 1 and replicate 2.
    2. Join them on (chrom, position, strand, site type).
    3. Measure how well the replicates agree (Pearson and Spearman on usage,
       for donors and acceptors separately, and for different read-support
       thresholds). This agreement is roughly the best score a model can
       reach, so it is saved to ``results/`` and reported next to every model
       result.
    4. Only after looking at that agreement, decide whether and how to
       combine the replicates (e.g. mean of the two, pooled-read estimate,
       or keeping only sites that agree).
    5. Apply the read-support filter and report how many sites pass it,
       by site type.
    6. Write the site table to ``data/processed/``.

TODO:
    * Depends on the output of ``parse_mcsplicer.py``.
    * Choose the read-support threshold (``labels.min_read_support`` in
      config.yaml) from the data, not in advance.
    * Decide on the output format (e.g. TSV or Parquet) and its columns.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented yet."""
    raise NotImplementedError("Requires parsed McSplicer tables.")


if __name__ == "__main__":
    main()
