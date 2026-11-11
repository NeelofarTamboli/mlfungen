"""Exploratory analysis of the labelled splice-site table.

Planned:
    * Usage distribution for donors and acceptors (expect many sites near 1).
    * Replicate 1 vs replicate 2 usage scatter, with Pearson and Spearman
      correlations at different read-support thresholds.
    * Number of sites passing the read-support filter, by site type and split.
    * Fraction of "alternative" sites (e.g. usage between 0.1 and 0.9).
    * Sanity check: the dinucleotide at each site (GT/AG for most donors and
      acceptors). If it is off, the coordinate convention is probably wrong.

TODO:
    * Depends on build_dataset.py and extract_sequences.py.
    * Save figures to ``figures/`` and summary tables to ``results/``.
"""

from __future__ import annotations
