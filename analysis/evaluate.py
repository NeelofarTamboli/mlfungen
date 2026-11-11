"""Evaluate model predictions against measured splice-site usage.

Planned metrics, reported separately for donors and acceptors, and also for
alternative sites only (e.g. usage between 0.1 and 0.9), because sites with
usage near 1 inflate correlations:
    * Pearson and Spearman correlation, R^2, MAE.
    * Replicate agreement on the same sites, as an approximate ceiling.
    * Mean and spread over random seeds for each context length.

Later analysis: per-site error change between the 400 bp and 10 kb models,
grouped by site properties (constitutive vs alternative, nearby competing
sites, exon/intron length), to see which sites benefit from long-range
sequence.

TODO:
    * Decide the prediction file format once models exist.
    * Use only the held-out test chromosomes for final numbers.
"""

from __future__ import annotations
