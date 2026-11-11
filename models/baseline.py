"""Simple baselines for predicting splice-site usage.

Planned baselines, simplest first:
    * Constant / mean predictor, split by site type (a floor any model
      must beat).
    * Local regression (linear or gradient boosting) on local features:
      MaxEntScan score (9-mer for donors, 23-mer for acceptors), k-mer counts
      in about +/-50 bp, GC content, exon/intron lengths.
    * Later, if practical: pretrained MaxEntScan, SpliceAI and Pangolin scores
      used directly as predictors. Pangolin has no K562 output, so that would
      be a cross-tissue comparison and must be described that way.

TODO:
    * Depends on the processed site table and the chromosome splits.
    * Choose the MaxEntScan implementation (e.g. ``maxentpy``) and check that
      its expected input windows match our coordinate convention.
    * Fit on train, tune on validation, and report on test only at the end.
"""

from __future__ import annotations
