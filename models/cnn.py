"""CNN models that predict splice-site usage from one-hot sequence.

NOT IMPLEMENTED YET, on purpose: the CNN comes after the data pipeline and
the baselines are working.

Plan (from implementation_plan.md):
    * One model per context length (80, 400, 2000, 10000 bp), from the same
      architecture family so they can be compared fairly, each trained with
      several random seeds.
    * A plain CNN for short windows; SpliceAI-style dilated residual blocks
      for 2 kb and 10 kb.
    * One output per site (the prediction at the centre position).
    * Loss: binary cross-entropy on the usage value (as in SpliceAI and
      Pangolin), optionally weighted by read support.

TODO:
    * Add PyTorch to requirements when this is implemented.
    * Confirm GPU access (Colab / Google Cloud / course cluster) before the
      2 kb and 10 kb models.
"""

from __future__ import annotations
