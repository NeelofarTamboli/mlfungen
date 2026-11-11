"""Assign each splice site to train / validation / test by chromosome.

The split is by chromosome so that nearby, correlated sites do not end up
on both sides. The planned test chromosomes (chr1, chr3, chr5, chr7, chr9)
are the ones SpliceAI and Pangolin hold out, so the pretrained baselines
have not been trained on our test sites. They are set in config.yaml.

TODO:
    * Choose the validation chromosomes (``splits.validation_chromosomes``).
    * Decide how to handle chrX, chrY, chrM and unplaced contigs.
    * Check that contig naming in the site table matches the config
      (``chr1`` vs ``1``).
    * Decide how to handle genes that SpliceAI/Pangolin training may have
      included through paralogs on other chromosomes, if at all.
    * Report how many sites fall in each split, by site type.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented yet."""
    raise NotImplementedError("Requires the processed site table.")


if __name__ == "__main__":
    main()
