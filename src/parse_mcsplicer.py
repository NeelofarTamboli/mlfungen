"""Parse McSplicer output into a tidy per-site usage table.

According to the McSplicer README, the output is a CSV with one row per
(bootstrap step, splice site) holding: the bootstrap step, a splice-site
index, gene strand, chromosome, splice-site genome position, and the usage
estimate. Step 0 is the estimate from the observed counts; steps 1..n are
bootstrap estimates. Site indices look like ``s0, s1, ...`` and
``e0, e1, ...`` and are numbered in order within each gene.

We have not seen a real output file yet. Do not write the parser until we
have one.

TODO (verify on McSplicer's bundled examples, then on our data):
    * Exact column names, delimiter, and header line.
    * Whether a gene identifier is included in each row, or only in the file
      name or prefix.
    * How ``s``/``e`` map to donor/acceptor on each strand. The README calls
      them "3' start" and "5' end" sites; check how this is defined for
      minus-strand genes.
    * Coordinate convention of the position: 0- or 1-based, and whether it
      is the last exonic base, the first intronic base, or the junction
      boundary. This must be right, or the sequence windows will be
      shifted.
    * How sites with no read support are reported (missing, NaN, or a prior).
    * Whether the output includes coverage or read support per site. If not,
      it must be computed separately (e.g. from the sigcount output or from
      junction reads in the BAM) for the read-support filter.
    * Duplicate sites shared by overlapping genes.

Planned output (one row per site per replicate): chrom, position, strand,
site type (donor/acceptor), gene, usage, bootstrap spread, read support.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until the format is known."""
    raise NotImplementedError("McSplicer output format has not been verified yet.")


if __name__ == "__main__":
    main()
