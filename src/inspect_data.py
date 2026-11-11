"""Sanity-check downloaded inputs before running anything expensive.

Intended checks (once files exist):
    * BAM: indexed? sorted? reference sequence names and lengths in the
      header; read length; paired-end flags; read counts per chromosome.
    * The BAM header contigs agree with the FASTA contigs (naming convention,
      e.g. ``chr1`` vs ``1``, and lengths).
    * GTF: same contig naming as the BAM/FASTA; annotation source/version
      recorded in the header lines; number of genes/transcripts/exons.
    * md5 checksums match the values published by ENCODE.

TODO:
    * Decide on a tool to read the BAM header (``samtools view -H`` or
      ``pysam``). Neither is a project dependency yet.
    * Write a short report of the findings into ``results/`` so that
      collaborators can see what was verified.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until data is downloaded."""
    raise NotImplementedError("No input data has been downloaded yet.")


if __name__ == "__main__":
    main()
