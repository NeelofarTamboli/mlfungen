"""Sanity-check downloaded inputs before running anything expensive.

Already checked on headers and a 60 kb subset (docs/encode_inputs.md; code in
src/validation/): BAM contigs vs FASTA vs GTF, sorted order, paired-end flags, read length,
no XS tag. Once the full files are downloaded, repeat on them:
    * md5 of each BAM, the GTF and the FASTA against the values ENCODE lists;
    * the BAM header contig list equals the one recorded in docs/encode_inputs.md;
    * BAM index (ENCODE ships none; create with pysam/samtools) and per-chromosome read
      counts for both replicates.

TODO:
    * Implement these checks and write a short report into ``results/``.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until data is downloaded."""
    raise NotImplementedError("No input data has been downloaded yet.")


if __name__ == "__main__":
    main()
