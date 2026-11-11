"""Run McSplicer on each replicate's BAM to estimate splice-site usage.

McSplicer: Alqassem et al., Bioinformatics 2021,
https://doi.org/10.1093/bioinformatics/btab050
Code: https://github.com/canzarlab/McSplicer

Steps, according to the McSplicer README (not yet tried by us):
    1. ``exonRefine <annotation.gtf> --prefix OUTPUT_PREFIX``
       splits the annotated exons into non-overlapping sub-exonic segments
       ("signatures") and writes a refined GTF.
    2. ``sigcount <alignments.bam> <annotation_refined.gtf> <outfile-prefix>``
       counts reads on each signature and signature junction. The BAM must be
       sorted and indexed.
    3. ``python3 python_scripts/McSplicer.py --gtf REFINED_GTF
       --count_file SIGNATURE_COUNT_FILE --out_dir OUT --bootstraps N
       --read_len READ_LENGTH --prefix PREFIX``

The README says the input GTF may be a reference annotation or a GTF
assembled from the RNA-seq reads (e.g. StringTie). That choice determines
which candidate splice sites get a usage estimate.

TODO (verify before implementing):
    * Do the ``bin/exonRefine`` and ``bin/sigcount`` binaries in the repo
      run on our platform, or do they have to be built from source?
    * Does exonRefine accept a full GENCODE GTF directly, or does it need
      filtering (e.g. by gene type or transcript support level)?
    * How does sigcount handle paired-end and stranded reads, multi-mapped
      reads, and duplicates?
    * What ``--read_len`` should be for our BAMs (check the actual reads)?
    * Run time and memory for a whole genome; first try a single chromosome
      (subset the BAM and GTF).
    * Whether to pool replicates by merging BAMs, as well as running each
      replicate on its own.
    * Record the McSplicer commit hash we use, for reproducibility.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until McSplicer is tested."""
    raise NotImplementedError(
        "Install McSplicer and test it on the bundled examples and on one "
        "chromosome before scripting the full run."
    )


if __name__ == "__main__":
    main()
