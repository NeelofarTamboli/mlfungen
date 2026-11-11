"""Download the input data for the project into ``data/raw/``.

Planned inputs:
    * Processed RNA-seq alignments (BAM + index) for both biological
      replicates of ENCODE experiment ENCSR000AEO (K562).
    * The reference genome FASTA matching the assembly the BAMs were aligned to.
    * The gene annotation GTF matching the annotation version ENCODE used.

Nothing is downloaded yet. The BAMs are large, and the files have not been
chosen. Before writing the download logic:

TODO:
    * On the ENCODE portal, list the files for ENCSR000AEO and record, for
      each candidate BAM: accession, biological/technical replicate, assembly,
      annotation version, pipeline version, file status (released vs
      archived/revoked), size, and md5sum.
    * Pick one BAM per biological replicate from the same pipeline run.
    * Record the chosen accessions in ``config.yaml`` (dataset.replicates).
    * Confirm the read length and whether the library is paired-end and
      stranded; McSplicer takes ``--read_len`` as an argument.
    * Choose the matching FASTA and GTF and record them in config.yaml.
    * Implement the download with md5 verification, and skip files that are
      already present.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until the files are chosen."""
    raise NotImplementedError(
        "Choose ENCODE file accessions and reference files first "
        "(see the TODOs in this module and in config.yaml)."
    )


if __name__ == "__main__":
    main()
