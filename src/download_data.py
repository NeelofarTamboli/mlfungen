"""Download the input data for the project into ``data/raw/``.

Inputs (chosen and checked, see docs/encode_inputs.md and config.yaml):
    * ENCFF864OUS (replicate 1) and ENCFF696CSX (replicate 2): GRCh38, GENCODE V29
      paired-end alignments for ENCSR000AEO, 10.4 GB and 13.4 GB.
    * gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz (42.7 MB).
    * GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz (873 MB).

Nothing is downloaded yet; the full BAMs are deliberately held back for now.

TODO:
    * Implement the download from the encode-public S3 URLs listed in the ENCODE file
      records (``/files/<accession>/?format=json`` -> ``cloud_metadata.url``), with md5
      verification against the values in docs/encode_inputs.md, resumable, skipping
      files that are already present.
    * Do not commit the downloaded files (data/ is git-ignored).
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
