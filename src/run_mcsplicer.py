"""Run McSplicer on each replicate's BAM to estimate splice-site usage.

McSplicer: Alqassem et al., Bioinformatics 2021,
https://doi.org/10.1093/bioinformatics/btab050 (code: canzarlab/McSplicer).
Verified settings and findings: docs/mcsplicer_notes.md. Planned flow per replicate
(K562, strand-aware, see config.yaml `mcsplicer`):

    1. Split the BAM by inferred transcript strand
       (src/validation/split_bam_by_strand.py; read 1 reverse or read 2 forward = +).
    2. Split the original GENCODE v29 GTF by strand.
    3. For each strand: exonRefine -> bin/sigcount (not sigcount_pe, no -s)
       -> python_scripts/McSplicer.py --read_len 100 --bootstraps N.
    4. Concatenate the per-gene CSVs of both strands.

All of this was run only on a 60 kb subset; nothing has been run on a full BAM.

TODO:
    * Implement the flow above (the commands are in docs/mcsplicer_notes.md).
    * Measure runtime/memory of sigcount and McSplicer on a full 10-13 GB BAM;
      try one chromosome first.
    * Choose the number of bootstraps.
    * Record the McSplicer commit hash used (d6a22c3 was tested).
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until McSplicer is tested."""
    raise NotImplementedError("Not implemented yet; the planned flow is in the module docstring.")


if __name__ == "__main__":
    main()
