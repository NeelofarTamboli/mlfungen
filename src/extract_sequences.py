"""Extract a genomic sequence window around each splice site.

For each site in the processed table and each context length in
``config.yaml`` (80, 400, 2000, 10000 bp), cut a window from the reference
FASTA centred on the site, reverse-complement it for minus-strand sites, and
pad with ``N`` where the window runs past a chromosome end. Then one-hot
encode it for the models.

TODO:
    * Genome build and FASTA are not chosen yet (see config.yaml).
    * Define exactly where the site sits within the window. This depends on
      the McSplicer coordinate convention, which is not yet verified (see
      parse_mcsplicer.py).
    * Define how the context length is split into flanks (e.g. 10000 bp =
      5000 upstream + 5000 downstream) and whether it counts the site
      itself. Keep this consistent with how SpliceAI defines its context
      lengths, so the results can be compared.
    * Choose a FASTA reader (e.g. pyfaidx or pysam). Neither is a project
      dependency yet.
    * Choose a storage format. At 10 kb per site, one-hot arrays get large,
      so storing coordinates and encoding when loading a batch may be
      better than storing every window.
    * Write tests: plus/minus strand, chromosome-end padding, a known
      GT/AG dinucleotide at annotated donors/acceptors.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented yet."""
    raise NotImplementedError("Genome build, FASTA, and site coordinates are not verified yet.")


if __name__ == "__main__":
    main()
