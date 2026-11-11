"""Extract a genomic sequence window around each splice site.

For each site in the processed table and each context length in config.yaml (80, 400,
2000, 10000), cut a window centred on the exon/intron junction from the GRCh38 FASTA,
reverse-complement it for minus-strand sites, pad with N past a chromosome end, then
one-hot encode.

Verified coordinate conversion (docs/mcsplicer_notes.md; reference implementation and
checks in src/validation/check_coordinates.py): with L the McSplicer locus (1-based exonic
base) and j the 0-based index of the first base right of the junction,
    + acceptor: j = L-1     + donor: j = L     - acceptor: j = L     - donor: j = L-1
    window = seq[j - W//2 : j + W//2], reverse-complemented if strand is '-'
so [0, W//2) is upstream and [W//2, W) downstream of the junction on both strands.

TODO:
    * Decide whether the context length W counts the junction as centre (as above) and
      how it matches SpliceAI's definition of its 80/400/2000/10000 contexts.
    * The FASTA is a plain gzip; extract or index the contigs first. chr1/chr2 contain
      IUPAC codes besides N (M, R, K, W, Y): decide how to one-hot them.
    * Choose a storage format; at 10 kb per site store coordinates and encode per batch.
    * Tests: plus/minus strand, chromosome-end padding, AG/GT at the window centre.
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented yet."""
    raise NotImplementedError("Not implemented yet; the verified conversion is in the module docstring.")


if __name__ == "__main__":
    main()
