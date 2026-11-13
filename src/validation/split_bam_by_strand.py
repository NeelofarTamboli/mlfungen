"""Split a paired-end BAM into plus-strand and minus-strand fragments (the strand-aware workaround).

sigcount ignores the library strand, so antisense-overlapping genes can pick up each
other's reads. Splitting the BAM by `check_strand.fragment_strand`, and running
exonRefine/sigcount/McSplicer once per strand against a GTF that contains only that
strand's genes, avoids that. Alignment order is preserved, so the outputs are sorted;
each is indexed.

Usage:
    python -m src.validation.split_bam_by_strand IN.bam OUT_PREFIX
writes OUT_PREFIX_plus.bam[.bai] and OUT_PREFIX_minus.bam[.bai].

Tested on the 60 kb K562 subsets only. On a whole 10-13 GB BAM it is one sequential pass.
"""

from __future__ import annotations

import sys

from src.validation.check_strand import fragment_strand


def split(in_bam: str, prefix: str, contig: str | None = None, start: int | None = None, end: int | None = None) -> dict[str, int]:
    """Split in_bam by fragment strand. With `contig` (and optional 0-based start/end) only that
    chromosome/region is read, which needs a .bai; without it the whole file is read in order."""
    import pysam

    src = pysam.AlignmentFile(in_bam)
    outs = {"+": pysam.AlignmentFile(f"{prefix}_plus.bam", "wb", template=src), "-": pysam.AlignmentFile(f"{prefix}_minus.bam", "wb", template=src)}
    n = {"+": 0, "-": 0}
    reads = src.fetch(contig, start, end) if contig else src.fetch(until_eof=True)
    for read in reads:
        s = fragment_strand(read)
        outs[s].write(read)
        n[s] += 1
    for o in outs.values():
        o.close()
    for name in ("plus", "minus"):
        pysam.index(f"{prefix}_{name}.bam")
    return n


if __name__ == "__main__":
    print(split(sys.argv[1], sys.argv[2]))
