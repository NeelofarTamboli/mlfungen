"""Strand rule for the reverse-stranded (dUTP) K562 libraries, and a check of it.

Rule: for a fragment, read 1 reverse OR read 2 forward  -> `+` strand transcript;
read 1 forward OR read 2 reverse -> `-` strand transcript.

`fragment_strand` is the rule. Running this module measures, for chosen genes, what
fraction of uniquely mapped alignments lying fully inside exons that overlap no gene
on the other strand follow it (expect ~100% for the library to be as annotated),
and the strand mix in exon regions shared by antisense genes.

Usage:
    python -m src.validation.check_strand REGION.bam GTF GENE_PLUS GENE_MINUS
(gene names as in the GTF `gene_name`, e.g. ATP5PB WDR77)
"""

from __future__ import annotations

import collections
import re
import sys


def fragment_strand(read) -> str:
    """Transcript strand implied by a pysam AlignedSegment, for a reverse-stranded paired library."""
    if read.is_read1:
        return "+" if read.is_reverse else "-"
    return "-" if read.is_reverse else "+"


def _merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def exons_by_gene(gtf: str):
    ex = collections.defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "exon":
                attrs = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
                ex[(attrs["gene_name"], f[6])].append((int(f[3]), int(f[4])))
    return {k: _merge(v) for k, v in ex.items()}


def main() -> None:
    import pysam

    bam_path, gtf, gene_plus, gene_minus = sys.argv[1:5]
    ex = exons_by_gene(gtf)
    plus, minus = ex[(gene_plus, "+")], ex[(gene_minus, "-")]
    overlap = [(max(a, c), min(b, d)) for a, b in plus for c, d in minus if max(a, c) <= min(b, d)]
    only = lambda ivs, other: [x for x in ivs if not any(x[0] <= o[1] and o[0] <= x[1] for o in other)]
    zones = {
        f"{gene_plus} exons not shared with {gene_minus} (+ gene)": only(plus, minus),
        f"{gene_minus} exons not shared with {gene_plus} (- gene)": only(minus, plus),
        f"exon overlap {gene_plus}/{gene_minus}": overlap,
    }
    counts = {k: collections.Counter() for k in zones}
    for r in pysam.AlignmentFile(bam_path).fetch(until_eof=True):
        if r.is_secondary or r.get_tag("NH") != 1:
            continue
        s, e = r.reference_start + 1, r.reference_end
        for name, ivs in zones.items():
            if any(s >= a and e <= b for a, b in ivs):
                counts[name][fragment_strand(r)] += 1
    print(f"exon-overlap intervals: {overlap}")
    for name, c in counts.items():
        n = sum(c.values())
        print(f"  {name:55s} n={n:6d}  +: {c['+'] / n:.3f}  -: {c['-'] / n:.3f}" if n else f"  {name}: no reads")


if __name__ == "__main__":
    main()
