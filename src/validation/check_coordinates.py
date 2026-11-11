"""Check McSplicer's locus convention against the reference sequence.

Convention (see docs/mcsplicer_notes.md): the locus L of an `s` (start) site is
the 1-based first exonic base in transcription direction, and of an `e` (end)
site the last exonic base. For every annotated internal junction of a
chromosome this reads the intron bases next to the junction and reports how
often they are GT (donor) / AG (acceptor), with the locus shifted by -2..+2 as a
control. It also checks the centred-window helper used for sequence extraction.

Needs a GENCODE-style GTF and `<seq-dir>/<chrom>.seq` (see stream_fasta.py).

Usage:
    python -m src.validation.check_coordinates GTF[.gz] SEQ_DIR [--chrom chr1]
"""

from __future__ import annotations

import argparse
import collections
import gzip
import re
from itertools import pairwise
from pathlib import Path

_COMP = str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def load_sequence(seq_dir: str | Path, chrom: str) -> str:
    return (Path(seq_dir) / f"{chrom}.seq").read_bytes().upper().decode()


def load_transcripts(gtf: str | Path, chrom: str):
    """Return ({gene_id: {transcript_id: [(start, end), ...]}}, {gene_id: strand}) for one chromosome."""
    transcripts = collections.defaultdict(lambda: collections.defaultdict(list))
    strand: dict[str, str] = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom or f[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
            transcripts[attrs["gene_id"]][attrs["transcript_id"]].append((int(f[3]), int(f[4])))
            strand[attrs["gene_id"]] = f[6]
    return transcripts, strand


def junctions(exons, strand: str):
    """Yield (donor_locus, acceptor_locus) for each intron of one transcript, in McSplicer's convention."""
    ordered = sorted(exons, reverse=(strand == "-"))
    for up, down in pairwise(ordered):
        yield (up[1], down[0]) if strand == "+" else (up[0], down[1])


def intron_dinucleotide(seq: str, locus: int, kind: str, strand: str) -> str:
    """The two intron bases next to the junction, read 5'->3' on the transcript strand."""
    if (kind, strand) == ("acceptor", "+"):
        return seq[locus - 3 : locus - 1]
    if (kind, strand) == ("donor", "+"):
        return seq[locus : locus + 2]
    if (kind, strand) == ("acceptor", "-"):
        return revcomp(seq[locus : locus + 2])
    if (kind, strand) == ("donor", "-"):
        return revcomp(seq[locus - 3 : locus - 1])
    raise ValueError((kind, strand))


def junction_index(locus: int, kind: str, strand: str) -> int:
    """0-based index of the first base to the right (higher coordinate) of the exon/intron junction."""
    return {("acceptor", "+"): locus - 1, ("donor", "+"): locus, ("acceptor", "-"): locus, ("donor", "-"): locus - 1}[
        (kind, strand)
    ]


def junction_window(seq: str, locus: int, kind: str, strand: str, width: int) -> str | None:
    """`width` bases centred on the junction, oriented 5'->3' on the transcript strand.

    Characters [0, width/2) are upstream of the junction, [width/2, width) downstream.
    Returns None if the window would run off the chromosome.
    """
    j = junction_index(locus, kind, strand)
    a, b = j - width // 2, j + width // 2
    if a < 0 or b > len(seq):
        return None
    w = seq[a:b]
    return w if strand == "+" else revcomp(w)


def window_dinucleotide(window: str, kind: str) -> str:
    h = len(window) // 2
    return window[h : h + 2] if kind == "donor" else window[h - 2 : h]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gtf")
    ap.add_argument("seq_dir")
    ap.add_argument("--chrom", default="chr1")
    a = ap.parse_args()

    seq = load_sequence(a.seq_dir, a.chrom)
    transcripts, strand = load_transcripts(a.gtf, a.chrom)
    donors, acceptors = set(), set()
    for gene, by_tx in transcripts.items():
        for exons in by_tx.values():
            for d, acc in junctions(exons, strand[gene]):
                donors.add((d, strand[gene]))
                acceptors.add((acc, strand[gene]))
    print(f"{a.chrom}: {len(donors)} distinct annotated donors, {len(acceptors)} acceptors")

    for kind, sites, canon in (("donor", donors, "GT"), ("acceptor", acceptors, "AG")):
        for st in "+-":
            loci = [l for l, s in sites if s == st]
            by_shift = {k: collections.Counter(intron_dinucleotide(seq, l + k, kind, st) for l in loci) for k in range(-2, 3)}
            n = len(loci)
            top = ", ".join(f"{d} {c / n:.3f}" for d, c in by_shift[0].most_common(3))
            shifted = ", ".join(f"{k:+d}: {by_shift[k][canon] / n:.3f}" for k in (-2, -1, 1, 2))
            print(f"  {kind:8s} strand {st}: n={n:6d}  at locus -> {top} | {canon} when shifted {shifted}")

    for width in (80, 400, 2000, 10000):
        ok = total = skipped = 0
        for kind, sites, canon in (("donor", donors, "GT"), ("acceptor", acceptors, "AG")):
            for locus, st in sites:
                w = junction_window(seq, locus, kind, st, width)
                if w is None:
                    skipped += 1
                    continue
                total += 1
                ok += window_dinucleotide(w, kind) == canon
        print(f"  window {width:5d}: GT/AG at the centre for {ok}/{total} ({ok / total:.3f}); skipped {skipped}")


if __name__ == "__main__":
    main()
