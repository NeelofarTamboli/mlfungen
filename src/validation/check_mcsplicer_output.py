"""Check McSplicer output against the annotation and the genome, and recover read support.

For every `s`/`e` row of a McSplicer run (step 0 plus bootstrap spread) this
  * classifies the site from the ORIGINAL GENCODE transcripts of the same gene:
    internal acceptor / internal donor / transcript start (TSS) / transcript end (TES)
    / both / neither (a boundary pulled in from an overlapping gene by exonRefine);
  * reads the intron dinucleotide at the locus (acceptor expects AG, donor GT);
  * recovers read support from the sigcount .cnt file:
      donors:    spliced (red) and contiguous (black) counts, via McSplicer's own function;
      acceptors: counts of subpaths entering the acceptor's subexon from a
                 non-adjacent node (spliced-in) or from the adjacent one (contiguous).
Writes a site table (CSV) and prints a summary. See docs/mcsplicer_notes.md.

Usage:
    python -m src.validation.check_mcsplicer_output \\
        --gtf ORIGINAL.gtf --refined-gtf REFINED.gtf --cnt X.cnt \\
        --mc-out OUT_DIR --seq-dir SEQ_DIR --chrom chr1 --read-len 100 --table sites.csv
"""

from __future__ import annotations

import argparse
import collections
import contextlib
import csv
import glob
import io
import os
import statistics
import sys
from itertools import pairwise
from pathlib import Path

from src.config import REPO_ROOT
from src.validation import check_coordinates as cc


def boundary_sets(transcripts, strand):
    """Per gene: internal acceptor/donor loci and transcript start/end loci, in transcription direction."""
    acc, don, tss, tes = (collections.defaultdict(set) for _ in range(4))
    for gene, by_tx in transcripts.items():
        st = strand[gene]
        for exons in by_tx.values():
            ordered = sorted(exons, reverse=(st == "-"))
            tss[gene].add(ordered[0][0] if st == "+" else ordered[0][1])
            tes[gene].add(ordered[-1][1] if st == "+" else ordered[-1][0])
            for d, a in cc.junctions(exons, st):
                don[gene].add(d)
                acc[gene].add(a)
    return acc, don, tss, tes


def classify(sets, gene: str, index: str, locus: int) -> str:
    acc, don, tss, tes = sets
    if index[0] == "s":
        internal, boundary, names = locus in acc[gene], locus in tss[gene], ("acceptor", "TSS")
    else:
        internal, boundary, names = locus in don[gene], locus in tes[gene], ("donor", "TES")
    if internal and boundary:
        return f"{names[0]}+{names[1]}"
    return names[0] if internal else names[1] if boundary else "neither"


def read_outputs(mc_out: str):
    """{(gene, index, locus, strand): {step: usage}} from <mc_out>/<chrom>/<gene>[_prefix].csv"""
    out = collections.defaultdict(dict)
    for path in glob.glob(os.path.join(mc_out, "*", "*.csv")):
        gene = os.path.basename(path)[:-4]
        with open(path) as fh:
            for r in csv.DictReader(fh):
                out[(gene, r["index"], int(r["splice site locus"]), r["strand"])][int(r["step"])] = float(r["usage estimate"])
    return out


def donor_support(M, cnt: str, refined_gtf: str, read_len: int):
    """{(gene, locus): (spliced, contiguous)} using McSplicer's own junction-read counts."""
    genes, _ = M.get_all_genes_dict(refined_gtf, False)
    orig, found = M.get_q_vector_from_junction_reads, {}
    for gene in genes:
        rows: list = []

        def hook(*args, _rows=rows, **kw):
            _rows.append(orig(*args, **kw))
            return _rows[-1]

        M.get_q_vector_from_junction_reads = hook
        with contextlib.redirect_stdout(io.StringIO()):  # McSplicer prints per-iteration warnings
            M.run_EM_bootstrap(cnt, gene, 0, read_len, genes[gene], True, True)
        for site, _key, _prob, red, black in (rows[0] if rows else []):
            found[(gene, site)] = (red, black)
    M.get_q_vector_from_junction_reads = orig
    return found


def acceptor_support(M, genes, subpaths, freqs, gene: str, locus: int, strand: str):
    """(spliced_in, contiguous_in) counts for the subexon whose transcription-direction start is `locus`."""
    _sd, _l, _s, _e, subexons = M.get_gene_data(gene, genes[gene])
    sub, fr = M.filter_subpath_list_by_subexon_ids(subpaths, freqs, list(subexons))
    nodes = {n for n, (start, _end) in subexons.items() if start == locus}
    spliced = contiguous = 0.0
    for path, f in zip(sub, fr):
        for m, n in pairwise(path):
            if n in nodes:
                gap = subexons[n][0] - subexons[m][1]
                if gap in ((0, 1) if strand == "+" else (0, -1)):
                    contiguous += f
                else:
                    spliced += f
    return spliced, contiguous


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="original (unrefined) GENCODE-style GTF")
    ap.add_argument("--refined-gtf", required=True)
    ap.add_argument("--cnt", required=True)
    ap.add_argument("--mc-out", required=True, help="McSplicer --out_dir")
    ap.add_argument("--seq-dir", required=True)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--mcsplicer-dir", default=str(REPO_ROOT / "external" / "McSplicer"))
    ap.add_argument("--table", help="write the per-site table here (CSV)")
    a = ap.parse_args()

    sys.path.insert(0, str(Path(a.mcsplicer_dir) / "python_scripts"))
    import McSplicer as M  # guarded by __main__ in that file, so importing is safe

    seq = cc.load_sequence(a.seq_dir, a.chrom)
    sets = boundary_sets(*cc.load_transcripts(a.gtf, a.chrom))
    outputs = read_outputs(a.mc_out)
    donors = donor_support(M, a.cnt, a.refined_gtf, a.read_len)
    genes, _ = M.get_all_genes_dict(a.refined_gtf, False)
    subpaths, freqs = M.parse_cnt_file(a.cnt)

    rows = []
    for (gene, index, locus, strand), steps in sorted(outputs.items(), key=lambda x: (x[0][0], x[0][3], x[0][1][0], int(x[0][1][1:]))):
        kind = "acceptor" if index[0] == "s" else "donor"
        boot = [v for k, v in steps.items() if k > 0]
        row = {
            "gene_id": gene, "index": index, "strand": strand, "locus": locus,
            "site_class": classify(sets, gene, index, locus),
            "dinucleotide": cc.intron_dinucleotide(seq, locus, kind, strand), "usage": steps[0],
            "bootstrap_sd": statistics.pstdev(boot) if len(boot) > 1 else "",
            "donor_spliced": "", "donor_contiguous": "", "acc_spliced_in": "", "acc_contiguous_in": "",
        }
        if kind == "donor":
            row["donor_spliced"], row["donor_contiguous"] = donors.get((gene, locus), ("", ""))
        else:
            row["acc_spliced_in"], row["acc_contiguous_in"] = acceptor_support(M, genes, subpaths, freqs, gene, locus, strand)
        rows.append(row)

    if a.table:
        with open(a.table, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    print("rows per (strand, class):", dict(collections.Counter((r["strand"], r["site_class"]) for r in rows)))
    for cls in sorted({r["site_class"] for r in rows}):
        rs = [r for r in rows if r["site_class"] == cls]
        ok = sum(r["dinucleotide"] in ("GT", "AG") for r in rs)
        print(f"  {cls:16s} n={len(rs):3d}  GT/AG at locus: {ok}")
    don = [r for r in rows if r["site_class"] == "donor"]
    acc = [r for r in rows if r["site_class"] == "acceptor"]
    if don:
        print(f"internal donors {len(don)}: usage>=0.999998: {sum(r['usage'] >= 0.999998 for r in don)}, "
              f"<=1.1e-6: {sum(r['usage'] <= 1.1e-6 for r in don)} (of which no reads: "
              f"{sum(r['usage'] <= 1.1e-6 and r['donor_spliced'] == 0 and r['donor_contiguous'] == 0 for r in don)})")
    if acc:
        print(f"internal acceptors {len(acc)}: >=0.99: {sum(r['usage'] >= 0.99 for r in acc)}, "
              f"0.05-0.95: {sum(0.05 < r['usage'] < 0.95 for r in acc)}, <0.05: {sum(r['usage'] <= 0.05 for r in acc)}")


if __name__ == "__main__":
    main()
