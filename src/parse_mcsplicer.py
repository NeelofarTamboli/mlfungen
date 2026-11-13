"""Turn the raw McSplicer output of one replicate into per-site tables.

Input: the directory written by src/run_mcsplicer.py for one replicate. Output:
    <rep_dir>/sites_all.tsv      every s/e row McSplicer reported, with its class (kept for QC/reporting)
    data/processed/acceptors_<rep>.tsv   only internal acceptors (the regression target)

Acceptors that are also the first-exon start of another transcript (class `acceptor+TSS`) stay in the
acceptor table with `also_transcript_start = True`: their McSplicer value mixes a transcription start and a
splice acceptor, so the primary QC leaves them out, but they are not removed from the data.

Columns (one row per site, all from the verified logic in docs/mcsplicer_notes.md):
    site_id        chrom:locus:strand:gene_id
    gene_id, chrom, locus (1-based exonic base), strand
    index          McSplicer's s<k>/e<k>
    site_class     acceptor / donor / TSS / TES / acceptor+TSS / donor+TES / neither (in the acceptor table: acceptor)
    also_transcript_start   True for an acceptor+TSS site (internal acceptor in one transcript, first exon start in another)
                   (from the original GENCODE transcripts: an `s` row is an acceptor if the locus
                   starts a non-first exon of some transcript of the gene; `e` likewise for donors)
    usage          McSplicer step-0 estimate
    bootstrap_sd   sample SD (n-1) of the bootstrap steps 1..N; empty if there are fewer than 2
    n_bootstraps
    spliced_in, contiguous_in   (acceptors) .cnt counts entering the acceptor's subexon from a
                   non-adjacent / the adjacent subexon. Counts are per mate (a read pair adds to
                   both mates) and include fractional weights for multi-mapped reads.
    read_support   spliced_in + contiguous_in = reads that cross the acceptor's left boundary
    n_incompat_warnings   McSplicer "Incompatible read" warnings for this gene (a gene-level flag)
    raw_csv, chunk        where the row came from, relative to the replicate directory

Usage:
    python -m src.parse_mcsplicer --rep-dir data/interim/mcsplicer/rep1 [--jobs 4]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import statistics
import sys
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

from src.config import REPO_ROOT, load_config, project_path
from src.validation import check_coordinates as cc
from src.validation.check_mcsplicer_output import boundary_sets, classify

COLUMNS = ["site_id", "gene_id", "chrom", "locus", "strand", "index", "site_class", "also_transcript_start", "usage", "bootstrap_sd", "n_bootstraps",
           "spliced_in", "contiguous_in", "read_support", "n_incompat_warnings", "raw_csv", "chunk"]


def import_mcsplicer():
    sys.path.insert(0, str(REPO_ROOT / load_config()["mcsplicer"]["dir"] / "python_scripts"))
    import McSplicer as M  # guarded by __main__ in that file, so importing is safe

    return M


def acceptor_support(M, genes, subpaths, freqs, gene: str, strand: str) -> dict[int, tuple[float, float]]:
    """{acceptor locus: (spliced_in, contiguous_in)} for one gene, in one pass over its subpaths.

    Same definition as check_mcsplicer_output.acceptor_support: every consecutive node pair (m, n) in a
    subpath enters the subexon n; it is contiguous if m ends right before n starts, else spliced.
    """
    _strand, _loc, _starts, _ends, subexons = M.get_gene_data(gene, genes[gene])
    sub, fr = M.filter_subpath_list_by_subexon_ids(subpaths, freqs, list(subexons))
    adjacent = (0, 1) if strand == "+" else (0, -1)
    out: dict[int, list[float]] = defaultdict(lambda: [0.0, 0.0])
    for path, f in zip(sub, fr):
        for m, n in pairwise(path):
            gap = subexons[n][0] - subexons[m][1]
            out[subexons[n][0]][1 if gap in adjacent else 0] += f
    return {locus: (v[0], v[1]) for locus, v in out.items()}


def parse_chrom(args: tuple[str, str, str]) -> list[dict]:
    rep_dir, chrom, gtf = Path(args[0]), args[1], args[2]
    M = import_mcsplicer()
    transcripts, gene_strand = cc.load_transcripts(gtf, chrom)
    sets = boundary_sets(transcripts, gene_strand)
    rows = []
    for sdir in sorted((rep_dir / chrom).glob("*")):
        if not sdir.is_dir():
            continue
        for chunk in sorted(sdir.glob("chunk_*")):
            genes, _ = M.get_all_genes_dict(str(chunk / "refined.gtf"), False)
            subpaths, freqs = M.parse_cnt_file(str(chunk / "counts.cnt"))
            warn = json.loads((chunk / "warnings.json").read_text())
            for path in sorted((chunk / "out").glob("*/*.csv")):
                gene = path.stem
                steps = defaultdict(dict)  # (index, locus, strand) -> {step: usage}
                with open(path) as fh:
                    for r in csv.DictReader(fh):
                        steps[(r["index"], int(r["splice site locus"]), r["strand"])][int(r["step"])] = float(r["usage estimate"])
                strand = next(iter(steps))[2]
                support = acceptor_support(M, genes, subpaths, freqs, gene, strand)
                for (index, locus, strand), by_step in steps.items():
                    boot = [v for k, v in by_step.items() if k > 0]
                    row = {
                        "site_id": f"{chrom}:{locus}:{strand}:{gene}", "gene_id": gene, "chrom": chrom, "locus": locus, "strand": strand,
                        "index": index, "site_class": (cls := classify(sets, gene, index, locus)), "also_transcript_start": cls == "acceptor+TSS", "usage": by_step[0],
                        "bootstrap_sd": statistics.stdev(boot) if len(boot) > 1 else "", "n_bootstraps": len(boot),
                        "n_incompat_warnings": warn.get(gene, 0), "raw_csv": str(path.relative_to(rep_dir)), "chunk": str(chunk.relative_to(rep_dir)),
                        "spliced_in": "", "contiguous_in": "", "read_support": "",
                    }
                    if index[0] == "s":
                        sp, co = support.get(locus, (0.0, 0.0))
                        row.update(spliced_in=sp, contiguous_in=co, read_support=sp + co)
                    rows.append(row)
    return rows


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep-dir", required=True)
    ap.add_argument("--gtf", default=str(REPO_ROOT / cfg["reference"]["annotation_gtf"]), help="original GENCODE GTF")
    ap.add_argument("--rep-name", help="name used in the output file (default: the directory name)")
    ap.add_argument("--out-dir", default=str(project_path(cfg, "data_processed")))
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    rep_dir = Path(a.rep_dir)
    name = a.rep_name or rep_dir.name
    chroms = sorted(p.name for p in rep_dir.iterdir() if p.is_dir() and (p / "done.json").exists())
    with cf.ProcessPoolExecutor(max_workers=min(a.jobs, len(chroms))) as pool:
        rows = [r for part in pool.map(parse_chrom, [(str(rep_dir), c, a.gtf) for c in chroms]) for r in part]
    rows.sort(key=lambda r: (r["chrom"], r["locus"], r["strand"], r["gene_id"]))
    acceptors = [dict(r, site_class="acceptor") for r in rows if r["site_class"] in ("acceptor", "acceptor+TSS")]
    write_tsv(rows, rep_dir / "sites_all.tsv")
    write_tsv(acceptors, Path(a.out_dir) / f"acceptors_{name}.tsv")
    by_class = defaultdict(int)
    for r in rows:
        by_class[r["site_class"]] += 1
    print(f"{name}: {len(rows)} s/e rows from {len(chroms)} chromosomes: {dict(by_class)}")
    flagged = sum(r["also_transcript_start"] for r in acceptors)
    print(f"  wrote {rep_dir / 'sites_all.tsv'} and {Path(a.out_dir) / f'acceptors_{name}.tsv'} ({len(acceptors)} internal acceptors, {flagged} flagged also_transcript_start)")
    failed = rep_dir / "failed_genes.json"
    n_failed = len(json.loads(failed.read_text())) if failed.exists() else 0
    print(f"  genes that crashed McSplicer and have no rows: {n_failed}" + (f" (see {failed})" if n_failed else ""))


if __name__ == "__main__":
    main()
