"""Run McSplicer on one replicate's BAM, strand-aware (see docs/mcsplicer_notes.md).

For each chromosome (one job each, jobs run in parallel):
    1. split the chromosome's reads by fragment strand (src/validation/split_bam_by_strand.py;
       read 1 reverse or read 2 forward = plus strand transcript);
    2. for each strand, take only that strand's genes from the GENCODE v29 GTF;
    3. exonRefine -> bin/sigcount (no -s) -> McSplicer.py --read_len 100.

McSplicer re-reads the whole .cnt file for every gene, which is far too slow for a whole
chromosome, so genes are run in chunks of `genes_per_chunk` with a .cnt trimmed to the lines
that gene can use (a count line is used by a gene if one of its mates, i.e. one side of
`left^right`, only touches that gene's nodes - the same test McSplicer applies itself).
Results per gene are the same as in one big run, except that bootstrap steps 1..N use the
random numbers of their own chunk (step 0 is deterministic).

Output layout (nothing is deleted except the temporary split BAMs):
    <out>/runs.jsonl                one line per invocation (settings, McSplicer commit)
    <out>/<chrom>/<plus|minus>/refined.gtf, counts.cnt.gz, annotation.gtf
    <out>/<chrom>/<plus|minus>/chunk_NNN/{refined.gtf, counts.cnt, mcsplicer.log, out/<chrom>/<gene>.csv}
    <out>/<chrom>/done.json

Usage:
    python -m src.run_mcsplicer --bam data/raw/ENCFF864OUS.bam --out data/interim/mcsplicer/rep1
    python -m src.run_mcsplicer --bam region.bam --region chr1:111420000-111480000 --out DIR
The BAM must be coordinate-sorted and indexed (ENCODE ships no .bai: `pysam.index(bam)`).
McSplicer is not in the repo: git clone https://github.com/canzarlab/McSplicer external/McSplicer
(tested at commit d6a22c3). The bundled binaries are Linux x86_64; bin_osx is untested.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from src.config import REPO_ROOT, load_config
from src.validation.split_bam_by_strand import split

STRAND_NAME = {"+": "plus", "-": "minus"}


def mcsplicer_tools(mc_dir: Path, sigcount: str = "sigcount") -> dict[str, Path]:
    bins = mc_dir / ("bin_osx" if platform.system() == "Darwin" else "bin")
    tools = {"exonRefine": bins / "exonRefine", "sigcount": bins / sigcount, "McSplicer": mc_dir / "python_scripts" / "McSplicer.py"}
    missing = [str(p) for p in tools.values() if not p.exists()]
    if missing:
        sys.exit(f"missing McSplicer files {missing}\nclone it: git clone https://github.com/canzarlab/McSplicer {mc_dir}")
    return tools


def run(cmd: list, log=None, cwd: Path | None = None) -> None:
    """Run a command; stdout/stderr go to `log` (an open file) if given. Raises on failure."""
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT, cwd=cwd, env=env, check=True)


def read_gtf_chrom(gtf: Path, chrom: str, region: tuple[int, int] | None):
    """GTF lines of one chromosome grouped by strand. With a region, keep only genes fully inside it."""
    lines = {"+": [], "-": []}
    span: dict[str, list[int]] = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.split("\t", 8)
            gene = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            s = span.setdefault(gene, [int(f[3]), int(f[4])])
            s[0], s[1] = min(s[0], int(f[3])), max(s[1], int(f[4]))
            lines[f[6]].append((gene, line))

    def inside(g: str) -> bool:
        return region is None or (span[g][0] >= region[0] and span[g][1] <= region[1])

    return {st: [line for g, line in ls if inside(g)] for st, ls in lines.items()}


def _nodes(side: str) -> frozenset:
    return frozenset(int(x) for x in side.split("-"))


def make_chunks(refined_gtf: Path, cnt: Path, out_dir: Path, per_chunk: int) -> list[Path]:
    """Split a refined GTF + .cnt into chunk directories of `per_chunk` genes that have reads."""
    subexons = defaultdict(list)  # gene -> subexon lines
    gene_ids = defaultdict(set)  # gene -> NodeIds
    gene_start: dict[str, int] = {}
    with open(refined_gtf) as fh:
        for line in fh:
            f = line.split("\t", 8)
            if f[2] != "subexon":
                continue
            gene = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            subexons[gene].append(line)
            gene_ids[gene].add(int(re.search(r"(?<![A-Za-z])NodeId (\d+)", f[8]).group(1)))
            gene_start[gene] = min(gene_start.get(gene, int(f[3])), int(f[3]))
    by_node = defaultdict(list)
    for g, ids in gene_ids.items():
        for n in ids:
            by_node[n].append(g)

    header, data = [], []
    with open(cnt) as fh:
        for line in fh:
            (header if line.startswith("#") else data).append(line)
    lines_of_gene = defaultdict(list)
    for i, line in enumerate(data):
        used = set()
        for side in line.split("\t", 1)[0].split("^"):
            nodes = _nodes(side)
            used.update(g for g in by_node[next(iter(nodes))] if nodes <= gene_ids[g])
        for g in used:
            lines_of_gene[g].append(i)

    genes = sorted(lines_of_gene, key=lambda g: (gene_start[g], g))  # genes with reads, in genome order
    chunks = []
    for k in range(0, len(genes), per_chunk):
        group = genes[k : k + per_chunk]
        d = out_dir / f"chunk_{k // per_chunk:03d}"
        d.mkdir(parents=True, exist_ok=True)
        keep = sorted({i for g in group for i in lines_of_gene[g]})
        (d / "refined.gtf").write_text("".join(line for g in group for line in subexons[g]))
        (d / "counts.cnt").write_text("".join(header) + "".join(data[i] for i in keep))
        chunks.append(d)
    return chunks


def chunk_genes(refined_gtf: Path) -> dict[str, list[str]]:
    """gene -> its subexon lines, in file order (McSplicer.py runs genes in this order)."""
    genes: dict[str, list[str]] = {}
    with open(refined_gtf) as fh:
        for line in fh:
            if line.split("\t", 3)[2] == "subexon":
                genes.setdefault(re.search(r'gene_id "([^"]+)"', line).group(1), []).append(line)
    return genes


def run_chunk(d: Path, tools: dict, a: dict, chrom: str, strand: str) -> dict:
    """Run McSplicer.py on one chunk; a crash on one gene does not lose the others.

    McSplicer.py has no error handling: one gene that raises kills the whole process, and with it every
    gene after it in the chunk. Here the failing gene is read off the log (it is the last "Running McSplicer
    on gene id" line before the traceback), recorded, and McSplicer is started again on the genes that are
    neither finished nor failed. Failures go to <chunk>/failed_genes.json (also printed to the run log);
    nothing is dropped silently.

    Its console output is mostly 'Incompatible read' warnings (one per read path and EM iteration). They are
    counted per gene into warnings.json instead of being stored; every other line goes to mcsplicer.log.

    Resumable: a chunk with warnings.json (written last) is skipped; an interrupted chunk keeps its finished
    genes (their CSVs) and its failure list in progress.json and only runs the rest.
    """
    if (d / "warnings.json").exists():
        return {"skipped": True}
    progress = json.loads((d / "progress.json").read_text()) if (d / "progress.json").exists() else {"warnings": {}, "failed": [], "no_output": []}
    genes = chunk_genes(d / "refined.gtf")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    while True:
        finished = {p.stem for p in (d / "out").glob("*/*.csv")}
        skip = finished | {f["gene_id"] for f in progress["failed"]} | set(progress["no_output"])
        todo = [g for g in genes if g not in skip]
        if not todo:
            break
        gtf = d / "refined.gtf" if len(todo) == len(genes) else d / "todo.gtf"
        if gtf.name == "todo.gtf":
            gtf.write_text("".join(line for g in todo for line in genes[g]))
        cmd = [sys.executable, "-u", tools["McSplicer"], "--gtf", gtf, "--count_file", d / "counts.cnt", "--out_dir", d / "out",
               "--read_len", a["read_len"], "--bootstraps", a["bootstraps"]]
        gene, tail, traceback_seen = None, [], False
        with open(d / "mcsplicer.log", "a") as log:
            proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
            for line in proc.stdout:
                m = re.match(r"Running McSplicer on gene id (\S+)", line)
                gene = m.group(1) if m else gene
                if "Incompatible read" in line:
                    progress["warnings"][gene] = progress["warnings"].get(gene, 0) + 1
                    continue
                if line.startswith("Warning: file"):
                    continue
                if (m2 := re.search(r"No counts are found for gene (\S+)", line)):
                    progress["no_output"].append(m2.group(1))
                traceback_seen = traceback_seen or line.startswith("Traceback")
                tail = (tail + [line.rstrip()])[-12:]
                log.write(line)
            rc = proc.wait()
        if rc == 0:
            break
        if not traceback_seen or gene is None:  # killed or broken environment, not a gene problem: stop here
            raise RuntimeError(f"McSplicer.py exited with {rc} without a Python traceback in {d}; last lines: {tail[-3:]}")
        failure = {"gene_id": gene, "chrom": chrom, "strand": strand, "chunk": d.name, "error": tail[-1], "traceback": tail}
        progress["failed"].append(failure)
        print(f"FAILED gene {gene} ({chrom} {strand}, {d.name}): {tail[-1]}", file=sys.stderr, flush=True)
        (d / "progress.json").write_text(json.dumps(progress))
    (d / "failed_genes.json").write_text(json.dumps(progress["failed"], indent=1))
    (d / "warnings.json").write_text(json.dumps(progress["warnings"]))  # last: marks the chunk as finished
    (d / "progress.json").unlink(missing_ok=True)
    return {"failed": len(progress["failed"]), "genes": len(genes)}


def run_strand(chrom: str, strand: str, gtf_lines: list[str], bam: Path, out: Path, tools: dict, a: dict) -> dict:
    sdir = out / STRAND_NAME[strand]
    sdir.mkdir(parents=True, exist_ok=True)
    if (sdir / "strand.json").exists():
        return json.loads((sdir / "strand.json").read_text())
    t = {}
    (sdir / "annotation.gtf").write_text("".join(gtf_lines))
    t0 = time.time()
    with open(sdir / "exonRefine.log", "w") as log:
        run([tools["exonRefine"], sdir / "annotation.gtf", "-p", sdir / "refined"], log)
    t["exonRefine"] = time.time() - t0
    t0 = time.time()
    with open(sdir / "sigcount.log", "w") as log:
        run([tools["sigcount"], bam, sdir / "refined.gtf", sdir / "counts"], log)
    t["sigcount"] = time.time() - t0
    t0 = time.time()
    chunks = make_chunks(sdir / "refined.gtf", sdir / "counts.cnt", sdir, a["genes_per_chunk"])
    t["chunking"] = time.time() - t0
    t0 = time.time()
    results = [run_chunk(d, tools, a, chrom, strand) for d in chunks]
    t["mcsplicer"] = time.time() - t0
    with open(sdir / "counts.cnt", "rb") as src, gzip.open(sdir / "counts.cnt.gz", "wb") as dst:
        dst.write(src.read())
    (sdir / "counts.cnt").unlink()
    info = {"chunks": len(chunks), "failed_genes": sum(r.get("failed", 0) for r in results), "seconds": {k: round(v, 1) for k, v in t.items()}}
    (sdir / "strand.json").write_text(json.dumps(info))
    return info


def run_chrom(job: dict) -> dict:
    chrom, out = job["chrom"], Path(job["out"]) / job["chrom"]
    if job["force"] and out.exists():
        shutil.rmtree(out)
    if (out / "done.json").exists():
        return {"chrom": chrom, "skipped": True}
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    region = job["region"]
    gtf = read_gtf_chrom(Path(job["gtf"]), chrom, region)
    has_genes = {s: any("\texon\t" in line for line in gtf[s]) for s in "+-"}
    pending = [s for s in "+-" if has_genes[s] and not (out / STRAND_NAME[s] / "strand.json").exists()]
    n = split(job["bam"], str(out / "reads"), chrom, region[0] - 1 if region else None, region[1] if region else None) if pending else "from an earlier run"
    info = {"chrom": chrom, "bootstraps": job["bootstraps"], "read_len": job["read_len"], "alignments": n,
            "seconds_split_and_gtf": round(time.time() - t0, 1), "strands": {}}
    tools = {k: Path(v) for k, v in job["tools"].items()}
    for strand in "+-":
        if not has_genes[strand]:
            info["strands"][strand] = "no genes"
            continue
        info["strands"][strand] = run_strand(chrom, strand, gtf[strand], out / f"reads_{STRAND_NAME[strand]}.bam", out, tools, job)
    if not job["keep_bams"]:
        for p in out.glob("reads_*.bam*"):
            p.unlink()
    info["failed_genes"] = sum(v["failed_genes"] for v in info["strands"].values() if isinstance(v, dict))
    info["seconds_total"] = round(time.time() - t0, 1)
    (out / "done.json").write_text(json.dumps(info, indent=1))
    return info


def collect_failed_genes(out: Path) -> list[dict]:
    """All gene-level failures recorded under `out` (from this and earlier invocations)."""
    return [f for p in sorted(out.glob("*/*/chunk_*/failed_genes.json")) for f in json.loads(p.read_text())]


def main() -> None:
    cfg = load_config()
    mc = cfg["mcsplicer"]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True, help="output directory for this replicate")
    ap.add_argument("--gtf", default=str(REPO_ROOT / cfg["reference"]["annotation_gtf"]))
    ap.add_argument("--chrom", nargs="+", help="chromosomes to run (default: config mcsplicer.chromosomes)")
    ap.add_argument("--region", help="chrN:START-END (1-based); only genes fully inside are used")
    ap.add_argument("--bootstraps", type=int, default=mc["bootstraps"])
    ap.add_argument("--read-len", type=int, default=mc["read_len"])
    ap.add_argument("--genes-per-chunk", type=int, default=mc["genes_per_chunk"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--keep-bams", action="store_true")
    ap.add_argument("--force", action="store_true", help="delete and redo the chromosomes asked for (default: resume where they stopped)")
    a = ap.parse_args()

    region, chroms = None, a.chrom or mc["chromosomes"]
    if a.region:
        m = re.fullmatch(r"(\w+):(\d+)-(\d+)", a.region)
        chroms, region = [m.group(1)], (int(m.group(2)), int(m.group(3)))
    mc_dir = REPO_ROOT / mc["dir"]
    tools = mcsplicer_tools(mc_dir, mc["sigcount_binary"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(["git", "-C", str(mc_dir), "rev-parse", "HEAD"], capture_output=True, text=True, check=False).stdout.strip()
    with open(out / "runs.jsonl", "a") as fh:  # one line per invocation; per-chromosome settings are also in <chrom>/done.json
        fh.write(json.dumps({
            "bam": a.bam, "gtf": a.gtf, "region": a.region, "read_len": a.read_len, "bootstraps": a.bootstraps,
            "genes_per_chunk": a.genes_per_chunk, "mcsplicer_commit": commit, "chromosomes": chroms,
            "strand_rule": mc["strand_rule"], "started": time.strftime("%Y-%m-%d %H:%M:%S"),
        }) + "\n")

    jobs = [{"chrom": c, "out": str(out), "bam": a.bam, "gtf": a.gtf, "region": region, "tools": {k: str(v) for k, v in tools.items()},
             "read_len": a.read_len, "bootstraps": a.bootstraps, "genes_per_chunk": a.genes_per_chunk,
             "keep_bams": a.keep_bams, "force": a.force} for c in chroms]
    t0, bad_chroms = time.time(), []
    with cf.ProcessPoolExecutor(max_workers=min(a.jobs, len(jobs))) as pool:
        futs = {pool.submit(run_chrom, j): j["chrom"] for j in jobs}
        for fut in cf.as_completed(futs):
            try:
                r = fut.result()
            except Exception as e:  # noqa: BLE001  one failed chromosome must not stop the others; it resumes on the next run
                bad_chroms.append(futs[fut])
                (out / futs[fut]).mkdir(parents=True, exist_ok=True)
                (out / futs[fut] / "failed.txt").write_text(f"{type(e).__name__}: {e}\n")
                print(f"[{time.time() - t0:7.0f}s] {futs[fut]}: CHROMOSOME FAILED, {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                continue
            print(f"[{time.time() - t0:7.0f}s] {futs[fut]}: " + ("skipped (done)" if r.get("skipped") else
                  f"{r['seconds_total']} s, {r['failed_genes']} failed genes, strands {r['strands']}"), flush=True)
    failed = collect_failed_genes(out)
    (out / "failed_genes.json").write_text(json.dumps(failed, indent=1))
    print(f"failed genes under {out}: {len(failed)}" + "".join(f"\n  {f['gene_id']} {f['chrom']}{f['strand']} {f['chunk']}: {f['error']}" for f in failed), flush=True)
    if bad_chroms:
        sys.exit(f"chromosomes that failed (rerun to resume): {bad_chroms}")


if __name__ == "__main__":
    main()
