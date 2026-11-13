"""Download the project inputs into data/raw/ (resumable, md5-checked).

    python -m src.download_data gtf          # GENCODE v29 annotation (43 MB)
    python -m src.download_data bam rep1     # ENCFF864OUS (10.4 GB)
    python -m src.download_data bam rep2     # ENCFF696CSX (13.4 GB)
    python -m src.download_data fasta        # GRCh38 no_alt analysis set (873 MB)

Accessions and md5s come from config.yaml. The BAMs are large: check disk first.
data/ is git-ignored, so nothing here is committed.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import subprocess
import time
import urllib.request
from pathlib import Path

from src.config import load_config, project_path

ENCODE = "https://www.encodeproject.org/files"


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, md5: str | None) -> Path:
    """Download url to dest (resuming a .part file), then check md5 if given."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and (md5 is None or md5sum(dest) == md5):
        print(f"already have {dest}")
        return dest
    part = dest.with_name(dest.name + ".part")
    subprocess.run(["curl", "-fL", "-C", "-", "--retry", "5", "--retry-delay", "5", "-o", str(part), url], check=True)
    got = md5sum(part)
    if md5 and got != md5:
        raise RuntimeError(f"md5 mismatch for {dest.name}: got {got}, expected {md5}")
    part.rename(dest)
    print(f"downloaded {dest} (md5 {got})")
    return dest


def cloud_url(accession: str) -> str:
    """Public S3 URL of an ENCODE file (supports Range requests, no signature needed)."""
    req = urllib.request.Request(f"{ENCODE}/{accession}/?format=json", headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.load(r)["cloud_metadata"]["url"]


def fetch_ranges(url: str, dest: Path, md5: str, workers: int = 8, chunk: int = 64 << 20) -> Path:
    """Download a big file as parallel byte ranges (resumable), then check md5.

    One connection to S3 is sometimes much slower than the others, so a single curl can crawl at
    2 MB/s while several ranges together run at 50 MB/s. A chunk that stays below 1 MB/s for 15 s is
    dropped and requested again. Progress is kept in <dest>.part.done so a restart skips finished chunks.
    """
    if dest.exists() and md5sum(dest) == md5:
        print(f"already have {dest}")
        return dest
    head = subprocess.run(["curl", "-sSI", url], capture_output=True, text=True, check=True).stdout
    size = int(next(ln.split(":")[1] for ln in head.splitlines() if ln.lower().startswith("content-length")))
    part, state = dest.with_name(dest.name + ".part"), dest.with_name(dest.name + ".part.done")
    n = -(-size // chunk)
    done = set(json.loads(state.read_text())) if state.exists() else set()
    if part.exists() and not state.exists():  # prefix written by a plain `curl -C -`: contiguous from byte 0
        done = set(range(part.stat().st_size // chunk))
    fd = os.open(part, os.O_RDWR | os.O_CREAT)
    os.ftruncate(fd, size)

    def get(i: int) -> int:
        a, b = i * chunk, min((i + 1) * chunk, size) - 1
        for attempt in range(10):
            r = subprocess.run(["curl", "-sS", "--speed-limit", "1000000", "--speed-time", "15", "-r", f"{a}-{b}", url], capture_output=True, check=False)
            if r.returncode == 0 and len(r.stdout) == b - a + 1:
                os.pwrite(fd, r.stdout, a)
                return i
            time.sleep(2 + attempt)
        raise RuntimeError(f"chunk {i} failed repeatedly")

    t0, finished = time.time(), len(done)
    print(f"{dest.name}: {size / 1e9:.1f} GB in {n} chunks, {len(done)} already done", flush=True)
    with cf.ThreadPoolExecutor(workers) as pool:
        for i in pool.map(get, [i for i in range(n) if i not in done]):
            done.add(i)
            finished += 1
            state.write_text(json.dumps(sorted(done)))
            if finished % 16 == 0 or finished == n:
                print(f"  {finished}/{n} chunks, {time.time() - t0:.0f} s", flush=True)
    os.close(fd)
    got = md5sum(part)
    if got != md5:
        raise RuntimeError(f"md5 mismatch for {dest.name}: got {got}, expected {md5}")
    part.rename(dest)
    state.unlink()
    print(f"downloaded {dest} (md5 {got})")
    return dest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["gtf", "bam", "fasta"])
    ap.add_argument("rep", nargs="?", choices=["rep1", "rep2"])
    a = ap.parse_args()
    cfg = load_config()
    raw = project_path(cfg, "data_raw")
    ref = cfg["reference"]
    if a.what == "gtf":
        name = Path(ref["annotation_gtf"]).name
        fetch(f"{ENCODE}/{ref['annotation_id']}/@@download/{name}", raw / name, ref["annotation_md5"])
    elif a.what == "fasta":
        name = Path(ref["genome_fasta"]).name
        fetch(f"{ENCODE}/{ref['genome_fasta_id']}/@@download/{name}", raw / name, ref["genome_fasta_md5"])
    else:
        if not a.rep:
            ap.error("bam needs rep1 or rep2")
        rep = cfg["dataset"]["replicates"][a.rep]
        acc = rep["encode_file_accession"]
        fetch_ranges(cloud_url(acc), raw / f"{acc}.bam", rep["md5"])


if __name__ == "__main__":
    main()
