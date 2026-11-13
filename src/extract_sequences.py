"""Prepare GRCh38 contigs and extract strand-oriented acceptor windows.

The compact on-disk representation stores only the largest requested window as
``uint8`` nucleotide codes (A=0, C=1, G=2, T=3, N/other=4). Shorter contexts
are centred views of that array, and one-hot encoding is produced lazily for a
batch. This avoids writing four redundant sequence datasets or a very large
float one-hot tensor.

For the McSplicer acceptor locus L, the 0-based junction index is L-1 on the
plus strand and L on the minus strand. A window is cut symmetrically around
that junction and reverse-complemented on the minus strand. Consequently the
canonical acceptor ``AG`` is always at ``window[half-2:half]``.

Outputs for ``--out-prefix PREFIX``:
    PREFIX.sequences.npy       uint8 array, shape (sites, largest_context)
    PREFIX.metadata.tsv.gz     input rows plus extraction QC and row index
    PREFIX.manifest.json       contexts, codebook, shapes and source paths

Example:
    python -m src.extract_sequences \
      --sites results/person_b_acceptors_chr1-4.tsv \
      --prepare-reference
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import BinaryIO, Iterable

import numpy as np
import pandas as pd

from src.config import REPO_ROOT, load_config, project_path

BASE_TO_CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
CODE_TO_BASE = np.array(list("ACGTN"), dtype="U1")
_ENCODE_LUT = np.full(256, 4, dtype=np.uint8)
for _base, _code in BASE_TO_CODE.items():
    _ENCODE_LUT[ord(_base)] = _code

# Preserve IUPAC complements while orienting the window; all non-ACGT symbols
# are subsequently encoded as N so models never receive a guessed nucleotide.
_IUPAC = b"ACGTRYMKWSBDHVN"
_IUPAC_COMP = b"TGCAYRKMWSVHDBN"
_REVCOMP_TABLE = bytes.maketrans(_IUPAC, _IUPAC_COMP)


def portable_path(path: str | Path) -> str:
    """Use a repository-relative path in reports when possible."""
    path = Path(path)
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def junction_index(locus: int, strand: str) -> int:
    """Return the 0-based first base to the genomic right of an acceptor junction."""
    if strand == "+":
        return locus - 1
    if strand == "-":
        return locus
    raise ValueError(f"strand must be '+' or '-', got {strand!r}")


def reverse_complement(sequence: bytes | str) -> bytes | str:
    """Reverse-complement a DNA sequence, including standard IUPAC symbols."""
    if isinstance(sequence, bytes):
        return sequence.upper().translate(_REVCOMP_TABLE)[::-1]
    table = str.maketrans(_IUPAC.decode(), _IUPAC_COMP.decode())
    return sequence.upper().translate(table)[::-1]


def acceptor_window(sequence: bytes | str, locus: int, strand: str, width: int) -> bytes | str:
    """Extract an oriented, N-padded acceptor window from a chromosome sequence."""
    if width <= 0 or width % 2:
        raise ValueError(f"width must be a positive even integer, got {width}")
    j = junction_index(int(locus), strand)
    start, end = j - width // 2, j + width // 2
    left_pad, right_pad = max(0, -start), max(0, end - len(sequence))
    n = b"N" if isinstance(sequence, bytes) else "N"
    window = n * left_pad + sequence[max(0, start) : min(len(sequence), end)] + n * right_pad
    if strand == "-":
        window = reverse_complement(window)
    if len(window) != width:
        raise AssertionError(f"extracted {len(window)} bases, expected {width}")
    return window


def encode_sequence(sequence: bytes | str) -> np.ndarray:
    """Encode A/C/G/T as 0..3 and N or any ambiguous IUPAC base as 4."""
    raw = sequence.encode() if isinstance(sequence, str) else sequence
    return _ENCODE_LUT[np.frombuffer(raw.upper(), dtype=np.uint8)]


def decode_codes(codes: np.ndarray) -> np.ndarray | str:
    """Decode one row to a string, or a matrix to an equally shaped string array."""
    decoded = CODE_TO_BASE[np.asarray(codes)]
    if decoded.ndim == 1:
        return "".join(decoded.tolist())
    return decoded


def one_hot_encode(codes: np.ndarray, dtype: np.dtype | type = np.uint8) -> np.ndarray:
    """Lazily one-hot encode codes; N/ambiguous bases become four zeros."""
    codes = np.asarray(codes)
    out = np.zeros(codes.shape + (4,), dtype=dtype)
    valid = codes < 4
    positions = np.nonzero(valid)
    out[positions + (codes[valid],)] = 1
    return out


def load_context(
    sequence_path: str | Path,
    width: int,
    rows: slice | Iterable[int] | None = None,
) -> np.ndarray:
    """Load a centred context from the stored maximum-width sequence array."""
    store = np.load(sequence_path, mmap_mode="r")
    if width <= 0 or width % 2 or width > store.shape[1]:
        raise ValueError(f"width must be positive, even, and <= {store.shape[1]}; got {width}")
    start = (store.shape[1] - width) // 2
    selected = store if rows is None else store[rows]
    return selected[..., start : start + width]


def _open_fasta(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def prepare_reference(fasta: str | Path, out_dir: str | Path, chromosomes: Iterable[str]) -> dict[str, int]:
    """Stream a FASTA once and write selected chromosomes as headerless ``.seq`` files."""
    fasta, out_dir = Path(fasta), Path(out_dir)
    wanted = set(chromosomes)
    if not fasta.exists():
        raise FileNotFoundError(f"reference FASTA not found: {fasta}")
    out_dir.mkdir(parents=True, exist_ok=True)
    handles: dict[str, BinaryIO] = {}
    lengths = {chrom: 0 for chrom in wanted}
    current: str | None = None
    try:
        with _open_fasta(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    if current in wanted:
                        handles[current] = open(out_dir / f"{current}.seq", "wb")
                    continue
                if current in wanted:
                    block = line.strip().upper().encode()
                    handles[current].write(block)
                    lengths[current] += len(block)
    finally:
        for fh in handles.values():
            fh.close()
    missing = sorted(chrom for chrom, length in lengths.items() if length == 0)
    if missing:
        raise ValueError(f"FASTA did not contain requested contigs: {missing}")
    (out_dir / "contigs.json").write_text(json.dumps(lengths, indent=2) + "\n")
    return lengths


def extract_sequences(
    sites: pd.DataFrame,
    sequence_dir: str | Path,
    out_prefix: str | Path,
    context_lengths: Iterable[int],
) -> tuple[Path, Path, Path, dict]:
    """Extract all sites into a compact maximum-context store plus QC metadata."""
    required = {"site_id", "chrom", "locus", "strand"}
    missing = required - set(sites.columns)
    if missing:
        raise ValueError(f"site table is missing columns: {sorted(missing)}")
    widths = sorted(set(int(width) for width in context_lengths))
    if not widths or any(width <= 0 or width % 2 for width in widths):
        raise ValueError(f"all context lengths must be positive and even: {widths}")

    sequence_dir, out_prefix = Path(sequence_dir), Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    sequence_path = Path(f"{out_prefix}.sequences.npy")
    metadata_path = Path(f"{out_prefix}.metadata.tsv.gz")
    manifest_path = Path(f"{out_prefix}.manifest.json")
    max_width = max(widths)
    store = np.lib.format.open_memmap(sequence_path, mode="w+", dtype=np.uint8, shape=(len(sites), max_width))

    metadata = sites.reset_index(drop=True).copy()
    metadata["sequence_row"] = np.arange(len(metadata), dtype=np.int64)
    metadata["junction_index_0based"] = -1
    metadata["genomic_start_0based"] = -1
    metadata["genomic_end_0based_exclusive"] = -1
    metadata["upstream_n_padding"] = 0
    metadata["downstream_n_padding"] = 0
    metadata["n_non_acgt_in_max_context"] = 0
    metadata["acceptor_23mer"] = ""
    metadata["center_dinucleotide"] = ""
    metadata["canonical_acceptor"] = False

    for chrom, group in metadata.groupby("chrom", sort=False):
        chrom_path = sequence_dir / f"{chrom}.seq"
        if not chrom_path.exists():
            raise FileNotFoundError(f"prepared chromosome sequence not found: {chrom_path}")
        chromosome = chrom_path.read_bytes().strip().upper()
        chrom_len = len(chromosome)
        for row in group.itertuples():
            j = junction_index(int(row.locus), row.strand)
            raw_start, raw_end = j - max_width // 2, j + max_width // 2
            left_pad, right_pad = max(0, -raw_start), max(0, raw_end - chrom_len)
            window = acceptor_window(chromosome, int(row.locus), row.strand, max_width)
            codes = encode_sequence(window)
            store[row.Index] = codes
            half = max_width // 2
            upstream_pad, downstream_pad = (left_pad, right_pad) if row.strand == "+" else (right_pad, left_pad)
            metadata.at[row.Index, "junction_index_0based"] = j
            metadata.at[row.Index, "genomic_start_0based"] = max(0, raw_start)
            metadata.at[row.Index, "genomic_end_0based_exclusive"] = min(chrom_len, raw_end)
            metadata.at[row.Index, "upstream_n_padding"] = upstream_pad
            metadata.at[row.Index, "downstream_n_padding"] = downstream_pad
            metadata.at[row.Index, "n_non_acgt_in_max_context"] = int((codes == 4).sum())
            metadata.at[row.Index, "acceptor_23mer"] = window[half - 20 : half + 3].decode()
            center = window[half - 2 : half].decode()
            metadata.at[row.Index, "center_dinucleotide"] = center
            metadata.at[row.Index, "canonical_acceptor"] = center == "AG"
        del chromosome

    store.flush()
    del store
    metadata.to_csv(metadata_path, sep="\t", index=False, compression="gzip")
    by_strand = {}
    for strand, group in metadata.groupby("strand"):
        canonical = int(group["canonical_acceptor"].sum())
        by_strand[strand] = {
            "n_sites": len(group),
            "canonical_acceptors": canonical,
            "canonical_acceptor_fraction": canonical / len(group),
        }
    validation_examples = []
    for strand in ("+", "-"):
        for row in metadata[metadata["strand"] == strand].head(3).itertuples():
            validation_examples.append({
                "site_id": row.site_id,
                "strand": strand,
                "acceptor_23mer": row.acceptor_23mer,
                "center_dinucleotide": row.center_dinucleotide,
                "canonical_acceptor": bool(row.canonical_acceptor),
            })
    summary = {
        "n_sites": len(metadata),
        "context_lengths": widths,
        "stored_context_length": max_width,
        "sequence_shape": [len(metadata), max_width],
        "sequence_dtype": "uint8",
        "codebook": BASE_TO_CODE,
        "ambiguous_base_policy": "All non-ACGT symbols, including IUPAC ambiguity codes, are encoded as N=4.",
        "one_hot_policy": "Created on demand; N=4 maps to [0,0,0,0].",
        "canonical_acceptors": int(metadata["canonical_acceptor"].sum()),
        "canonical_acceptor_fraction": float(metadata["canonical_acceptor"].mean()) if len(metadata) else None,
        "canonical_acceptors_by_strand": by_strand,
        "sites_with_padding": int(((metadata["upstream_n_padding"] + metadata["downstream_n_padding"]) > 0).sum()),
        "sites_with_non_acgt_in_max_context": int((metadata["n_non_acgt_in_max_context"] > 0).sum()),
        "validation_examples": validation_examples,
        "sites_source": None,
        "prepared_sequence_dir": portable_path(sequence_dir),
        "sequence_file": portable_path(sequence_path),
        "metadata_file": portable_path(metadata_path),
    }
    manifest_path.write_text(json.dumps(summary, indent=2) + "\n")
    return sequence_path, metadata_path, manifest_path, summary


def main() -> None:
    cfg = load_config()
    default_sites = REPO_ROOT / "results" / "person_b_acceptors_chr1-4.tsv"
    default_seq_dir = project_path(cfg, "data_processed") / "grch38_contigs"
    default_prefix = project_path(cfg, "data_processed") / "person_b_acceptors_chr1-4"
    default_fasta = REPO_ROOT / cfg["reference"]["genome_fasta"]

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default=str(default_sites))
    ap.add_argument("--fasta", default=str(default_fasta))
    ap.add_argument("--sequence-dir", default=str(default_seq_dir))
    ap.add_argument("--out-prefix", default=str(default_prefix))
    ap.add_argument("--qc-out", default=str(project_path(cfg, "results") / "sequence_qc_chr1-4.json"))
    ap.add_argument("--contexts", nargs="+", type=int, default=cfg["sequence"]["context_lengths"])
    ap.add_argument("--prepare-reference", action="store_true", help="extract the required contigs from the FASTA first")
    ap.add_argument("--prepare-only", action="store_true", help="prepare contigs but do not extract site windows")
    ap.add_argument("--limit", type=int, help="extract only the first N rows (for a smoke test)")
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str, "strand": str})
    if args.limit is not None:
        sites = sites.head(args.limit).copy()
    chromosomes = sites["chrom"].drop_duplicates().tolist()
    seq_dir = Path(args.sequence_dir)
    need_prepare = args.prepare_reference or any(not (seq_dir / f"{chrom}.seq").exists() for chrom in chromosomes)
    if need_prepare:
        lengths = prepare_reference(args.fasta, seq_dir, chromosomes)
        print("prepared contigs:", ", ".join(f"{chrom}={lengths[chrom]:,}" for chrom in chromosomes))
    if args.prepare_only:
        return
    sequence_path, metadata_path, manifest_path, summary = extract_sequences(sites, seq_dir, args.out_prefix, args.contexts)
    summary["sites_source"] = portable_path(args.sites)
    manifest_path.write_text(json.dumps(summary, indent=2) + "\n")
    qc_path = Path(args.qc_out)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"{summary['n_sites']:,} sites -> {sequence_path}")
    print(f"metadata -> {metadata_path}")
    print(f"canonical AG at centre: {summary['canonical_acceptors']:,}/{summary['n_sites']:,} "
          f"({summary['canonical_acceptor_fraction']:.3%})")
    print(f"sequence QC -> {qc_path}")


if __name__ == "__main__":
    main()
