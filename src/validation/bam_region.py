"""Pull one genomic region out of a remote coordinate-sorted BAM, without an index.

ENCODE BAMs come without a .bai. This binary-searches the BGZF block stream with
HTTP range requests to find the region, downloads only the blocks that cover it
(a few MB), and writes them as a new small BAM that keeps the full original header.

Reads are kept by their start position, so a pair whose mate starts outside the
region is incomplete; choose a region with margin around the genes you care about.

Usage:
    python -m src.validation.bam_region URL CHROM START END OUT.bam [--size BYTES] [--index]

START and END are 1-based inclusive. The URL must support HTTP Range requests
(the encode-public S3 URLs listed in the ENCODE file records do).
"""

from __future__ import annotations

import argparse
import struct
import subprocess
import sys
import zlib

BGZF_HEADER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

_downloaded = 0


def _fetch(url: str, first: int, last: int) -> bytes:
    global _downloaded
    out = subprocess.run(
        ["curl", "-sS", "-m", "120", "-r", f"{first}-{last}", url], capture_output=True, check=True
    ).stdout
    _downloaded += len(out)
    return out


def _file_size(url: str) -> int:
    head = subprocess.run(["curl", "-sSIL", url], capture_output=True, check=True).stdout.decode()
    sizes = [int(line.split(":")[1]) for line in head.splitlines() if line.lower().startswith("content-length")]
    return sizes[-1]


def _blocks(chunk: bytes):
    """Yield (start, end, uncompressed) for each BGZF block fully inside chunk."""
    p = chunk.find(BGZF_HEADER)
    while p >= 0 and p + 18 <= len(chunk):
        end = p + struct.unpack("<H", chunk[p + 16 : p + 18])[0] + 1
        if end > len(chunk):
            return
        # the next block must also start with the BGZF magic, else this was a false hit
        if end + len(BGZF_HEADER) <= len(chunk) and chunk[end : end + len(BGZF_HEADER)] != BGZF_HEADER:
            p = chunk.find(BGZF_HEADER, p + 1)
            continue
        yield p, end, zlib.decompress(chunk[p + 18 : end - 8], -15)
        p = end


def _parse_header(buf: bytes):
    """Return (header_bytes, reference_names, reference_lengths) from the start of a BAM stream."""
    assert buf[:4] == b"BAM\x01", "not a BAM stream"
    text_len = struct.unpack("<i", buf[4:8])[0]
    p = 8 + text_len
    n_ref = struct.unpack("<i", buf[p : p + 4])[0]
    p += 4
    names, lengths = [], []
    for _ in range(n_ref):
        l_name = struct.unpack("<i", buf[p : p + 4])[0]
        names.append(buf[p + 4 : p + 4 + l_name - 1].decode())
        lengths.append(struct.unpack("<i", buf[p + 4 + l_name : p + 8 + l_name])[0])
        p += 8 + l_name
    return buf[:p], names, lengths


def _record_at(buf: bytes, i: int, n_ref: int):
    """If a plausible BAM record starts at i return (refid, pos, next_index), else None."""
    if i + 36 > len(buf):
        return None
    size = struct.unpack("<i", buf[i : i + 4])[0]
    if not 32 <= size <= 200000:
        return None
    refid, pos, l_name, _mapq, _bin, n_cigar, _flag, l_seq, mate_ref, _mpos, _tlen = struct.unpack(
        "<iiBBHHHiiii", buf[i + 4 : i + 36]
    )
    if not (-1 <= refid < n_ref and -1 <= pos < 300_000_000 and 2 <= l_name <= 255 and -1 <= mate_ref < n_ref):
        return None
    if i + 4 + size > len(buf) or buf[i + 35 + l_name] != 0:
        return None
    if not all(32 < c < 127 for c in buf[i + 36 : i + 35 + l_name]):
        return None
    if 32 + l_name + 4 * n_cigar + (l_seq + 1) // 2 + l_seq > size:
        return None
    return refid, pos, i + 4 + size


def _first_record(buf: bytes, n_ref: int, chain: int = 4):
    """Find where the first complete record starts by requiring `chain` consecutive valid records."""
    for i in range(min(len(buf), 400000)):
        j, first, ok = i, None, True
        for _ in range(chain):
            r = _record_at(buf, j, n_ref)
            if r is None:
                ok = False
                break
            first = first or r
            j = r[2]
        if ok:
            return i, first[0], first[1]
    return None


def _probe(url: str, offset: int, n_ref: int, size: int):
    chunk = _fetch(url, offset, min(offset + 400_000, size - 1))
    buf, first = b"", None
    for start, _end, data in _blocks(chunk):
        first = start if first is None else first
        buf += data
        if len(buf) > 250_000:
            break
    rec = _first_record(buf, n_ref)
    return (offset + first, rec[1], rec[2]) if rec else None


def _bgzf(data: bytes) -> bytes:
    out = b""
    for i in range(0, len(data), 0xFF00):
        part = data[i : i + 0xFF00]
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = comp.compress(part) + comp.flush()
        out += BGZF_HEADER[:16] + struct.pack("<H", len(body) + 25) + body
        out += struct.pack("<II", zlib.crc32(part) & 0xFFFFFFFF, len(part))
    return out


def extract_region(url: str, chrom: str, start: int, end: int, out_path: str, size: int | None = None) -> int:
    """Write records starting in chrom:start-end to out_path. Returns the number of records kept."""
    size = size or _file_size(url)
    head = _fetch(url, 0, 2_097_151)
    header_bytes, names, _lengths = _parse_header(b"".join(d for _, _, d in _blocks(head)))
    n_ref, tid = len(names), names.index(chrom)

    target = (tid, start - 1 - 30000)  # land a bit before the region
    lo, hi = 3_000_000, size - 1_000_000
    while hi - lo > 300_000:
        mid = (lo + hi) // 2
        probe = _probe(url, mid, n_ref, size)
        if probe is None or (probe[1], probe[2]) < target:
            lo = mid
        else:
            hi = mid

    records, buf, offset, idx, begun, done = [], b"", lo, 0, False, False
    while not done:
        chunk = _fetch(url, offset, min(offset + 4_000_000 - 1, size - 1))
        consumed = None
        for _s, e, data in _blocks(chunk):
            buf += data
            consumed = e
        if consumed is None:
            raise RuntimeError("no complete BGZF block in downloaded chunk")
        offset += consumed
        if not begun:
            found = _first_record(buf, n_ref)
            if not found:
                raise RuntimeError("could not find a record boundary")
            idx, begun = found[0], True
        while True:
            rec = _record_at(buf, idx, n_ref)
            if rec is None:
                break
            refid, pos, nxt = rec
            records.append((refid, pos, buf[idx:nxt]))
            idx = nxt
            if refid > tid or (refid == tid and pos > end):
                done = True
                break
        buf, idx = buf[idx:], 0

    kept = [r for r in records if r[0] == tid and start - 1 <= r[1] <= end]
    with open(out_path, "wb") as fh:
        fh.write(_bgzf(header_bytes + b"".join(r[2] for r in kept)) + BGZF_EOF)
    return len(kept)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("url")
    ap.add_argument("chrom")
    ap.add_argument("start", type=int)
    ap.add_argument("end", type=int)
    ap.add_argument("out")
    ap.add_argument("--size", type=int, help="file size in bytes (default: ask the server)")
    ap.add_argument("--index", action="store_true", help="also write OUT.bai (needs pysam)")
    a = ap.parse_args()
    n = extract_region(a.url, a.chrom, a.start, a.end, a.out, a.size)
    if a.index:
        import pysam

        pysam.index(a.out)
    print(f"{a.chrom}:{a.start}-{a.end}: kept {n} records, downloaded {_downloaded / 1e6:.1f} MB -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
