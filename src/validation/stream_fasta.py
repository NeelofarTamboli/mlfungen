"""Stream a gzipped FASTA once: list every contig and save chosen contigs as plain text.

The ENCODE GRCh38 FASTA is a plain gzip (not bgzip), so it cannot be indexed or
range-read. This reads it sequentially (873 MB compressed), writes
`contigs.json` (name, full header line, length for every contig) and one
`<name>.seq` file (sequence on a single line, no header) for each requested contig.

Usage:
    python -m src.validation.stream_fasta URL OUT_DIR chr1,chr2
"""

from __future__ import annotations

import json
import subprocess
import sys
import zlib
from pathlib import Path


def stream_fasta(url: str, out_dir: Path, keep: set[str]) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    contigs: list[dict] = []
    state = {"cur": None, "fh": None}

    def feed(block: bytes) -> None:
        for line in block.split(b"\n"):
            if not line:
                continue
            if line.startswith(b">"):
                if state["fh"]:
                    state["fh"].close()
                    state["fh"] = None
                header = line[1:].decode()
                state["cur"] = {"name": header.split()[0], "header": header, "length": 0}
                contigs.append(state["cur"])
                if state["cur"]["name"] in keep:
                    state["fh"] = open(out_dir / f"{state['cur']['name']}.seq", "wb")  # noqa: SIM115 (stays open across chunks)
            else:
                state["cur"]["length"] += len(line)
                if state["fh"]:
                    state["fh"].write(line)

    proc = subprocess.Popen(["curl", "-sS", url], stdout=subprocess.PIPE)
    dec, carry = zlib.decompressobj(31), b""
    while chunk := proc.stdout.read(1 << 20):
        data = carry + dec.decompress(chunk)
        cut = data.rfind(b"\n")
        if cut < 0:
            carry = data
            continue
        feed(data[:cut])
        carry = data[cut + 1 :]
    feed(carry)
    if state["fh"]:
        state["fh"].close()
    if proc.wait() != 0:
        raise RuntimeError("curl failed")
    (out_dir / "contigs.json").write_text(json.dumps(contigs))
    return contigs


if __name__ == "__main__":
    url, out, keep = sys.argv[1], Path(sys.argv[2]), set(sys.argv[3].split(","))
    result = stream_fasta(url, out, keep)
    print(f"{len(result)} contigs written to {out / 'contigs.json'}", file=sys.stderr)
