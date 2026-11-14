"""Score acceptor sites with pretrained SpliceAI and Pangolin (no training on our labels).

Both models take 5,000 bases on each side of a position and predict splicing at that position.
For each site we build the oriented sequence (transcript strand, N-padded at chromosome ends),
centred on the first exonic base (`src.extract_sequences.acceptor_window`), and keep the
predictions for offsets -K..+K around it, so the position convention of each model can be checked
instead of assumed.

    SpliceAI (Jaganathan et al. 2019): Keras models spliceai1-5.h5 from the `spliceai` package,
        output (neither, acceptor, donor); we keep the acceptor probability, averaged over the 5 models.
    Pangolin (Zeng & Li 2022): PyTorch weights `final.{1..3}.{m}.3.v2` from external/Pangolin (the
        weights the current pangolin CLI loads), applied as in the authors' scripts/custom_usage.py.
        m = 0..7 is heart / liver / brain / testis x (P(splice), usage); each output is averaged over
        its 3 models.

Both run in bfloat16 on CPU (about 4x faster for Pangolin and 1.7x for SpliceAI here). On 64 sites the
bf16 outputs matched float32 to within 0.011 (correlation > 0.9999). `--fp32` turns it off.

SpliceAI and Pangolin were both trained with chromosomes 1, 3, 5, 7 and 9 held out, so only those
chromosomes are scored here (our current test set is chr1 and chr3). Predictions are cached per
chromosome under data/processed/pretrained/ (git-ignored); metrics, together with the existing
local baselines on the same test rows, go to results/pretrained_metrics_test_chr1_chr3.tsv.

    python -m models.pretrained_baselines score --model spliceai --chrom chr1 chr3
    python -m models.pretrained_baselines score --model pangolin --chrom chr1 chr3
    python -m models.pretrained_baselines evaluate
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from src.config import REPO_ROOT, load_config, project_path
from src.extract_sequences import acceptor_window, encode_sequence

FLANK = 5000
K = 2  # offsets -K..+K around the first exonic base
TISSUES = ["heart", "liver", "brain", "testis"]
PANGOLIN_HEADS = [f"{t}_{kind}" for t in TISSUES for kind in ("p_splice", "usage")]  # model numbers 0..7
PANGOLIN_CHANNEL = {0: 1, 1: 2, 2: 4, 3: 5, 4: 7, 5: 8, 6: 10, 7: 11}  # from custom_usage.py INDEX_MAP
HELD_OUT = {"chr1", "chr3", "chr5", "chr7", "chr9"}


def site_codes(sites: pd.DataFrame, chrom_seq: bytes, k: int = K) -> np.ndarray:
    """(n, 2*FLANK + 2k + 1) uint8 codes; column FLANK + k is the first exonic base."""
    width = 2 * FLANK + 2 * k + 2
    rows = [encode_sequence(acceptor_window(chrom_seq, loc, s, width))[1:] for loc, s in zip(sites["locus"], sites["strand"])]
    return np.stack(rows)


def one_hot(codes: np.ndarray) -> np.ndarray:
    """A,C,G,T -> unit vectors, N -> zeros (both models use this encoding). Shape (n, L, 4)."""
    return np.eye(5, 4, dtype=np.float32)[codes]


class SpliceAI:
    def __init__(self, bf16: bool = True) -> None:
        import keras
        import spliceai

        if bf16:
            keras.mixed_precision.set_global_policy("mixed_bfloat16")

        d = Path(spliceai.__file__).parent / "models"
        self.models = [keras.models.load_model(d / f"spliceai{i}.h5", compile=False) for i in range(1, 6)]

    def __call__(self, codes: np.ndarray) -> dict[str, np.ndarray]:
        x = one_hot(codes)
        acc = np.mean([np.asarray(m(x, training=False), dtype=np.float32)[:, :, 1] for m in self.models], axis=0)
        return {"acceptor": acc}


class Pangolin:
    def __init__(self, bf16: bool = True) -> None:
        import sys

        import torch

        sys.path.insert(0, str(REPO_ROOT / "external" / "Pangolin"))
        from pangolin.model import AR, L, W
        from pangolin.model import Pangolin as Net

        self.torch, self.bf16 = torch, bf16
        wdir = REPO_ROOT / "external" / "Pangolin" / "pangolin" / "models"
        self.models: dict[int, list] = {}
        for m in range(8):
            self.models[m] = []
            for j in range(1, 4):
                net = Net(L, W, AR)
                net.load_state_dict(torch.load(wdir / f"final.{j}.{m}.3.v2", map_location="cpu"))
                net.eval()
                self.models[m].append(net)

    def __call__(self, codes: np.ndarray) -> dict[str, np.ndarray]:
        x = self.torch.from_numpy(one_hot(codes).transpose(0, 2, 1).copy())
        out = {}
        with self.torch.inference_mode(), self.torch.autocast("cpu", dtype=self.torch.bfloat16, enabled=self.bf16):
            for m, nets in self.models.items():
                out[PANGOLIN_HEADS[m]] = np.mean([net(x)[:, PANGOLIN_CHANNEL[m], :].float().numpy() for net in nets], axis=0)
        return out


def score_chrom(model_name: str, scorer, sites: pd.DataFrame, seq_dir: Path, out_dir: Path, batch: int, chunk: int = 512) -> Path:
    """Score one chromosome in chunks (resumable: finished chunks are kept in a .part directory)."""
    chrom = sites["chrom"].iloc[0]
    out = out_dir / f"{model_name}_{chrom}.npz"
    if out.exists():
        print(f"already scored {out}")
        return out
    part = out_dir / f"{model_name}_{chrom}.part"
    part.mkdir(parents=True, exist_ok=True)
    seq = (seq_dir / f"{chrom}.seq").read_bytes()
    t0 = time.time()
    n_chunks = -(-len(sites) // chunk)
    for c in range(n_chunks):
        f = part / f"{c:05d}.npz"
        if f.exists():
            continue
        block = sites.iloc[c * chunk : (c + 1) * chunk]
        codes = site_codes(block, seq)
        res: dict[str, list] = {}
        for b in range(0, len(codes), batch):
            for key, v in scorer(codes[b : b + batch]).items():
                res.setdefault(key, []).append(v)
        np.savez(f, **{key: np.concatenate(v) for key, v in res.items()})
        done = c + 1
        print(f"  {model_name} {chrom}: {min(done * chunk, len(sites))}/{len(sites)} sites, {time.time() - t0:.0f} s", flush=True)
    parts = [np.load(part / f"{c:05d}.npz") for c in range(n_chunks)]
    np.savez(out, site_id=sites["site_id"].to_numpy(), **{key: np.concatenate([p[key] for p in parts]) for key in parts[0].files})
    for f in part.glob("*.npz"):
        f.unlink()
    part.rmdir()
    return out


def load_predictions(model_name: str, chroms: list[str], out_dir: Path) -> pd.DataFrame:
    """Long table: site_id + one column per output and offset (e.g. acceptor_+0)."""
    frames = []
    for chrom in chroms:
        z = np.load(out_dir / f"{model_name}_{chrom}.npz", allow_pickle=True)
        cols = {"site_id": z["site_id"]}
        for key in z.files:
            if key == "site_id":
                continue
            for i, off in enumerate(range(-K, K + 1)):
                cols[f"{model_name}_{key}_{off:+d}"] = z[key][:, i]
        frames.append(pd.DataFrame(cols))
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    cfg = load_config()
    processed = project_path(cfg, "data_processed")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score")
    s.add_argument("--model", required=True, choices=["spliceai", "pangolin"])
    s.add_argument("--chrom", nargs="+", default=["chr1", "chr3"])
    s.add_argument("--limit", type=int, help="score only the first N sites per chromosome (benchmark)")
    s.add_argument("--batch", type=int, default=16)
    s.add_argument("--threads", type=int, default=os.cpu_count())
    s.add_argument("--fp32", action="store_true", help="run in float32 instead of bfloat16")
    e = sub.add_parser("evaluate")
    e.add_argument("--chrom", nargs="+", default=["chr1", "chr3"])
    for p in (s, e):
        p.add_argument("--sites", default=str(REPO_ROOT / "results" / "person_b_acceptors_chr1-4.tsv"))
        p.add_argument("--seq-dir", default=str(processed / "grch38_contigs"))
        p.add_argument("--out-dir", default=str(processed / "pretrained"))
    e.add_argument("--metrics-out", default=str(project_path(cfg, "results") / "pretrained_metrics_test_chr1_chr3.tsv"))
    e.add_argument("--offsets-out", default=str(project_path(cfg, "results") / "pretrained_offsets_test_chr1_chr3.tsv"))
    e.add_argument("--baseline-metrics", default=str(project_path(cfg, "results") / "baseline_metrics_chr1-4.tsv"))
    a = ap.parse_args()

    leaked = sorted(set(a.chrom) - HELD_OUT)
    if leaked:
        ap.error(f"{leaked} were training chromosomes for SpliceAI/Pangolin; score only {sorted(HELD_OUT)}")
    sites = pd.read_csv(a.sites, sep="\t", dtype={"chrom": str, "strand": str})
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if a.cmd == "score":
        if a.limit:
            out_dir = out_dir / f"bench_{a.limit}"
            out_dir.mkdir(exist_ok=True)
        if a.model == "pangolin":
            import torch

            torch.set_num_threads(a.threads)
        scorer = (SpliceAI if a.model == "spliceai" else Pangolin)(bf16=not a.fp32)
        for chrom in a.chrom:
            block = sites[sites["chrom"] == chrom].reset_index(drop=True)
            if a.limit:
                block = block.head(a.limit)
            t0 = time.time()
            path = score_chrom(a.model, scorer, block, Path(a.seq_dir), out_dir, a.batch)
            print(f"{a.model} {chrom}: {len(block)} sites in {time.time() - t0:.0f} s -> {path}")
        return

    evaluate(sites, a.chrom, out_dir, Path(a.metrics_out), Path(a.offsets_out), Path(a.baseline_metrics))


def evaluate(sites: pd.DataFrame, chroms: list[str], out_dir: Path, metrics_out: Path, offsets_out: Path, baseline_metrics: Path) -> None:
    from analysis.evaluate import evaluate_predictions

    data = sites[sites["chrom"].isin(chroms)].copy()
    data["target"] = data[["usage_rep1", "usage_rep2"]].mean(axis=1)  # same preliminary target as models/baseline.py --target mean
    data["split"] = "test"
    for name in ("spliceai", "pangolin"):
        if all((out_dir / f"{name}_{c}.npz").exists() for c in chroms):
            data = data.merge(load_predictions(name, chroms, out_dir), on="site_id", how="left", validate="one_to_one")
    score_cols = [c for c in data.columns if c.startswith(("spliceai_", "pangolin_"))]
    if not score_cols:
        raise SystemExit(f"no complete predictions in {out_dir}")

    # Offset convention: for each output, the offset with the highest mean score over all sites.
    # This uses the model scores only (no labels), so it does not tune anything to our targets.
    off = pd.DataFrame([{"output": c.rsplit("_", 1)[0], "offset": int(c.rsplit("_", 1)[1]), "mean_score": data[c].mean()} for c in score_cols])
    best = off.loc[off.groupby("output")["mean_score"].idxmax()].set_index("output")["offset"]
    offsets_out.parent.mkdir(parents=True, exist_ok=True)
    off.pivot(index="output", columns="offset", values="mean_score").to_csv(offsets_out, sep="\t", float_format="%.4f")

    preds = pd.DataFrame({"site_id": data["site_id"]})
    for output, o in best.items():
        preds[output] = data[f"{output}_{o:+d}"].to_numpy()
    usage = [f"pangolin_{t}_usage" for t in TISSUES]
    if set(usage) <= set(preds.columns):
        preds["pangolin_mean_usage"] = preds[usage].mean(axis=1)
    pred_cols = [c for c in preds.columns if c != "site_id"]
    data = data.drop(columns=score_cols).join(preds.drop(columns="site_id"))
    data[["site_id", "chrom", "target", "in_primary_qc", *pred_cols]].to_csv(out_dir / "predictions_test.tsv.gz", sep="\t", index=False)

    # Reported outputs. SpliceAI gives a splice-site probability, not a usage value, so only rank/linear
    # agreement (Pearson, Spearman) is meaningful; R2 and MAE are left empty for it.
    report = {"spliceai_acceptor": ("SpliceAI", "splice-site probability")}
    report.update({f"pangolin_{t}_usage": (f"Pangolin {t}", "usage") for t in TISSUES})
    report["pangolin_mean_usage"] = ("Pangolin mean of 4 tissues", "usage")
    report = {k: v for k, v in report.items() if k in data}
    m = evaluate_predictions(data, list(report), target_column="target", split="test")
    m["score_type"] = m["model"].map(lambda c: report[c][1])
    m["model"] = m["model"].map(lambda c: report[c][0])
    m.loc[m["score_type"] == "splice-site probability", ["r2", "mae"]] = np.nan

    # Existing local baselines (models/baseline.py --target mean) on the same test rows.
    base = pd.read_csv(baseline_metrics, sep="\t")
    base = base[base["split"] == "test"].copy()
    base["model"] = base["model"].map({"mean": "Training mean", "ridge": "Ridge regression", "gradient_boosting": "Gradient boosting"})
    base["score_type"] = "usage (trained on chr2, chr4 val)"
    for subset, n in m.groupby("subset")["n_sites"].max().items():
        nb = set(base.loc[base["subset"] == subset, "n_sites"])
        if nb != {n}:
            raise SystemExit(f"site counts differ from the baseline table for {subset}: {n} vs {sorted(nb)}")
    cols = ["subset", "model", "score_type", "n_sites", "pearson", "spearman", "r2", "mae"]
    metrics = pd.concat([base[cols], m[cols]], ignore_index=True)
    metrics["subset"] = pd.Categorical(metrics["subset"], ["all", "alternative_0.1_to_0.9"])
    metrics = metrics.sort_values("subset", kind="stable")
    metrics.to_csv(metrics_out, sep="\t", index=False, float_format="%.4f")
    print("offsets (highest mean score):", json.dumps({k: int(v) for k, v in best.items()}))
    print(metrics.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"metrics -> {metrics_out}")

if __name__ == "__main__":
    main()
