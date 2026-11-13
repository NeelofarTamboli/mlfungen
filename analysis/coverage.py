"""Coverage analysis: how does replicate disagreement change with read support?

This is evidence for a later read-support filter. It does NOT choose one, and a threshold must not be picked from
tables made on held-out test chromosomes: use training/validation chromosomes only. Read support of a site is
taken as the smaller of its two replicate supports (a site is only as reliable as its weaker replicate).
Bins are powers of two (neutral); each row says how big the bin is and how well the replicates agree there.

Used on matched sites that are not also_transcript_start. The "0" bin is sites with no reads in at least one
replicate; every other bin is nonzero support in both replicates, i.e. the primary QC set.

For every bin ("by_bin") and for every cumulative cut "support >= t" ("cumulative"):
    n sites, median / mean / 90th percentile of |rep1 - rep2|, share with |diff| > 0.1,
    Pearson and Spearman (n >= 10), median bootstrap SD, and the median |diff| expected from bootstrap noise
    alone (0.6745 * sqrt(sd1^2 + sd2^2) per site). Observed well above expected means the replicates
    differ by more than McSplicer's own sampling noise (replicate-specific effects), so more reads alone
    would not help there.

Outputs:
    results/coverage_by_bin.tsv, results/coverage_cumulative.tsv
    figures/coverage_disagreement.png      per-bin median and 90th percentile of |diff| with the bin sizes
    figures/coverage_cumulative.png        sites kept and median |diff| / noise floor if a cut were applied

Usage:
    python -m analysis.coverage [--labels data/processed/acceptor_labels.tsv] [--tag NAME]
"""

from __future__ import annotations

import argparse
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analysis.common import (
    ACCENT,
    INK_2,
    SUPPORT_EDGES,
    load_labels,
    new_figure,
    support_bin_label,
)
from src.config import REPO_ROOT, load_config, project_path


def summarize(d: pd.DataFrame) -> dict:
    out = {"n_sites": len(d)}
    if len(d) == 0:
        return out
    out.update(median_abs_diff=d.abs_diff.median(), mean_abs_diff=d.abs_diff.mean(), p90_abs_diff=d.abs_diff.quantile(0.9),
               share_diff_over_0_1=(d.abs_diff > 0.1).mean(), median_bootstrap_sd=d.sd_mean.median(),
               median_expected_from_bootstrap=d.expected_abs_diff.median())
    if len(d) >= 10 and d.usage_rep1.nunique() > 1 and d.usage_rep2.nunique() > 1:
        out.update(pearson=stats.pearsonr(d.usage_rep1, d.usage_rep2)[0], spearman=stats.spearmanr(d.usage_rep1, d.usage_rep2)[0])
    return out


def tables(both: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_bin = [{"support_min": support_bin_label(lo, hi), "lo": lo, "hi": hi, **summarize(both[(both.support_min >= lo) & (both.support_min < hi)])}
              for lo, hi in pairwise(SUPPORT_EDGES)]
    cumulative = [{"support_min_at_least": lo, **summarize(both[both.support_min >= lo]), "share_of_matched_kept": (both.support_min >= lo).mean()}
                  for lo in SUPPORT_EDGES[:-1]]
    return pd.DataFrame(by_bin), pd.DataFrame(cumulative)


def plots(by_bin: pd.DataFrame, cum: pd.DataFrame, fig_dir: Path, tag: str) -> list[Path]:
    paths = []
    b = by_bin[by_bin.n_sites >= 5]
    fig, (ax, ax2) = new_figure(2, 1, figsize=(7.2, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1.6]})
    x = np.arange(len(b))
    ax.plot(x, b.p90_abs_diff, color=INK_2, lw=1.4, marker="s", ms=4, ls="--")
    ax.plot(x, b.median_abs_diff, color=ACCENT, lw=2, marker="o", ms=5)
    ax.set_ylabel("|usage rep 1 - usage rep 2|")
    ax.set_title("Replicate disagreement by read-support bin", fontsize=11, loc="left")
    ax.text(0.99, 0.97, "solid: median\ndashed: 90th percentile\n(bins with >= 5 matched sites)", transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK_2)
    ax2.bar(x, b.n_sites, color=ACCENT, width=0.7)
    ax2.set_yscale("log"); ax2.set_ylabel("matched sites"); ax2.set_xlabel("read support (the smaller of the two replicates)")
    ax2.set_xticks(x); ax2.set_xticklabels(b.support_min, rotation=45, ha="right", fontsize=8.5)
    p = fig_dir / f"coverage_disagreement{tag}.png"; fig.tight_layout(); fig.savefig(p, dpi=200); paths.append(p)

    c = cum[cum.n_sites >= 5]
    fig, (ax, ax2) = new_figure(2, 1, figsize=(7.2, 6.0), sharex=True, gridspec_kw={"height_ratios": [1.6, 3]})
    ax.plot(c.support_min_at_least.clip(lower=0.5), c.n_sites, color=ACCENT, lw=2, marker="o", ms=4)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel("matched sites kept")
    ax.set_title("What a minimum-support cut would keep, and how well the kept sites agree", fontsize=11, loc="left")
    ax2.plot(c.support_min_at_least.clip(lower=0.5), c.median_abs_diff, color=ACCENT, lw=2, marker="o", ms=4)
    ax2.plot(c.support_min_at_least.clip(lower=0.5), c.median_expected_from_bootstrap, color=INK_2, lw=1.4, ls="--", marker="s", ms=3.5)
    ax2.set_ylabel("median |rep 1 - rep 2| of kept sites")
    ax2.set_xlabel("minimum read support (both replicates)")
    ax2.set_xticks([0.5, 1, 10, 100, 1000]); ax2.set_xticklabels(["0", "1", "10", "100", "1k"]); ax2.minorticks_off()
    ax2.text(0.99, 0.97, "solid: observed\ndashed: expected from bootstrap noise alone", transform=ax2.transAxes, ha="right", va="top", fontsize=8.5, color=INK_2)
    p = fig_dir / f"coverage_cumulative{tag}.png"; fig.tight_layout(); fig.savefig(p, dpi=200); paths.append(p)
    return paths


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", default=str(project_path(cfg, "data_processed") / "acceptor_labels.tsv"))
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    tag = f"_{a.tag.lstrip('_')}" if a.tag else ""
    res_dir, fig_dir = project_path(cfg, "results"), project_path(cfg, "figures")
    res_dir.mkdir(exist_ok=True); fig_dir.mkdir(exist_ok=True)
    df = load_labels(a.labels)
    both = df[df.in_both & ~df.also_transcript_start]
    by_bin, cum = tables(both)
    by_bin.to_csv(res_dir / f"coverage_by_bin{tag}.tsv", sep="\t", index=False, float_format="%.4g")
    cum.to_csv(res_dir / f"coverage_cumulative{tag}.tsv", sep="\t", index=False, float_format="%.4g")
    figs = plots(by_bin, cum, fig_dir, tag)
    show = ["support_min", "n_sites", "median_abs_diff", "p90_abs_diff", "share_diff_over_0_1", "median_expected_from_bootstrap", "spearman"]
    print(by_bin[[c for c in show if c in by_bin]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("figures:", [str(p.relative_to(REPO_ROOT)) for p in figs])


if __name__ == "__main__":
    main()
