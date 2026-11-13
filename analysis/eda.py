"""Interim replicate QC of the internal-acceptor label table (K562, McSplicer).

Sets used (all from the label table built by src/build_dataset.py):
    matched        site present in both replicates
    nonzero_both   matched, with read support > 0 in both replicates
    primary        nonzero_both and not also_transcript_start   <- the main QC set (column in_primary_qc)
No other read-support threshold is applied, and replicates are not averaged.

Outputs (TAG, e.g. chr1-4, goes in the file names):
    results/replicate_qc_TAG.json, results/replicate_qc_TAG.txt
    results/support_bins_TAG.tsv                 agreement by read-support bin (matched, not also_transcript_start)
    figures/qc_A_replicate_scatter_TAG.png       usage rep 1 vs rep 2 (primary set)
    figures/qc_B_usage_distribution_TAG.png      usage distribution per replicate (primary set)
    figures/qc_C_support_vs_disagreement_TAG.png read support vs |rep1 - rep2| (primary set; zero-support sites in grey)
    figures/qc_D_sites_by_support_bin_TAG.png    number of matched sites per read-support bin
    figures/qc_E_bootstrap_sd_vs_support_TAG.png bootstrap SD vs support, only where bootstraps exist (see --bootstrap-table);
                                                 exploratory: only a few dozen sites, not a main result

Usage:
    python -m analysis.eda --labels data/processed/acceptor_labels_chr1-4.tsv --tag chr1-4 \\
        [--bootstrap-table data/interim/region_test_processed/acceptor_labels.tsv]
"""

from __future__ import annotations

import argparse
import json
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analysis.common import (
    ACCENT,
    INK,
    INK_2,
    SUPPORT_EDGES,
    as_plot_support,
    load_labels,
    log_support_axis,
    new_figure,
    support_bin_label,
)
from analysis.coverage import tables
from src.config import REPO_ROOT, load_config, project_path

USAGE_BINS = [0, 0.05, 0.2, 0.5, 0.8, 0.95, 1.0000001]
USAGE_BIN_NAMES = ["<0.05", "0.05-0.2", "0.2-0.5", "0.5-0.8", "0.8-0.95", ">=0.95"]
GREY = "#9a9993"  # the zero-support group


def describe(s: pd.Series) -> dict:
    s = s.dropna()
    if s.empty:
        return {"n": 0}
    q = s.quantile([0.1, 0.25, 0.5, 0.75, 0.9])
    return {"n": len(s), "min": float(s.min()), "p10": float(q[0.1]), "p25": float(q[0.25]), "median": float(q[0.5]),
            "p75": float(q[0.75]), "p90": float(q[0.9]), "max": float(s.max())}


def usage_counts(s: pd.Series) -> dict:
    c = pd.cut(s.dropna(), USAGE_BINS, right=False, labels=USAGE_BIN_NAMES).value_counts().reindex(USAGE_BIN_NAMES)
    return {k: int(v) for k, v in c.items()}


def agreement(d: pd.DataFrame) -> dict:
    """Replicate agreement on a set of matched sites."""
    if len(d) < 3:
        return {"n": len(d)}
    diff = (d.usage_rep1 - d.usage_rep2).abs()
    return {"n": len(d), "pearson": float(stats.pearsonr(d.usage_rep1, d.usage_rep2)[0]), "spearman": float(stats.spearmanr(d.usage_rep1, d.usage_rep2)[0]),
            "median_abs_diff": float(diff.median()), "mean_abs_diff": float(diff.mean()), "p90_abs_diff": float(diff.quantile(0.9)),
            "share_abs_diff_over_0.1": float((diff > 0.1).mean())}


def group_profile(d: pd.DataFrame) -> dict:
    """How a group of matched sites behaves: usage, support and replicate agreement."""
    u = pd.concat([d.usage_rep1, d.usage_rep2])
    return {**agreement(d), "usage_median": float(u.median()), "usage_share_below_0.05": float((u < 0.05).mean()), "usage_share_at_least_0.95": float((u >= 0.95).mean()),
            "read_support_min_median": float(d.support_min.median())}


def metrics(df: pd.DataFrame) -> dict:
    matched = df[df.in_both]
    nonzero = matched[matched.has_read_support_rep1 & matched.has_read_support_rep2]
    primary = df[df.in_primary_qc]
    flagged = df[df.also_transcript_start]
    flagged_nz = nonzero[nonzero.also_transcript_start]
    clean_nz = nonzero[~nonzero.also_transcript_start]
    return {
        "counts": {
            "total_internal_acceptors_union": len(df), "in_rep1": int(df.in_rep1.sum()), "in_rep2": int(df.in_rep2.sum()),
            "only_rep1": int((df.in_rep1 & ~df.in_rep2).sum()), "only_rep2": int((df.in_rep2 & ~df.in_rep1).sum()),
            "matched": len(matched),
            "zero_support_rep1": int((df.read_support_rep1 == 0).sum()), "zero_support_rep2": int((df.read_support_rep2 == 0).sum()),
            "zero_support_in_matched": {"rep1": int((~matched.has_read_support_rep1).sum()), "rep2": int((~matched.has_read_support_rep2).sum()),
                                        "either": int((~(matched.has_read_support_rep1 & matched.has_read_support_rep2)).sum())},
            "nonzero_support_in_both": len(nonzero),
            "also_transcript_start_flagged": {"all_sites": len(flagged), "matched": int(matched.also_transcript_start.sum()), "nonzero_in_both": len(flagged_nz)},
            "clean_primary_qc_sites": len(primary),
            "per_chrom": {k: {"total": int(v), "matched": int((matched.chrom == k).sum()), "primary": int((primary.chrom == k).sum())}
                          for k, v in df.chrom.value_counts().sort_index().items()},
        },
        "agreement": {"matched_all": agreement(matched), "matched_nonzero_both": agreement(nonzero), "primary": agreement(primary)},
        "ambiguous_vs_clean_nonzero_both": {"also_transcript_start": group_profile(flagged_nz), "not_flagged": group_profile(clean_nz)},
        "usage_distribution_primary": {"rep1": usage_counts(primary.usage_rep1), "rep2": usage_counts(primary.usage_rep2)},
        "read_support_primary_min_of_both": describe(primary.support_min),
        "primary_sites_in_genes_with_incompatible_read_warnings": int(((primary.n_incompat_warnings_rep1 > 0) | (primary.n_incompat_warnings_rep2 > 0)).sum()),
        "bootstrap_sd_available": {"rep1": int(df.bootstrap_sd_rep1.notna().sum()), "rep2": int(df.bootstrap_sd_rep2.notna().sum())},
    }


def _save(fig, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    return path


def plots(df: pd.DataFrame, fig_dir: Path, tag: str, boot: pd.DataFrame | None) -> list[Path]:
    primary = df[df.in_primary_qc]
    matched_clean = df[df.in_both & ~df.also_transcript_start]
    zero_grp = matched_clean[matched_clean.support_min == 0]
    out = []
    a = agreement(primary)

    # A. replicate scatter
    fig, ax = new_figure(figsize=(5.6, 5.6))
    ax.plot([0, 1], [0, 1], color=INK_2, lw=0.8, ls="--", zorder=1)
    ax.scatter(primary.usage_rep1, primary.usage_rep2, s=3, color=ACCENT, alpha=0.12, linewidths=0, zorder=2, rasterized=True)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.set_xlabel("usage, replicate 1")
    ax.set_ylabel("usage, replicate 2")
    ax.set_title(f"Internal acceptor usage, K562 replicates ({tag})", fontsize=11, loc="left")
    ax.text(0.03, 0.97, f"n = {a['n']:,} sites (primary set)\nPearson r = {a['pearson']:.3f}\nSpearman rho = {a['spearman']:.3f}\n"
            f"median |diff| = {a['median_abs_diff']:.3f}\ndashed: y = x", transform=ax.transAxes, va="top", fontsize=9, color=INK)
    out.append(_save(fig, fig_dir / f"qc_A_replicate_scatter_{tag}.png"))

    # B. usage distribution
    fig, axes = new_figure(1, 2, figsize=(8.6, 3.6), sharey=True)
    for ax, rep in zip(axes, ("rep1", "rep2")):
        ax.hist(primary[f"usage_{rep}"], bins=40, range=(0, 1), color=ACCENT, edgecolor="#fcfcfb", linewidth=0.4)
        ax.set_xlabel("usage")
        ax.set_title(f"Replicate {rep[-1]}", fontsize=10, loc="left")
    axes[0].set_ylabel("number of sites")
    fig.suptitle(f"Usage of the primary-set acceptors, K562 ({tag}, n = {len(primary):,})", fontsize=11, x=0.01, ha="left", color=INK)
    out.append(_save(fig, fig_dir / f"qc_B_usage_distribution_{tag}.png"))

    # C. support vs disagreement
    fig, ax = new_figure(figsize=(7.0, 4.6))
    ax.scatter(as_plot_support(zero_grp.support_min), zero_grp.abs_diff, s=3, color=GREY, alpha=0.12, linewidths=0, rasterized=True)
    ax.scatter(as_plot_support(primary.support_min), primary.abs_diff, s=3, color=ACCENT, alpha=0.12, linewidths=0, rasterized=True)
    mids, med, p90 = [], [], []
    for lo, hi in pairwise(SUPPORT_EDGES):
        d = primary[(primary.support_min >= lo) & (primary.support_min < hi)]
        if len(d) >= 20:
            mids.append(np.sqrt(max(lo, 0.5) * min(hi, 4096)))
            med.append(d.abs_diff.median())
            p90.append(d.abs_diff.quantile(0.9))
    ax.plot(mids, p90, color=INK_2, lw=1.3, ls="--", marker="s", ms=3.5)
    ax.plot(mids, med, color=INK, lw=2, marker="o", ms=4.5)
    log_support_axis(ax)
    ax.set_xlabel("read support (the smaller of the two replicates)")
    ax.set_ylabel("|usage rep 1 - usage rep 2|")
    ax.set_title(f"Replicate disagreement vs read support ({tag})", fontsize=11, loc="left")
    ax.text(0.99, 0.97, "blue: primary set\ngrey at 0: no reads in at least one replicate\nsolid / dashed: median / 90th percentile\n(support bins with >= 20 sites)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK_2)
    out.append(_save(fig, fig_dir / f"qc_C_support_vs_disagreement_{tag}.png"))

    # D. counts by support bin
    labels = [support_bin_label(lo, hi) for lo, hi in pairwise(SUPPORT_EDGES)]
    counts = [int(((matched_clean.support_min >= lo) & (matched_clean.support_min < hi)).sum()) for lo, hi in pairwise(SUPPORT_EDGES)]
    fig, ax = new_figure(figsize=(7.6, 4.0))
    x = np.arange(len(labels))
    ax.bar(x, counts, color=[GREY] + [ACCENT] * (len(labels) - 1), width=0.72)
    for xi, c in zip(x, counts):
        ax.text(xi, c, f"{c:,}", ha="center", va="bottom", fontsize=7.5, color=INK_2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8.5)
    ax.set_xlabel("read support (the smaller of the two replicates)")
    ax.set_ylabel("matched sites")
    ax.set_title(f"Matched acceptors by read-support bin ({tag})", fontsize=11, loc="left")
    ax.text(0.99, 0.97, "grey: no reads in at least one replicate (removed)\nblue: primary set", transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK_2)
    ax.set_ylim(top=max(counts) * 1.14)
    out.append(_save(fig, fig_dir / f"qc_D_sites_by_support_bin_{tag}.png"))

    # E. bootstrap SD vs support, only where bootstraps exist
    if boot is not None:
        b = boot[boot.in_primary_qc].dropna(subset=["sd_mean"])
        if len(b):
            fig, ax = new_figure(figsize=(6.4, 4.4))
            ax.scatter(as_plot_support(b.support_min), b.sd_mean, s=22, color=ACCENT, alpha=0.8, linewidths=0)
            log_support_axis(ax)
            ax.set_xlabel("read support (the smaller of the two replicates)")
            ax.set_ylabel("bootstrap SD (mean of the two replicates)")
            ax.set_title("Bootstrap uncertainty vs read support (exploratory)", fontsize=11, loc="left")
            ax.text(0.98, 0.97, f"only the 60 kb test region on chr1,\nwhere 10 bootstraps were run\n(n = {len(b)} primary-set sites)",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK_2)
            out.append(_save(fig, fig_dir / f"qc_E_bootstrap_sd_vs_support_{tag}.png"))
    return out


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", default=str(project_path(cfg, "data_processed") / "acceptor_labels.tsv"))
    ap.add_argument("--tag", required=True, help="used in output file names, e.g. chr1-4")
    ap.add_argument("--bootstrap-table", help="a label table that has bootstrap SDs (the 60 kb region run), for plot E")
    a = ap.parse_args()
    fig_dir, res_dir = project_path(cfg, "figures"), project_path(cfg, "results")
    fig_dir.mkdir(exist_ok=True)
    res_dir.mkdir(exist_ok=True)

    df = load_labels(a.labels)
    m = metrics(df)
    m["labels_file"] = a.labels
    by_bin, _ = tables(df[df.in_both & ~df.also_transcript_start])
    by_bin.to_csv(res_dir / f"support_bins_{a.tag}.tsv", sep="\t", index=False, float_format="%.4g")
    cols = [c for c in ["support_min", "n_sites", "median_abs_diff", "p90_abs_diff", "share_diff_over_0_1", "pearson", "spearman"] if c in by_bin]
    m["agreement_by_support_bin"] = json.loads(by_bin[cols].to_json(orient="records"))
    figs = plots(df, fig_dir, a.tag, load_labels(a.bootstrap_table) if a.bootstrap_table else None)
    (res_dir / f"replicate_qc_{a.tag}.json").write_text(json.dumps(m, indent=1))

    c, ag, amb = m["counts"], m["agreement"], m["ambiguous_vs_clean_nonzero_both"]

    def fmt(d: dict) -> str:
        return f"n={d['n']:,}, Pearson {d.get('pearson', float('nan')):.3f}, Spearman {d.get('spearman', float('nan')):.3f}, median |diff| {d.get('median_abs_diff', float('nan')):.4f}"

    def prof(d: dict) -> str:
        return (f"{fmt(d)}; usage median {d['usage_median']:.3f}, share<0.05 {d['usage_share_below_0.05']:.3f}, "
                f"share>=0.95 {d['usage_share_at_least_0.95']:.3f}, median support {d['read_support_min_median']:.0f}")

    lines = [
        f"total internal acceptors (union of both replicates): {c['total_internal_acceptors_union']:,}  (rep1 {c['in_rep1']:,}, rep2 {c['in_rep2']:,}; only rep1 {c['only_rep1']:,}, only rep2 {c['only_rep2']:,})",
        f"matched between replicates: {c['matched']:,}",
        f"zero-support sites: rep1 {c['zero_support_rep1']:,}, rep2 {c['zero_support_rep2']:,}; among matched: {c['zero_support_in_matched']}",
        f"matched with nonzero support in both: {c['nonzero_support_in_both']:,}",
        f"also_transcript_start flagged: {c['also_transcript_start_flagged']}",
        f"clean sites used for QC (primary): {c['clean_primary_qc_sites']:,}",
        f"agreement, all matched:                {fmt(ag['matched_all'])}",
        f"agreement, nonzero support in both:    {fmt(ag['matched_nonzero_both'])}",
        f"agreement, primary (also not flagged): {fmt(ag['primary'])}",
        f"flagged, nonzero in both:   {prof(amb['also_transcript_start'])}",
        f"unflagged, nonzero in both: {prof(amb['not_flagged'])}",
        f"figures: {[str(p.relative_to(REPO_ROOT)) for p in figs]}",
    ]
    (res_dir / f"replicate_qc_{a.tag}.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(by_bin[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
