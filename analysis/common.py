"""Shared helpers for the label QC: loading the label table, plot style, support bins."""

from __future__ import annotations

import numpy as np
import pandas as pd

# Plot colours: one accent hue, neutral greys (light surface).
ACCENT = "#2a78d6"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

# Support bins: exactly zero, then 0<s<1 (multi-mapped reads give fractional counts), then powers of two.
# Powers of two are neutral on purpose: they are not a pre-chosen cutoff.
ZERO = 1e-9
SUPPORT_EDGES = [0, ZERO, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, np.inf]


def load_labels(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "strand": str})
    df["support_min"] = df[["read_support_rep1", "read_support_rep2"]].min(axis=1)  # a site is as reliable as its weaker replicate
    df["abs_diff"] = (df["usage_rep1"] - df["usage_rep2"]).abs()
    df["sd_mean"] = df[["bootstrap_sd_rep1", "bootstrap_sd_rep2"]].mean(axis=1)
    # what bootstrap noise alone would give as the median |difference| between replicates
    df["expected_abs_diff"] = 0.6745 * np.sqrt(df["bootstrap_sd_rep1"] ** 2 + df["bootstrap_sd_rep2"] ** 2)
    return df


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)


def new_figure(*args, **kw):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(*args, **kw)
    fig.patch.set_facecolor(SURFACE)
    for ax in np.atleast_1d(axes).ravel():
        style_axes(ax)
    return fig, axes


def log_support_axis(ax) -> None:
    """Log x axis for read support, with 0 drawn at 0.5."""
    ax.set_xscale("log")
    ticks = [0.5, 1, 10, 100, 1000, 10000]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0", "1", "10", "100", "1k", "10k"])
    ax.minorticks_off()


def as_plot_support(s: pd.Series) -> pd.Series:
    return s.clip(lower=0.5)


def support_bin_label(lo: float, hi: float) -> str:
    if hi <= ZERO:
        return "0"
    if lo <= ZERO:
        return "(0,1)"
    if np.isinf(hi):
        return f">={int(lo)}"
    return f"{int(lo)}" if hi - lo == 1 else f"{int(lo)}-{int(hi) - 1}"
