"""Unified publication figure style for the SMVT manuscript.

Single source of truth for fonts, sizes, and a colorblind-safe palette.
Import and call apply_style() at the top of every figure script:
    from fig_style import apply_style, PALETTE
    apply_style()
"""
import matplotlib
import matplotlib.pyplot as plt

# Colorblind-safe semantic palette (Brewer-inspired, blue/red/green/neutral)
PALETTE = {
    "blue":   "#2166AC",   # main / high
    "red":    "#B2182B",   # contrast / low / negative control
    "green":  "#4DAF4A",   # positive / enrichment
    "purple": "#984EA3",
    "orange": "#FF7F00",
    "gray":   "#999999",
    "navy":   "#0F4D92",
    "crimson":"#B64342",
}

# Sequential gradient for heatmaps (colorblind-safe blue->gray->red)
CMAP = "RdBu_r"


def apply_style():
    matplotlib.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "legend.title_fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "lines.linewidth": 1.2,
        "legend.frameon": False,
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,   # embed fonts as TrueType
    })


def style_axis(ax, xlabel=None, ylabel=None):
    """Apply consistent spine/linewidth, remove top-right spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(width=0.8)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def add_panel_label(ax, label, x=-0.12, y=1.08, size=11):
    """Bold panel label (a/b/c/d) in upper-left, publication style."""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="top", ha="right")
