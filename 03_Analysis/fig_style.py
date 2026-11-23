"""Unified publication figure style for the SMVT manuscript.

Single source of truth for fonts, sizes, palette, and export. Import and call
apply_style() at the top of every figure script, then save with save_pub():

    from fig_style import apply_style, PALETTE, save_pub, add_panel_label
    apply_style()
    ...
    save_pub(fig, "Fig1_SMVT_TCGA_pan_cancer")

Target: npj Systems Biology and Applications (Nature Portfolio) — 89/127/183 mm
widths, Arial/Helvetica, nothing below 5 pt, line art >=600 dpi, TIFF compressed,
vector PDF with embedded TrueType text. Semantics (data, colours' meaning, chart
types) are deliberately unchanged; this module only governs presentation.
"""
import os
import matplotlib
import matplotlib.pyplot as plt

MM = 1 / 25.4  # millimetre -> inch

# Nature Portfolio content widths
SINGLE = 89 * MM      # 89 mm
ONE_HALF = 127 * MM   # 127 mm
DOUBLE = 183 * MM     # 183 mm

# Colorblind-safe semantic palette (Brewer-inspired, blue/red/green/neutral).
# Meanings are fixed manuscript-wide; do not re-assign per figure.
PALETTE = {
    "blue":   "#2166AC",   # main / high
    "red":    "#B2182B",   # contrast / low / negative control
    "green":  "#4DAF4A",   # positive / enrichment
    "purple": "#984EA3",
    "orange": "#FF7F00",
    "gray":   "#999999",
    "navy":   "#0F4D92",
    "crimson": "#B64342",
}

# Sequential gradient for heatmaps (colorblind-safe blue->gray->red)
CMAP = "RdBu_r"

# Redundant (non-colour) encodings so panels survive greyscale printing and
# colour-vision deficiency: pair hatch/markers with colour, never colour alone.
HATCH = ["", "///", "...", "xxx", "\\\\\\", "|||"]
MARKERS = ["o", "s", "^", "D", "v", "P"]

MIN_FONT = 5.0  # journal floor: no rendered text below this


def apply_style():
    """rcParams for a full-width-or-single Nature Portfolio figure page."""
    matplotlib.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "legend.title_fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.4,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "svg.fonttype": "none",   # editable text in SVG
        "pdf.fonttype": 42,       # embedded TrueType text in PDF
        "hatch.linewidth": 0.5,
        "axes.formatter.use_mathtext": True,
    })


def style_axis(ax, xlabel=None, ylabel=None):
    """Apply consistent spine/linewidth, remove top-right spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(width=0.6, length=2.5)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def add_panel_label(ax, label, x=-0.08, y=1.02, size=8):
    """Bold lowercase panel letter at the panel's upper-left (Nature style)."""
    ax.text(x, y, label.lower(), transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="top", ha="right")


def save_pub(fig, stem, outdir=None, formats=("pdf", "svg", "png", "tiff"),
             preview_dpi=300, line_art_dpi=600):
    """Write the vector masters + raster required by Nature Portfolio venues.

    pdf  : vector, embedded TrueType (text stays editable/selectable)
    svg  : vector, live text
    png  : 300 dpi preview for the manuscript draft
    tiff : 600 dpi, LZW-compressed (line-art requirement)
    Returns the list of files written.
    """
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        stem = os.path.join(outdir, os.path.basename(stem))
    written = []
    for ext in formats:
        path = f"{stem}.{ext}"
        # pdf.fonttype / svg.fonttype come from rcParams (set by apply_style);
        # savefig does not take rcParam keys as kwargs.
        if ext == "tiff":
            fig.savefig(path, dpi=line_art_dpi, bbox_inches="tight",
                        pil_kwargs={"compression": "tiff_lzw"})
        elif ext == "svg":
            fig.savefig(path, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight",
                        dpi=line_art_dpi if ext == "pdf" else preview_dpi)
        written.append(path)
    return written


def check_text_floor(fig, floor=MIN_FONT, print_mm=183):
    """Return [(text, set_size, printed_size)] below the journal's 5 pt floor.

    A figure drawn wider than it will be printed shrinks its text by the same
    factor, so the printed size is what gets judged, not the one on screen.
    """
    width_in = fig.get_size_inches()[0]
    scale = (print_mm * MM) / width_in if width_in > 0 else 1.0
    bad = []
    for obj in fig.findobj(matplotlib.text.Text):
        txt = obj.get_text().strip()
        if not txt or not obj.get_visible():
            continue
        s = obj.get_fontsize() * scale
        if s < floor:
            bad.append((txt[:28], obj.get_fontsize(), round(s, 2)))
    return bad
