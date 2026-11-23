"""FigS25 — UMAP of the GSE178341 CRC atlas showing SLC5A6 tumour-epithelial localisation.

Conclusion this figure must defend:
  SLC5A6 expression on the atlas-wide UMAP is confined to the epithelial compartment and
  is denser in the tumour-epithelial region than in normal epithelium, so the pan-cancer
  bulk RNA up-regulation is tumour-cell-intrinsic rather than stromal or immune.

Panels (one independent question each; all 370,115 cells, no subsampling):
  a  lineage colouring — establishes the atlas architecture.
  b  SLC5A6 log1p CPM over all cells — where the transcript actually sits.
  c  epithelial zoom split into tumour vs normal, positivity marked — the claim itself.

Coordinates come from umap_gse178341_compute.py (scanpy, HVG 2000, PCA 50,
n_neighbors 15, all cells). SLC5A6 values are read from the same per-cell table used for
FigS23, and this script re-asserts the manuscript's positivity numbers before drawing.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, r"D:\Researching\SMVT\03_Analysis")
from fig_style import (apply_style, PALETTE, save_pub, add_panel_label,
                       check_text_floor, DOUBLE)

COORD = r"D:\Researching\SMVT\03_Analysis\outputs\umap_gse178341\umap_coordinates.csv.gz"
OUTDIR = r"D:\Researching\SMVT\04_Manuscript\figures"
GROUP_COLOR = {
    "Tumor epithelial": PALETTE["red"], "Normal epithelial": PALETTE["orange"],
    "Strom": PALETTE["purple"], "Myeloid": PALETTE["green"], "TNKILC": PALETTE["blue"],
    "B": "#66C2A5", "Plasma": "#8DA0CB", "Mast": "#E78AC3",
}
apply_style()
df = pd.read_csv(COORD)
assert len(df) == 370115, f"expected 370,115 cells, got {len(df)}"
te = 100 * df.loc[df.group == "Tumor epithelial", "slc5a6_pos"].mean()
ne = 100 * df.loc[df.group == "Normal epithelial", "slc5a6_pos"].mean()
epi = 100 * df.loc[df.clTopLevel == "Epi", "slc5a6_pos"].mean()
assert abs(te - 28.25) < 0.05 and abs(ne - 6.95) < 0.05 and abs(epi - 21.28) < 0.05, \
    f"UMAP table does not reproduce FigS23: {te:.2f}/{ne:.2f}/{epi:.2f}"
print(f"consistency OK: tumor epi {te:.2f}% | normal epi {ne:.2f}% | epi {epi:.2f}%")

expr_cmap = LinearSegmentedColormap.from_list("slc", ["#F7F7F7", "#C6DBEF",
                                                      "#6BAED6", PALETTE["blue"], PALETTE["navy"]])
rng = np.random.default_rng(0)
order = rng.permutation(len(df))          # draw order only; every cell is plotted

fig = plt.figure(figsize=(DOUBLE, 68 / 25.4))
gs = fig.add_gridspec(1, 3, wspace=0.30)

axa = fig.add_subplot(gs[0, 0])
for g, c in GROUP_COLOR.items():
    s = df.group == g
    axa.scatter(df.UMAP_1[s], df.UMAP_2[s], s=0.35, c=c, linewidths=0, alpha=0.5,
                rasterized=True, label=f"{g} ({int(s.sum()):,})")
axa.set_title("a  lineage", loc="left")
axa.set_xlabel("UMAP_1"); axa.set_ylabel("UMAP_2")
# Direct labels at cluster centroids instead of a legend: the legend collided with
# panel b, and journal style prefers no eye-travel for spatially fixed categories.
SHORT = {"Tumor epithelial": "Tumour epi.", "Normal epithelial": "Normal epi.",
         "Strom": "Stroma", "Myeloid": "Myeloid", "TNKILC": "T/NK/ILC",
         "B": "B cells", "Plasma": "Plasma", "Mast": "Mast"}
cent = df.groupby("group")[["UMAP_1", "UMAP_2"]].mean()
for g, row in cent.iterrows():
    if g not in GROUP_COLOR:
        continue
    axa.annotate(SHORT[g], (row.UMAP_1, row.UMAP_2), fontsize=5.6,
                 ha="center", va="center", color="black",
                 bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="#CCCCCC", lw=0.4),
                 annotation_clip=False)

axb = fig.add_subplot(gs[0, 1])
sc = axb.scatter(df.UMAP_1[order], df.UMAP_2[order], s=0.35,
                 c=df.slc5a6_log1p[order], cmap=expr_cmap, linewidths=0,
                 vmin=0, vmax=0.6, rasterized=True)
axb.set_title("b  SLC5A6 log$_1$CPM", loc="left")
axb.set_xlabel("UMAP_1")
axb.text(0.03, 0.03, f"epithelial {epi:.1f}% positive", transform=axb.transAxes,
         fontsize=5.6, va="bottom",
         bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#CCCCCC", lw=0.4))
cb = fig.colorbar(sc, ax=axb, fraction=0.046, pad=0.03)
cb.set_ticks([0, 0.2, 0.4, 0.6]); cb.ax.tick_params(labelsize=5.5)

axc = fig.add_subplot(gs[0, 2])
e = df[df.clTopLevel == "Epi"]
for g, c in (("Normal epithelial", GROUP_COLOR["Normal epithelial"]),
             ("Tumor epithelial", GROUP_COLOR["Tumor epithelial"])):
    s = e.group == g
    axc.scatter(e.UMAP_1[s], e.UMAP_2[s], s=0.5, c=c, linewidths=0, alpha=0.45,
                rasterized=True, label=f"{g} (n={int(s.sum()):,})")
pos = e[e.slc5a6_pos]
# Unfilled rings, so the tumour/normal colour underneath stays visible: solid black
# markers swallowed exactly the contrast this panel exists to show.
axc.scatter(pos.UMAP_1, pos.UMAP_2, s=3.2, facecolors="none", edgecolors="k",
            linewidths=0.35, alpha=0.55, rasterized=True,
            label=f"SLC5A6-positive (n={len(pos):,})")
axc.text(0.03, 0.03, f"positive: tumour epi {te:.1f}%  vs  normal epi {ne:.1f}%",
         transform=axc.transAxes, fontsize=5.6, va="bottom",
         bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#CCCCCC", lw=0.4))
axc.set_title("c  epithelium, tumour vs normal", loc="left")
axc.set_xlabel("UMAP_1")
axc.legend(fontsize=5.6, loc="upper left", frameon=True, framealpha=0.85,
                   edgecolor="#DDDDDD", handletextpad=0.2,
           labelspacing=0.3)
for ax in fig.axes:
    ax.set_xticks(ax.get_xticks()[:4]); ax.tick_params(labelsize=6)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

low = check_text_floor(fig, print_mm=183)
stem = os.path.join(OUTDIR, "FigS25_SLC5A6_umap")
files = save_pub(fig, stem)
print("wrote:", [os.path.basename(f) for f in files])
print("cells plotted:", len(df), "| epithelium panel:", len(e),
      "| positive in epi:", int(e.slc5a6_pos.sum()))
print("sub-5pt text at 183 mm:", low if low else "none")
plt.close(fig)
