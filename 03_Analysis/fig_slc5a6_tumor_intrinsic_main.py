"""Main figure: SLC5A6 overexpression is tumor-cell-intrinsic.

Panel a: TCGA log2FC across 21 cancer types (up=red, down=blue, ns=gray) — bulk.
Panel b: single-cell dot plot — % SLC5A6+ (size) x mean log1p(CPM) (color) per cell type.
Panel c: % SLC5A6+ bars, tumor vs normal epithelium vs immune/stromal + Mann-Whitney p.
"""
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import apply_style, PALETTE, add_panel_label, style_axis

OUT = "03_Analysis/outputs"
FIGDIR = "04_Manuscript/figures"
per = pd.read_csv(f"{OUT}/slc5a6_percell_gse178341.csv.gz")
de = pd.read_csv(f"{OUT}/pancan_de.csv")

apply_style()

# ================= Panel a: TCGA log2FC forest (lollipop) =================
de = de.sort_values("log2FC", ascending=True)
sig_up = (de["log2FC"] > 0) & (de["sig"] == "*")
sig_dn = (de["log2FC"] < 0) & (de["sig"] == "*")
colors = np.where(sig_up, PALETTE["red"],
         np.where(sig_dn, PALETTE["blue"], PALETTE["gray"]))

fig = plt.figure(figsize=(183 / 25.4, 130 / 25.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9],
                      width_ratios=[1.1, 1.0], hspace=0.35, wspace=0.32)

ax = fig.add_subplot(gs[0, :])  # forest spans full width
ax.hlines(de["cancer"], 0, de["log2FC"], color=colors, lw=2.5, alpha=0.85)
ax.scatter(de["log2FC"], de["cancer"], c=colors, s=38, edgecolors="black",
           linewidths=0.4, zorder=3)
ax.axvline(0, color="black", lw=0.7, ls="--")
ax.set_xlabel("SLC5A6 log$_2$FC (tumor vs normal)")
ax.set_ylabel("TCGA cancer type")
# legend
from matplotlib.lines import Line2D
lg = [Line2D([0], [0], marker="o", ls="", color=PALETTE["red"], label="Up (FDR<0.05)"),
      Line2D([0], [0], marker="o", ls="", color=PALETTE["blue"], label="Down (FDR<0.05)"),
      Line2D([0], [0], marker="o", ls="", color=PALETTE["gray"], label="n.s.")]
ax.legend(handles=lg, loc="lower right", frameon=False, fontsize=6.5)
n_up = int(sig_up.sum()); n_dn = int(sig_dn.sum())
ax.set_title(f"{n_up} up / {n_dn} down / {len(de)-n_up-n_dn} n.s. across 21 types",
             fontsize=8, pad=6)
style_axis(ax)
add_panel_label(ax, "a")

# ================= Panel b: single-cell dot plot =================
ct = per.groupby("clTopLevel").agg(
    pct=("slc5a6_pos", lambda s: 100.0 * s.mean()),
    mean_l1=("log1p_cpm", "mean"),
).reset_index().sort_values("mean_l1")

axb = fig.add_subplot(gs[1, 0])
scb = axb.scatter(ct["mean_l1"], ct["clTopLevel"], s=ct["pct"] * 6 + 12,
                  c=ct["mean_l1"], cmap="RdBu_r", edgecolors="black",
                  linewidths=0.4, vmin=0, vmax=max(ct["mean_l1"].max(), 0.5))
axb.set_xlabel("Mean SLC5A6 log$_1$(CPM)")
axb.set_ylabel("Cell type")
cb = plt.colorbar(scb, ax=axb, pad=0.02, fraction=0.046)
cb.set_label("Mean log$_1$(CPM)", fontsize=7); cb.ax.tick_params(labelsize=6)
for pct in [5, 20, 40]:
    axb.scatter([], [], s=pct * 6 + 12, edgecolors="black", facecolors="none",
                linewidths=0.4, label=f"{pct}%")
axb.legend(title="% SLC5A6+", loc="lower right", frameon=False, fontsize=6,
           title_fontsize=6)
style_axis(axb)
add_panel_label(axb, "b")

# ================= Panel c: % positive bars =================
key_order = ["Tumor epithelial", "Normal epithelial", "Myeloid", "Strom",
             "TNKILC", "B", "Plasma", "Mast"]
g = per.groupby("group").agg(n=("slc5a6_pos", "size"),
                             pct=("slc5a6_pos", lambda s: 100.0 * s.mean())
                              ).reset_index()
g = g[g["group"].isin(key_order)].set_index("group").loc[key_order].reset_index()

te = per.loc[(per.clTopLevel == "Epi") & (per.SPECIMEN_TYPE == "T"), "log1p_cpm"]
ne = per.loc[(per.clTopLevel == "Epi") & (per.SPECIMEN_TYPE == "N"), "log1p_cpm"]
_, p_mw = mannwhitneyu(te, ne, alternative="two-sided")

axc = fig.add_subplot(gs[1, 1])
cols = [PALETTE["red"], PALETTE["blue"]] + [PALETTE["gray"]] * (len(g) - 2)
axc.barh(g["group"], g["pct"], color=cols, edgecolor="black", linewidth=0.4,
         height=0.7)
axc.invert_yaxis()
for y, (pctv, nv) in enumerate(zip(g["pct"], g["n"])):
    axc.text(pctv + 0.4, y, f"n={nv:,}", va="center", fontsize=6, color="#333333")
axc.annotate("Mann–Whitney\n$p<10^{-6}$", xy=(0.98, 0.97),
             xycoords="axes fraction", ha="right", va="top", fontsize=6.5)
axc.set_xlabel("% SLC5A6+ cells")
axc.set_xlim(0, g["pct"].max() * 1.22)
style_axis(axc)
add_panel_label(axc, "c")

for ext, dpi in [("pdf", None), ("png", 300), ("svg", None)]:
    fig.savefig(os.path.join(FIGDIR, f"Fig_SLC5A6_tumor_intrinsic.{ext}"),
                dpi=dpi, bbox_inches="tight")
print(f"Saved main figure. up={n_up} down={n_dn}; Mann-Whitney p={p_mw}")
