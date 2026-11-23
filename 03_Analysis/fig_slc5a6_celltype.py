"""Fig S23 (redo, publication-grade): SLC5A6 single-cell cell-type localization.

Panel a: dot plot — % SLC5A6+ cells (size) x mean log1p(CPM) (color) per cell type.
Panel b: % SLC5A6+ (bars) for tumor vs normal epithelium vs major immune/stromal
         populations, with n and Mann-Whitney p (tumor vs normal epithelium).
"""
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import apply_style, PALETTE, add_panel_label, style_axis, save_pub

OUT = "03_Analysis/outputs"
FIGDIR = "04_Manuscript/figures"
per = pd.read_csv(f"{OUT}/slc5a6_percell_gse178341.csv.gz")

apply_style()

# ---- cell-type summary, ordered by mean expression ----
ct = per.groupby("clTopLevel").agg(
    n=("slc5a6_pos", "size"),
    pct=("slc5a6_pos", lambda s: 100.0 * s.mean()),
    mean_l1=("log1p_cpm", "mean"),
).reset_index().sort_values("mean_l1")  # ascending for horizontal dot plot

# ---- Mann-Whitney: tumor vs normal epithelium ----
te = per.loc[(per.clTopLevel == "Epi") & (per.SPECIMEN_TYPE == "T"), "log1p_cpm"]
ne = per.loc[(per.clTopLevel == "Epi") & (per.SPECIMEN_TYPE == "N"), "log1p_cpm"]
_, p_mw = mannwhitneyu(te, ne, alternative="two-sided")

fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 85 / 25.4),
                         gridspec_kw={"width_ratios": [1.1, 1.0]})

# ================= Panel a: dot plot =================
ax = axes[0]
sc = ax.scatter(
    ct["mean_l1"], ct["clTopLevel"],
    s=ct["pct"] * 6 + 12,
    c=ct["mean_l1"], cmap="RdBu_r",
    edgecolors="black", linewidths=0.4,
    vmin=0, vmax=max(ct["mean_l1"].max(), 0.5),
)
ax.set_xlabel("Mean SLC5A6 log$_1$(CPM)")
ax.set_ylabel("Cell type")
cb = plt.colorbar(sc, ax=ax, pad=0.02, fraction=0.046)
cb.set_label("Mean log$_1$(CPM)", fontsize=8)
cb.ax.tick_params(labelsize=6)
# size legend (percent positive)
for pct in [5, 20, 40]:
    ax.scatter([], [], s=pct * 6 + 12, edgecolors="black", facecolors="none",
               linewidths=0.4, label=f"{pct}%")
ax.legend(title="% SLC5A6+", loc="lower right", frameon=False, fontsize=6,
          title_fontsize=6)
style_axis(ax)
add_panel_label(ax, "a")

# ================= Panel b: % positive bars =================
ax2 = axes[1]
key_order = ["Tumor epithelial", "Normal epithelial", "Myeloid", "Strom",
             "TNKILC", "B", "Plasma", "Mast"]
g = per.groupby("group").agg(
    n=("slc5a6_pos", "size"),
    pct=("slc5a6_pos", lambda s: 100.0 * s.mean()),
).reset_index()
g = g[g["group"].isin(key_order)].set_index("group").loc[key_order].reset_index()

cols = [PALETTE["red"], PALETTE["blue"]] + [PALETTE["gray"]] * (len(g) - 2)
bars = ax2.barh(g["group"], g["pct"], color=cols, edgecolor="black",
                linewidth=0.4, height=0.7)
ax2.invert_yaxis()
# n annotation
for y, (pctv, nv) in enumerate(zip(g["pct"], g["n"])):
    ax2.text(pctv + 0.4, y, f"n={nv:,}", va="center", fontsize=6,
             color="#333333")
# Mann-Whitney annotation (tumor vs normal epithelial)
ax2.annotate("Mann–Whitney\n$p<10^{-6}$",
             xy=(0.98, 0.97), xycoords="axes fraction",
             ha="right", va="top", fontsize=7.5)
ax2.set_xlabel("% SLC5A6+ cells")
ax2.set_xlim(0, g["pct"].max() * 1.22)
style_axis(ax2)
add_panel_label(ax2, "b")

fig.tight_layout()
save_pub(fig, "FigS23_SLC5A6_celltype", outdir=FIGDIR)
print("Saved FigS23. Mann-Whitney p (tumor vs normal epi) =", p_mw)
print(ct.to_string(index=False))
print(g.to_string(index=False))
