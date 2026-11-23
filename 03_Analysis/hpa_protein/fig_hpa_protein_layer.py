"""FigS24 — protein-level evidence for SLC5A6 and the biotin-ACC-FASN axis.

Conclusion this figure must defend (one sentence):
  Independent mass-spectrometry proteomics confirms SLC5A6 protein up-regulation in
  tumour versus adjacent normal colon and lung squamous epithelium, while its
  normal-tissue protein signal is sparse because the transporter sits at the MS/IHC
  detection floor — so the pan-cancer RNA signal survives at the protein level in
  tumour, and normal-tissue silence is a detection limit, not a contradiction.

Panels (one independent question each):
  a  CPTAC tumour-vs-normal protein log2FC for SLC5A6 across 10 cancers; FDR<0.05
     marked. The positive replication.
  b  HPA normal-tissue IHC best level per tissue (45 tissues), as a colour-coded tile
     strip with rotated labels, so every tissue is readable at print size.
  c  MS quantification completeness: SLC5A6 0/20 tissues against axis enzymes that
     are quantified 12-20/20 — the detection-floor evidence behind panel b.

Source: Human Protein Atlas v26.0 bulk TSV (Ensembl 109), downloaded 2026-11-23 and
subset by 02_Data/HPA/extract/*.tsv. HPA's ordinal levels are used verbatim; no
thresholds are invented here.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, r"D:\Researching\SMVT\03_Analysis")
from fig_style import (apply_style, PALETTE, save_pub, add_panel_label,
                       check_text_floor, DOUBLE)

HPA = r"D:\Researching\SMVT\02_Data\HPA\extract"
OUTDIR = r"D:\Researching\SMVT\04_Manuscript\figures"
SLC5A6 = "ENSG00000138074"
LEVEL_ORDER = ["Not detected", "Low", "Medium", "High"]
LEVEL_COLOR = {"Not detected": "#F4F4F4", "Low": "#9ECAE1",
               "Medium": PALETTE["blue"], "High": PALETTE["navy"]}
apply_style()

cp = pd.read_csv(os.path.join(HPA, "cancer_cptac_axis.tsv"), sep="\t")
cp["p"] = pd.to_numeric(cp["pval_adj"], errors="coerce")
cp["logFC"] = pd.to_numeric(cp["logFC"], errors="coerce")
ihc = pd.read_csv(os.path.join(HPA, "normal_ihc_data_axis.tsv"), sep="\t")
ms = pd.read_csv(os.path.join(HPA, "ms_tissue_axis.tsv"), sep="\t")
ms["det"] = ms["Intensity"].notna()

fig = plt.figure(figsize=(DOUBLE, 158 / 25.4))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.05], hspace=0.30, wspace=0.85)

# ---- a: CPTAC forest, SLC5A6 across cancers ----------------------------------
axa = fig.add_subplot(gs[0, 0:4])
slc = cp[cp.Gene == SLC5A6].sort_values("logFC")
sig = slc.p < 0.05
y = np.arange(len(slc))
axa.axvline(0, color="k", lw=0.8)
axa.scatter(slc.logFC[~sig], y[~sig], s=26, facecolor="white", edgecolor=PALETTE["gray"],
            linewidth=0.9, zorder=3, label="not significant")
axa.scatter(slc.logFC[sig], y[sig], s=34, color=PALETTE["red"], edgecolor="k",
            linewidth=0.5, zorder=4, label="FDR < 0.05")
axa.set_yticks(y)
axa.set_yticklabels(slc.Cancer, fontsize=6.5)
axa.set_ylim(-0.7, len(slc) - 0.3)
lo, hi = slc.logFC.min(), slc.logFC.max()
axa.set_xlim(lo - 0.35, hi + 0.42)
axa.set_xlabel(r"protein log$_2$(tumour / adjacent normal), CPTAC")
axa.set_title("a  SLC5A6 protein in tumour vs normal", loc="left", fontsize=8)
axa.legend(fontsize=5.5, loc="lower right", ncol=1, handlelength=1.0)
for i, r in zip(y, slc.itertuples()):
    star = r.p < 0.05
    axa.text(r.logFC + (0.05 if r.logFC >= 0 else -0.05), i, f"{r.logFC:+.2f}",
             va="center", ha="left" if r.logFC >= 0 else "right", fontsize=5.5,
             color=PALETTE["red"] if star else PALETTE["gray"],
             fontweight="bold" if star else "normal")

# ---- c: MS detection floor ---------------------------------------------------
axc = fig.add_subplot(gs[0, 4:6])
tab = ms.groupby("Gene name")["det"].agg(["sum", "size"]).sort_values("sum")
pct = tab["sum"] / tab["size"] * 100
colors = [PALETTE["red"] if g == "SLC5A6" else PALETTE["gray"] for g in tab.index]
axc.barh(np.arange(len(tab)), pct, color=colors, edgecolor="k", linewidth=0.4, height=0.72)
axc.set_yticks(np.arange(len(tab)))
axc.set_yticklabels(tab.index, fontsize=6.5)
axc.set_xlim(0, 118)
axc.set_xticks([0, 50, 100])
axc.set_xlabel("% of 20 tissues quantified")
axc.set_title("c  MS detection floor", loc="left", fontsize=8)
for i, r in zip(np.arange(len(tab)), tab.itertuples()):
    inside = r.sum / r.size > 0.15
    axc.text((r.sum / r.size * 100) - 4 if inside else 2, i, f"{r.sum}/{r.size}",
             va="center", ha="right" if inside else "left", fontsize=5.5,
             color="white" if inside else "black",
             fontweight="bold" if r.Index == "SLC5A6" else "normal")

# ---- b: IHC tile strip, 45 tissues ------------------------------------------
axb = fig.add_subplot(gs[1, :])
s = ihc[ihc.Gene == SLC5A6].copy()
s["rank"] = s.Level.map({l: i for i, l in enumerate(LEVEL_ORDER)})
best = s.loc[s.groupby("Tissue")["rank"].idxmax(), ["Tissue", "Level"]].copy()
best = best.sort_values(["Level", "Tissue"],
                        key=lambda c: c.map({l: i for i, l in enumerate(LEVEL_ORDER)})
                        if c.name == "Level" else c.str.lower(),
                        ascending=[False, True]).reset_index(drop=True)
x = np.arange(len(best))
for i, r in enumerate(best.itertuples()):
    axb.bar(i, 1, color=LEVEL_COLOR[r.Level], edgecolor="#7A7A7A", linewidth=0.35, width=0.92)
axb.set_xlim(-0.6, len(best) - 0.4)
axb.set_ylim(0, 1)
axb.set_yticks([0.5])
axb.set_yticklabels(["HPA level"], fontsize=6.5)
axb.set_xticks(x)
axb.set_xticklabels(best.Tissue.str.title(), rotation=90, fontsize=5.5, ha="center")
axb.tick_params(axis="x", length=0, pad=1.5)
axb.spines["top"].set_visible(False)
axb.spines["right"].set_visible(False)
axb.spines["left"].set_visible(False)
axb.set_title(f"b  Normal-tissue IHC, best level per tissue (n={len(best)}; "
              f"High in 0, antibody HPA036958, Reliability Approved)",
              loc="left", fontsize=8)
axb.legend(handles=[Patch(facecolor=LEVEL_COLOR[l], edgecolor="#7A7A7A",
                          label=f"{l}: {(best.Level == l).sum()}") for l in LEVEL_ORDER],
           fontsize=5.5, loc="upper right", ncol=4, handlelength=1.1,
           columnspacing=0.9, borderpad=0.2)

low = check_text_floor(fig, print_mm=183)
stem = os.path.join(OUTDIR, "FigS24_HPA_protein_layer")
files = save_pub(fig, stem)
print("wrote:", [os.path.basename(f) for f in files])
counts = {l: int((best.Level == l).sum()) for l in LEVEL_ORDER}
print("IHC tissue-level counts:", counts, "| total tissues:", len(best))
print("CPTAC SLC5A6:", len(slc), "cancers,", int(sig.sum()), "FDR<0.05",
      "| colon logFC", float(slc.loc[slc.Cancer == "Colon AC", "logFC"].iloc[0]),
      "p", float(slc.loc[slc.Cancer == "Colon AC", "p"].iloc[0]))
print("MS:", {r.Index: f"{r.sum}/{r.size}" for r in tab.itertuples()})
print("sub-5pt text at 183 mm:", low if low else "none")
plt.close(fig)
