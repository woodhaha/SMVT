"""Regenerate paper figures (fig1 discovery/cross-validation, fig2 transport refutation,
fig3 chemotype) in nature-figure style. Reads real CSVs. Run with base python (matplotlib only)."""
import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "paper", "figures")
EXP = os.path.join(BASE, "experiments")

# nature-figure style
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.6, "axes.labelsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})
BLUE = "#0F4D92"; RED = "#B64342"; TEAL = "#42949E"; PURPLE = "#9A4D8E"; GRAY = "#767676"


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


# ---------- Fig 1: three-state cross-validation (discovery) ----------
ligands = ["COOH-bisCF3", "COOH-CF3", "COOH-OCF3", "COOH-Ac", "biotin"]
af2  = [-9.68, -9.04, -8.94, -8.84, -6.76]
va26 = [-8.55, -8.52, -8.13, -7.84, -6.54]
vd26 = [-8.37, -7.51, -7.32, -6.92, -6.01]

x = np.arange(len(ligands)); w = 0.26
fig, ax = plt.subplots(figsize=(5.4, 2.9))
b1 = ax.bar(x - w, af2,  w, color=BLUE,  label="AlphaFold (AF2)")
b2 = ax.bar(x,      va26, w, color=TEAL,  label="cryo-EM 26va (occluded)")
b3 = ax.bar(x + w,  vd26, w, color=RED,   label="cryo-EM 26vd (outward-open)")
for bars in (b1, b2, b3):
    for r in bars:
        ax.text(r.get_x() + r.get_width()/2, r.get_y() - 0.35, f"{r.get_height():.2f}",
                ha="center", va="top", fontsize=6.5)
ax.axhline(-6.76, color=GRAY, ls=":", lw=0.8)  # biotin substrate ref
ax.set_xticks(x); ax.set_xticklabels(ligands, rotation=18, ha="right")
ax.set_ylabel("Vina score (kcal/mol)"); ax.set_ylim(-10.6, -5.4)
ax.set_title("COOH-bisCF3 is top-ranked on all three receptor states", fontsize=8)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout(); save(fig, "fig1_three_state")


# ---------- Fig 2: transport refutation (linker size-intolerance) ----------
rows = {}
with open(os.path.join(EXP, "rq5_linker_design", "rq5_results.csv")) as f:
    for r in csv.DictReader(f):
        if r.get("ok") == "True":
            rows[r["name"]] = float(r["best_affinity"])

order = ["7-CC#C", "7-CH3", "base bisCF3 (ref -9.68)", "7-CH2OH", "7-CH2NH2",
         "7-CH2SH", "7-CH2N3", "7-CH2CH2NH2", "7-CH2CH2OH",
         "7-CH2OCH2CH2NH2", "7-CH2CH2CH2NH2"]
labels = ["7-ethynyl", "7-methyl", "7-H (parent)", "7-CH2OH", "7-CH2NH2",
          "7-CH2SH", "7-CH2N3", "7-(CH2)2NH2", "7-(CH2)2OH",
          "7-CH2O(CH2)2NH2", "7-(CH2)3NH2"]
vals = [rows[k] for k in order]
colors = [TEAL if k in ("7-CC#C", "7-CH3", "base bisCF3 (ref -9.68)")
          else (RED if v > -8.3 else GRAY) for k, v in zip(order, vals)]

fig, ax = plt.subplots(figsize=(5.6, 2.9))
y = np.arange(len(order))[::-1]
ax.barh(y, vals, color=colors, height=0.62)
ax.axvline(-9.68, color=GRAY, ls="--", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
for yi, v in zip(y, vals):
    ax.text(v + 0.05, yi, f"{v:.2f}", va="center", fontsize=6.5)
ax.set_xlabel("Vina score (AF2, kcal/mol)"); ax.set_xlim(-10.4, -7.0)
ax.set_title("7-position is size-intolerant: small groups hold, longer linkers collapse", fontsize=8)
ax.annotate("full SN-38 conjugate (763.7 Da, 22 Å)\n→ rejected on all three structures\n(AF2 +59.0, occluded +23.7, outward-open +3.2)",
            xy=(0.5, -0.62), xycoords="axes fraction", ha="center", va="top",
            fontsize=7, color=RED, fontweight="bold")
fig.subplots_adjust(bottom=0.28); save(fig, "fig2_transport_refutation")


# ---------- Fig 3: chemotype (Table 1) ----------
chem = [("Adamantane + CF3  ·  COOH-bisCF3", -9.68, RED),
        ("Adamantane + Ac (prior best)", -8.84, BLUE),
        ("Adamantane diCOOH", -8.48, BLUE),
        ("Aromatic amino-acid (1-naphthylalanine)", -7.74, BLUE),
        ("Alternative cage (cubane diCOOH)", -7.57, BLUE),
        ("biotin (substrate)", -6.76, GRAY)]
chem = chem[::-1]  # worst at top for horizontal bar

fig, ax = plt.subplots(figsize=(5.4, 2.5))
names = [c[0] for c in chem]; cv = [c[1] for c in chem]; cc = [c[2] for c in chem]
yy = np.arange(len(chem))
ax.barh(yy, cv, color=cc, height=0.6)
for yi, v in zip(yy, cv):
    ax.text(v - 0.12, yi, f"{v:.2f}", va="center", ha="right", fontsize=7, color="white")
ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=7)
ax.axvline(-8.84, color=TEAL, ls="--", lw=0.9)
ax.text(-8.84, len(chem)-0.25, "prior best −8.84", fontsize=6.5, color=TEAL, ha="right")
ax.set_xlabel("Vina score (AF2, kcal/mol)"); ax.set_xlim(-10.3, 0)
ax.set_title("The adamantane cage, not aromatics or alternative cages, reaches the pocket", fontsize=8)
fig.tight_layout(); save(fig, "fig3_chemotype")

# clean up the old mislabeled figure
old = os.path.join(FIG, "fig2_linker_overfit")
for ext in (".pdf", ".png"):
    if os.path.exists(old + ext):
        os.remove(old + ext)
        print("removed", old + ext)
