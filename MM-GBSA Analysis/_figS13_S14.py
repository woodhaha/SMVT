"""
FigS13 + FigS14 — Reproducible vector PDF generation for SMVT manuscript.

FigS13: GB model comparison (OBC1 vs OBC2) — grouped bar chart
FigS14: GB model cross-validation (OBC1 vs OBC2) — scatter + Spearman

Usage:
    conda run -n smvt-md python3 _figS13_S14.py

Requires:
    - mmgbsa_results.json per compound (OBC2, from _run_mmgbsa.py)
    - obc1_results.json (OBC1, from _run_obc1.py)
"""
import json, os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.titlesize": 11,
                     "figure.dpi": 300, "savefig.dpi": 300, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial"]})

ANALYSIS = r"D:\Researching\SMVT\MM-GBSA Analysis"
FIGURES = r"D:\Researching\SMVT\04_Manuscript\figures"
COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
PALETTE = ["#e67e22","#3498db","#1abc9c","#e74c3c",
           "#2ecc71","#9b59b6","#f39c12","#95a5a6"]
LABELS = {"NAFTAZONE":"Naftazone","BIOTIN":"Biotin","ESKETAMINE":"Esketamine",
          "FUROSEMIDE":"Furosemide","GABAPENTIN_ENACARBIL":"Gabapentin Enacarbil",
          "HYDROMORPHONE":"Hydromorphone","PHENOBARBITAL":"Phenobarbital","RIBOFLAVIN":"Riboflavin"}
SHORT = {"NAFTAZONE":"NAF","BIOTIN":"BIO","ESKETAMINE":"ESK","FUROSEMIDE":"FUR",
         "GABAPENTIN_ENACARBIL":"GAB","HYDROMORPHONE":"HYD","PHENOBARBITAL":"PHE","RIBOFLAVIN":"RIB"}

# ── Load OBC2 data ──
obc2 = {}
for c in COMPS:
    fn = os.path.join(ANALYSIS, c, "mmgbsa_results.json")
    if os.path.exists(fn):
        obc2[c] = json.load(open(fn))["dG_kcal"]

# ── Load OBC1 data ──
obc1_fn = os.path.join(ANALYSIS, "obc1_results.json")
obc1 = json.load(open(obc1_fn)) if os.path.exists(obc1_fn) else {}

valid = [c for c in COMPS if c in obc1 and c in obc2]
print(f"Loaded OBC1={len(obc1)} OBC2={len(obc2)} valid={len(valid)}")

# ════════════════════════════════════════════════════════════════
# FIGURE S13: GB Model Comparison — OBC1 vs OBC2 Grouped Bar
# ════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle("GB Model Comparison: OBC1 vs OBC2 (SMVT-8 Compounds, 100ns MD)",
             fontsize=12, fontweight="bold", y=0.98)

# Panel A: Grouped bar chart
ax = ax1
x = np.arange(len(valid))
w = 0.32
colors = [PALETTE[COMPS.index(c)] for c in valid]
b1 = ax.bar(x - w/2, [obc1[c] for c in valid], w, color=colors, alpha=0.8, edgecolor="white", linewidth=0.5, label="GB(OBC1)")
b2 = ax.bar(x + w/2, [obc2[c] for c in valid], w, color=colors, alpha=0.4, edgecolor="white", linewidth=0.5, hatch="///", label="GB(OBC2)")
ax.set_xticks(x)
ax.set_xticklabels([SHORT[c] for c in valid], fontsize=8, rotation=20)
ax.set_ylabel("ΔG (kcal/mol)")
ax.set_title("a. Absolute binding free energies", fontsize=10, loc="left")
ax.axhline(y=0, color="#888", linewidth=0.5)
ax.legend(fontsize=7, ncol=2)
ax.grid(axis="y", alpha=0.2)

# Annotate bars
for i, c in enumerate(valid):
    v1, v2 = obc1[c], obc2[c]
    ax.text(i - w/2, v1 + (0.3 if v1 > 0 else -1.5), f"{v1:.0f}", ha="center", fontsize=6, color="black")
    ax.text(i + w/2, v2 + (0.3 if v2 > 0 else -1.5), f"{v2:.0f}", ha="center", fontsize=6, color="black")

# Panel B: Rank correlation table/plot
ax = ax2
obc1_rank = {c: i+1 for i, (c,_) in enumerate(sorted(obc1.items(), key=lambda x: x[1]))}
obc2_rank = {c: i+1 for i, (c,_) in enumerate(sorted(obc2.items(), key=lambda x: x[1]))}
r_vals = np.array([obc1_rank[c] for c in valid])
s_vals = np.array([obc2_rank[c] for c in valid])
rho, pval = spearmanr(r_vals, s_vals)

ax.scatter(r_vals, s_vals, s=120, color=[PALETTE[COMPS.index(c)] for c in valid],
           edgecolors="k", linewidths=0.5, zorder=5)
for c in valid:
    ax.annotate(SHORT[c], (obc1_rank[c], obc2_rank[c]), fontsize=8,
                ha="center", va="bottom", xytext=(0, 5), textcoords="offset points")
ax.plot([0, 9], [0, 9], "--", color="#888", linewidth=0.8, alpha=0.5)
ax.set_xlabel("OBC1 Rank")
ax.set_ylabel("OBC2 Rank")
ax.set_title(f"b. Rank correlation (Spearman ρ={rho:.3f}, p={pval:.2e})", fontsize=10, loc="left")
ax.set_xlim(0.5, 8.5)
ax.set_ylim(0.5, 8.5)
ax.set_xticks(range(1, 9))
ax.set_yticks(range(1, 9))
ax.grid(alpha=0.2)

plt.tight_layout()
fig_s13_png = os.path.join(ANALYSIS, "FIG_gb_model_comparison.png")
fig_s13_pdf = os.path.join(ANALYSIS, "FIG_gb_model_comparison.pdf")
plt.savefig(fig_s13_png, dpi=300, bbox_inches="tight")
plt.savefig(fig_s13_pdf, dpi=300, bbox_inches="tight")
plt.close()
print(f"1. FigS13: {fig_s13_png}, {fig_s13_pdf}")

# Copy to manuscript figures
for ext in [".png", ".pdf"]:
    src = os.path.join(ANALYSIS, f"FIG_gb_model_comparison{ext}")
    dst = os.path.join(FIGURES, f"FigS13_gb_model_comparison{ext}")
    if os.path.exists(src):
        with open(src, "rb") as s, open(dst, "wb") as d:
            d.write(s.read())
        print(f"   Copied: {dst}")

# ════════════════════════════════════════════════════════════════
# FIGURE S14: GB Cross-Validation — OBC1 vs OBC2 Correlation
# ════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle("GB Model Cross-Validation: OBC1 vs OBC2",
             fontsize=12, fontweight="bold", y=0.98)

# Panel A: Scatter with regression
ax = ax1
x_vals = np.array([obc1[c] for c in valid])
y_vals = np.array([obc2[c] for c in valid])
for c in valid:
    ax.scatter(obc1[c], obc2[c], s=140, color=PALETTE[COMPS.index(c)],
               edgecolors="k", linewidths=0.5, zorder=5)
    ax.annotate(SHORT[c], (obc1[c], obc2[c]), fontsize=8, ha="center",
                va="bottom", xytext=(0, 5), textcoords="offset points")

z = np.polyfit(x_vals, y_vals, 1)
p = np.poly1d(z)
xl = np.linspace(x_vals.min() - 2, x_vals.max() + 2, 50)
ax.plot(xl, p(xl), "--", color="#888", linewidth=0.8, alpha=0.5)
r2 = np.corrcoef(x_vals, y_vals)[0, 1] ** 2
rho2, pval2 = spearmanr(x_vals, y_vals)

ax.set_xlabel("GB(OBC1) ΔG (kcal/mol)")
ax.set_ylabel("GB(OBC2) ΔG (kcal/mol)")
ax.set_title(f"a. Correlation (R²={r2:.3f}, Spearman ρ={rho2:.3f})", fontsize=10, loc="left")
ax.axhline(y=0, color="#bbb", linewidth=0.5)
ax.axvline(x=0, color="#bbb", linewidth=0.5)
ax.text(0.05, 0.95, f"y = {z[0]:.2f}x + {z[1]:.1f}", transform=ax.transAxes,
        fontsize=8, va="top", ha="left", bbox=dict(boxstyle="round", fc="white", alpha=0.7))
ax.grid(alpha=0.2)

# Panel B: Absolute difference (OBC2 - OBC1)
ax = ax2
diffs = [obc2[c] - obc1[c] for c in valid]
colors_b = [PALETTE[COMPS.index(c)] for c in valid]
bars = ax.bar(range(len(valid)), diffs, color=colors_b, edgecolor="white", linewidth=0.6, width=0.6)
ax.set_xticks(range(len(valid)))
ax.set_xticklabels([SHORT[c] for c in valid], fontsize=8, rotation=20)
ax.set_ylabel("ΔG difference: OBC2 − OBC1 (kcal/mol)")
ax.set_title(f"b. Model deviation (RMSD={np.std(diffs):.1f} kcal/mol)", fontsize=10, loc="left")
ax.axhline(y=0, color="#888", linewidth=0.5)
ax.grid(axis="y", alpha=0.2)

# Annotate
for i, d in enumerate(diffs):
    ax.text(i, d + (0.3 if d >= 0 else -1.2), f"{d:+.1f}", ha="center", fontsize=7)

plt.tight_layout()
fig_s14_png = os.path.join(ANALYSIS, "FIG_gb_vs_pb_comparison.png")
fig_s14_pdf = os.path.join(ANALYSIS, "FIG_gb_vs_pb_comparison.pdf")
plt.savefig(fig_s14_png, dpi=300, bbox_inches="tight")
plt.savefig(fig_s14_pdf, dpi=300, bbox_inches="tight")
plt.close()
print(f"2. FigS14: {fig_s14_png}, {fig_s14_pdf}")

# Copy to manuscript figures
for ext in [".png", ".pdf"]:
    src = os.path.join(ANALYSIS, f"FIG_gb_vs_pb_comparison{ext}")
    dst = os.path.join(FIGURES, f"FigS14_gb_vs_pb{ext}")
    if os.path.exists(src):
        with open(src, "rb") as s, open(dst, "wb") as d:
            d.write(s.read())
        print(f"   Copied: {dst}")

# The OBC1/OBC2 comparison is used because PB(PBSA) data requires APBS
# which is unavailable in the current environment. True MM-PBSA results
# can be regenerated by installing APBS in WSL and running _wsl_run_apbs.sh.
print("\n=== ALL COMPLETE ===")
