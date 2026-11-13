#!/usr/bin/env python
"""MM-GBSA ranking + decomposition + docking cross-check + BISCF3 MD stability
for the SMVT membrane-MD reinforcement — restyled to match the existing
"MM-GBSA Analysis" publication figures (flat palette, white bar edges, grid,
value labels). See _pub_figures.py / _figS13_S14.py in MM-GBSA Analysis/.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, r"D:\Researching\SMVT\03_Analysis")
from fig_style import apply_style  # noqa: E402
apply_style()
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.titlesize": 11,
                     "figure.dpi": 300, "savefig.dpi": 300})

# --- flat palette matching the existing 8-compound MM-GBSA figures ---
PAL = {"BIOTIN": "#3498db",      # blue  (matches existing)
       "BISCF3": "#9b59b6",      # purple (candidate targeting head)
       "5HTP": "#1abc9c",        # teal
       "PHENOBARBITAL": "#f39c12"}  # amber (matches existing)
SHORT = {"BIOTIN": "BIO", "BISCF3": "BIS", "5HTP": "5HT", "PHENOBARBITAL": "PHE"}
LBL = {"BIOTIN": "Biotin (COO-)", "BISCF3": "BISCF3 (COO-)",
       "5HTP": "5-HTP (zwitterion)", "PHENOBARBITAL": "Phenobarbital (neutral)"}

BASE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(BASE), "04_Manuscript", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# --- MM-GBSA dG_bind (kcal/mol) ---
MMGBSA = {"BIOTIN": -56.6550, "BISCF3": -54.6334,
          "5HTP": -36.8197, "PHENOBARBITAL": -22.8419}
MMGBSA_ERR = {"BIOTIN": 4.8474, "BISCF3": 2.9956,
              "5HTP": 5.4146, "PHENOBARBITAL": 2.8063}
# --- per-component (kcal/mol): EEL + VDWAALS = "Coulomb+LJ" (favorable), EGB = desolvation penalty ---
DECOMP = {"BIOTIN":       {"VDWAALS": -29.5402, "EEL": -206.1785, "EGB": 184.2337, "ESURF": -5.1700},
          "BISCF3":       {"VDWAALS": -31.3695, "EEL": -137.4018, "EGB": 118.7044, "ESURF": -4.5664},
          "5HTP":         {"VDWAALS": -25.9145, "EEL": -78.2909,  "EGB": 71.6545,  "ESURF": -4.2689},
          "PHENOBARBITAL": {"VDWAALS": -31.0001, "EEL": -14.6292,  "EGB": 27.0580,  "ESURF": -4.2706}}
# --- vina re-dock (zhen2026 pocket, meeko+vina, Na+-free) ---
DOCK = {"BIOTIN": -7.34, "BISCF3": -8.95, "5HTP": -7.34, "PHENOBARBITAL": -7.90}

ORDER = ["BIOTIN", "BISCF3", "5HTP", "PHENOBARBITAL"]  # strongest -> weakest


def load(name, f):
    return np.loadtxt(os.path.join(BASE, name, "analysis", f), delimiter=",", skiprows=1)


def save(fig, stem):
    for ext, kw in [(".pdf", {}), (".svg", {}), (".png", {"dpi": 300})]:
        fig.savefig(os.path.join(FIGDIR, stem + ext), **kw)


def fig_binding():
    """Panel a style: horizontal barh, sorted, white edges, value labels."""
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    vals = [MMGBSA[n] for n in ORDER]
    errs = [MMGBSA_ERR[n] for n in ORDER]
    y = np.arange(len(ORDER))
    bars = ax.barh(y, vals, xerr=errs, color=[PAL[n] for n in ORDER],
                   edgecolor="white", height=0.6, capsize=3,
                   error_kw=dict(lw=0.8, capsize=3))
    ax.set_yticks(y)
    ax.set_yticklabels([LBL[n] for n in ORDER], fontsize=9)
    ax.set_xlabel("MM-GBSA ΔG$_{bind}$ (kcal/mol)")
    ax.axvline(x=0, color="#888", linewidth=0.5)
    ax.grid(axis="x", alpha=0.2)
    for bar, v in zip(bars, vals):
        ax.text(v - 0.8, bar.get_y() + bar.get_height() / 2, f"{v:.1f}",
                ha="right", va="center", fontsize=8, fontweight="bold")
    ax.set_xlim(-72, 4)
    fig.tight_layout()
    save(fig, "Fig_MMGBSA_binding")
    plt.close(fig)


def fig_decomposition():
    """3-term decomposition matching _plot_decomposition.py: dG_total / Coulomb+LJ / GB polar."""
    x = np.arange(len(ORDER))
    dg = [MMGBSA[n] for n in ORDER]
    clj = [DECOMP[n]["EEL"] + DECOMP[n]["VDWAALS"] for n in ORDER]
    gb = [DECOMP[n]["EGB"] for n in ORDER]
    w = 0.26
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    b1 = ax.bar(x - w, dg, w, color=[PAL[n] for n in ORDER], edgecolor="white",
                linewidth=0.5, label="ΔG total", zorder=3)
    b2 = ax.bar(x, clj, w, color="#2c3e50", edgecolor="white", linewidth=0.5,
                label="Coulomb + LJ (vdW)", alpha=0.85, zorder=3)
    b3 = ax.bar(x + w, gb, w, color="#e74c3c", edgecolor="white", linewidth=0.5,
                label="GB polar (desolvation)", alpha=0.85, zorder=3)
    ax.axhline(y=0, color="#555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{LBL[n].split(' ')[0]}\n({SHORT[n]})" for n in ORDER], fontsize=8)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(axis="y", alpha=0.2)
    for bar, v in zip(b1, dg):
        ax.text(bar.get_x() + bar.get_width() / 2, v - 8, f"{v:.0f}", ha="center",
                fontsize=7, fontweight="bold", color="black")
    for bar, v in zip(b2, clj):
        ax.text(bar.get_x() + bar.get_width() / 2, v - 10, f"{v:.0f}", ha="center",
                fontsize=6, color="white", fontweight="bold")
    for bar, v in zip(b3, gb):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 4, f"+{v:.0f}", ha="center",
                fontsize=6, color="#c0392b", fontweight="bold")
    fig.tight_layout()
    save(fig, "Fig_MMGBSA_decomposition")
    plt.close(fig)


def fig_scatter():
    """Vina vs MM-GBSA cross-check, matching panel c of _pub_figures.py."""
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for n in ORDER:
        ax.scatter(DOCK[n], MMGBSA[n], s=110, color=PAL[n], edgecolors="k",
                   linewidths=0.5, zorder=5)
        ax.annotate(SHORT[n], (DOCK[n], MMGBSA[n]), fontsize=8, ha="center",
                    va="bottom", xytext=(0, 5), textcoords="offset points")
    ax.set_xlabel("Vina docking score (kcal/mol)")
    ax.set_ylabel("MM-GBSA ΔG (kcal/mol)")
    ax.grid(alpha=0.2)
    ax.axhline(y=0, color="#bbb", linewidth=0.5)
    ax.set_xlim(-9.3, -6.9)
    fig.tight_layout()
    save(fig, "Fig_docking_vs_MMGBSA")
    plt.close(fig)


def fig_stability():
    """BISCF3 added to the 3-ligand MD stability panels (4-ligand overlay)."""
    ORDER4 = ["BIOTIN", "PHENOBARBITAL", "5HTP", "BISCF3"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0))
    panels = [("rmsd_ligand.csv", "Ligand RMSD (Å)", "a", None),
              ("rmsd_protein.csv", "Protein Cα RMSD (Å)", "b", None),
              ("membrane_thickness.csv", "Membrane thickness (Å)", "c", None),
              ("na_ligand_dist.csv", "Min Na+-ligand distance (Å)", "d", "na")]
    for ax, (f, ylab, letter, kind) in zip(axes.flat, panels):
        for n in ORDER4:
            d = load(n, f)
            y = d[:, 1] if kind != "na" else d[:, 1:].min(axis=1)
            ax.plot(d[:, 0], y, color=PAL[n], lw=1.2)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, None)
        ax.grid(alpha=0.2)
        ax.set_title(letter, loc="left", fontsize=11, fontweight="bold")
    handles = [plt.Line2D([], [], color=PAL[n], lw=1.6, label=LBL[n]) for n in ORDER4]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, "Fig_BISCF3_MD_stability")
    plt.close(fig)


def fig_summary():
    """2x2 combined summary (a: ranking, b: decomposition, c: Vina vs MM-GBSA,
    d: ligand stability), matching FIG_comprehensive_summary layout."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle("SMVT membrane MD — MM-GBSA & docking summary (4 ligands)",
                 fontsize=13, fontweight="bold", y=0.99)

    # (a) MM-GBSA ranking barh
    ax = axes[0, 0]
    vals = [MMGBSA[n] for n in ORDER]
    errs = [MMGBSA_ERR[n] for n in ORDER]
    y = np.arange(len(ORDER))
    bars = ax.barh(y, vals, xerr=errs, color=[PAL[n] for n in ORDER],
                   edgecolor="white", height=0.6, capsize=3,
                   error_kw=dict(lw=0.8, capsize=3))
    ax.set_yticks(y)
    ax.set_yticklabels([LBL[n] for n in ORDER], fontsize=8)
    ax.set_xlabel("MM-GBSA ΔG (kcal/mol)")
    ax.set_title("a. Binding free energy", loc="left", fontsize=10, fontweight="bold")
    ax.axvline(x=0, color="#888", lw=0.5)
    ax.grid(axis="x", alpha=0.2)
    for bar, v in zip(bars, vals):
        ax.text(v - 0.6, bar.get_y() + bar.get_height() / 2, f"{v:.1f}",
                ha="right", va="center", fontsize=7, fontweight="bold")
    ax.set_xlim(-72, 4)

    # (b) energy decomposition (3-term)
    ax = axes[0, 1]
    x = np.arange(len(ORDER))
    dg = [MMGBSA[n] for n in ORDER]
    clj = [DECOMP[n]["EEL"] + DECOMP[n]["VDWAALS"] for n in ORDER]
    gb = [DECOMP[n]["EGB"] for n in ORDER]
    w = 0.26
    ax.bar(x - w, dg, w, color=[PAL[n] for n in ORDER], edgecolor="white",
           lw=0.5, label="ΔG total", zorder=3)
    ax.bar(x, clj, w, color="#2c3e50", edgecolor="white", lw=0.5,
           label="Coulomb + LJ", alpha=0.85, zorder=3)
    ax.bar(x + w, gb, w, color="#e74c3c", edgecolor="white", lw=0.5,
           label="GB desolvation", alpha=0.85, zorder=3)
    ax.axhline(0, color="#555", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[n] for n in ORDER], fontsize=8)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("b. Energy decomposition", loc="left", fontsize=10, fontweight="bold")
    ax.legend(fontsize=6.5, loc="upper right")
    ax.grid(axis="y", alpha=0.2)

    # (c) Vina vs MM-GBSA scatter
    ax = axes[1, 0]
    for n in ORDER:
        ax.scatter(DOCK[n], MMGBSA[n], s=100, color=PAL[n], edgecolors="k",
                   linewidths=0.5, zorder=5)
        ax.annotate(SHORT[n], (DOCK[n], MMGBSA[n]), fontsize=8, ha="center",
                    va="bottom", xytext=(0, 5), textcoords="offset points")
    ax.set_xlabel("Vina docking score (kcal/mol)")
    ax.set_ylabel("MM-GBSA ΔG (kcal/mol)")
    ax.set_title("c. Vina vs MM-GBSA", loc="left", fontsize=10, fontweight="bold")
    ax.grid(alpha=0.2)
    ax.axhline(0, color="#bbb", lw=0.5)
    ax.set_xlim(-9.3, -6.9)

    # (d) ligand RMSD stability (4-ligand overlay)
    ax = axes[1, 1]
    for n in ORDER:
        d = load(n, "rmsd_ligand.csv")
        ax.plot(d[:, 0], d[:, 1], color=PAL[n], lw=1.2)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("d. Ligand stability", loc="left", fontsize=10, fontweight="bold")
    ax.set_ylim(0, None)
    ax.grid(alpha=0.2)
    handles = [plt.Line2D([], [], color=PAL[n], lw=1.6, label=LBL[n]) for n in ORDER]
    ax.legend(handles=handles, fontsize=6.5, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, "Fig_summary_MMGBSA")
    plt.close(fig)


if __name__ == "__main__":
    fig_binding()
    fig_decomposition()
    fig_scatter()
    fig_stability()
    fig_summary()
    print("done: 4 MM-GBSA/docking/BISCF3 figures + 1 combined summary ->", FIGDIR)
