#!/usr/bin/env python
"""Regenerate Nature-style Fig S17 (RDF) and Fig S18 (MD stability) for the
SMVT membrane-embedded Na+-bound MD reinforcement (PLOS CB revision).

Inputs : SMVT_membrane_MD/{BIOTIN,PHENOBARBITAL,5HTP}/analysis/*.csv
Outputs: 04_Manuscript/figures/FigS17_Na_ligand_RDF.{pdf,svg,png}
         04_Manuscript/figures/FigS18_membrane_MD_stability.{pdf,svg,png}

Nature-figure conventions (Arial/Helvetica 8pt, editable SVG/PDF text,
muted colorblind-reasonable palette, bold lowercase panel labels).
"""
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- mandatory Nature rcParams ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"   # editable SVG text
plt.rcParams["pdf.fonttype"] = 42       # editable TrueType text in PDF
plt.rcParams["font.size"] = 8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.titlesize"] = 8
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["legend.frameon"] = False

# --- palette (muted) ---
PAL = {"BIOTIN": "#0F4D92",        # blue_main  — natural substrate (reference)
       "PHENOBARBITAL": "#B64342",  # red_strong — barbiturate hit
       "5HTP": "#42949E"}           # teal       — amino-acid candidate
LBL = {"BIOTIN": "Biotin (COO⁻)",
       "PHENOBARBITAL": "Phenobarbital",
       "5HTP": "5-HTP"}
ORDER = ["BIOTIN", "PHENOBARBITAL", "5HTP"]

BASE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(BASE), "04_Manuscript", "figures")


def load(name, f):
    return np.loadtxt(os.path.join(BASE, name, "analysis", f), delimiter=",", skiprows=1)


def panel(ax, letter):
    ax.text(-0.10, 1.05, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="bottom")


def shared_legend(fig, x=0.5, y=1.02):
    h = [plt.Line2D([], [], color=PAL[n], lw=1.4, label=LBL[n]) for n in ORDER]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(x, y))


def save(fig, stem):
    for ext, kw in [(".pdf", {}), (".svg", {}), (".png", {"dpi": 300})]:
        fig.savefig(os.path.join(FIGDIR, stem + ext), **kw)


def fig_s17():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    for n in ORDER:
        r, g = load(n, "rdf_Na_ligandHeavy.csv").T
        ax.plot(r, g, color=PAL[n], lw=1.2)
    ax.axvline(3.2, color="#767676", ls="--", lw=0.8)
    ax.text(3.35, ax.get_ylim()[1] * 0.88, "3.2 Å", fontsize=7, color="#767676", va="top")
    ax.set_xlabel("Na⁺–ligand distance (Å)")
    ax.set_ylabel("g(r)")
    ax.set_xlim(0, 16); ax.set_ylim(0, None)
    panel(ax, "a")

    ax = axes[1]
    for n in ORDER:
        rp, hp, cp = load(n, "rdf_funcNa_proteinON.csv").T
        ax.plot(rp, cp, color=PAL[n], lw=1.2)
    ax.set_xlabel("Na⁺–protein O/N distance (Å)")
    ax.set_ylabel("Coordination number n(r)")
    ax.set_xlim(0, 6); ax.set_ylim(0, None)
    panel(ax, "b")
    shared_legend(fig)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, "FigS17_Na_ligand_RDF")
    plt.close(fig)


def fig_s18():
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    panels = [("rmsd_ligand.csv", "Ligand RMSD (Å)", "a", None),
              ("rmsd_protein.csv", "Protein Cα RMSD (Å)", "b", None),
              ("membrane_thickness.csv", "Membrane thickness (Å)", "c", None),
              ("na_ligand_dist.csv", "Min Na⁺–ligand distance (Å)", "d", "na")]
    for ax, (f, ylab, letter, kind) in zip(axes.flat, panels):
        for n in ORDER:
            d = load(n, f)
            y = d[:, 1] if kind != "na" else d[:, 1:].min(axis=1)
            ax.plot(d[:, 0], y, color=PAL[n], lw=1.0)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, None)
        panel(ax, letter)
    shared_legend(fig)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, "FigS18_membrane_MD_stability")
    plt.close(fig)


if __name__ == "__main__":
    fig_s17()
    fig_s18()
    print("done: FigS17 + FigS18 regenerated")
