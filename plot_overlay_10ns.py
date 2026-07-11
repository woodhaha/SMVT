#!/usr/bin/env python3
"""Overlay plots — all 8 compounds, x-axis capped at 10ns for zoomed view."""
import mdtraj as md, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11, 'figure.dpi': 200, 'savefig.dpi': 300})

BASE = "SMVT_MD/trajectories"
OUT = "out"
os.makedirs(OUT, exist_ok=True)

COMPOUNDS = ["BIOTIN","ESKETAMINE","FUROSEMIDE","GABAPENTIN_ENACARBIL",
             "HYDROMORPHONE","NAFTAZONE","PHENOBARBITAL","RIBOFLAVIN"]
COLORS    = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f']
LABELS    = {"BIOTIN":"Biotin(ctrl)","ESKETAMINE":"Esketamine","FUROSEMIDE":"Furosemide",
             "GABAPENTIN_ENACARBIL":"Gabapentin enac.(ctrl)","HYDROMORPHONE":"Hydromorphone",
             "NAFTAZONE":"Naftazone","PHENOBARBITAL":"Phenobarbital","RIBOFLAVIN":"Riboflavin(ctrl)"}

XMAX = 1  # ns

def load(name):
    d = f"{BASE}/{name}"
    t = md.load(f"{d}/{name}_100ns.dcd", top=f"{d}/{name}_final.pdb")
    return t

def mask_10ns(t):
    """Return slice indices for frames within first 10ns."""
    idx = np.where(t.time / 1000 <= XMAX)[0]
    return idx

# ─── 1. RMSD ───
print("[1/3] RMSD overlay (0-10ns)...")
fig, ax = plt.subplots(figsize=(12, 5))
for name, clr in zip(COMPOUNDS, COLORS):
    t = load(name); idx = mask_10ns(t)
    bb = t.topology.select("backbone")
    rmsd = md.rmsd(t, t, 0, atom_indices=bb) * 10  # nm→Å
    ax.plot(t.time[idx]/1000, rmsd[idx], color=clr, lw=0.7, label=LABELS[name])
ax.set(xlabel="Time (ns)", ylabel="Backbone RMSD (Å)",
       title=f"SMVT — Backbone RMSD Overlay (0–{XMAX}ns, 8 Compounds)")
ax.legend(fontsize=7, ncol=2); ax.set_xlim(0, XMAX); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/rmsd_overlay.png", dpi=300); plt.close()

# ─── 2. Rg ───
print("[2/3] Rg overlay (0-10ns)...")
fig, ax = plt.subplots(figsize=(12, 5))
for name, clr in zip(COMPOUNDS, COLORS):
    t = load(name); idx = mask_10ns(t)
    prot = t.topology.select("protein")
    rg = md.compute_rg(t.atom_slice(prot)) * 10  # nm→Å
    ax.plot(t.time[idx]/1000, rg[idx], color=clr, lw=0.7, label=LABELS[name])
ax.set(xlabel="Time (ns)", ylabel="Rg (Å)",
       title=f"SMVT — Radius of Gyration Overlay (0–{XMAX}ns)")
ax.legend(fontsize=7, ncol=2); ax.set_xlim(0, XMAX); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/rg_overlay.png", dpi=300); plt.close()

# ─── 3. RMSF ───
print("[3/3] RMSF overlay...")
fig, ax = plt.subplots(figsize=(14, 5))
for name, clr in zip(COMPOUNDS, COLORS):
    t = load(name)
    bb = t.topology.select("backbone")
    ca = t.topology.select("name CA")
    t_align = t.superpose(t, 0, atom_indices=bb)
    rmsf = md.rmsf(t_align, t_align, 0, atom_indices=ca) * 10
    resids = [t.topology.atom(i).residue.resSeq for i in ca]
    ax.plot(resids, rmsf, color=clr, lw=0.5, label=LABELS[name])
ax.set(xlabel="Residue", ylabel="RMSF (Å)",
       title="SMVT — Per-Residue RMSF Overlay")
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/rmsf_overlay.png", dpi=300); plt.close()

print(f"✅ Done → {OUT}/{rmsd,rg,rmsf}_overlay.png (0–{XMAX}ns zoom)")
