#!/usr/bin/env python3
"""Compute RMSD/Rg/RMSF overlay from 8 DCDs using mdtraj in batch mode."""
import os, json
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/root/autodl-tmp/SMVT_AutoDL_Package/trajectories'
COMPOUNDS = [
    'BIOTIN', 'FUROSEMIDE', 'GABAPENTIN_ENACARBIL',
    'HYDROMORPHONE', 'NAFTAZONE', 'PHENOBARBITAL',
    'RIBOFLAVIN', 'ESKETAMINE'
]
COLORS = ['#1f77b4','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#ff7f0e']
OUT = os.path.join(BASE, 'overlay')
os.makedirs(OUT, exist_ok=True)

all_data = {}

for i, name in enumerate(COMPOUNDS):
    d = os.path.join(BASE, name)
    dcdf = os.path.join(d, f'{name}_100ns.dcd')
    pdbf = os.path.join(d, f'{name}_final.pdb')
    print(f"[{i+1}/8] Loading {name}...", flush=True)

    traj = md.load(dcdf, top=pdbf)
    time_ns = traj.time / 1000

    # Use first frame as reference
    ref = traj[0]

    # ── Backbone RMSD ── (vectorized)
    bb = traj.topology.select('backbone')
    rmsd = md.rmsd(traj, ref, atom_indices=bb) * 10  # nm→Å
    print(f"  RMSD: {rmsd.shape} frames", flush=True)

    # ── Protein Rg ── (vectorized, whole trajectory)
    prot = traj.topology.select('protein')
    # subsample 1/10 for speed — 100 frames is plenty for Rg trend
    step = max(1, traj.n_frames // 100)
    # compute_rg returns shape (n_frames, 1) — need flatten
    rg_all = np.array([md.compute_rg(traj[j]).item() for j in range(0, traj.n_frames, step)]) * 10
    # Expand back via interpolation
    rg = np.interp(np.arange(traj.n_frames), np.arange(0, traj.n_frames, step), rg_all)
    print(f"  Rg: {rg.shape} frames (sampled every {step})", flush=True)

    # ── CA RMSF ── (vectorized)
    ca = traj.topology.select('name CA')
    rmsf = md.rmsf(traj, ref, atom_indices=ca) * 10
    ca_idx = np.arange(1, len(rmsf) + 1)
    print(f"  RMSF: {len(rmsf)} residues", flush=True)

    all_data[name] = {
        'color': COLORS[i],
        'time': time_ns, 'rmsd': rmsd, 'rg': rg,
        'rmsf': rmsf, 'ca_idx': ca_idx,
    }

# Save numeric data as JSON
out_json = {}
for name, d in all_data.items():
    out_json[name] = {
        'time_ns': d['time'].tolist(),
        'rmsd': d['rmsd'].tolist(),
        'rg': d['rg'].tolist(),
        'rmsf': d['rmsf'].tolist(),
    }
with open(os.path.join(OUT, 'overlay_metrics.json'), 'w') as f:
    json.dump(out_json, f)

# ── RMSD Overlay ──
fig, ax = plt.subplots(figsize=(12, 5))
for name, d in all_data.items():
    ax.plot(d['time'], d['rmsd'], color=d['color'], label=name, lw=0.8)
ax.set(xlabel='Time (ns)', ylabel='Backbone RMSD (Å)',
       title='SMVT — Backbone RMSD Overlay (8 Compounds, 100ns)')
ax.legend(fontsize=8, ncol=2); ax.set_xlim(0, 100); ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'rmsd_overlay.png'), dpi=200); plt.close(fig)
print("RMSD overlay saved.", flush=True)

# ── Rg Overlay ──
fig, ax = plt.subplots(figsize=(12, 5))
for name, d in all_data.items():
    ax.plot(d['time'], d['rg'], color=d['color'], label=name, lw=0.8)
ax.set(xlabel='Time (ns)', ylabel='Rg (Å)',
       title='SMVT — Radius of Gyration Overlay (8 Compounds, 100ns)')
ax.legend(fontsize=8, ncol=2); ax.set_xlim(0, 100); ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'rg_overlay.png'), dpi=200); plt.close(fig)
print("Rg overlay saved.", flush=True)

# ── RMSF Overlay ──
fig, ax = plt.subplots(figsize=(14, 5))
for name, d in all_data.items():
    ax.plot(d['ca_idx'], d['rmsf'], color=d['color'], label=name, lw=0.6)
ax.set(xlabel='Residue Index (CA)', ylabel='RMSF (Å)',
       title='SMVT — Per-Residue RMSF Overlay (8 Compounds, 100ns)')
ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'rmsf_overlay.png'), dpi=200); plt.close(fig)
print("RMSF overlay saved.", flush=True)

# ── Summary Table ──
print(f"\n{'Compound':25s} {'Mean RMSD':>9s} {'SD RMSD':>8s} {'Mean Rg':>8s} {'SD Rg':>6s}")
print("-"*56)
for name in COMPOUNDS:
    d = all_data[name]
    print(f"{name:25s} {np.mean(d['rmsd']):>8.3f} {np.std(d['rmsd']):>8.3f} "
          f"{np.mean(d['rg']):>8.3f} {np.std(d['rg']):>6.3f}")

print("\n✅ All overlays generated.", flush=True)
