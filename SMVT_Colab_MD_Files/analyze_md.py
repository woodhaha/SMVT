#!/usr/bin/env python3
"""
Analyze SMVT MD trajectory — 5-panel publication-quality output.
RMSD | RMSF | H-bond occupancy | SASA | Radius of Gyration
"""
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, hydrogenbonds, sasa
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 11, 'figure.dpi': 150, 'savefig.dpi': 300})
import numpy as np
import os, sys, json, argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--dir', required=True, help='MD output directory')
parser.add_argument('--name', default='PHENOBARBITAL', help='Compound name')
parser.add_argument('--prod-ns', type=int, default=100, help='Production length in ns')
args = parser.parse_args()

OUT_DIR = args.dir
NAME = args.name
PROD_NS = args.prod_ns

dcd = f'{OUT_DIR}/{NAME}_{PROD_NS}ns.dcd'
pdb = f'{OUT_DIR}/{NAME}_final.pdb'

if not os.path.exists(dcd):
    print(f'ERROR: No trajectory at {dcd}')
    sys.exit(1)
if not os.path.exists(pdb):
    print(f'WARNING: No final PDB at {pdb}, using receptor PDB')
    sys.exit(1)

print(f'Loading: {dcd}')
u = mda.Universe(pdb, dcd)
protein = u.select_atoms('protein')
backbone = u.select_atoms('protein and backbone')
ca = u.select_atoms('protein and name CA')
ligand = u.select_atoms('not protein and not (resname HOH or resname WAT or resname SOL or resname NA or resname CL)')
n_frames = len(u.trajectory)
dt_ns = u.trajectory.dt / 1000  # timestep in ns

print(f'Protein: {len(protein)} atoms | Backbone: {len(backbone)} | CA: {len(ca)} | Ligand: {len(ligand)}')
print(f'Frames: {n_frames} | Time: {n_frames * u.trajectory.dt / 1000:.0f} ns')

results = {}

# ═══ Panel 1: Backbone RMSD ═══
print('\n[1/5] Computing backbone RMSD...')
R = rms.RMSD(u, backbone, ref_frame=0).run()
time_ns = R.results['time'] / 1000
rmsd_vals = R.results.rmsd[:, 2]
final_20ns = rmsd_vals[-int(len(rmsd_vals)*0.2):]
rmsd_drift = abs(final_20ns.mean() - rmsd_vals.mean())
results['rmsd_mean'] = float(rmsd_vals.mean())
results['rmsd_std'] = float(rmsd_vals.std())
results['rmsd_drift'] = float(rmsd_drift)
results['rmsd_stable'] = rmsd_drift < 0.5 and rmsd_vals.mean() < 3.0

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(time_ns, rmsd_vals, color='#2171b5', linewidth=0.8)
ax.axhline(y=rmsd_vals.mean(), color='#cb181d', linestyle='--', label=f'Mean = {rmsd_vals.mean():.2f} Å')
ax.axhline(y=3.0, color='gray', linestyle=':', alpha=0.5, label='3.0 Å')
ax.set_xlabel('Time (ns)'); ax.set_ylabel('RMSD (Å)')
ax.set_title(f'{NAME} — SMVT Backbone RMSD ({PROD_NS}ns)')
ax.legend(loc='upper right')
plt.tight_layout(); plt.savefig(f'{OUT_DIR}/rmsd_plot.png', dpi=300, bbox_inches='tight'); plt.close()
print(f'  Mean RMSD: {rmsd_vals.mean():.2f} ± {rmsd_vals.std():.2f} Å | Drift: {rmsd_drift:.2f} Å')

# ═══ Panel 2: Per-Residue RMSF ═══
print('[2/5] Computing per-residue RMSF...')
# Align to average first
average = align.AverageStructure(u, u, select='protein and backbone', ref_frame=0).run()
ref = average.results.universe
align.AlignTraj(u, ref, select='protein and backbone', in_memory=True).run()

# Per-CA RMSF
ca_positions = np.zeros((n_frames, len(ca), 3))
for i, ts in enumerate(u.trajectory):
    ca_positions[i] = ca.positions
ca_mean = ca_positions.mean(axis=0)
rmsf_vals = np.sqrt(3 * ((ca_positions - ca_mean)**2).mean(axis=0).sum(axis=1))
resids = ca.residues.resids
results['rmsf_mean'] = float(rmsf_vals.mean())
results['rmsf_std'] = float(rmsf_vals.std())

threshold = rmsf_vals.mean() + 2*rmsf_vals.std()
high_flex = [(resids[i], float(rmsf_vals[i])) for i in range(len(rmsf_vals)) if rmsf_vals[i] > threshold]
results['flexible_residues'] = [int(r[0]) for r in high_flex]

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(resids, rmsf_vals, color='#2171b5', linewidth=0.8)
ax.fill_between(resids, 0, rmsf_vals, alpha=0.15, color='#2171b5')
ax.set_xlabel('Residue number'); ax.set_ylabel('RMSF (Å)')
ax.set_title(f'{NAME} — SMVT Per-Residue RMSF')
for rid, val in high_flex[:10]:
    ax.annotate(f'R{rid}', (rid, val), fontsize=7, alpha=0.7, xytext=(0, 5), textcoords='offset points')
plt.tight_layout(); plt.savefig(f'{OUT_DIR}/rmsf_plot.png', dpi=300, bbox_inches='tight'); plt.close()
print(f'  Mean RMSF: {rmsf_vals.mean():.2f} Å | Flexible (>2σ): {len(high_flex)} residues')

# ═══ Panel 3: Hydrogen Bond Occupancy ═══
print('[3/5] Computing hydrogen bond occupancy...')
if len(ligand) > 0:
    hba = hydrogenbonds.HydrogenBondAnalysis(u, between={
        'protein', f'index {ligand.indices[0]}:{ligand.indices[-1]}'
    })
    hba.run(verbose=False)
    hbond_count = hba.count_by_time()
    occupancy = (hbond_count > 0).sum() / n_frames * 100
    results['hbond_occupancy_pct'] = float(occupancy)
    results['hbond_mean_count'] = float(hbond_count.mean())
else:
    hbond_count = np.zeros(n_frames)
    occupancy = 0
    results['hbond_occupancy_pct'] = 0
    results['hbond_mean_count'] = 0
    print('  WARNING: No ligand detected for H-bond analysis')

fig, ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(time_ns[1:], hbond_count, color='#238b45', linewidth=0.5)
ax1.set_xlabel('Time (ns)'); ax1.set_ylabel('H-bond count')
ax1.set_title(f'{NAME} — Ligand-SMVT H-bonds (occupancy={occupancy:.0f}%)')
ax1.axhline(y=hbond_count.mean(), color='#cb181d', linestyle='--')
plt.tight_layout(); plt.savefig(f'{OUT_DIR}/hbond_occupancy_plot.png', dpi=300, bbox_inches='tight'); plt.close()
print(f'  H-bond occupancy: {occupancy:.1f}% | Mean: {hbond_count.mean():.1f}/frame')

# ═══ Panel 4: SASA ═══
print('[4/5] Computing SASA...')
sasa_calc = sasa.SASA(u, select='protein').run()
sasa_vals = sasa_calc.results.sasa / 100  # Å² → nm²
results['sasa_mean_nm2'] = float(np.mean(sasa_vals))
results['sasa_std_nm2'] = float(np.std(sasa_vals))

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(time_ns[::max(1, len(sasa_vals)//len(time_ns))], sasa_vals, color='#6a51a3', linewidth=0.8)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('SASA (nm²)')
ax.set_title(f'{NAME} — SMVT Solvent Accessible Surface Area')
ax.axhline(y=np.mean(sasa_vals), color='#cb181d', linestyle='--', label=f'Mean = {np.mean(sasa_vals):.1f} nm²')
ax.legend()
plt.tight_layout(); plt.savefig(f'{OUT_DIR}/sasa_plot.png', dpi=300, bbox_inches='tight'); plt.close()
print(f'  Mean SASA: {np.mean(sasa_vals):.1f} ± {np.std(sasa_vals):.1f} nm²')

# ═══ Panel 5: Radius of Gyration ═══
print('[5/5] Computing radius of gyration...')
rg_vals = np.array([protein.radius_of_gyration() for ts in u.trajectory])
results['rg_mean_A'] = float(np.mean(rg_vals))
results['rg_std_A'] = float(np.std(rg_vals))
results['rg_drift_A'] = float(abs(rg_vals[-1] - rg_vals[0]))

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(time_ns, rg_vals, color='#d94801', linewidth=0.8)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('Rg (Å)')
ax.set_title(f'{NAME} — SMVT Radius of Gyration')
ax.axhline(y=np.mean(rg_vals), color='#cb181d', linestyle='--', label=f'Mean = {np.mean(rg_vals):.1f} Å')
ax.legend()
plt.tight_layout(); plt.savefig(f'{OUT_DIR}/rg_plot.png', dpi=300, bbox_inches='tight'); plt.close()
print(f'  Mean Rg: {np.mean(rg_vals):.1f} Å | Drift: {abs(rg_vals[-1]-rg_vals[0]):.2f} Å')

# ═══ Generate Report ═══
print('\n' + '='*60)
report_path = f'{OUT_DIR}/analysis_report.md'
with open(report_path, 'w') as f:
    f.write(f"""# {NAME} — SMVT {PROD_NS}ns MD Analysis Report

## Simulation Summary
| Parameter | Value |
|-----------|-------|
| Protein | SLC5A6/SMVT (AlphaFold, OpenMM minimized) |
| Ligand | {NAME} (barbiturate, docking -8.30 kcal/mol) |
| Force Field | AMBER ff14SB + GAFF2 + TIP3P |
| Ensemble | NPT, 310K, 1 atm |
| Production | {PROD_NS} ns, 2 fs timestep |
| Water model | TIP3P, 0.15M NaCl |

## Key Metrics
| Metric | Value | Status |
|--------|-------|--------|
| Backbone RMSD | {rmsd_vals.mean():.2f} ± {rmsd_vals.std():.2f} Å | {'✅ PASS < 3.0' if rmsd_vals.mean() < 3.0 else '⚠ WARN'} |
| RMSD drift (final 20ns) | {rmsd_drift:.2f} Å | {'✅ Stable' if rmsd_drift < 0.5 else '⚠ Drifting'} |
| Mean RMSF | {rmsf_vals.mean():.2f} Å | {'—'} |
| H-bond occupancy | {occupancy:.1f}% | {'✅ > 50%' if occupancy > 50 else '⚠ < 50%'} |
| Mean H-bonds/frame | {hbond_count.mean():.1f} | {'—'} |
| Mean SASA | {np.mean(sasa_vals):.1f} nm² | {'—'} |
| Mean Rg | {np.mean(rg_vals):.1f} Å | {'✅ Stable' if abs(rg_vals[-1]-rg_vals[0]) < 2.0 else '⚠ Compacting/Expanding'} |
| Flexible regions | {len(high_flex)} residues | Residues: {[r[0] for r in high_flex[:5]]}{'...' if len(high_flex) > 5 else ''} |

## Interpretation
- The {NAME}-SMVT complex {'is **stable**' if rmsd_vals.mean() < 3.0 and rmsd_drift < 0.5 else 'shows **instability**'} over {PROD_NS}ns of MD simulation.
- RMSD {'plateaued' if rmsd_drift < 0.5 else 'continued to drift'} in the final 20ns, suggesting {'equilibration was achieved' if rmsd_drift < 0.5 else 'longer simulation may be needed'}.
- {NAME} maintained H-bonds with SMVT in {occupancy:.0f}% of trajectory frames, {'supporting' if occupancy > 50 else 'weakly supporting'} the barbiturate-ureido ring mimicry mechanism.
- The major flexibility hotspots are at residues {[r[0] for r in high_flex[:5]]}, which {'may correspond to extracellular loops' if high_flex else 'are within expected range for membrane proteins'}.

## Conclusions
1. {'✅' if rmsd_vals.mean() < 3.0 and rmsd_drift < 0.5 else '⚠'} The complex is {'stable under physiological MD' if rmsd_vals.mean() < 3.0 and rmsd_drift < 0.5 else 'unstable — consider extending simulation or re-docking'}.
2. {'✅' if occupancy > 50 else '⚠'} The carboxyl-independent binding mode is {'supported' if occupancy > 50 else 'weakly supported'} by H-bond data.
3. {'Proceed to production MD on additional barbiturates.' if rmsd_vals.mean() < 3.0 else 'Debug and optimize before expanding to other compounds.'}

## Output Files
| File | Path |
|------|------|
| RMSD plot | {OUT_DIR}/rmsd_plot.png |
| RMSF plot | {OUT_DIR}/rmsf_plot.png |
| H-bond plot | {OUT_DIR}/hbond_occupancy_plot.png |
| SASA plot | {OUT_DIR}/sasa_plot.png |
| Rg plot | {OUT_DIR}/rg_plot.png |
| Metrics JSON | {OUT_DIR}/analysis_metrics.json |
""")

# Save metrics
with open(f'{OUT_DIR}/analysis_metrics.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f'Report: {report_path}')
print(f'Metrics: {OUT_DIR}/analysis_metrics.json')
print('\nAll 5 plots generated:')
for p in ['rmsd', 'rmsf', 'hbond_occupancy', 'sasa', 'rg']:
    print(f'  ✅ {OUT_DIR}/{p}_plot.png')
