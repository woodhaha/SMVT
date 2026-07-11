"""Publication-quality figures for manuscript — all remaining analyses."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from math import pi
plt.rcParams.update({'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
                     'figure.dpi': 300, 'savefig.dpi': 300, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial']})

ANALYSIS = r'D:\Researching\SMVT\MM-GBSA Analysis'
SMVT_MD = r'D:\Researching\SMVT\SMVT-MD-Analysis'
COMPS = ['BIOTIN','HYDROMORPHONE','GABAPENTIN_ENACARBIL','NAFTAZONE',
         'ESKETAMINE','FUROSEMIDE','PHENOBARBITAL','RIBOFLAVIN']
PAL = {'BIOTIN':'#3498db','HYDROMORPHONE':'#9b59b6','GABAPENTIN_ENACARBIL':'#2ecc71',
       'NAFTAZONE':'#e67e22','ESKETAMINE':'#1abc9c','FUROSEMIDE':'#e74c3c',
       'PHENOBARBITAL':'#f39c12','RIBOFLAVIN':'#95a5a6'}
LAB = {'BIOTIN':'Biotin (REF)','HYDROMORPHONE':'Hydromorphone',
       'GABAPENTIN_ENACARBIL':'Gabapentin Enacarbil','NAFTAZONE':'Nafazone',
       'ESKETAMINE':'Esketamine','FUROSEMIDE':'Furosemide',
       'PHENOBARBITAL':'Phenobarbital','RIBOFLAVIN':'Riboflavin'}

with open(os.path.join(ANALYSIS, 'md_master_metrics.json')) as f:
    MET = json.load(f)

with open(os.path.join(ANALYSIS, 'compounds.json')) as f:
    CMP = json.load(f)

# Get SMILES
SMILES = {c['id']: c['smiles'] for c in CMP['compounds']}

# ================================================================
# FIGURE 1: PC1×PC2 FEL — 2×4 grid (publication quality)
# ================================================================
KT = 0.596
fig, axes = plt.subplots(2, 4, figsize=(20, 10), constrained_layout=True)
fig.suptitle('Free Energy Landscape (PC1 × PC2) — 100 ns MD Trajectories',
             fontsize=14, fontweight='bold')

for i, c in enumerate(COMPS):
    d = np.load(os.path.join(ANALYSIS, f'{c}_pca.npz'))
    pc1, pc2 = d['pc1'], d['pc2']
    ve = d['var_exp']
    ax = axes[i//4, i%4]

    H, xe, ye = np.histogram2d(pc1, pc2, bins=(50, 50), density=True)
    H = gaussian_filter(H, sigma=1.0)
    H[H < 1e-12] = 1e-12
    FE = -KT * np.log(H)
    FE -= FE.min()
    X, Y = np.meshgrid((xe[:-1]+xe[1:])/2, (ye[:-1]+ye[1:])/2, indexing='ij')

    cnt = ax.contourf(X, Y, FE.T, levels=np.arange(0, 7.5, 0.5),
                      cmap='viridis', extend='max')
    ax.plot(pc1[0], pc2[0], '*', color='cyan', ms=10, mew=0.5, label='Start')
    ax.plot(pc1[-1], pc2[-1], 'D', color='red', ms=6, mew=0.5, label='End')
    ax.set_xlabel(f'PC1 ({ve[0]*100:.0f}%)')
    ax.set_ylabel(f'PC2 ({ve[1]*100:.0f}%)')
    ax.set_title(LAB[c], fontsize=10)
    if i == 0:
        ax.legend(fontsize=7, loc='upper right')

plt.savefig(os.path.join(ANALYSIS, 'FIG_pc12_fel_pub.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_pc12_fel_pub.pdf'), dpi=300, bbox_inches='tight')
plt.close()
print('1. FIG_pc12_fel_pub.png + .pdf')

# ================================================================
# FIGURE 2: Comprehensive summary — 6-panel (manuscript ready)
# ================================================================
DG = {}
for c in COMPS:
    fn = os.path.join(ANALYSIS, c, 'mmgbsa_results.json')
    if os.path.exists(fn):
        with open(fn) as f:
            DG[c] = json.load(f)

VINA = {}
for c in COMPS:
    fn = os.path.join(ANALYSIS, f'{c}_vina_log.txt')
    if os.path.exists(fn):
        with open(fn) as f:
            for line in f:
                if line.strip().startswith('   1'):
                    VINA[c] = float(line.strip().split()[1])

fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle('SMVT Computational Validation — Integrated Summary',
             fontsize=15, fontweight='bold', y=1.01)
axs = axes.flatten()

# (a) MM-GBSA ranking
ax = axs[0]
sorted_c = sorted(DG.keys(), key=lambda x: DG[x]['dG_kcal'])
vals = [DG[c]['dG_kcal'] for c in sorted_c]
colors = [PAL[c] for c in sorted_c]
bars = ax.barh(range(len(sorted_c)), vals, color=colors, edgecolor='white', height=0.6)
ax.set_yticks(range(len(sorted_c)))
ax.set_yticklabels([sorted_c[i] for i in range(len(sorted_c))], fontsize=8)
ax.set_xlabel('ΔG (kcal/mol)')
ax.set_title('a. MM-GBSA Binding Free Energy', fontsize=11, loc='left')
ax.axvline(x=0, color='#888', linewidth=0.5)
ax.grid(axis='x', alpha=0.2)
# Annotate
for bar, v in zip(bars, vals):
    ax.text(v - 0.3 if v < 0 else v + 0.3, bar.get_y()+bar.get_height()/2,
            f'{v:.1f}', ha='right' if v < 0 else 'left', va='center', fontsize=7)

# (b) Energy decomposition
ax = axs[1]
x = np.arange(len(COMPS))
clj = [DG[c]['E_complex_kJ'] - DG[c]['E_receptor_kJ'] - DG[c]['E_ligand_kJ'] for c in COMPS]
ax.bar(x, clj, color=[PAL[c] for c in COMPS], edgecolor='white', width=0.6)
ax.set_xticks(x)
ax.set_xticklabels([c[:4] for c in COMPS], fontsize=7, rotation=30)
ax.set_ylabel('ΔG (kJ/mol)')
ax.set_title('b. MM-GBSA by Compound', fontsize=11, loc='left')
ax.axhline(y=0, color='#888', linewidth=0.5)
ax.grid(axis='y', alpha=0.2)

# (c) Vina vs MM-GBSA comparison
ax = axs[2]
for c in COMPS:
    if c in VINA and c in DG:
        ax.scatter(VINA[c], DG[c]['dG_kcal'], s=100, color=PAL[c], edgecolors='k', linewidths=0.5, zorder=5)
        ax.annotate(c[:4], (VINA[c], DG[c]['dG_kcal']), fontsize=7, ha='center', va='bottom', xytext=(0,5), textcoords='offset points')
ax.set_xlabel('Vina score (kcal/mol)')
ax.set_ylabel('MM-GBSA ΔG (kcal/mol)')
ax.set_title('c. Cross-Validation: Vina vs MM-GBSA', fontsize=11, loc='left')
ax.grid(alpha=0.2)
ax.axhline(y=0, color='#bbb', linewidth=0.5)

# (d) Ligand efficiency
ax = axs[3]
le_vals = [DG[c]['dG_kcal']/DG[c]['n_lig_atoms'] for c in COMPS]
sorted_le = sorted(zip(le_vals, COMPS), key=lambda x: x[0])
bars = ax.barh(range(len(sorted_le)), [v for v,c in sorted_le],
               color=[PAL[c] for v,c in sorted_le], edgecolor='white', height=0.6)
ax.set_yticks(range(len(sorted_le)))
ax.set_yticklabels([c for v,c in sorted_le], fontsize=8)
ax.set_xlabel('LE (kcal/mol/atom)')
ax.set_title('d. Ligand Efficiency (ΔG/atom)', fontsize=11, loc='left')
ax.grid(axis='x', alpha=0.2)

# (e) H-bond occupancy
ax = axs[4]
hb = [MET[c]['hbond_occupancy_pct'] for c in COMPS]
bars = ax.barh(range(len(COMPS)), hb, color=[PAL[c] for c in COMPS], edgecolor='white', height=0.6)
ax.set_yticks(range(len(COMPS)))
ax.set_yticklabels(COMPS, fontsize=8)
ax.set_xlabel('H-bond Occupancy (%)')
ax.set_title('e. Binding Stability (H-bond%)', fontsize=11, loc='left')
ax.set_xlim(0, 100)
ax.grid(axis='x', alpha=0.2)

# (f) RMSD stability
ax = axs[5]
rmsd = [MET[c]['rmsd_mean'] for c in COMPS]
drift = [MET[c]['rmsd_drift'] for c in COMPS]
x = np.arange(len(COMPS))
w = 0.3
ax.bar(x - w/2, rmsd, w, label='RMSD mean', color='steelblue', alpha=0.8)
ax.bar(x + w/2, drift, w, label='Drift (20ns)', color='coral', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(COMPS, fontsize=7, rotation=30)
ax.set_ylabel('RMSD (Å)')
ax.set_title('f. Structural Stability (RMSD)', fontsize=11, loc='left')
ax.legend(fontsize=7)
ax.grid(axis='y', alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'FIG_comprehensive_summary.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_comprehensive_summary.pdf'), dpi=300, bbox_inches='tight')
plt.close()
print('2. FIG_comprehensive_summary.png + .pdf')

# ================================================================
# FIGURE 3: DCCM grid at 300 DPI
# ================================================================
KEY = {"NAFTAZONE": "Most compact FEL",
       "ESKETAMINE": "Highest PC1",
       "RIBOFLAVIN": "Lowest RMSD",
       "HYDROMORPHONE": "Highest RMSD"}

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Dynamic Cross-Correlation (DCCM) — Domain-Blocked View',
             fontsize=14, fontweight='bold')
for i, (c, note) in enumerate(KEY.items()):
    corr = np.load(os.path.join(SMVT_MD, f'{c}_dccm.npy'))
    ax = axes[i//2, i%2]
    block = 50
    n_blocks = corr.shape[0] // block
    blocked = np.zeros((n_blocks, n_blocks))
    for bi in range(n_blocks):
        for bj in range(n_blocks):
            blocked[bi,bj] = corr[bi*block:(bi+1)*block, bj*block:(bj+1)*block].mean()
    im = ax.imshow(blocked, cmap='RdBu_r', vmin=-0.4, vmax=0.4, aspect='auto')
    ax.set_xlabel('Residue block (×50)'); ax.set_ylabel('Residue block (×50)')
    ax.set_title(f'{c}: {note}', fontsize=10)
    plt.colorbar(im, ax=ax, label='Pearson R', shrink=0.75)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(ANALYSIS, 'FIG_dccm_grid_pub.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_dccm_grid_pub.pdf'), dpi=300, bbox_inches='tight')
plt.close()
print('3. FIG_dccm_grid_pub.png + .pdf')

print('\n=== ALL COMPLETE ===')
