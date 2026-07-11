"""Paper upgrade: per-residue comparison + H-bond occupancy chart."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
plt.rcParams.update({'font.size': 10, 'figure.dpi': 200})

ANALYSIS = r'D:\Researching\SMVT\MM-GBSA Analysis'
COMPS = ['BIOTIN','HYDROMORPHONE','GABAPENTIN_ENACARBIL','NAFTAZONE',
         'ESKETAMINE','FUROSEMIDE','PHENOBARBITAL','RIBOFLAVIN']
PALETTE = ['#3498db','#9b59b6','#2ecc71','#e67e22','#1abc9c','#e74c3c','#f39c12','#95a5a6']
SHORT = {'NAFTAZONE':'NAF','BIOTIN':'BIO','ESKETAMINE':'ESK','FUROSEMIDE':'FUR',
         'GABAPENTIN_ENACARBIL':'GAB','HYDROMORPHONE':'HYD','PHENOBARBITAL':'PHE','RIBOFLAVIN':'RIB'}

with open(os.path.join(ANALYSIS, 'md_master_metrics.json')) as f:
    METRICS = json.load(f)

# ================================================================
# FIGURE 1: Per-Residue Decomposition — BIOTIN vs ESKETAMINE
# ================================================================
with open(os.path.join(ANALYSIS, 'per_residue_decomp.json')) as f:
    DECOMP = json.load(f)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Per-Residue Energy Decomposition — Top Binding Pocket Hotspots',
             fontsize=13, fontweight='bold', y=1.02)

for idx, c in enumerate(['BIOTIN', 'ESKETAMINE']):
    ax = axes[idx]
    residues = DECOMP[c]['pocket_residues']
    names = [r['residue'] for r in residues][:12]
    totals = [r['dE_total'] for r in residues][:12]
    elecs = [r['dE_elec'] for r in residues][:12]
    vdws = [r['dE_vdw'] for r in residues][:12]

    x = np.arange(len(names))
    w = 0.25
    ax.bar(x - w, elecs, w, label='Elec.', color='#e74c3c', alpha=0.85)
    ax.bar(x, vdws, w, label='VDW', color='#3498db', alpha=0.85)
    ax.bar(x + w, totals, w, label='Total', color='#2c3e50', alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('Energy (kJ/mol)')
    ax.set_title(f'{c} ({"REF" if c=="BIOTIN" else "Test compound"})', fontsize=11)
    ax.legend(fontsize=7)
    ax.axhline(y=0, color='#888', linewidth=0.5)
    ax.grid(axis='y', alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'FIG_per_residue_decomp.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_per_residue_decomp.pdf'), dpi=200, bbox_inches='tight')
plt.close()
print('1. FIG_per_residue_decomp.png + .pdf done')

# ================================================================
# FIGURE 2: H-bond Occupancy Comparison (all 8 compounds)
# ================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Hydrogen Bond Analysis — 100ns MD Trajectories',
             fontsize=13, fontweight='bold', y=1.02)

x = np.arange(len(COMPS))
occ = [METRICS[c]['hbond_occupancy_pct'] for c in COMPS]
cnt = [METRICS[c]['hbond_mean_count'] for c in COMPS]

bars1 = ax1.bar(x, occ, color=PALETTE, edgecolor='white', linewidth=0.5)
ax1.set_xticks(x)
ax1.set_xticklabels([SHORT[c] for c in COMPS], fontsize=9)
ax1.set_ylabel('H-bond Occupancy (%)')
ax1.set_title('Occupancy (fraction of frames with contacts)')
ax1.grid(axis='y', alpha=0.2)
for bar, v in zip(bars1, occ):
    ax1.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
             f'{v:.0f}%', ha='center', va='bottom', fontsize=8)

bars2 = ax2.bar(x, cnt, color=PALETTE, edgecolor='white', linewidth=0.5)
ax2.set_xticks(x)
ax2.set_xticklabels([SHORT[c] for c in COMPS], fontsize=9)
ax2.set_ylabel('Mean H-bond Count')
ax2.set_title('Average number of concurrent H-bonds')
ax2.grid(axis='y', alpha=0.2)
for bar, v in zip(bars2, cnt):
    ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
             f'{v:.1f}', ha='center', va='bottom', fontsize=8)

# Color bar legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=PALETTE[i], label=f'{COMPS[i]}') for i in range(len(COMPS))]
ax2.legend(handles=legend_elements, fontsize=6, loc='upper left', ncol=2)

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'FIG_hbond_occupancy.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_hbond_occupancy.pdf'), dpi=200, bbox_inches='tight')
plt.close()
print('2. FIG_hbond_occupancy.png + .pdf done')

# ================================================================
# FIGURE 3: Pocket Contact Distance Heatmap (all 8 compounds)
# ================================================================
with open(os.path.join(ANALYSIS, 'pocket_contacts.json')) as f:
    CONTACTS = json.load(f)

# Get union of top pocket residues
all_res = []
for c in COMPS:
    for r in CONTACTS[c]['pocket_residues']:
        if r['residue'] not in all_res:
            all_res.append(r['residue'])

res_to_idx = {r:i for i,r in enumerate(all_res)}
n_res = len(all_res)
n_comp = len(COMPS)

dist_mat = np.full((n_comp, n_res), np.nan)
hbond_mat = np.zeros((n_comp, n_res), dtype=bool)

for ci, c in enumerate(COMPS):
    for r in CONTACTS[c]['pocket_residues']:
        ri = res_to_idx[r['residue']]
        dist_mat[ci, ri] = r['dist_A']
        if r.get('hbond'):
            hbond_mat[ci, ri] = True

# Filter to residues present in at least 4 compounds
keep = np.sum(~np.isnan(dist_mat), axis=0) >= 4
dist_mat = dist_mat[:, keep]
hbond_mat = hbond_mat[:, keep]
labels = [all_res[i] for i in range(n_res) if keep[i]]

fig, ax = plt.subplots(figsize=(14, 6))
cmap = plt.cm.YlOrRd_r
cmap.set_bad('#f0f0f0')
im = ax.imshow(dist_mat, cmap=cmap, aspect='auto', vmin=1.5, vmax=6.0)

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
ax.set_yticks(range(n_comp))
ax.set_yticklabels(COMPS, fontsize=8)

# Overlay H-bond markers
for ci in range(n_comp):
    for ri in range(len(labels)):
        if hbond_mat[ci, ri]:
            ax.plot(ri, ci, '*', color='gold', markersize=8, markeredgecolor='#333', markeredgewidth=0.3)

ax.set_title('Ligand-Protein Contact Distances (Å) — ★ = H-bond', fontsize=11, fontweight='bold')
fig.colorbar(im, ax=ax, shrink=0.8, label='Distance (Å)')

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'FIG_contact_heatmap.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'FIG_contact_heatmap.pdf'), dpi=200, bbox_inches='tight')
plt.close()
print('3. FIG_contact_heatmap.png + .pdf done')

print('\n=== ALL PAPER UPGRADE FIGURES DONE ===')
