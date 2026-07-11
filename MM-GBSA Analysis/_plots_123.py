"""Figures 1-3: RMSF flexibility map, 100ns time series, PCA scatter."""
import json, os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))

comps = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
         'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']

palette = {'GABAPENTIN_ENACARBIL':'#2ecc71','RIBOFLAVIN':'#95a5a6',
           'BIOTIN':'#3498db','HYDROMORPHONE':'#9b59b6',
           'ESKETAMINE':'#1abc9c','PHENOBARBITAL':'#f39c12',
           'NAFTAZONE':'#e67e22','FUROSEMIDE':'#e74c3c'}
short = {'NAFTAZONE':'NAF','BIOTIN':'BIO','ESKETAMINE':'ESK','FUROSEMIDE':'FUR',
         'GABAPENTIN_ENACARBIL':'GAB','HYDROMORPHONE':'HYD','PHENOBARBITAL':'PHE','RIBOFLAVIN':'RIB'}
metrics = json.load(open(os.path.join(ANALYSIS, 'md_master_metrics.json')))

# ============================================================
# FIGURE 1: Per-Residue Flexibility Consensus Map
# ============================================================
all_flex = {}
max_res = 0
for c in comps:
    flex = set(metrics[c]['flexible_residues'])
    all_flex[c] = flex
    if flex: max_res = max(max_res, max(flex))

n_res = max_res
flex_matrix = np.zeros((len(comps), n_res), dtype=float)
for i, c in enumerate(comps):
    for r in all_flex[c]:
        if r <= n_res: flex_matrix[i, r-1] = 1.0

consensus = flex_matrix.sum(axis=0)

fig = plt.figure(figsize=(16, 5))
fig.suptitle('SMVT — Per-Residue Flexibility (RMSF-driven, 8 compounds)',
             fontsize=14, fontweight='bold', y=0.98)

gs = fig.add_gridspec(2, 1, height_ratios=[1, 2.5], hspace=0.08)

# Top consensus
ax0 = fig.add_subplot(gs[0])
ax0.fill_between(range(1, n_res+1), consensus, alpha=0.6, color='#e74c3c', step='mid')
ax0.axhline(y=4, color='#e74c3c', linestyle='--', linewidth=0.8, alpha=0.5)
ax0.set_ylabel('Flexible\ncompounds', fontsize=8)
ax0.set_ylim(0, 8.5)
ax0.set_xlim(0, n_res)
ax0.set_xticks([])
ax0.tick_params(labelsize=7)
ax0.grid(axis='y', alpha=0.2)

# Bottom heatmap
ax1 = fig.add_subplot(gs[1])
ax1.imshow(flex_matrix, aspect='auto', cmap='Reds', vmin=0, vmax=1,
           extent=[1, n_res, len(comps)-0.5, -0.5])
ax1.set_yticks(range(len(comps)))
ax1.set_yticklabels([f'{c} ({metrics[c]["rmsf_mean"]:.1f}A)' for c in comps], fontsize=7)
ax1.set_xlabel('Residue Index', fontsize=9)
ax1.set_xlim(0, n_res)
ax1.tick_params(labelsize=7)

# Hotspot regions (>=6 compounds)
hotspots = [(i+1, int(consensus[i])) for i in range(n_res) if consensus[i] >= 6]
if hotspots:
    regions = [[hotspots[0][0], hotspots[0][0]]]
    for r, _ in hotspots[1:]:
        if r == regions[-1][1] + 1: regions[-1][1] = r
        else: regions.append([r, r])
    for s, e in regions:
        ax0.axvspan(s, e, alpha=0.15, color='red')
        ax1.axvspan(s, e, alpha=0.15, color='red')
        ax1.text((s+e)/2, len(comps)-0.3, f'{s}-{e}', ha='center', fontsize=5.5,
                 color='#c0392b', fontweight='bold')

# Domain markers
for start, end, label in [(1, 40, 'N-term'), (470, n_res, 'C-term')]:
    for ax_i in [ax0, ax1]:
        ax_i.axvspan(start, end, alpha=0.06, color='#2ecc71')
    ax1.text((start+end)/2, -0.5, label, ha='center', fontsize=6, color='#2ecc71',
             transform=ax1.get_xaxis_transform())

plt.savefig(os.path.join(ANALYSIS, 'ALL_flexibility_map.png'), dpi=150, bbox_inches='tight')
plt.close()
print('1. ALL_flexibility_map.png done')

# ============================================================
# FIGURE 2: 100ns Time Series — RMSD + Potential Energy
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(14, 7.5), sharex=True)
fig.suptitle('SMVT 100ns MD — Time Series', fontsize=14, fontweight='bold')

ax = axes[0]
for c in comps:
    fn = os.path.join(ANALYSIS, f'{c}_100ns_full.csv')
    if not os.path.exists(fn): continue
    try:
        data = np.loadtxt(fn, delimiter=',', skiprows=1)
        ax.plot(data[:,0], data[:,1], color=palette[c], linewidth=0.5, alpha=0.8)
        ax.annotate(short[c], (data[-1,0], data[-1,1]), fontsize=6, color=palette[c],
                    ha='left', va='center', xytext=(2, 0), textcoords='offset points')
    except: pass
ax.set_ylabel(r'C$\alpha$ RMSD ($\AA$)', fontsize=9)
ax.legend(fontsize=6, ncol=4, loc='upper right')
ax.grid(alpha=0.2)

ax = axes[1]
for c in comps:
    fn = os.path.join(ANALYSIS, f'{c}_100ns.csv')
    if not os.path.exists(fn): continue
    try:
        with open(fn) as f:
            reader = csv.reader(f)
            next(reader)
            rows = list(reader)
        times = np.array([float(r[1])/1000 for r in rows if r])
        potential = np.array([float(r[2])/1000 for r in rows if r])
        ax.plot(times, potential, color=palette[c], linewidth=0.3, alpha=0.7)
    except: pass

ax.set_xlabel('Time (ns)', fontsize=9)
ax.set_ylabel('Potential E (1000 kJ/mol)', fontsize=9)
ax.legend(fontsize=5.5, ncol=4, loc='lower right')
ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_timeseries.png'), dpi=150, bbox_inches='tight')
plt.close()
print('2. ALL_timeseries.png done')

# ============================================================
# FIGURE 3: PCA Projection — PC1 vs PC2
# ============================================================
fig, ax = plt.subplots(figsize=(11, 8))
fig.suptitle('SMVT 100ns MD — PCA Projection (Calpha, protein only)',
             fontsize=14, fontweight='bold')

for c in comps:
    fn = os.path.join(ANALYSIS, f'{c}_pca.npz')
    if not os.path.exists(fn): continue
    d = np.load(fn)
    ax.scatter(d['pc1'], d['pc2'], s=4, c=palette[c], alpha=0.4,
               label=f'{short[c]} ({c})', zorder=2)

    # Mark median as pseudo-centroid
    cx, cy = np.median(d['pc1']), np.median(d['pc2'])
    ax.scatter(cx, cy, s=100, c=palette[c], edgecolors='k', linewidths=1.2,
               marker='X', zorder=5)
    ax.annotate(short[c], (cx, cy), fontsize=9, fontweight='bold', color=palette[c],
                ha='center', va='bottom', xytext=(0, 6), textcoords='offset points')

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.legend(fontsize=6, ncol=2, loc='upper right', markerscale=3)
ax.grid(alpha=0.2)

var = np.load(os.path.join(ANALYSIS, 'BIOTIN_pca.npz'))['var_exp']
ax.text(0.98, 0.02, 'PC1={:.0f}%, PC2={:.0f}% (BIOTIN ref)'.format(var[0], var[1]),
        transform=ax.transAxes, fontsize=8, ha='right', color='#555')

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_pca_scatter.png'), dpi=150, bbox_inches='tight')
plt.close()
print('3. ALL_pca_scatter.png done')

print('\nALL 3 FIGURES DONE')
