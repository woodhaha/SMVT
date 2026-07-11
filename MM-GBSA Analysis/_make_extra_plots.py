import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import pi
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

comps = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
         'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']

dg = {}
metrics = json.load(open('md_master_metrics.json'))
for c in comps:
    fn = f'{c}/mmgbsa_results.json'
    if os.path.exists(fn):
        dg[c] = json.load(open(fn))

palette = {'GABAPENTIN_ENACARBIL':'#2ecc71','RIBOFLAVIN':'#95a5a6',
           'BIOTIN':'#3498db','HYDROMORPHONE':'#9b59b6',
           'ESKETAMINE':'#1abc9c','PHENOBARBITAL':'#f39c12',
           'NAFTAZONE':'#e67e22','FUROSEMIDE':'#e74c3c'}
short = {'NAFTAZONE':'NAF','BIOTIN':'BIO','ESKETAMINE':'ESK','FUROSEMIDE':'FUR',
         'GABAPENTIN_ENACARBIL':'GAB','HYDROMORPHONE':'HYD','PHENOBARBITAL':'PHE','RIBOFLAVIN':'RIB'}

# ============================================================
# FIGURE 1: Multi-metric correlations with dG (3 panels)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
fig.suptitle('SMVT -- Binding Free Energy Correlations with MD Metrics',
             fontsize=14, fontweight='bold', y=1.02)

x_pairs = [
    ('rmsd_mean', 'RMSD_mean (A)', False),
    ('hbond_occupancy_pct', 'H-bond Occupancy (%)', False),
    ('sasa_mean_nm2', 'SASA_mean (nm2)', False),
]

for idx, (key, xlabel, logx) in enumerate(x_pairs):
    ax = axes[idx]
    for c in comps:
        dg_v = dg[c]['dG_kcal']
        mv = metrics[c][key]
        ax.scatter(mv, dg_v, s=120, c=palette[c], edgecolors='k', linewidths=0.5, zorder=5)
        ax.annotate(short[c], (mv, dg_v), fontsize=7, ha='center', va='bottom',
                    xytext=(0, 5), textcoords='offset points')

    xs = np.array([metrics[c][key] for c in comps])
    ys = np.array([dg[c]['dG_kcal'] for c in comps])
    r = np.corrcoef(xs, ys)[0,1]

    if abs(r) > 0.2:
        z = np.polyfit(xs, ys, 1)
        p = np.poly1d(z)
        x_line = np.linspace(xs.min(), xs.max(), 50)
        ax.plot(x_line, p(x_line), '--', color='#888', linewidth=0.8, alpha=0.5)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('dG_bind (kcal/mol)')
    ax.set_title(f'Pearson r = {r:.3f}', fontsize=10, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.axhline(y=0, color='#bbb', linewidth=0.5)

plt.tight_layout()
plt.savefig('ALL_correlations.png', dpi=150, bbox_inches='tight')
plt.close()
print('1. ALL_correlations.png done')

# ============================================================
# FIGURE 2: Ligand Efficiency bar + Multi-Metric Radar
# ============================================================
dg_all = [dg[c]['dG_kcal'] for c in comps]
lig_atoms = [dg[c]['n_lig_atoms'] for c in comps]
le = [dg_all[i]/lig_atoms[i] for i in range(len(comps))]
hbond_occ = [metrics[c]['hbond_occupancy_pct'] for c in comps]
hbond_cnt = [metrics[c]['hbond_mean_count'] for c in comps]
rmsd_vals = [metrics[c]['rmsd_mean'] for c in comps]
sasa_vals = [metrics[c]['sasa_mean_nm2'] for c in comps]

def normalize(vals, higher_better=True):
    arr = np.array(vals, dtype=float)
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.zeros_like(arr)
    if higher_better:
        return (arr - mn) / (mx - mn)
    else:
        return 1 - (arr - mn) / (mx - mn)

fig = plt.figure(figsize=(14, 6))
fig.suptitle('Ligand Efficiency + Normalized Metrics Profile', fontsize=14, fontweight='bold')

# Left: LE bar
ax1 = plt.subplot(1,2,1)
colors_le = [palette[c] for c in comps]
x = np.arange(len(comps))
bars = ax1.barh(x, le, color=colors_le, edgecolor='white', height=0.6)
ax1.set_yticks(x)
ax1.set_yticklabels([f'{c} ({dg[c]["dG_kcal"]:.1f})' for c in comps], fontsize=7)
ax1.set_xlabel('Ligand Efficiency (dG / heavy atom)')
ax1.set_title('Ligand Efficiency (lower = better per-atom binding)', fontsize=10)
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.2)
for bar, v in zip(bars, le):
    ax1.text(max(v-0.02, 0.01), bar.get_y()+bar.get_height()/2, f'{v:.3f}',
             ha='right' if v > 0.01 else 'left', va='center', fontsize=7,
             color='white' if v > 0.01 else 'black', fontweight='bold')

# Right: Radar
ax2 = plt.subplot(1,2,2, projection='polar')
norm_dg = normalize(dg_all, higher_better=True)
norm_rmsd = normalize(rmsd_vals, higher_better=False)
norm_hb = normalize(hbond_occ, higher_better=True)
norm_hbc = normalize(hbond_cnt, higher_better=True)
norm_sasa = normalize(sasa_vals, higher_better=False)

categories = ['dG_bind', 'RMSD (lower)', 'H-bond occ.', 'H-bond count', 'SASA (lower)']
N = len(categories)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

for i, c in enumerate(comps):
    vals = [norm_dg[i], norm_rmsd[i], norm_hb[i], norm_hbc[i], norm_sasa[i]]
    vals += vals[:1]
    ax2.plot(angles, vals, 'o-', linewidth=1.5, label=short[c], color=palette[c], alpha=0.8)
    ax2.fill(angles, vals, alpha=0.05, color=palette[c])

ax2.set_xticks(angles[:-1])
ax2.set_xticklabels(categories, fontsize=7)
ax2.set_ylim(0, 1.15)
ax2.set_title('Multi-Metric Radar (1=best within set)', fontsize=10, pad=10)
ax2.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), fontsize=7, ncol=1)

plt.tight_layout()
plt.savefig('ALL_ligand_efficiency_radar.png', dpi=150, bbox_inches='tight')
plt.close()
print('2. ALL_ligand_efficiency_radar.png done')

# ============================================================
# FIGURE 3: Multi-metric heatmap (z-score)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 5.5))

metric_names = ['dG_bind', 'LE', 'RMSD', 'H-bond%', 'H-bond#', 'SASA', 'Rg']

raw = np.zeros((len(comps), len(metric_names)))
for i, c in enumerate(comps):
    raw[i,0] = dg[c]['dG_kcal']
    raw[i,1] = dg[c]['dG_kcal'] / dg[c]['n_lig_atoms']
    raw[i,2] = metrics[c]['rmsd_mean']
    raw[i,3] = metrics[c]['hbond_occupancy_pct']
    raw[i,4] = metrics[c]['hbond_mean_count']
    raw[i,5] = metrics[c]['sasa_mean_nm2']
    raw[i,6] = metrics[c]['rg_mean_A']

z = (raw - raw.mean(axis=0)) / raw.std(axis=0)

# Flip: lower = bad, so negate RMSD, SASA, Rg
z[:,2] *= -1
z[:,5] *= -1
z[:,6] *= -1

im = ax.imshow(z, cmap='RdYlGn', aspect='auto', vmin=-2.5, vmax=2.5)

ax.set_xticks(range(len(metric_names)))
ax.set_xticklabels(metric_names, fontsize=9)
ax.set_yticks(range(len(comps)))
ax.set_yticklabels([f'{c} ({dg[c]["dG_kcal"]:.1f})' for c in comps], fontsize=8)

for i in range(len(comps)):
    for j in range(len(metric_names)):
        v = z[i,j]
        color = 'white' if abs(v) > 1.2 else 'black'
        ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=6.5, color=color)

ax.set_title('Multi-Metric Profile (z-score, green=better, red=worse)', fontsize=11, fontweight='bold')
fig.colorbar(im, ax=ax, shrink=0.8, label='z-score (sigma from mean)')

plt.tight_layout()
plt.savefig('ALL_metrics_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print('3. ALL_metrics_heatmap.png done')

print('\nALL 3 EXTRA FIGURES DONE')
print('Files: ALL_correlations.png, ALL_ligand_efficiency_radar.png, ALL_metrics_heatmap.png')
