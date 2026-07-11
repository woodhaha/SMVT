"""Figure 4: Energy decomposition bar chart."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Load merged decomposition results
ed = json.load(open(os.path.join(ANALYSIS, 'energy_decomposition.json')))

comps = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
         'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']

palette = {'GABAPENTIN_ENACARBIL':'#2ecc71','RIBOFLAVIN':'#95a5a6',
           'BIOTIN':'#3498db','HYDROMORPHONE':'#9b59b6',
           'ESKETAMINE':'#1abc9c','PHENOBARBITAL':'#f39c12',
           'NAFTAZONE':'#e67e22','FUROSEMIDE':'#e74c3c'}
short = {'NAFTAZONE':'NAF','BIOTIN':'BIO','ESKETAMINE':'ESK','FUROSEMIDE':'FUR',
         'GABAPENTIN_ENACARBIL':'GAB','HYDROMORPHONE':'HYD','PHENOBARBITAL':'PHE','RIBOFLAVIN':'RIB'}

dg = [ed[c]['dG_total'] for c in comps]
clj = [ed[c]['dE_coulomb_LJ'] for c in comps]
gb = [ed[c]['dG_GB_polar'] for c in comps]

fig, ax = plt.subplots(figsize=(13, 6))
fig.suptitle('SMVT — MM-GBSA Energy Decomposition (kJ/mol)', fontsize=14, fontweight='bold')

x = np.arange(len(comps))
w = 0.25

bars_dg = ax.bar(x - w, dg, w, color=[palette[c] for c in comps], edgecolor='white', linewidth=0.5, label='Delta G total', zorder=3)
bars_clj = ax.bar(x, clj, w, color='#2c3e50', edgecolor='white', linewidth=0.5, label='Coulomb + LJ (vdW)', alpha=0.8, zorder=3)
bars_gb = ax.bar(x + w, gb, w, color='#e74c3c', edgecolor='white', linewidth=0.5, label='GB polar (desolvation)', alpha=0.8, zorder=3)

# dG value labels
for bar, v in zip(bars_dg, dg):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height() + (3 if v>=0 else -8),
            f'{v:.0f}', ha='center', fontsize=8, fontweight='bold', color='black')

# CLJ labels
for bar, v in zip(bars_clj, clj):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height() - 8,
            f'{v:.0f}', ha='center', fontsize=6.5, color='white', fontweight='bold')

# GB labels
for bar, v in zip(bars_gb, gb):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height() + 3,
            f'+{v:.0f}', ha='center', fontsize=6.5, color='#c0392b', fontweight='bold')

ax.axhline(y=0, color='#555', linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels([f'{c}\n({short[c]})' for c in comps], fontsize=8)
ax.set_ylabel('Energy (kJ/mol)')
ax.legend(fontsize=8, loc='upper right')
ax.grid(axis='y', alpha=0.2)

# Annotations for key insights
# Hydromorphone best CLJ/GB ratio
ax.annotate('Best CLJ/GB\nbalance', xy=(4, clj[4]),
            xytext=(5.5, clj[4]-60), fontsize=8, color='#9b59b6',
            arrowprops=dict(arrowstyle='->', color='#9b59b6', lw=1.2))

# Furosemide massive GB penalty
ax.annotate('GB penalty exceeds\nCLJ gain', xy=(2, gb[2]),
            xytext=(1.5, gb[2]+40), fontsize=8, color='#e74c3c',
            arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=1.2))

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_energy_decomposition.png'), dpi=150, bbox_inches='tight')
plt.close()
print('ALL_energy_decomposition.png done')

# ============================================================
# FIGURE 5: Energy vs dG scatter with trend
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle('SMVT — Energy Component Correlations', fontsize=13, fontweight='bold')

for idx, (comp_name, comp_vals, xlab) in enumerate([
    ('dE_Coulomb+LJ', clj, 'Coulomb+LJ (kJ/mol)'),
    ('dG_GB_polar', gb, 'GB polar desolvation (kJ/mol)'),
]):
    ax = axes[idx]
    for i, c in enumerate(comps):
        ax.scatter(comp_vals[i], dg[i], s=100, c=palette[c], edgecolors='k', linewidths=0.5, zorder=5)
        ax.annotate(short[c], (comp_vals[i], dg[i]), fontsize=7, ha='center', va='bottom',
                    xytext=(0, 4), textcoords='offset points')

    xs = np.array(comp_vals)
    ys = np.array(dg)
    r = np.corrcoef(xs, ys)[0,1]
    z = np.polyfit(xs, ys, 1)
    p = np.poly1d(z)
    x_line = np.linspace(xs.min()-10, xs.max()+10, 50)
    ax.plot(x_line, p(x_line), '--', color='#888', linewidth=0.8, alpha=0.5)

    ax.set_xlabel(xlab)
    ax.set_ylabel('dG_total (kJ/mol)')
    ax.set_title(f'Pearson r = {r:.3f}', fontsize=10, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.axhline(y=0, color='#bbb', linewidth=0.5)

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_energy_correlations.png'), dpi=150, bbox_inches='tight')
plt.close()
print('ALL_energy_correlations.png done')

print('\nDONE')
