"""Final comparison: Vina vs MM-GBSA rankings + summary figure."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Load data
vina = {}
for name in ['NAFTAZONE','BIOTIN','ESKETAMINE','FUROSEMIDE',
             'GABAPENTIN_ENACARBIL','HYDROMORPHONE','PHENOBARBITAL','RIBOFLAVIN']:
    log = os.path.join(ANALYSIS, f'{name}_vina_log.txt')
    if os.path.exists(log):
        with open(log) as f:
            for ln in f:
                parts = ln.strip().split()
                if len(parts) >= 2 and parts[0] == '1':
                    vina[name] = float(parts[1])
                    break

mg = {}
for name in vina:
    fn = os.path.join(ANALYSIS, name, 'mmgbsa_results.json')
    if os.path.exists(fn):
        mg[name] = json.load(open(fn))['dG_kcal']

# Print comparison table
print('=' * 72)
print('  {:28s}  {:>8s}  {:>8s}  {:>9s}  {:>8s}'.format('Compound', 'Vina', 'MM-GBSA', 'Vina Rank', 'MG Rank'))
print('-' * 72)
v_sorted = sorted(vina.items(), key=lambda x: x[1])
m_sorted = sorted(mg.items(), key=lambda x: x[1])
for name in sorted(vina.keys()):
    vs = vina.get(name)
    ms = mg.get(name)
    vr = next(i+1 for i,(k,_) in enumerate(v_sorted) if k == name)
    mr = next(i+1 for i,(k,_) in enumerate(m_sorted) if k == name)
    print('  {:28s}  {:>+7.2f}  {:>+7.1f}  #{:<7d}  #{:<7d}'.format(name, vs, ms, vr, mr))

# Save full results
full = {k: {'vina': vina[k], 'mmgbsa': mg.get(k)} for k in vina}
json.dump(full, open(os.path.join(ANALYSIS, 'docking_comparison.json'), 'w'), indent=2)

# ============================================================
# Comparison figure: side-by-side bar chart
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('SMVT-8: Vina Docking vs MM-GBSA — Binding Affinity Comparison',
             fontsize=13, fontweight='bold')

palette = {'GABAPENTIN_ENACARBIL':'#2ecc71','RIBOFLAVIN':'#95a5a6',
           'BIOTIN':'#3498db','HYDROMORPHONE':'#9b59b6',
           'ESKETAMINE':'#1abc9c','PHENOBARBITAL':'#f39c12',
           'NAFTAZONE':'#e67e22','FUROSEMIDE':'#e74c3c'}

for idx, (dataset, title, unit) in enumerate([
        (vina, 'AutoDock Vina (rigid receptor)', 'kcal/mol'),
        (mg, 'MM-GBSA (amber14/GB(OBC2))', 'kcal/mol')]):
    ax = axes[idx]
    sorted_data = sorted(dataset.items(), key=lambda x: x[1])
    names = [x[0] for x in sorted_data]
    values = [x[1] for x in sorted_data]
    colors = [palette[n] for n in names]

    bars = ax.bar(range(len(names)), values, color=colors, edgecolor='white', linewidth=0.8, width=0.6)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=7, rotation=30, ha='right')
    ax.set_ylabel(f'Score ({unit})')
    ax.set_title(title, fontsize=10)
    ax.grid(axis='y', alpha=0.2)
    ax.axhline(y=0, color='#666', linewidth=0.5)

    for bar, v in zip(bars, values):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+(0.2 if v<0 else -0.8),
                f'{v:.2f}' if abs(v) < 10 else f'{v:.1f}',
                ha='center', fontsize=7, fontweight='bold', color='black')

    # Spearman rank correlation
    from scipy.stats import spearmanr
    # Compare only overlapping compounds
    common = [n for n in names if n in vina and n in mg]
    if common and idx == 1:
        x_vals = [vina[n] for n in common]
        y_vals = [mg[n] for n in common]
        r, p = spearmanr(x_vals, y_vals)
        ax.text(0.98, 0.02, f'Spearman rho = {r:.3f} (p={p:.3f})',
                transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
                bbox=dict(boxstyle='round', fc='white', alpha=0.7))

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_vina_vs_mmgbsa.png'), dpi=150, bbox_inches='tight')
plt.close()
print('\nALL_vina_vs_mmgbsa.png done')
