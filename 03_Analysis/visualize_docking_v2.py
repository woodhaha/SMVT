"""
SMVT Docking Composite Figure — Redrawn for clarity.
Improved: spacing, font sizes, color contrast, hit annotations.
"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import os; os.chdir("D:/Researching/SMVT")

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 9,
    'axes.titlesize': 12, 'axes.labelsize': 10,
    'figure.dpi': 200, 'savefig.dpi': 300, 'savefig.bbox': 'tight'
})

df = pd.read_csv("03_Analysis/docking/docking_expanded_results.csv")
df = df.sort_values('best_affinity')

# Color scheme
C_CTRL = '#27AE60'    # natural substrate / control
C_NSAID = '#E74C3C'   # NSAID
C_FDA = '#2980B9'     # FDA drug
C_VIT = '#F39C12'     # vitamin
C_BARB = '#8E44AD'    # barbiturate (highlight)

def classify_color(row):
    name = str(row['name']).lower()
    t = str(row['type']).lower()
    # Detect barbiturates
    barb_names = ['phenobarbital', 'pentobarbital', 'secobarbital',
                  'butabarbital', 'amobarbital', 'barbital', 'metharbital']
    if any(b in name for b in barb_names):
        return C_BARB
    if 'barbit' in name:
        return C_BARB
    if t == 'substrate': return C_CTRL
    if t == 'nsaid': return C_NSAID
    if t == 'vitamin': return C_VIT
    return C_FDA

df['color'] = df.apply(classify_color, axis=1)

# ═══ FIGURE ═══
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(2, 2, height_ratios=[2.5, 1], hspace=0.35, wspace=0.25)

# ── AX1: Ranking bar chart ──
ax1 = fig.add_subplot(gs[0, :])
df_plot = df.iloc[::-1]
y_pos = range(len(df_plot))
vals = -df_plot['best_affinity'].values
colors = df_plot['color'].values

bars = ax1.barh(y_pos, vals, color=colors, edgecolor='white',
                linewidth=0.6, height=0.68, alpha=0.92)

# Thresholds
ax1.axvline(x=7.0, color='#E74C3C', linestyle='--', linewidth=1.5, alpha=0.5)
ax1.axvline(x=8.0, color='#C0392B', linestyle=':', linewidth=1, alpha=0.35)
ax1.axvline(x=6.0, color='#F39C12', linestyle='--', linewidth=1, alpha=0.35)

# Annotations
for i, (_, r) in enumerate(df_plot.iterrows()):
    score = -r['best_affinity']
    name = str(r['name']).replace('_', ' ')
    # Hit labels
    if score > 8.0:
        ax1.text(score + 0.15, len(df_plot)-1-i, '★ ELITE',
                fontsize=7, fontweight='bold', color='#C0392B', va='center')
    elif score > 7.0:
        ax1.text(score + 0.15, len(df_plot)-1-i, 'HIT',
                fontsize=6.5, fontweight='bold', color='#E74C3C', va='center')
    # Substrate marker
    if r['type'] == 'substrate':
        ax1.text(score + 0.15, len(df_plot)-1-i, '  (ctrl)',
                fontsize=6, color='#7F8C8D', va='center')

# Threshold labels
ax1.text(7.05, len(df_plot)-1.5, 'Hit threshold\nΔG < −7.0', fontsize=7.5,
         color='#E74C3C', alpha=0.8)
ax1.text(8.05, len(df_plot)-3, 'Elite\nΔG < −8.0', fontsize=7,
         color='#C0392B', alpha=0.7)

# Axis
ax1.set_yticks(y_pos)
ax1.set_yticklabels(df_plot['name'].str.replace('_', ' ').values,
                    fontfamily='monospace', fontsize=7.5)
ax1.set_xlabel('Binding Affinity −ΔG (kcal/mol)', fontsize=11)
ax1.set_title('SMVT (SLC5A6) — AutoDock Vina Docking Results\n49 Compounds | Exhaustiveness=16 | AlphaFold Receptor',
              fontsize=13, fontweight='bold')

# Legend
leg = [
    mpatches.Patch(color=C_BARB, label='Barbiturate (100% hit)'),
    mpatches.Patch(color=C_NSAID, label='NSAID'),
    mpatches.Patch(color=C_FDA, label='FDA Drug'),
    mpatches.Patch(color=C_VIT, label='Vitamin/Cofactor'),
    mpatches.Patch(color=C_CTRL, label='Natural Substrate'),
]
ax1.legend(handles=leg, loc='lower right', fontsize=8, framealpha=0.9, ncol=3)

# ── AX2: Chemical family comparison ──
ax2 = fig.add_subplot(gs[1, 0])

# Classify compounds into families
families = {}
for _, r in df.iterrows():
    name = str(r['name']).lower()
    score = -r['best_affinity']
    t = str(r['type']).lower()

    if any(b in name for b in ['phenobarbital', 'pentobarbital', 'butabarbital',
                                'secobarbital', 'amobarbital', 'barbital']):
        fam = 'Barbiturate'
    elif t == 'nsaid':
        fam = 'NSAID'
    elif 'opioid' in name or 'morph' in name or 'ketamine' in name:
        fam = 'Opioid/Ketamine'
    elif t == 'substrate':
        fam = 'Substrate'
    elif t == 'vitamin':
        fam = 'Vitamin'
    elif 'furosemide' in name or 'sulfonamide' in name or 'sulfa' in name:
        fam = 'Sulfonamide'
    else:
        fam = 'Other FDA'
    families.setdefault(fam, []).append(score)

# Order by median
fam_order = sorted(families.keys(), key=lambda k: np.median(families[k]), reverse=True)
fam_data = [families[k] for k in fam_order]
fam_colors = ['#8E44AD', '#E74C3C', '#9B59B6', '#3498DB', '#2ECC71', '#F39C12', '#95A5A6']

bp = ax2.boxplot(fam_data, patch_artist=True, widths=0.55, showfliers=True,
                 flierprops=dict(marker='o', markersize=5, alpha=0.5))

for patch, color in zip(bp['boxes'], fam_colors[:len(fam_order)]):
    patch.set_facecolor(color); patch.set_alpha(0.3)

# Jittered points
for i, vals in enumerate(fam_data):
    jitter = np.random.normal(0, 0.05, len(vals))
    ax2.scatter(np.ones(len(vals))*(i+1)+jitter, vals, s=30, alpha=0.75,
                color=fam_colors[i], edgecolors='white', linewidth=0.5, zorder=5)

# Count labels
for i, vals in enumerate(fam_data):
    n = len(vals)
    n_hit = sum(1 for v in vals if v > 7.0)
    best = max(vals)
    ax2.text(i+1, max(vals)+0.3, f'n={n}\n{n_hit} hit{"" if n_hit==1 else "s"}',
             ha='center', fontsize=7, fontweight='bold', color='#2C3E50')

ax2.set_xticklabels(fam_order, rotation=20, ha='right', fontsize=9)
ax2.set_ylabel('−ΔG (kcal/mol)', fontsize=10)
ax2.set_title('Affinity by Chemical Family', fontweight='bold', loc='left', fontsize=12)
ax2.axhline(y=7.0, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.4)
ax2.text(0.5, 7.1, 'hit', fontsize=7, color='#E74C3C', alpha=0.5)

# ── AX3: Summary stats ──
ax3 = fig.add_subplot(gs[1, 1])
ax3.axis('off')
ax3.set_title('Key Findings', fontweight='bold', loc='left', fontsize=12)

n_total = len(df)
n_hits = sum(-df['best_affinity'] > 7.0)
n_elite = sum(-df['best_affinity'] > 8.0)
n_barb = sum(1 for _, r in df.iterrows()
             if any(b in str(r['name']).lower() for b in
                    ['phenobarbital', 'pentobarbital', 'butabarbital',
                     'secobarbital', 'amobarbital', 'barbital']))
n_barb_hit = sum(1 for _, r in df.iterrows()
                 if any(b in str(r['name']).lower() for b in
                        ['phenobarbital', 'pentobarbital', 'butabarbital',
                         'secobarbital', 'amobarbital', 'barbital'])
                 and -r['best_affinity'] > 7.0)

findings = [
    f'Compounds Screened: {n_total}',
    f'Hits (ΔG < −7.0): {n_hits}',
    f'Elite Hits (ΔG < −8.0): {n_elite}',
    '',
    '★ Barbiturate Class',
    f'   {n_barb_hit}/{n_barb} hits — 100% hit rate',
    '   Barbituric acid = Biotin ureido mimic',
    '',
    '★ Novel SMVT Ligand Classes',
    '   • Opioids (Hydromorphone −8.58)',
    '   • Sulfonamides (Furosemide −8.36)',
    '   • Arylcyclohexylamines (Esketamine −7.58)',
    '',
    '★ Pharmacophore',
    '   Cyclic ureide/C=O core > Carboxyl',
    '   Aromatic ring + ≥2 HBA = optimal',
    '',
    'Biotin baseline: −6.76 kcal/mol',
]

y = 0.95
for line in findings:
    if line.startswith('★'):
        ax3.text(0.05, y, line, fontsize=9.5, fontweight='bold', color='#C0392B',
                 transform=ax3.transAxes, va='top')
        y -= 0.07
    elif line.startswith('  '):
        ax3.text(0.08, y, line, fontsize=8.5, color='#34495E',
                 transform=ax3.transAxes, va='top')
        y -= 0.05
    elif line == '':
        y -= 0.02
    else:
        ax3.text(0.05, y, line, fontsize=9, fontweight='bold', color='#2C3E50',
                 transform=ax3.transAxes, va='top')
        y -= 0.07

# Main title
fig.suptitle('SMVT (SLC5A6) Virtual Screening — Docking Results & SAR Analysis',
             fontsize=15, fontweight='bold', y=1.01)

# Save
os.makedirs("03_Analysis/figures", exist_ok=True)
fig.savefig("03_Analysis/figures/Fig_Docking_composite_v2.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Docking_composite_v2.pdf", facecolor='white')
plt.close()
print("✓ Fig_Docking_composite_v2 saved")
