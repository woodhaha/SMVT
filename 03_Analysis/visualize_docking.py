#!/usr/bin/env python3
"""SMVT Docking Results Visualization"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import os; os.chdir(os.path.join(os.path.dirname(__file__), ".."))

plt.rcParams.update({'font.family':'sans-serif','font.size':10,'axes.titlesize':13,
    'axes.labelsize':11,'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight'})

df = pd.read_csv("03_Analysis/docking/docking_expanded_results.csv")
df = df.sort_values('best_affinity')

# Colors
c_sub = '#2ECC71'; c_nsaid = '#E74C3C'; c_fda = '#3498DB'; c_vit = '#F39C12'
def color_for(row):
    if row['type'] == 'substrate': return c_sub
    if row['type'] == 'nsaid': return c_nsaid
    if row['type'] == 'vitamin': return c_vit
    return c_fda
df['color'] = df.apply(color_for, axis=1)

# ═══ FIGURE: Docking affinity bar chart ═══
fig, ax = plt.subplots(figsize=(12, 8))
df_plot = df.iloc[::-1]  # bottom-to-top
bars = ax.barh(range(len(df_plot)), -df_plot['best_affinity'],
               color=df_plot['color'].values, edgecolor='white', linewidth=0.8, height=0.7)

# Highlight hits: affinity < -7
for i, (_, r) in enumerate(df_plot.iterrows()):
    if r['best_affinity'] < -7:
        ax.text(-r['best_affinity'] + 0.1, len(df_plot) - 1 - i, 'HIT',
                fontsize=9, fontweight='bold', color='#E74C3C', va='center')
    if r['type'] == 'substrate':
        ax.text(-r['best_affinity'] + 0.1, len(df_plot) - 1 - i, '*(ctrl)',
                fontsize=8, color='gray', va='center' if r['best_affinity'] > -7 else 'bottom')

# Threshold lines
ax.axvline(x=7, color='#E74C3C', linestyle='--', alpha=0.4, linewidth=0.8)
ax.text(7.1, len(df_plot)-0.5, '< -7 kcal/mol', fontsize=8, color='#E74C3C', alpha=0.7)
ax.axvline(x=6, color='#F39C12', linestyle='--', alpha=0.4, linewidth=0.8)

ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot['name'].values, fontfamily='monospace')
ax.set_xlabel('Binding Affinity — more negative = stronger binding (kcal/mol)', fontsize=11)
ax.set_title('SMVT (SLC5A6) Molecular Docking — Expanded FDA Drug Screening\nAutoDock Vina | 49 compounds | exhaustiveness=16',
             fontsize=13, fontweight='bold')

# Legend
leg = [mpatches.Patch(color=c_sub, label='Natural Substrate'),
       mpatches.Patch(color=c_nsaid, label='NSAID Inhibitor'),
       mpatches.Patch(color=c_fda, label='FDA Drug'),
       mpatches.Patch(color=c_vit, label='Vitamin / Cofactor')]
ax.legend(handles=leg, loc='lower right', fontsize=9, framealpha=0.9)
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Docking_ranking.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Docking_ranking.pdf", facecolor='white')
plt.close()
print("Fig_Docking_ranking OK")

# ═══ FIGURE: Type comparison boxplot ═══
fig, ax = plt.subplots(figsize=(7, 5))
types_data = {
    'Natural\nSubstrates': -df[df['type']=='substrate']['best_affinity'].values,
    'NSAID\nInhibitors': -df[df['type']=='nsaid']['best_affinity'].values,
    'FDA\nDrugs': -df[df['type']=='fda']['best_affinity'].values,
    'Vitamin\nCofactors': -df[df['type']=='vitamin']['best_affinity'].values,
}
bp = ax.boxplot(types_data.values(), patch_artist=True, widths=0.5)
for patch, color in zip(bp['boxes'], [c_sub, c_nsaid, c_fda, c_vit]):
    patch.set_facecolor(color); patch.set_alpha(0.3)
for i, (label, vals) in enumerate(types_data.items()):
    colors_4 = [c_sub, c_nsaid, c_fda, c_vit]
    jitter = np.random.normal(0, 0.04, len(vals))
    ax.scatter(np.ones(len(vals))*(i+1)+jitter, vals, s=40, alpha=0.7,
               color=colors_4[i], edgecolors='white', linewidth=0.5, zorder=5)

ax.set_xticklabels(types_data.keys())
ax.set_ylabel('Binding Affinity — more negative = stronger (kcal/mol)', fontsize=10)
ax.set_title('Docking Affinity by Compound Type (49 compounds)', fontsize=12, fontweight='bold')
ax.axhline(y=7, color='#E74C3C', linestyle='--', alpha=0.3)
ax.text(0.6, 7.1, 'hit threshold', fontsize=7, color='#E74C3C', alpha=0.5)
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Docking_type_comparison.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Docking_type_comparison.pdf", facecolor='white')
plt.close()
print("Fig_Docking_type_comparison OK")

# ═══ FIGURE: Combined panel ═══
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Ranking
df_plot2 = df.iloc[::-1]
ax1.barh(range(len(df_plot2)), -df_plot2['best_affinity'],
         color=df_plot2['color'].values, edgecolor='white', linewidth=0.6, height=0.65)
ax1.set_yticks(range(len(df_plot2)))
ax1.set_yticklabels(df_plot2['name'].values, fontfamily='monospace', fontsize=8)
ax1.set_xlabel('Binding Affinity (kcal/mol)', fontsize=10)
ax1.axvline(x=7, color='#E74C3C', linestyle='--', alpha=0.4)
ax1.set_title('Compound Ranking', fontweight='bold', loc='left')
leg = [mpatches.Patch(color=c_sub, label='Substrate'), mpatches.Patch(color=c_nsaid, label='NSAID'), mpatches.Patch(color=c_fda, label='FDA'), mpatches.Patch(color=c_vit, label='Vitamin')]
ax1.legend(handles=leg, loc='lower right', fontsize=7, framealpha=0.9)

# Right: Type comparison
bp2 = ax2.boxplot(types_data.values(), patch_artist=True, widths=0.45)
for patch, color in zip(bp2['boxes'], [c_sub, c_nsaid, c_fda, c_vit]):
    patch.set_facecolor(color); patch.set_alpha(0.25)
for i, (label, vals) in enumerate(types_data.items()):
    colors_4 = [c_sub, c_nsaid, c_fda, c_vit]
    jitter = np.random.normal(0, 0.04, len(vals))
    ax2.scatter(np.ones(len(vals))*(i+1)+jitter, vals, s=50, alpha=0.8,
                color=colors_4[i], edgecolors='white', linewidth=0.5, zorder=5)
    ax2.text(i+1, np.max(vals)+0.5, f'best: {np.max(vals):.1f}', ha='center', fontsize=8, fontweight='bold')
ax2.set_xticklabels(types_data.keys())
ax2.set_ylabel('Binding Affinity (kcal/mol)', fontsize=10)
ax2.set_title('Affinity by Compound Type', fontweight='bold', loc='left')
ax2.axhline(y=7, color='#E74C3C', linestyle='--', alpha=0.3)

# Annotation
ax2.annotate(f'Diclofenac\n(NSAID, -7.15)', xy=(2, 7.15), xytext=(2.4, 8.5),
             fontsize=8, fontweight='bold', ha='center', color='#C0392B',
             arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=1))
ax2.annotate(f'Aspirin\n(FDA, -6.68)', xy=(3, 6.68), xytext=(3.4, 7.8),
             fontsize=8, fontweight='bold', ha='center',
             arrowprops=dict(arrowstyle='->', color='#2980B9', lw=1))

fig.suptitle('SMVT (SLC5A6) Molecular Docking — Expanded Virtual Screening\nAutoDock Vina | 49 compounds | exhaustiveness=16 | SMVT AlphaFold structure',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Docking_composite.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Docking_composite.pdf", facecolor='white')
plt.close()
print("Fig_Docking_composite OK")

print("\n=== 3 DOCKING FIGURES GENERATED ===")
print("  Fig_Docking_ranking.png/pdf")
print("  Fig_Docking_type_comparison.png/pdf")
print("  Fig_Docking_composite.png/pdf")
