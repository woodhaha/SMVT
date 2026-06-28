#!/usr/bin/env python3
"""SMVT — Final recommended visualizations for publication"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os; os.chdir("D:/Researching/SMVT")

plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.titlesize':11,
    'axes.labelsize':9,'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight'})

# Color scheme
c_tca = '#E74C3C'; c_bio = '#2ECC71'; c_slc = '#3498DB'
c_mem = '#9B59B6'; c_lip = '#F39C12'; c_gray = '#95A5A6'
smvt_red = '#C0392B'

# Load data
pearson = pd.read_csv("03_Analysis/outputs/scTenifoldKnk_DRGs.csv")
pc_val = pd.read_csv("03_Analysis/outputs/PC_regression_vs_Pearson_validation.csv")
surv = pd.read_csv("03_Analysis/outputs/survival_results.csv")
coad = pd.read_csv("03_Analysis/outputs/coad_survival_results.csv")

# ═══════════════════════════════════════════════════════════════
# FIGURE X: Method comparison scatter (Pearson vs PC regression)
# ═══════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 7))
merged = pc_val[pc_val['gene'] != 'SLC5A6']

# Color by module
def module(g):
    tca={'CS','PDHA1','PDHB','PDHX','PC','ACLY','DLAT'}
    bio={'HLCS','ACACB','PCCA','MCCC1','BTD','ACACA','FASN','SCD'}
    slc={'SLC19A2','SLC26A4','SLC22A12','SLC5A3','SLC23A1','SLC5A7'}
    mem={'PDZD11','CFTR','CHAT'}
    if g in tca: return c_tca
    if g in bio: return c_bio
    if g in slc: return c_slc
    if g in mem: return c_mem
    return c_lip if g == 'SREBF1' else c_gray

colors = [module(g) for g in merged['gene']]

# Main scatter
ax.scatter(merged['pearson_impact'], merged['pc_regression_impact'],
           c=colors, s=120, edgecolors='white', linewidth=1, zorder=5, alpha=0.85)

# Zoomed inset for crowded consensus region
consensus_genes = ['CS','SLC19A2','DPH2','SLC26A4','PDHA1','ACACB']
# Highlight the crowded zone with a rectangle
from matplotlib.patches import Rectangle
rect = Rectangle((0.375, 0.035), 0.055, 0.065, linewidth=1, edgecolor='#2C3E50',
                 facecolor='none', linestyle='-', zorder=4)
ax.add_patch(rect)
ax.annotate('', xy=(0.50, 0.10), xytext=(0.43, 0.10),
            arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=0.8))

# Inset axes: zoom into the crowded region
ax_inset = ax.inset_axes([0.48, 0.55, 0.48, 0.40])
ax_inset.set_xlim(0.375, 0.430)
ax_inset.set_ylim(0.035, 0.100)

# Plot all genes in inset (lighter)
for _, r in merged.iterrows():
    ax_inset.scatter(r['pearson_impact'], r['pc_regression_impact'],
                     s=30, c='lightgray', edgecolors='white', linewidth=0.5, alpha=0.5)

# Consensus genes highlighted in inset
offsets = {
    'CS':       (0.002, 0.002),
    'SLC19A2':  (-0.006, 0.003),
    'DPH2':     (0.003, -0.003),
    'SLC26A4':  (0.001, -0.003),
    'PDHA1':    (-0.004, -0.001),
    'ACACB':    (-0.004, 0.003),
}
for g in consensus_genes:
    r = merged[merged['gene'] == g].iloc[0]
    px, py = r['pearson_impact'], r['pc_regression_impact']
    ax_inset.scatter(px, py, s=80, c='#2C3E50', edgecolors='white', linewidth=1.5, zorder=10)
    ox, oy = offsets[g]
    ax_inset.annotate(g, (px, py), (px+ox, py+oy), fontsize=7.5, fontweight='bold',
                      ha='center', va='center',
                      bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.9,
                                edgecolor='#2C3E50', linewidth=0.5))

ax_inset.set_title('Zoom: Consensus DRGs', fontsize=8, fontweight='bold', color='#2C3E50')
ax_inset.tick_params(labelsize=6)
ax.indicate_inset_zoom(ax_inset, edgecolor='#2C3E50', alpha=0.8)

# Diagonal
mn = min(merged['pearson_impact'].min(), merged['pc_regression_impact'].min()) - 0.01
mx = max(merged['pearson_impact'].max(), merged['pc_regression_impact'].max()) + 0.01
ax.plot([mn, mx], [mn, mx], '--', color='gray', alpha=0.5, linewidth=0.8)

# Correlation annotation
from scipy.stats import spearmanr
rho, p = spearmanr(merged['pearson_impact'], merged['pc_regression_impact'])
ax.text(0.95, 0.95, f"Spearman rho = {rho:.3f}\np = {p:.1e}",
        transform=ax.transAxes, ha='right', va='top', fontsize=10, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_xlabel('Pearson Co-expression Impact', fontsize=11)
ax.set_ylabel('PC Regression Impact\n(scTenifoldKnk pcNet equivalent)', fontsize=11)
ax.set_title('Method Validation: Pearson vs PC Regression Virtual KO', fontsize=13, fontweight='bold')

# Legend
from matplotlib.lines import Line2D
leg = [Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markersize=10,label=l)
       for c,l in [(c_tca,'TCA Cycle'),(c_bio,'Biotin Carboxylase'),(c_slc,'SLC Transporter'),
                    (c_mem,'Membrane Polarity'),(c_lip,'Lipid Metabolism')]]
ax.legend(handles=leg, loc='upper left', fontsize=8, framealpha=0.9,
          bbox_to_anchor=(0.01, 0.99))
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Method_Validation_scatter.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Method_Validation_scatter.pdf", facecolor='white')
plt.close()
print("Fig_Method_Validation_scatter OK")

# ═══════════════════════════════════════════════════════════════
# FIGURE Y: Multi-omics integration summary (radar-style bar)
# ═══════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Evidence strength by analysis type
evidence = {
    'TCGA Expression\n(6 cancer types)': 5,
    'Survival Analysis\n(4 significant)': 4,
    'STRING PPI\n(PDZD11 0.969)': 4,
    'Pathway Enrichment\n(GO/KEGG/Reactome)': 4,
    'Virtual KO\n(qualitative + quantitative)': 5,
    'PC Regression\nValidation': 4,
    'Mutation Analysis\n(zero TCGA mutations)': 5,
    'scRNA Localization\n(HPA + GSE178341)': 3,
}
labels = list(evidence.keys())
values = list(evidence.values())
colors_bar = ['#27AE60' if v >= 4 else '#F39C12' for v in values]

ax1.barh(range(len(labels)), values, color=colors_bar, edgecolor='white', height=0.6)
ax1.set_yticks(range(len(labels)))
ax1.set_yticklabels(labels, fontsize=9)
ax1.set_xlim(0, 6)
ax1.set_xlabel('Evidence Strength (1-5)', fontsize=10)
ax1.set_title('Multi-Omics Evidence Summary', fontweight='bold', fontsize=12)
for i, v in enumerate(values):
    ax1.text(v + 0.1, i, '★★★★★'[:v] + '☆'*(5-v), fontsize=8, va='center')

# Right: Consensus DRG ranking with validation
consensus_genes = ['CS','SLC19A2','DPH2','SLC26A4','PDHA1','ACACB',
                   'PC','HLCS','PDHB','PDZD11','PDHX','ACLY','PCCA','BTD','MCCC1']
consensus_data = merged[merged['gene'].isin(consensus_genes)].set_index('gene')
mean_impact = consensus_data[['pearson_impact','pc_regression_impact']].mean(axis=1)
mean_impact = mean_impact.sort_values(ascending=True)

colors_c = [module(g) for g in mean_impact.index]
ax2.barh(range(len(mean_impact)), mean_impact.values, color=colors_c,
         edgecolor='white', height=0.7)
ax2.set_yticks(range(len(mean_impact)))
ax2.set_yticklabels(mean_impact.index, fontfamily='monospace', fontsize=9)
ax2.set_xlabel('Mean Impact (Pearson + PC Regression)', fontsize=10)
ax2.set_title('Consensus DRG Ranking', fontweight='bold', fontsize=12)

fig.suptitle('SMVT (SLC5A6) — Multi-Omics Integration & Cross-Method Validation',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Evidence_Summary.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Evidence_Summary.pdf", facecolor='white')
plt.close()
print("Fig_Evidence_Summary OK")

# ═══════════════════════════════════════════════════════════════
# FIGURE Z: Pan-cancer survival forest plot (publication ready)
# ═══════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 6))
surv_plot = surv.sort_values('cox_hr', ascending=True)

y_positions = range(len(surv_plot))
colors_s = ['#E74C3C' if p < 0.05 else '#95A5A6' for p in surv_plot['cox_p']]
for i, (_, r) in enumerate(surv_plot.iterrows()):
    ax.errorbar(r['cox_hr'], i, xerr=[[r['cox_hr'] - r['cox_hr_lower']], [r['cox_hr_upper'] - r['cox_hr']]],
                fmt='o', color=colors_s[i], capsize=3, markersize=8, linewidth=1.5)
ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
ax.set_yticks(y_positions)
ax.set_yticklabels([f"{r['cancer_type']} (n={int(r['n_total'])})" for _, r in surv_plot.iterrows()], fontsize=9)
ax.set_xlabel('Hazard Ratio (High vs Low SMVT)', fontsize=10)
ax.set_title('Pan-Cancer Survival — SMVT Expression Stratification', fontsize=12, fontweight='bold')

for i, (_, r) in enumerate(surv_plot.iterrows()):
    p = r['cox_p']
    if p < 0.001: sig = '***'
    elif p < 0.01: sig = '**'
    elif p < 0.05: sig = '*'
    else: sig = 'ns'
    ax.text(r['cox_hr_upper'] + 0.05, i, sig, fontsize=9, fontweight='bold', va='center',
            color='#E74C3C' if sig != 'ns' else '#95A5A6')

plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_PanCancer_Forest.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_PanCancer_Forest.pdf", facecolor='white')
plt.close()
print("Fig_PanCancer_Forest OK")

print("\n=== 3 NEW FIGURES GENERATED ===")
for f in ['Fig_Method_Validation_scatter','Fig_Evidence_Summary','Fig_PanCancer_Forest']:
    print(f"  {f}.png/pdf")
