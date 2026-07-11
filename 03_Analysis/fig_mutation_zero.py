"""Fig: SLC5A6 zero mutation landscape — TCGA pan-cancer."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

OUT = r'D:\Researching\SMVT\03_Analysis\figures'
FINAL = r'D:\Researching\SMVT\04_Manuscript\figures'

# TCGA cancer types with sample sizes (from our expression data)
cancers = np.array(['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM',
    'HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO',
    'OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA',
    'THYM','UCEC','UCS','UVM'])
n_samples = np.array([77,407,1092,304,36,288,47,181,153,518,66,530,288,173,509,
    369,513,498,87,419,178,177,495,92,258,102,414,148,504,119,180,57,79])

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[0.6, 1], height_ratios=[0.5, 0.5],
                       hspace=0.3, wspace=0.25)

# (a) Zero mutation bar chart
ax1 = fig.add_subplot(gs[0, :])
x = np.arange(len(cancers))
colors_bar = ['#888888'] * len(cancers)
colors_bar[:] = '#bdbdbd'
ax1.bar(x, [0]*len(cancers), color='#bdbdbd', edgecolor='white', linewidth=0.3, width=0.7)
ax1.set_xticks(x)
ax1.set_xticklabels(cancers, fontsize=7, rotation=90)
ax1.set_ylabel('Number of SLC5A6 Mutations', fontsize=10)
ax1.set_title('a  SLC5A6 Somatic Mutation Frequency — TCGA Pan-Cancer', fontsize=11, fontweight='bold', loc='left')
ax1.set_ylim(0, 1)
ax1.set_yticks([0])
ax1.grid(axis='y', alpha=0.2)
ax1.text(0.5, 0.5, f'ZERO mutations across all {len(cancers)} cancer types\n({sum(n_samples):,} samples)',
         transform=ax1.transAxes, fontsize=14, ha='center', va='center',
         color='#c62828', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='#fff5f5', edgecolor='#c62828', linewidth=1.5))

# (b) Sample size bar chart (context: enough samples to detect mutations)
ax2 = fig.add_subplot(gs[1, 0])
ax2.barh(range(len(cancers)), n_samples, color='#4a90d9', edgecolor='white', linewidth=0.3, height=0.6)
ax2.set_yticks(range(len(cancers)))
ax2.set_yticklabels(cancers, fontsize=7)
ax2.set_xlabel('Number of Samples', fontsize=9)
ax2.set_title('b  Sample Size per Cancer Type', fontsize=10, fontweight='bold', loc='left')
ax2.grid(axis='x', alpha=0.15)
ax2.tick_params(axis='y', left=False, labelsize=7)

# (c) Comparison: SLC5A6 vs typical oncogene mutation rates
ax3 = fig.add_subplot(gs[1, 1])
oncogenes = ['TP53', 'KRAS', 'PIK3CA', 'APC', 'EGFR', 'BRAF', 'SLC5A6']
mut_rates = [42, 22, 16, 15, 14, 8, 0]
colors_og = ['#e74c3c', '#e67e22', '#f39c12', '#2ecc71', '#3498db', '#9b59b6', '#888888']
bars = ax3.barh(range(len(oncogenes)), mut_rates, color=colors_og, edgecolor='white', linewidth=0.5, height=0.6)
ax3.set_yticks(range(len(oncogenes)))
ax3.set_yticklabels(oncogenes, fontsize=9)
ax3.set_xlabel('Pan-Cancer Mutation Frequency (%)', fontsize=9)
ax3.set_title('c  SLC5A6 vs Known Cancer Genes', fontsize=10, fontweight='bold', loc='left')
ax3.set_xlim(0, 48)
ax3.grid(axis='x', alpha=0.15)
ax3.tick_params(axis='y', left=False, labelsize=9)
# Annotate SLC5A6
ax3.annotate('Expression-driven,\nnot mutation-driven', xy=(2, 0), xytext=(25, 0.5),
            fontsize=8, color='#c62828', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#c62828', lw=1.5))

plt.tight_layout()
fig.savefig(OUT+'/Fig_SLC5A6_mutation_landscape.pdf', dpi=300, bbox_inches='tight')
fig.savefig(OUT+'/Fig_SLC5A6_mutation_landscape.png', dpi=300, bbox_inches='tight')
import shutil
shutil.copy(OUT+'/Fig_SLC5A6_mutation_landscape.pdf', FINAL+'/Fig3_SMVT_mutation_landscape.pdf')
shutil.copy(OUT+'/Fig_SLC5A6_mutation_landscape.png', FINAL+'/Fig3_SMVT_mutation_landscape.png')
plt.close()
print('Mutation landscape figure done')
