"""Figure 1: Pan-cancer SLC5A6 differential expression."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = r'D:\Researching\SMVT\04_Manuscript\figures'
de = pd.read_csv(r'D:\Researching\SMVT\03_Analysis\outputs\pancan_de.csv')
de = de.sort_values('log2FC', ascending=True)

colors = []
for _, r in de.iterrows():
    if r['sig'] == '*' and r['log2FC'] > 0:
        colors.append('#d62728')
    elif r['sig'] == '*' and r['log2FC'] < 0:
        colors.append('#1f77b4')
    else:
        colors.append('#cccccc')

fig, ax = plt.subplots(figsize=(8, 7))
ax.barh(range(len(de)), de['log2FC'].values, color=colors,
        edgecolor='white', linewidth=0.4, height=0.65)
ax.axvline(0, color='black', linewidth=0.6)
ax.axvline(0.2, color='gray', linestyle='--', linewidth=0.4, alpha=0.4)

yticklabels = []
for _, r in de.iterrows():
    label = r['cancer'] + ' (' + str(r['n_t']) + '/' + str(r['n_n']) + ')'
    if r['sig'] == '*':
        if r['pval'] < 1e-20:
            label += ' ***'
        elif r['pval'] < 1e-10:
            label += ' **'
        else:
            label += ' *'
    yticklabels.append(label)

ax.set_yticks(range(len(de)))
ax.set_yticklabels(yticklabels, fontsize=8)
ax.set_xlabel('Log$_2$ Fold Change (Tumor vs Normal)', fontsize=10)
ax.set_title('SLC5A6 Pan-Cancer Differential Expression (TCGA)', fontsize=11, fontweight='bold', loc='left')

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#d62728', label='Upregulated (P < 0.05)'),
    Patch(facecolor='#1f77b4', label='Downregulated (P < 0.05)'),
    Patch(facecolor='#cccccc', label='Not significant'),
]
ax.legend(handles=legend_elements, fontsize=7, loc='lower right')
ax.set_xlim(de['log2FC'].min() - 0.15, de['log2FC'].max() + 0.35)
ax.grid(axis='x', alpha=0.15)
plt.tight_layout()
fig.savefig(OUT + '/Fig1_SMVT_TCGA_pan_cancer.pdf', dpi=300, bbox_inches='tight')
fig.savefig(OUT + '/Fig1_SMVT_TCGA_pan_cancer.png', dpi=300, bbox_inches='tight')
plt.close()
print('Fig1 saved')
sig_up = sum(1 for _, r in de.iterrows() if r['sig'] == '*' and r['log2FC'] > 0)
sig_down = sum(1 for _, r in de.iterrows() if r['sig'] == '*' and r['log2FC'] < 0)
print('UP: ' + str(sig_up) + ', DOWN: ' + str(sig_down))
