"""Fig4: SLC5A6 mutation frequency (MC3) vs expression log2FC (recomputed DE)."""
import pandas as pd, numpy as np, re, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fig_style import apply_style, PALETTE
apply_style()

MAF = 'D:/Researching/SMVT/02_Data/cleaned/slc5a6_mc3_maf.csv'
DE = 'D:/Researching/SMVT/03_Analysis/outputs/de_table_fdr.json'
OUT = 'D:/Researching/SMVT/03_Analysis/figures/Fig4_SMVT_mutation_vs_expression.pdf'
FINAL = 'D:/Researching/SMVT/04_Manuscript/figures/Fig4_SMVT_mutation_vs_expression.pdf'

# code_map + samples from fig_mutation_final.py
src = open('D:/Researching/SMVT/03_Analysis/fig_mutation_final.py').read()
code_map = eval(re.search(r'code_map = (\{.*?\n\})', src, re.S).group(1))
samples = {'BRCA':1092,'LUAD':585,'LUSC':504,'BLCA':412,'COAD':458,'LIHC':373,'STAD':443,
           'HNSC':528,'KIRC':537,'ESCA':185,'SARC':261,'LGG':530,'PAAD':185,'CESC':307,
           'UCEC':548,'OV':419,'KIRP':291,'SKCM':470,'GBM':153,'PRAD':500}

df = pd.read_csv(MAF, sep='\t', dtype=str)
df['Cancer'] = df['Tumor_Sample_Barcode'].apply(lambda b: code_map.get(b.split('-')[1], 'OTHER'))

# mutation rate per cancer (unique patients / TCGA N)
mut = []
for c, n in samples.items():
    sub = df[df['Cancer'] == c]
    np_ = sub['Tumor_Sample_Barcode'].str[:12].nunique() if len(sub) else 0
    mut.append({'cancer': c, 'rate': np_ / n * 100, 'n_mut': np_})
mut = pd.DataFrame(mut)

# expression log2FC from recomputed DE
de = pd.DataFrame(json.load(open(DE)))
de = de[['cancer', 'log2FC']]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# (a) mutation frequency
m = mut.sort_values('rate', ascending=False)
cols = ['#e74c3c' if r > 2 else '#3498db' for r in m['rate']]
axes[0].bar(m['cancer'], m['rate'], color=cols, edgecolor='white', width=0.6)
axes[0].axhline(2, color='#888', ls='--', lw=0.8, label='2% threshold')
for i, (c, r) in enumerate(zip(m['cancer'], m['rate'])):
    axes[0].text(i, r + 0.08, f'{r:.1f}%', ha='center', fontsize=7)
axes[0].set_ylabel('Mutation frequency (% patients)', fontsize=10)
axes[0].set_title('a  SLC5A6 somatic mutation rate (MC3)', fontsize=11, fontweight='bold')
axes[0].legend(fontsize=8); axes[0].tick_params(axis='x', rotation=90, labelsize=8)
axes[0].grid(axis='y', alpha=0.15)

# (b) expression log2FC
de = de.sort_values('log2FC', ascending=False)
cols2 = ['#c62828' if x > 0 else '#1565c0' for x in de['log2FC']]
axes[1].bar(de['cancer'], de['log2FC'], color=cols2, edgecolor='white', width=0.6)
axes[1].axhline(0, color='#333', lw=0.8)
axes[1].set_ylabel('Log$_2$ fold change (tumor vs normal)', fontsize=10)
axes[1].set_title('b  SLC5A6 expression (recomputed, BH-FDR)', fontsize=11, fontweight='bold')
axes[1].tick_params(axis='x', rotation=90, labelsize=8)
axes[1].grid(axis='y', alpha=0.15)

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight')
plt.savefig(OUT.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
import shutil
shutil.copy(OUT, FINAL)
shutil.copy(OUT.replace('.pdf', '.png'), FINAL.replace('.pdf', '.png'))
print('Fig4 saved + synced')
