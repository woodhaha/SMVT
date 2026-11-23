"""Fig4: SLC5A6 mutation frequency (MC3) vs expression log2FC (recomputed DE)."""
import os
import pandas as pd, numpy as np, re, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fig_style import apply_style, PALETTE, MM, DOUBLE, save_pub
apply_style()

MAF = 'D:/Researching/SMVT/02_Data/cleaned/slc5a6_mc3_maf.csv'
DE = 'D:/Researching/SMVT/03_Analysis/outputs/de_table_fdr.json'
STEM = 'Fig4_SMVT_mutation_vs_expression'
OUT = 'D:/Researching/SMVT/03_Analysis/figures'
FINAL = 'D:/Researching/SMVT/04_Manuscript/figures'

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

fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 82 * MM))

# (a) mutation frequency
m = mut.sort_values('rate', ascending=False)
cols = ['#e74c3c' if r > 2 else '#3498db' for r in m['rate']]
axes[0].bar(m['cancer'], m['rate'], color=cols, edgecolor='white', width=0.6)
axes[0].axhline(2, color='#888', ls='--', lw=0.8, label='2% threshold')
for i, (c, r) in enumerate(zip(m['cancer'], m['rate'])):
    # stagger alternating labels: bars are ~4 mm apart at 183 mm, one label width
    axes[0].text(i, r + (0.08 if i % 2 else 0.55), f'{r:.1f}%', ha='center', fontsize=5)
axes[0].set_ylabel('Mutation frequency (% patients)')
axes[0].set_title('a  SLC5A6 somatic mutation rate (MC3)')
axes[0].legend(fontsize=8); axes[0].tick_params(axis='x', rotation=90, labelsize=8)
axes[0].grid(axis='y', alpha=0.15)

# (b) expression log2FC
de = de.sort_values('log2FC', ascending=False)
cols2 = ['#c62828' if x > 0 else '#1565c0' for x in de['log2FC']]
axes[1].bar(de['cancer'], de['log2FC'], color=cols2, edgecolor='white', width=0.6)
axes[1].axhline(0, color='#333', lw=0.8)
axes[1].set_ylabel('Log$_2$ fold change (tumor vs normal)')
axes[1].set_title('b  SLC5A6 expression (recomputed, BH-FDR)')
axes[1].tick_params(axis='x', rotation=90, labelsize=8)
axes[1].grid(axis='y', alpha=0.15)

plt.tight_layout()
import shutil
for path in save_pub(fig, STEM, outdir=OUT):
    shutil.copy(path, FINAL + '/' + os.path.basename(path))
print('Fig4 saved + synced')
