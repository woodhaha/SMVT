"""Fig1 multi-panel: bar + volcano + heatmap."""
import pandas as pd, numpy as np, json, csv, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

OUT = r'D:\Researching\SMVT\04_Manuscript\figures'
de = pd.read_csv(r'D:\Researching\SMVT\03_Analysis\outputs\pancan_de.csv')
de = de.sort_values('log2FC', ascending=True)
n = len(de)

colors = [
    '#c62828' if r['sig'] == '*' and r['log2FC'] > 0 else '#1565c0' if r['sig'] == '*' and r['log2FC'] < 0 else '#bdbdbd'
    for _, r in de.iterrows()
]

# Expression cache
cache = r'D:\Researching\SMVT\03_Analysis\data\slc5a6_by_cancer.json'
if os.path.exists(cache):
    by_cancer = json.load(open(cache))
else:
    CRC = r'D:\Researching\CRCproject\02_Data\raw'
    pheno = pd.read_csv(CRC + '/TcgaTargetGTEX_phenotype.txt', sep='\t', encoding='ISO-8859-1')
    s2c = dict(zip(pheno['sample'], pheno['primary disease or tissue']))
    t2a = {'Breast Invasive Carcinoma':'BRCA','Kidney Clear Cell Carcinoma':'KIRC','Kidney Papillary Cell Carcinoma':'KIRP',
        'Lung Adenocarcinoma':'LUAD','Thyroid Carcinoma':'THCA','Head & Neck Squamous Cell Carcinoma':'HNSC',
        'Prostate Adenocarcinoma':'PRAD','Lung Squamous Cell Carcinoma':'LUSC','Stomach Adenocarcinoma':'STAD',
        'Bladder Urothelial Carcinoma':'BLCA','Liver Hepatocellular Carcinoma':'LIHC',
        'Uterine Corpus Endometrioid Carcinoma':'UCEC','Colon Adenocarcinoma':'COAD','Sarcoma':'SARC',
        'Kidney Chromophobe':'KICH','Cervical Squamous Cell Carcinoma':'CESC','Pancreatic Adenocarcinoma':'PAAD',
        'Esophageal Carcinoma':'ESCA','Pheochromocytoma & Paraganglioma':'PCPG',
        'Rectal Adenocarcinoma':'READ','Uveal Melanoma':'UVM','Thymoma':'THYM',
        'Testicular Germ Cell Tumor':'TGCT','Cholangiocarcinoma':'CHOL','Uterine Carcinosarcoma':'UCS',
        'Mesothelioma':'MESO','Diffuse Large B-Cell Lymphoma':'DLBC','Adrenocortical Carcinoma':'ACC',
        'Glioblastoma Multiforme':'GBM','Rectum Adenocarcinoma':'READ','Cervical & Endocervical Cancer':'CESC',
        'Adrenocortical Cancer':'ACC','Ovarian Serous Cystadenocarcinoma':'OV',
        'Skin Cutaneous Melanoma':'SKCM','Acute Myeloid Leukemia':'LAML','Brain Lower Grade Glioma':'LGG'}
    print('Loading TOIL...')
    with open(CRC + '/TcgaTargetGtex_rsem_gene_fpkm', 'r', encoding='ISO-8859-1') as f:
        reader = csv.reader(f, delimiter='\t')
        hdr = next(reader)
        for row in reader:
            if row[0].startswith('ENSG00000138074'):
                vals = {hdr[i]: float(v) for i, v in enumerate(row[1:], 1) if v and float(v) > 0}
                break
    by_cancer = {}
    for sid, v in vals.items():
        if sid.startswith('TCGA') and sid.endswith(('-01','-03')):
            cn = s2c.get(sid, '')
            if not cn or cn != cn: continue
            ab = t2a.get(cn, '')
            if not ab: continue
            by_cancer.setdefault(ab, []).append(float(np.log2(v + 1)))
    json.dump(by_cancer, open(cache, 'w'))

# Figure
fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.7], height_ratios=[0.5, 0.5], hspace=0.28, wspace=0.18)

# (a) Bar
ax1 = fig.add_subplot(gs[0, 0])
ax1.barh(range(n), de['log2FC'].values, color=colors, ec='white', lw=0.4, height=0.6)
ax1.axvline(0, color='#333', lw=0.8)
ax1.axvline(0.2, color='#888', ls='--', lw=0.5, alpha=0.5)
yl = []
for _, r in de.iterrows():
    l = r['cancer'] + ' (' + str(r['n_t']) + '/' + str(r['n_n']) + ')'
    if r['sig'] == '*':
        if r['pval'] < 1e-30: l += ' ***'
        elif r['pval'] < 1e-10: l += ' **'
        else: l += ' *'
    yl.append(l)
ax1.set_yticks(range(n)); ax1.set_yticklabels(yl, fontsize=7.5)
ax1.set_xlabel('Log$_2$ Fold Change', fontsize=9)
ax1.set_title('a  Differential Expression', fontsize=10, fontweight='bold', loc='left')
ax1.set_xlim(de['log2FC'].min()-0.03, de['log2FC'].max()+max(0.12, abs(de['log2FC'].max())*0.12))
ax1.grid(axis='x', alpha=0.1); ax1.set_axisbelow(True)
ax1.tick_params(axis='y', left=False, labelsize=7.5)
from matplotlib.patches import Patch
ax1.legend(handles=[Patch(color='#c62828',label='Up'),Patch(color='#1565c0',label='Down'),Patch(color='#bdbdbd',label='NS')], fontsize=6.5, loc='lower right')

# (b) Volcano
ax2 = fig.add_subplot(gs[0, 1])
for _, r in de.iterrows():
    ax2.scatter(r['log2FC'], -np.log10(r['pval']), s=15 if r['sig']=='*' else 8,
               c='#c62828' if r['sig']=='*' and r['log2FC']>0 else '#1565c0' if r['sig']=='*' else '#bdbdbd', alpha=0.8)
for _, r in de.iterrows():
    if r['pval'] < 1e-15:
        ax2.annotate(r['cancer'], (r['log2FC'], -np.log10(r['pval'])), fontsize=6, ha='center', va='bottom', xytext=(0,3), textcoords='offset points')
ax2.axhline(-np.log10(0.05), color='#888', ls='--', lw=0.5, alpha=0.5)
ax2.set_xlabel('Log$_2$ FC', fontsize=9); ax2.set_ylabel('$-$Log$_{10}$(P)', fontsize=9)
ax2.set_title('b  Volcano Plot', fontsize=10, fontweight='bold', loc='left')
ax2.grid(alpha=0.1); ax2.axvline(0, color='#ccc', lw=0.3)

# (c) Heatmap
ax3 = fig.add_subplot(gs[1, :])
de_h = de.sort_values('log2FC', ascending=False)
cancers = list(de_h['cancer'])
expr_mat = []
for c in cancers:
    v = by_cancer.get(c,[])
    s = np.random.choice(v, min(80,len(v)), replace=False) if v else [np.nan]
    expr_mat.append(s)
mx = max(len(v) for v in expr_mat)
p = np.full((len(cancers),mx), np.nan)
for i,v in enumerate(expr_mat): p[i,:len(v)]=v
p = p[:, np.argsort(np.nanmedian(p, axis=0))]

ax3.imshow(p, aspect='auto', cmap='YlOrRd', extent=[0,mx,-0.5,len(cancers)-0.5], vmin=1.5, vmax=3.0)
for i,c in enumerate(cancers):
    m = float(np.nanmean(by_cancer.get(c,[np.nan])))
    s = de_h[de_h['cancer']==c]['sig'].values[0]
    ax3.plot(mx*1.04, i, 's', color='#c62828' if s=='*' else '#bdbdbd', markersize=5, mec='w', mew=0.3)
ax3.set_yticks(range(len(cancers))); ax3.set_yticklabels(cancers, fontsize=7.5)
ax3.set_xticks([]); ax3.set_xlim(0, mx*1.12)
ax3.set_xlabel('Tumor samples', fontsize=9)
ax3.set_title('c  SLC5A6 Tumor Expression', fontsize=10, fontweight='bold', loc='left')
ax3.tick_params(axis='y', left=False, labelsize=7.5)
plt.colorbar(ax3.images[0], ax=ax3, shrink=0.5, pad=0.03).set_label('log$_2$(FPKM+1)', fontsize=8)
from matplotlib.lines import Line2D
ax3.legend(handles=[Line2D([0],[0],marker='s',color='w',mfc='#c62828',ms=5,label='UP'),Line2D([0],[0],marker='s',color='w',mfc='#bdbdbd',ms=5,label='NS')], fontsize=6, loc='upper right')

plt.tight_layout()
fig.savefig(OUT+'/Fig1_SMVT_TCGA_pan_cancer.pdf', dpi=300, bbox_inches='tight')
fig.savefig(OUT+'/Fig1_SMVT_TCGA_pan_cancer.png', dpi=300, bbox_inches='tight')
plt.close()
print('Fig1 done')
