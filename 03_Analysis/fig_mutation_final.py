"""Fig3: SLC5A6 mutation landscape from MC3 MAF — final version."""
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

OUT = r'D:\Researching\SMVT\03_Analysis\figures'
FINAL = r'D:\Researching\SMVT\04_Manuscript\figures'
MAF = r'D:\Researching\SMVT\02_Data\cleaned\slc5a6_mc3_maf.csv'

# Load from fresh for figure data
code_map = {'BR':'BRCA','AP':'BRCA','BH':'BRCA','AA':'COAD','D5':'COAD','AD':'COAD',
    '55':'LUAD','44':'LUAD','49':'LUAD','50':'LUAD','56':'LUAD','78':'LUAD','85':'LUAD','95':'LUAD',
    '43':'LUSC','66':'LUSC','77':'LUSC','18':'LUSC','21':'LUSC','22':'LUSC','34':'LUSC',
    '63':'LUSC','69':'LUSC','S9':'LGG','A5':'UCEC','B5':'UCEC','D1':'UCEC','E6':'UCEC','F1':'UCEC',
    'BL':'BLCA','EE':'ESCA','LN':'ESCA','HU':'HNSC','BQ':'HNSC','AX':'LIHC','CC':'LIHC',
    'AY':'STAD','CG':'STAD','CK':'STAD','CM':'STAD','Q1':'STAD',
    'DI':'KIRC','CJ':'SARC','DX':'SARC','Z6':'LGG','EB':'GBM','EC':'GBM',
    '24':'CESC','D8':'SKCM','ER':'SKCM','F5':'SKCM','G4':'SKCM',
    'EW':'OV','A3':'UCEC','B1':'PRAD','CH':'PRAD','HR':'PRAD','A2':'PAAD',
    'C4':'KIRP','IR':'KIRP','M5':'LAML','PK':'LUSC','VS':'LUSC','WE':'LUSC',
    '61':'BRCA','DK':'COAD','IG':'ESCA','VL':'ESCA','VR':'ESCA','2W':'BLCA','G3':'BLCA',
    'GD':'BLCA','S8':'SARC','RU':'SARC','TF':'SARC','CS':'LGG','GN':'STAD','K7':'STAD',
    'QC':'STAD','R1':'STAD','SS':'STAD','P3':'LIHC','RC':'LIHC','WQ':'LIHC',
    'C8':'PRAD','EJ':'PRAD','GH':'PRAD','Z3':'PRAD',
    'CU':'UCEC','DA':'UCEC',
    'L5':'CESC','R5':'CESC',
    'RP':'OV','S0':'OV','QD':'OV','OE':'OV',
    'F4':'KIRC','SX':'KIRC','KQ':'KIRP','AK':'KIRP',
    'PM':'BRCA','P0':'BRCA',
    'DB':'SARC','HD':'SARC','MD':'SARC','F2':'PAAD','B8':'PAAD','B2':'LUAD','N8':'LUAD',
    'N7':'LUSC','QC':'BLCA','UZ':'PRAD','YA':'PRAD','VJ':'PRAD',
    'ED':'LIHC','AL':'LIHC','BM':'LIHC','EP':'LIHC','ES':'LIHC',
    'CE':'LUAD','N5':'LUAD','N8':'LUAD','38':'LUAD','67':'LUAD','93':'LUAD',
    '37':'LUSC','39':'LUSC',
    'V7':'BRCA','E1':'PRAD','G9':'PRAD',
    # Rare / resolved via GDC API
    'D3':'SKCM','FR':'SKCM','FS':'SKCM','W3':'SKCM',
    'A8':'BRCA','AN':'BRCA',
    'DF':'UCEC','EO':'UCEC','AJ':'UCEC',
    'A4':'KIRP','DJ':'THCA','EI':'READ','GU':'BLCA',
    'EK':'CESC','LP':'CESC','ND':'UCS','US':'PAAD','VQ':'STAD',
}

samples = {'BRCA':1092,'LUAD':585,'LUSC':504,'BLCA':412,'COAD':458,'LIHC':373,'STAD':443,
           'HNSC':528,'KIRC':537,'ESCA':185,'SARC':261,'LGG':530,'PAAD':185,'CESC':307,
           'UCEC':201,'OV':419,'KIRP':291,'SKCM':102,'GBM':153,'PRAD':500,'LAML':173}

def mc(b):
    p = b.split('-')
    return code_map.get(p[1], 'OTHER')

df = pd.read_csv(MAF, sep='\t', dtype=str)
df['Cancer'] = df['Tumor_Sample_Barcode'].apply(mc)
df = df[df['Cancer']!='OTHER']

# Per-cancer mutation rate
cancer_order = ['UCEC','ESCA','SKCM','LUSC','LUAD','BRCA','STAD','GBM','PAAD',
    'COAD','SARC','BLCA','HNSC','LIHC','CESC','PRAD','OV','KIRP','KIRC','LAML']
mut_rates = {}
for c in cancer_order:
    n = samples.get(c, 100)
    m = len(df[df['Cancer']==c])
    mut_rates[c] = m/n*100

# Classify mutation types
df['Class'] = 'Silent/Other'
df.loc[df['Variant_Classification'].isin(['Missense_Mutation','In_Frame_Del']), 'Class'] = 'Missense'
df.loc[df['Variant_Classification'].isin(['Nonsense_Mutation','Frame_Shift_Del']), 'Class'] = 'Truncating'
df.loc[df['Variant_Classification']=='Splice_Site', 'Class'] = 'Splice'

fig = plt.figure(figsize=(14, 8))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.6], height_ratios=[0.4, 0.6], hspace=0.3, wspace=0.25)

# (a) Mutation frequency bar chart
ax1 = fig.add_subplot(gs[0, :])
c = [c for c in cancer_order if c in mut_rates]
v = [mut_rates[cc] for cc in c]
bars = ax1.bar(range(len(c)), v, color='#e74c3c', edgecolor='white', linewidth=0.4, width=0.6)
ax1.set_xticks(range(len(c)))
ax1.set_xticklabels(c, fontsize=8, rotation=90)
ax1.set_ylabel('SLC5A6 Mutation Frequency (%)', fontsize=10)
ax1.set_title('a  SLC5A6 Somatic Mutation Rate (TCGA MC3, n=193)', fontsize=11, fontweight='bold', loc='left')
ax1.axhline(2, color='#888', ls='--', lw=0.6, alpha=0.5, label='2% (passenger threshold)')
ax1.legend(fontsize=7)
ax1.grid(axis='y', alpha=0.15)

# (b) Mutation type pie
ax2 = fig.add_subplot(gs[1, 0])
class_counts = df['Class'].value_counts()
colors_pie = ['#e74c3c', '#f39c12', '#3498db', '#95a5a6']
ax2.pie(class_counts.values, labels=[f'{k}\n({v})' for k,v in class_counts.items()],
        colors=colors_pie[:len(class_counts)], autopct='%1.1f%%', startangle=90,
        textprops={'fontsize':8})
ax2.set_title('b  Mutation Type Distribution', fontsize=10, fontweight='bold')

# (c) Lollipop-like plot
ax3 = fig.add_subplot(gs[1, 1])
protein_len = 635
missense = df[df['Class']=='Missense'][['HGVSp_Short','Variant_Classification']].dropna()
trunc = df[df['Class']=='Truncating'][['HGVSp_Short','Variant_Classification']].dropna()

def parse_pos(hgvsp):
    if not hgvsp or hgvsp == '.' or hgvsp != hgvsp: return None
    for p in hgvsp.split('/'):
        m = re.search(r'p\.([A-Z]\d+)', str(p))
        if m:
            try: return int(re.search(r'\d+', m.group(1)).group())
            except: return None
    return None

# Count at each position
pos_counts = {}
for _, r in df[df['Class'].isin(['Missense','Truncating'])].iterrows():
    pos = parse_pos(r.get('HGVSp_Short',''))
    if pos and 1 <= pos <= protein_len:
        typ = r['Class']
        key = (pos, typ)
        pos_counts[key] = pos_counts.get(key, 0) + 1

# Plot missense
ms_pos = [(p, c) for (p, t), c in pos_counts.items() if t == 'Missense']
tr_pos = [(p, c) for (p, t), c in pos_counts.items() if t == 'Truncating']

ax3.scatter([p for p,_ in ms_pos], [c for _,c in ms_pos],
           color='#f39c12', s=30, alpha=0.6, label='Missense')
ax3.scatter([p for p,_ in tr_pos], [c for _,c in tr_pos],
           color='#e74c3c', s=50, marker='v', alpha=0.8, label='Truncating')

# Add labels — cascade down per cluster with generous headroom
labels_list = []
for (p, t), c in pos_counts.items():
    if c >= 3:
        label = ''
        for hgvsp in df[df['Class']==t]['HGVSp_Short'].dropna():
            if parse_pos(hgvsp) == p:
                label = hgvsp.split('/')[0] if '/' in str(hgvsp) else hgvsp
                break
        labels_list.append((p, c, label))
labels_list.sort()
max_y = max(pos_counts.values()) if pos_counts else 5
extra = int(1.5 * len(labels_list)) + 3
ax3.set_ylim(-0.5, max_y + extra)
cluster = []
for i, (p, c, label) in enumerate(labels_list):
    cluster.append((p, c, label))
    if i == len(labels_list)-1 or labels_list[i+1][0] - p > 50:
        for j, (cp, cc, clabel) in enumerate(cluster):
            dy = 2 + j * 3  # data-units offset above the point
            ax3.annotate(clabel, (cp, cc + dy), fontsize=5, ha='center', va='bottom', alpha=0.7)
        cluster = []

for start, end, name in [(1,40,'N-term'), (470,635,'C-term')]:
    ax3.axvspan(start, end, alpha=0.06, color='#2ecc71')
ax3.text(20, 0.02, 'N', fontsize=8, color='#2ecc71')
ax3.text(610, 0.02, 'C', fontsize=8, color='#2ecc71')

ax3.set_xlabel('SLC5A6 Protein Position', fontsize=9)
ax3.set_ylabel('Mutation Count', fontsize=9)
ax3.set_title('c  Mutation Distribution', fontsize=10, fontweight='bold', loc='left')
ax3.set_xlim(0, protein_len+5)
ax3.legend(fontsize=7)
ax3.grid(alpha=0.15)

fig.subplots_adjust(left=0.06, right=0.98, bottom=0.06, top=0.95)
fig.savefig(OUT+'/Fig3_SMVT_mutation_landscape.pdf', dpi=300, bbox_inches='tight')
fig.savefig(OUT+'/Fig3_SMVT_mutation_landscape.png', dpi=300, bbox_inches='tight')
import shutil
shutil.copy(OUT+'/Fig3_SMVT_mutation_landscape.pdf', FINAL+'/Fig3_SMVT_mutation_landscape.pdf')
shutil.copy(OUT+'/Fig3_SMVT_mutation_landscape.png', FINAL+'/Fig3_SMVT_mutation_landscape.png')
plt.close()

print('Fig3 done — MC3 MAF based')
print(f'Total mutations: {len(df)}')
ns = sum(1 for _, r in df.iterrows() if r['Class'] != 'Silent/Other')
print(f'Non-silent: {ns}')
pts = df["Tumor_Sample_Barcode"].nunique()
print(f"Patients: {pts}")
