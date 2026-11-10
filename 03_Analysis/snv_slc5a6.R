"""SLC5A6 pan-cancer mutation analysis — MC3 MAF."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN = 'D:/Researching/SMVT/02_Data/raw/mc3_maf.gz'
OUT = 'D:/Researching/SMVT/03_Analysis'
FIG = OUT + '/figures'

# Mapping: TCGA barcode TSS codes -> cancer type abbreviations
# Merged from snv_slc5a6.R + fig_mutation_final.py for max coverage
code_map = {
    # BRCA
    'BR':'BRCA','AP':'BRCA','BH':'BRCA','A7':'BRCA','AQ':'BRCA','AR':'BRCA',
    '61':'BRCA','PM':'BRCA','P0':'BRCA','V7':'BRCA',
    # LUAD
    'LU':'LUAD','55':'LUAD','44':'LUAD','49':'LUAD','50':'LUAD','56':'LUAD',
    '78':'LUAD','85':'LUAD','95':'LUAD','67':'LUAD','93':'LUAD',
    '12':'LUAD','17':'LUAD','27':'LUAD','29':'LUAD','33':'LUAD','38':'LUAD',
    'N5':'LUAD','N8':'LUAD','B2':'LUAD','CE':'LUAD',
    # LUSC
    '43':'LUSC','66':'LUSC','77':'LUSC','18':'LUSC','21':'LUSC','22':'LUSC',
    '34':'LUSC','37':'LUSC','39':'LUSC','63':'LUSC','69':'LUSC',
    'N7':'LUSC','PK':'LUSC','VS':'LUSC','WE':'LUSC',
    '23':'LUSC',
    # BLCA
    'BL':'BLCA','CQ':'BLCA','CF':'BLCA','BT':'BLCA',
    'AB':'BLCA','2W':'BLCA','G3':'BLCA','GD':'BLCA','EM':'BLCA',
    'QC':'BLCA','94':'BLCA',
    # COAD / COADREAD
    'AA':'COAD','D5':'COAD','AD':'COAD','DK':'COADREAD',
    # ESCA
    'EE':'ESCA','LN':'ESCA','IG':'ESCA','VL':'ESCA','VR':'ESCA',
    # HNSC
    'HU':'HNSC','CN':'HNSC','CX':'HNSC','H6':'HNSC','BQ':'HNSC','CR':'HNSC',
    # SARC
    'CJ':'SARC','DX':'SARC','DB':'SARC','HD':'SARC','MD':'SARC',
    'RU':'SARC','S8':'SARC','TF':'SARC','86':'SARC',
    # LGG
    'Z6':'LGG','W7':'LGG','S9':'LGG','KM':'LGG','CS':'LGG',
    # STAD
    'AY':'STAD','CD':'STAD','CG':'STAD','CK':'STAD','CM':'STAD',
    'GN':'STAD','K7':'STAD','Q1':'STAD','QC':'STAD','R1':'STAD','SS':'STAD',
    '3B':'STAD','5K':'STAD','3A':'STAD','4Z':'STAD','3C':'STAD','3D':'STAD',
    # LIHC
    'AX':'LIHC','BC':'LIHC','CC':'LIHC','DD':'LIHC','EP':'LIHC','ES':'LIHC',
    'BM':'LIHC','AL':'LIHC','ED':'LIHC','P3':'LIHC','RC':'LIHC','WQ':'LIHC',
    'R9':'LIHC','SA':'LIHC','T1':'LIHC',
    # PAAD
    'A2':'PAAD','FZ':'PAAD','HZ':'PAAD','IB':'PAAD','Q3':'PAAD',
    'B8':'PAAD','F2':'PAAD',
    # GBM
    'EB':'GBM','EC':'GBM','E2':'GBM','E9':'GBM','19':'GBM','28':'GBM','32':'GBM','74':'GBM',
    # PRAD
    'B1':'PRAD','CH':'PRAD','E1':'PRAD','G9':'PRAD','HR':'PRAD',
    'C8':'PRAD','EJ':'PRAD','GH':'PRAD','UZ':'PRAD','VJ':'PRAD','YA':'PRAD','Z3':'PRAD',
    # KIRC
    'KC':'KIRC','DI':'KIRC','DV':'KIRC','CW':'KIRC','C3':'KIRC','F4':'KIRC','SX':'KIRC',
    # KIRP
    'C4':'KIRP','DZ':'KIRP','G2':'KIRP','AK':'KIRP','IR':'KIRP','KQ':'KIRP',
    # UCEC
    'A3':'UCEC','A5':'UCEC','B5':'UCEC','D1':'UCEC','E6':'UCEC','F1':'UCEC',
    'CU':'UCEC','DA':'UCEC','GM':'UCEC','J1':'UCEC','M8':'UCEC','MA':'UCEC',
    # CESC
    '24':'CESC','EA':'CESC','DS':'CESC','FI':'CESC','L5':'CESC','R5':'CESC',
    # OV
    'A6':'OV','EW':'OV','QD':'OV','OE':'OV','RP':'OV','S0':'OV',
    # SKCM
    'D8':'SKCM','ER':'SKCM','F5':'SKCM','G4':'SKCM',
    # Other
    'M5':'LAML',
    # Rare / unmapped → resolved via GDC API
    'D3':'SKCM','FR':'SKCM','FS':'SKCM','W3':'SKCM',
    'A8':'BRCA','AN':'BRCA',
    'DF':'UCEC','EO':'UCEC','AJ':'UCEC',
    'A4':'KIRP','DJ':'THCA','EI':'READ','GU':'BLCA',
    'EK':'CESC','LP':'CESC','ND':'UCS','US':'PAAD','VQ':'STAD',
}

# Load MC3 MAF
print('Loading MC3 MAF...')
import csv, gzip
mutations = []
with gzip.open(IN, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    h = next(reader)
    for row in reader:
        if row[0] == 'SLC5A6':
            mutations.append(dict(zip(h, row)))

df = pd.DataFrame(mutations)
print(f'Loaded {len(df)} SLC5A6 mutations')

# Map cancer types
def map_cancer(barcode):
    parts = barcode.split('-')
    if len(parts) >= 2:
        code = parts[0] + '-' + parts[1]
        if code in code_map: return code_map[code]
        if parts[1] in code_map: return code_map[parts[1]]
    return 'OTHER'

df['Cancer'] = df['Tumor_Sample_Barcode'].apply(map_cancer)
other = df[df['Cancer']=='OTHER']
if len(other) > 0:
    for _, r in other.iterrows():
        b = r['Tumor_Sample_Barcode']
        print(f'  UNMAPPED: {b}')

print(f'\n=== Per Cancer Type ===')
ct = df.groupby('Cancer').size().sort_values(ascending=False)
print(ct.to_string())

# Mutation frequency per cancer type (approximate sample sizes)
n_samples = {'BRCA':1092,'LUAD':585,'LUSC':504,'BLCA':412,'COAD':458,'LIHC':373,
             'STAD':443,'HNSC':528,'KIRC':537,'ESCA':185,'SARC':261,'LGG':530,
             'PAAD':185,'CESC':307,'UCEC':201,'OV':419,'KIRP':291,'SKCM':102}
freq = {}
for k, v in ct.items():
    n = n_samples.get(k, 300)
    freq[k] = v / n * 100

print(f'\n=== Mutation Frequency ===')
for k, v in sorted(freq.items(), key=lambda x: -x[1]):
    print(f'{k:6s} {ct[k]:4d} / {n_samples.get(k,300):4d} = {v:.2f}%')

print(f'\n=== By Type ===')
print(df['Variant_Classification'].value_counts().to_string())

print(f'\n=== Unique patients: {df["Tumor_Sample_Barcode"].nunique()}')

# Save
df.to_csv(f'{OUT}/outputs/slc5a6_mutations.csv', index=False)
print(f'\nSaved to outputs/slc5a6_mutations.csv')
