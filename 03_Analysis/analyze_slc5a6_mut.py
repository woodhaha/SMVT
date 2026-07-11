"""SLC5A6 MC3 MAF mutation analysis."""
import pandas as pd, numpy as np, csv, gzip

MAF = 'D:/Researching/SMVT/02_Data/raw/mc3_maf.gz'
OUT = 'D:/Researching/SMVT/03_Analysis'

code_map = {'BR':'BRCA','AP':'BRCA','BH':'BRCA','A7':'BRCA','AQ':'BRCA','AR':'BRCA','PM':'BRCA','P0':'BRCA',
    'AA':'COAD','D5':'COAD','AD':'COAD','DK':'COAD',
    '55':'LUAD','44':'LUAD','49':'LUAD','50':'LUAD','56':'LUAD','78':'LUAD','85':'LUAD','95':'LUAD','12':'LUAD','17':'LUAD','27':'LUAD','29':'LUAD','33':'LUAD','38':'LUAD','67':'LUAD','93':'LUAD',
    '43':'LUSC','66':'LUSC','77':'LUSC','18':'LUSC','21':'LUSC','22':'LUSC','34':'LUSC','37':'LUSC','39':'LUSC','63':'LUSC','69':'LUSC','N7':'LUSC','PK':'LUSC','VS':'LUSC','WE':'LUSC',
    'BL':'BLCA','CQ':'BLCA','CF':'BLCA','BT':'BLCA','AB':'BLCA','2W':'BLCA','G3':'BLCA','EM':'BLCA','GD':'BLCA','DK':'COAD',
    'EE':'ESCA','LN':'ESCA','IG':'ESCA','VL':'ESCA','VR':'ESCA',
    'HU':'HNSC','CN':'HNSC','CX':'HNSC','H6':'HNSC','BQ':'HNSC','CR':'HNSC',
    'CJ':'SARC','DX':'SARC','DB':'SARC','HD':'SARC','MD':'SARC','TF':'SARC','RU':'SARC','S8':'SARC',
    'Z6':'LGG','W7':'LGG','S9':'LGG','CS':'LGG',
    'AY':'STAD','CD':'STAD','CG':'STAD','CK':'STAD','CM':'STAD','GN':'STAD','K7':'STAD','QC':'STAD','R1':'STAD','SS':'STAD',
    'AX':'LIHC','BC':'LIHC','CC':'LIHC','DD':'LIHC','EP':'LIHC','ES':'LIHC','BM':'LIHC','AL':'LIHC','ED':'LIHC','P3':'LIHC','RC':'LIHC','WQ':'LIHC',
    'A2':'PAAD','FZ':'PAAD','HZ':'PAAD','IB':'PAAD','Q3':'PAAD','B8':'PAAD','F2':'PAAD',
    'EB':'GBM','EC':'GBM','19':'GBM','28':'GBM','32':'GBM','74':'GBM','E2':'GBM','E9':'GBM',
    'B1':'PRAD','CH':'PRAD','E1':'PRAD','G9':'PRAD','HR':'PRAD','C8':'PRAD','EJ':'PRAD','GH':'PRAD','YA':'PRAD','Z3':'PRAD',
    'A3':'UCEC','A5':'UCEC','B5':'UCEC','D1':'UCEC','E6':'UCEC','F1':'UCEC','CU':'UCEC','DA':'UCEC',
    'DI':'KIRC','DV':'KIRC','CW':'KIRC','C3':'KIRC','KC':'KIRC','F4':'KIRC','SX':'KIRC',
    'C4':'KIRP','DZ':'KIRP','G2':'KIRP','AK':'KIRP','IR':'KIRP','KQ':'KIRP',
    '24':'CESC','EA':'CESC','DS':'CESC','FI':'CESC','L5':'CESC','R5':'CESC',
    'D8':'SKCM','ER':'SKCM','F5':'SKCM','G4':'SKCM',
    'EW':'OV','QD':'OV','OE':'OV','RP':'OV','S0':'OV',
    'M5':'LAML',
    '61':'BRCA',
}

samples = {'BRCA':1092,'LUAD':585,'LUSC':504,'BLCA':412,'COAD':458,'LIHC':373,'STAD':443,
           'HNSC':528,'KIRC':537,'ESCA':185,'SARC':261,'LGG':530,'PAAD':185,'CESC':307,
           'UCEC':201,'OV':419,'KIRP':291,'SKCM':102,'GBM':153,'PRAD':500,'LAML':173}

def map_cancer(b):
    p = b.split('-')
    return code_map.get(p[1], 'OTHER')

print('Scanning MC3 MAF...')
rows = []
with gzip.open(MAF, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    h = next(reader)
    for row in reader:
        if row[0] == 'SLC5A6':
            rows.append(dict(zip(h, row)))

df = pd.DataFrame(rows)
print(f'Total mutations: {len(df)}')
df['Cancer'] = df['Tumor_Sample_Barcode'].apply(map_cancer)
unm = df[df['Cancer']=='OTHER']
if len(unm) > 0:
    print(f'\nUnmapped ({len(unm)}):')
    for _, r in unm.iterrows():
        pass  # skip verbose

ct = df[df['Cancer']!='OTHER'].groupby('Cancer').size().sort_values(ascending=False)
print(f'\n{"Cancer":>6s}  N_mut  N_samples  Freq')
print('-'*35)
for k, v in ct.items():
    n = samples.get(k, 100)
    print(f'{k:>6s}  {v:5d}  {n:9d}  {v/n*100:5.2f}%')

print(f'\nMutation types:')
print(df['Variant_Classification'].value_counts().to_string())
ns = df[df['Variant_Classification'].isin(['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','In_Frame_Del'])]
print(f'\nNon-silent: {len(ns)} / {len(df)}')
print(f'Unique patients: {df["Tumor_Sample_Barcode"].nunique()}')

# Recurrent
hot = df.groupby('HGVSp_Short').size().sort_values(ascending=False)
rec = hot[hot > 1]
if len(rec) > 0:
    print(f'\nRecurrent (>1 patient):')
    for k, v in rec.items():
        if k and k != 'nan':
            print(f'  {k:15s} x{v}')

df.to_csv(OUT+'/outputs/slc5a6_mutations.csv', index=False)
print(f'\nSaved to outputs/slc5a6_mutations.csv')
