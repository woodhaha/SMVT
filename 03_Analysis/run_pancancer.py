"""Pan-cancer SLC5A6: DE + survival analysis using local TCGA data."""
import json, os, csv, sqlite3, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu
warnings.filterwarnings('ignore')

CRC = r'D:\Researching\CRCproject\02_Data\raw'
OUT = r'D:\Researching\SMVT\03_Analysis'
os.makedirs(f'{OUT}/data', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)

# Abbreviation mapping
t2a = {
    'Breast Invasive Carcinoma': 'BRCA', 'Kidney Clear Cell Carcinoma': 'KIRC',
    'Kidney Papillary Cell Carcinoma': 'KIRP', 'Lung Adenocarcinoma': 'LUAD',
    'Thyroid Carcinoma': 'THCA', 'Head & Neck Squamous Cell Carcinoma': 'HNSC',
    'Prostate Adenocarcinoma': 'PRAD', 'Lung Squamous Cell Carcinoma': 'LUSC',
    'Brain Lower Grade Glioma': 'LGG', 'Skin Cutaneous Melanoma': 'SKCM',
    'Stomach Adenocarcinoma': 'STAD', 'Ovarian Serous Cystadenocarcinoma': 'OV',
    'Bladder Urothelial Carcinoma': 'BLCA', 'Liver Hepatocellular Carcinoma': 'LIHC',
    'Acute Myeloid Leukemia': 'LAML', 'Uterine Corpus Endometrioid Carcinoma': 'UCEC',
    'Colon Adenocarcinoma': 'COAD', 'Sarcoma': 'SARC', 'Kidney Chromophobe': 'KICH',
    'Cervical Squamous Cell Carcinoma': 'CESC', 'Pancreatic Adenocarcinoma': 'PAAD',
    'Esophageal Carcinoma': 'ESCA', 'Pheochromocytoma & Paraganglioma': 'PCPG',
    'Rectal Adenocarcinoma': 'READ', 'Uveal Melanoma': 'UVM', 'Thymoma': 'THYM',
    'Testicular Germ Cell Tumor': 'TGCT', 'Cholangiocarcinoma': 'CHOL',
    'Uterine Carcinosarcoma': 'UCS', 'Mesothelioma': 'MESO',
    'Diffuse Large B-Cell Lymphoma': 'DLBC', 'Adrenocortical Carcinoma': 'ACC',
    'Glioblastoma Multiforme': 'GBM',
    'Rectum Adenocarcinoma': 'READ', 'Cervical & Endocervical Cancer': 'CESC',
    'Adrenocortical Cancer': 'ACC',
}

# 1. SLC5A6 expression
print('Loading SLC5A6...')
with open(f'{CRC}/TcgaTargetGtex_rsem_gene_fpkm', 'r', encoding='ISO-8859-1') as f:
    reader = csv.reader(f, delimiter='\t')
    hdr = next(reader)
    for row in reader:
        if row[0].startswith('ENSG00000138074'):
            slc = {}
            for i, v in enumerate(row[1:], 1):
                try:
                    x = float(v)
                    if x > 0: slc[hdr[i]] = np.log2(x + 1)
                except: pass
            break

# 2. Phenotype mapping
pheno = pd.read_csv(f'{CRC}/TcgaTargetGTEX_phenotype.txt', sep='\t', encoding='ISO-8859-1')
pheno['abbrev'] = pheno['primary disease or tissue'].map(t2a)
s2c = dict(zip(pheno['sample'], pheno['abbrev']))
print(f'Expression samples: {len(slc)}')

# 3. Build per-cancer TCGA datasets
cancer_expr = defaultdict(lambda: {'tumor': [], 'normal': []})
tcga_by_patient = {}  # patient_id -> {cancer, expr}

for sid, v in slc.items():
    if sid.startswith('TCGA'):
        cancer = s2c.get(sid)
        if not cancer: continue
        patient = sid[:12]
        expr = float(v)
        if sid.endswith(('-01','-03')):  # tumor
            cancer_expr[cancer]['tumor'].append(expr)
            tcga_by_patient.setdefault(patient, {})[sid] = {'cancer': cancer, 'expr': expr}
        elif sid.endswith(('-11','-10')):  # normal
            cancer_expr[cancer]['normal'].append(expr)

print(f'Tumor samples across {len(cancer_expr)} cancer types')

# 4. Differential expression
print('\n=== DIFFERENTIAL EXPRESSION ===')
de = []
for cancer, d in sorted(cancer_expr.items()):
    tv, nv = d['tumor'], d['normal']
    if len(tv) >= 10 and len(nv) >= 3:
        tm, nm = float(np.mean(tv)), float(np.mean(nv))
        l2fc = round(tm - nm, 3)
        _, p = mannwhitneyu(tv, nv, alternative='two-sided')
        sigma = '*' if p < 0.05 else ' '
        de.append({'cancer':cancer, 'n_t':len(tv), 'n_n':len(nv),
                   't_mean':tm, 'n_mean':nm, 'log2FC':l2fc, 'pval':p, 'sig':sigma.strip()})
        print(f'{cancer:6s} n_t={len(tv):4d} n_n={len(nv):3d} '
              f't_mean={tm:.3f} n_mean={nm:.3f} l2FC={l2fc:+7.3f} P={p:.2e} {sigma}')

pd.DataFrame(de).to_csv(f'{OUT}/outputs/pancan_de.csv', index=False)
sig_up = [r['cancer'] for r in de if r['log2FC'] > 0 and r['pval'] < 0.05]
sig_down = [r['cancer'] for r in de if r['log2FC'] < 0 and r['pval'] < 0.05]
print(f'\nSignificant UP: {sig_up}')
print(f'Significant DOWN: {sig_down}')

# 5. Per-cancer survival analysis
print('\n=== SURVIVAL ANALYSIS ===')
# Build patient → cancer + expr mapping (one entry per patient)
pat_to_expr = {}
for sid, v in slc.items():
    if sid.startswith('TCGA') and sid.endswith(('-01','-03')):
        cancer = s2c.get(sid)
        if cancer:
            pid = sid[:12]
            pat_to_expr[pid] = {'cancer': cancer, 'expr': float(v)}

conn = sqlite3.connect(f'{CRC}/TCGA.survival.sqlite')
tables = pd.read_sql('SELECT name FROM sqlite_master WHERE type="table"', conn)
surv_out = []

for tbl in tables['name']:
    abbr = tbl.split('.')[1]
    s = pd.read_sql(f'SELECT * FROM "{tbl}"', conn)
    s.columns = ['pid','OS','ev','rfs_d','rfs_e']

    # Match patients by cancer type via their expression records
    matched_pats = {pid: d for pid, d in pat_to_expr.items() if d['cancer'] == abbr}
    pids = set(matched_pats.keys()) & set(s['pid'].str[:12])
    if len(pids) < 30: continue

    # Build analysis dataframe
    rows = []
    for pid in pids:
        sr = s[s['pid'].str[:12] == pid].iloc[0]
        rows.append({'pid': pid, 'OS': sr['OS'], 'ev': sr['ev'],
                     'SLC': matched_pats[pid]['expr']})
    m = pd.DataFrame(rows).dropna(subset=['OS','ev','SLC'])
    n = len(m)
    if n < 30: continue

    m = m.copy()
    med = m['SLC'].median()
    m['grp'] = np.where(m['SLC'] >= med, 'High', 'Low')
    h, l = m[m.grp=='High'], m[m.grp=='Low']
    if len(h) < 10 or len(l) < 10: continue

    kh = KaplanMeierFitter().fit(h['OS'], h['ev'])
    kl = KaplanMeierFitter().fit(l['OS'], l['ev'])
    lr = logrank_test(h['OS'], l['OS'], event_observed_A=h['ev'], event_observed_B=l['ev'])

    cph = CoxPHFitter()
    cx = m[['OS','ev','SLC']].rename(columns={'OS':'t','ev':'e'})
    try:
        cph.fit(cx, duration_col='t', event_col='e')
        hr = float(np.exp(cph.params_.iloc[0]))
        hl = float(np.exp(cph.confidence_intervals_.iloc[0,0]))
        hu = float(np.exp(cph.confidence_intervals_.iloc[0,1]))
        cp = float(cph.summary.p.iloc[0])
    except:
        hr, hl, hu, cp = float('nan'), float('nan'), float('nan'), float('nan')

    sig = lr.p_value < 0.05
    surv_out.append({'cancer':abbr, 'n':n, 'ev':int(m['ev'].sum()),
                     'med_hi':kh.median_survival_time_, 'med_lo':kl.median_survival_time_,
                     'logrank':lr.p_value, 'hr':hr, 'hr_lo':hl, 'hr_hi':hu, 'cox_p':cp})
    mark = ' *' if sig else ''
    print(f'{abbr:6s} n={n:4d} ev={int(m["ev"].sum()):3d} '
          f'med_hi={kh.median_survival_time_:5.1f} med_lo={kl.median_survival_time_:5.1f} '
          f'P={lr.p_value:.4f} HR={hr:.4f}{mark}')

    if sig:
        fig, ax = plt.subplots(figsize=(7,5))
        kh.plot(ax=ax, color='#d62728', lw=1.5)
        kl.plot(ax=ax, color='#1f77b4', lw=1.5)
        ax.set_xlabel('Months'); ax.set_ylabel('OS')
        ax.set_title(f'{abbr}: SLC5A6 (n={n})\nHR={hr:.3f}, P={lr.p_value:.4f}')
        ax.legend(); fig.savefig(f'{OUT}/figures/KM_{abbr}_SMVT_survival.pdf', dpi=300, bbox_inches='tight')
        plt.close()

pd.DataFrame(surv_out).to_csv(f'{OUT}/outputs/pancan_survival.csv', index=False)
sig_cts = [r['cancer'] for r in surv_out if r['logrank'] < 0.05]
print(f'\nTotal: {len(surv_out)} cancers analyzed')
print(f'Significant (P<0.05): {sig_cts if sig_cts else "NONE"}')
