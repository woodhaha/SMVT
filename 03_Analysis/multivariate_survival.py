"""Corrected survival analysis for the 6 significant SMVT cancer types.

Diagnoses a bug in run_pancancer.py: the manuscript Table 2 HR used CONTINUOUS
SLC5A6 expression (per log2-unit), not the median-split binary claimed in Methods.
Computes three HRs per type:
  HR_cont  — continuous SLC5A6 (reproduces manuscript Table 2; inflated)
  HR_med   — median-split binary (correct)
  HR_adj   — median-split + AGE + GENDER + STAGE (multivariable)
"""
import json, sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

DATA = Path(r'D:\Researching\SMVT\03_Analysis\data')
PHENO = r'D:\Researching\CRC_Data_Mining\CRC_GEGES\TcgaTargetGTEX_phenotype.txt'
SQLITE = r'C:\Users\woodh\Documents\databases\TCGA.survival.sqlite'
GDC_CACHE = Path(r'D:\Researching\SMVT\03_Analysis\outputs\gdc_cache')

TYPES = {
    'KIRC': ('Kidney Clear Cell Carcinoma', 'TCGA-KIRC'),
    'LGG':  ('Brain Lower Grade Glioma', 'TCGA-LGG'),
    'ACC':  ('Adrenocortical Cancer', 'TCGA-ACC'),
    'UVM':  ('Uveal Melanoma', 'TCGA-UVM'),
    'MESO': ('Mesothelioma', 'TCGA-MESO'),
    'SARC': ('Sarcoma', 'TCGA-SARC'),
}

def barcode12(s):
    return s[:12] if s else ''

# 1. TOIL SLC5A6 -> {sample: log2(rsem+1)}
lines = open(DATA / 'slc5a6_toil_row.tsv', encoding='utf-8').read().splitlines()
hdr = lines[0].split('\t')[1:]
vals = lines[1].split('\t')[1:]
expr = {}
for s, v in zip(hdr, vals):
    try:
        x = float(v)
        if x >= 0:
            expr[s] = np.log2(x + 1)
    except ValueError:
        pass

# 2. Phenotype
pheno = pd.read_csv(PHENO, sep='\t',
                    usecols=['sample', 'primary disease or tissue', '_sample_type', '_gender'],
                    encoding='ISO-8859-1')
pheno['sample'] = pheno['sample'].astype(str)

# 3. Survival sqlite
conn = sqlite3.connect(SQLITE)

# 4. GDC clinical -> AGE + AJCC STAGE
def gdc_clinical(project):
    j = json.load(open(GDC_CACHE / f'{project}_clinical.json', encoding='utf-8'))
    recs = []
    for case in j.get('data', {}).get('hits', []):
        demo = case.get('demographic', {}) or {}
        diag = (case.get('diagnoses') or [{}])[0]
        rec = {'patient_id': case.get('submitter_id', ''), 'AGE': None, 'STAGE': None}
        dtb = demo.get('days_to_birth')
        if dtb is not None:
            rec['AGE'] = abs(float(dtb)) / 365.25
        st = diag.get('ajcc_pathologic_stage')
        if isinstance(st, str):
            for r in [('IV', 4), ('III', 3), ('II', 2), ('I', 1)]:
                if r[0] in st:
                    rec['STAGE'] = r[1]; break
        recs.append(rec)
    return pd.DataFrame(recs)

def fit_hr(df, duration, event):
    c = CoxPHFitter().fit(df, duration_col=duration, event_col=event)
    return float(np.exp(c.params_['expr'])), float(c.summary['p']['expr'])

print(f"{'CT':6s} {'N':>4s} {'HR_cont':>8s} {'HR_med':>8s} {'HR_adj':>8s} {'P_adj':>9s} {'N_adj':>5s} {'cov':>15s}")
print('-' * 82)

rows_out = []
for abbrev, (disease, project) in TYPES.items():
    tum = pheno[(pheno['primary disease or tissue'] == disease) &
                (pheno['_sample_type'].isin(['Primary Tumor', 'Primary Solid Tumor']))]
    gmap = dict(zip(tum['sample'], tum['_gender']))
    df = pd.DataFrame({'sample': tum['sample'].tolist()})
    df['pat'] = df['sample'].map(barcode12)
    df['expr'] = df['sample'].map(expr)
    df['gender'] = df['sample'].map(gmap)
    df = df.dropna(subset=['expr'])

    surv = pd.read_sql_query(f'SELECT row_names AS sample, OS, OS_status FROM "survival.{abbrev}"', conn)
    surv['pat'] = surv['sample'].map(barcode12)
    surv = surv.drop_duplicates(subset='pat').set_index('pat')[['OS', 'OS_status']]

    clin = gdc_clinical(project).set_index('patient_id')[['AGE', 'STAGE']]

    m = df.join(surv, on='pat').join(clin, on='pat')
    m = m.dropna(subset=['OS', 'OS_status']).copy()
    m['OS'] = m['OS'].astype(float)
    m['OS_status'] = m['OS_status'].astype(int)
    m['expr_hi'] = (m['expr'] >= m['expr'].median()).astype(int)
    m['gender'] = m['gender'].map(lambda g: 1 if str(g).upper().startswith('M') else 0)
    n = len(m)

    # continuous (reproduces manuscript bug)
    d_cont = m[['OS', 'OS_status', 'expr']].rename(columns={'OS': 'duration', 'OS_status': 'event'})
    hr_cont, _ = fit_hr(d_cont, 'duration', 'event')

    # median-split binary (correct)
    d_med = m[['OS', 'OS_status', 'expr_hi']].rename(
        columns={'OS': 'duration', 'OS_status': 'event', 'expr_hi': 'expr'})
    hr_med, _ = fit_hr(d_med, 'duration', 'event')

    # multivariable: median-split + AGE + GENDER + STAGE(if avail)
    covs = ['expr', 'AGE', 'gender']
    has_stage = m['STAGE'].notna().sum() > n * 0.5
    if has_stage:
        covs.append('STAGE')
    d_adj = m[['OS', 'OS_status', 'expr_hi'] + covs[1:]].dropna().rename(
        columns={'OS': 'duration', 'OS_status': 'event', 'expr_hi': 'expr'})
    n_adj = len(d_adj)
    hr_adj = p_adj = None
    if n_adj >= 30:
        try:
            hr_adj, p_adj = fit_hr(d_adj, 'duration', 'event')
        except Exception as e:
            hr_adj, p_adj = f'ERR:{e}', None

    covlab = '+'.join(covs[1:])
    hr_adj_s = f"{hr_adj:.2f}" if isinstance(hr_adj, float) else str(hr_adj)
    p_adj_s = f"{p_adj:.2e}" if isinstance(p_adj, float) else str(p_adj)
    print(f"{abbrev:6s} {n:>4d} {hr_cont:>8.2f} {hr_med:>8.2f} {hr_adj_s:>8s} {p_adj_s:>9s} {n_adj:>5d} {covlab:>15s}")
    rows_out.append({'cancer': abbrev, 'n': n, 'HR_cont': hr_cont, 'HR_med': hr_med,
                     'HR_adj': hr_adj, 'P_adj': p_adj, 'N_adj': n_adj, 'cov': covlab})

pd.DataFrame(rows_out).to_csv(DATA / 'multivariate_survival_results.csv', index=False)
conn.close()
print(f"\nSaved: {DATA / 'multivariate_survival_results.csv'}")
