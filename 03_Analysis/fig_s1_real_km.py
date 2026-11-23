"""FigS1 — REAL Kaplan-Meier curves for the 6 FDR-significant cancer types.

Replaces the Weibull-approximation in fig_s1_s2_real.py with actual per-patient
KM curves (lifelines KaplanMeierFitter), log-rank P, and corrected median-split HR.
"""
import sqlite3, csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fig_style import apply_style, PALETTE, MM, DOUBLE, save_pub
apply_style()
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA = Path(r'D:\Researching\SMVT\03_Analysis\data')
PHENO = r'D:\Researching\CRC_Data_Mining\CRC_GEGES\TcgaTargetGTEX_phenotype.txt'
SQLITE = r'C:\Users\woodh\Documents\databases\TCGA.survival.sqlite'
OUT = Path(r'D:\Researching\SMVT\04_Manuscript\figures')

TYPES = {
    'KIRC': 'Kidney Clear Cell Carcinoma',
    'LGG':  'Brain Lower Grade Glioma',
    'ACC':  'Adrenocortical Cancer',
    'UVM':  'Uveal Melanoma',
    'MESO': 'Mesothelioma',
    'SARC': 'Sarcoma',
}

# corrected HR + q (from median-split analysis)
hr_q = {}
for r in csv.DictReader(open(DATA / '../outputs/pancan_survival_corrected.csv')):
    hr_q[r['cancer']] = (float(r['hr']), float(r['q']))

# expression
lines = open(DATA / 'slc5a6_toil_row.tsv', encoding='utf-8').read().splitlines()
hdr = lines[0].split('\t')[1:]; vals = lines[1].split('\t')[1:]
expr = {}
for s, v in zip(hdr, vals):
    try:
        x = float(v)
        if x >= 0: expr[s] = np.log2(x + 1)
    except ValueError:
        pass

pheno = pd.read_csv(PHENO, sep='\t', usecols=['sample', 'primary disease or tissue', '_sample_type'],
                    encoding='ISO-8859-1')
pheno['sample'] = pheno['sample'].astype(str)
conn = sqlite3.connect(SQLITE)

# order by significance (q ascending)
order = sorted(TYPES, key=lambda a: hr_q.get(a, (1, 1))[1])

fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 110 * MM))
for ax, abbrev in zip(axes.flat, order):
    disease = TYPES[abbrev]
    tum = pheno[(pheno['primary disease or tissue'] == disease) &
                (pheno['_sample_type'].isin(['Primary Tumor', 'Primary Solid Tumor']))]
    df = pd.DataFrame({'sample': tum['sample'].tolist()})
    df['pat'] = df['sample'].str[:12]
    df['expr'] = df['sample'].map(expr)
    df = df.dropna(subset=['expr'])

    surv = pd.read_sql_query(f'SELECT row_names AS sample, OS, OS_status FROM "survival.{abbrev}"', conn)
    surv['pat'] = surv['sample'].str[:12]
    surv = surv.drop_duplicates(subset='pat').set_index('pat')[['OS', 'OS_status']]

    m = df.join(surv, on='pat').dropna(subset=['OS', 'OS_status']).copy()
    m['OS'] = m['OS'].astype(float); m['OS_status'] = m['OS_status'].astype(int)
    med = m['expr'].median()
    hi = m[m['expr'] >= med]; lo = m[m['expr'] < med]

    kmf_h = KaplanMeierFitter().fit(hi['OS'], hi['OS_status'], label=f'High (n={len(hi)})')
    kmf_l = KaplanMeierFitter().fit(lo['OS'], lo['OS_status'], label=f'Low (n={len(lo)})')
    lr = logrank_test(hi['OS'], lo['OS'], event_observed_A=hi['OS_status'], event_observed_B=lo['OS_status'])

    kmf_h.plot_survival_function(ax=ax, color='#d62728', lw=1.5)
    kmf_l.plot_survival_function(ax=ax, color='#1f77b4', lw=1.5)
    hr, q = hr_q[abbrev]
    ax.set_title(f'{abbrev} (n={len(m)})\nHR={hr:.2f}, q={q:.2e}, log-rank P={lr.p_value:.2g}', fontsize=7)
    ax.set_xlabel('Months'); ax.set_ylabel('Overall survival')
    ax.set_ylim(0, 1.05); ax.grid(alpha=0.2); ax.legend(fontsize=6)

fig.suptitle('Fig. S1 | Kaplan-Meier overall survival, SLC5A6 median split — 6 FDR-significant cancer types (real per-patient data)', fontsize=8)
plt.tight_layout()
save_pub(fig, 'FigS1_SMVT_survival_KM', outdir=OUT)
plt.close()
print('saved FigS1 (real KM, 6 types)')
conn.close()
