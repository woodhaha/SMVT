"""Pan-cancer SLC5A6: differential expression + survival analysis."""
import json, os, csv, sqlite3, warnings, sys
from collections import defaultdict
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

CRC = r'D:\Researching\CRCproject\02_Data\raw'
OUT = r'D:\Researching\SMVT\03_Analysis'
os.makedirs(f'{OUT}/data', exist_ok=True)
os.makedirs(f'{OUT}/figures', exist_ok=True)

# 1. Load SLC5A6 expression from TOIL
print('Loading SLC5A6 expression from TOIL...')
with open(f'{CRC}/TcgaTargetGtex_rsem_gene_fpkm', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if row[0].startswith('ENSG00000138074'):
            slc_vals = {}
            for i, v in enumerate(row[1:], 1):
                try:
                    val = float(v)
                    if val > 0:
                        slc_vals[header[i]] = np.log2(val + 1)
                except: pass
            break
print(f'SLC5A6 expression: {len(slc_vals)} samples')

# 2. Load phenotype (sample → cancer type mapping)
pheno = pd.read_csv(f'{CRC}/TcgaTargetGTEX_phenotype.txt', sep='\t')
print(f'Phenotype: {len(pheno)} samples')
# Map sample to cancer type
sample_to_cancer = dict(zip(pheno['sample'], pheno['_primary_sic']))

# 3. Classify samples as TCGA tumor/normal
tcga_tumor, tcga_normal = [], []
all_tcga = {}
for sid, expr in slc_vals.items():
    if sid.startswith('TCGA'):
        cancer = sample_to_cancer.get(sid, 'UNKNOWN')
        all_tcga[sid] = {'expr': expr, 'cancer': cancer}
        if sid.endswith('-01') or sid.endswith('-03'):
            tcga_tumor.append(sid)
        elif sid.endswith('-11') or sid.endswith('-10'):
            tcga_normal.append(sid)

print(f'TCGA tumor: {len(tcga_tumor)}, normal: {len(tcga_normal)}')

# 4. Differential expression: tumor vs normal per cancer type
print('\n=== DIFFERENTIAL EXPRESSION ===')
cancer_expr = defaultdict(lambda: {'tumor': [], 'normal': []})
for sid in tcga_tumor:
    d = all_tcga[sid]
    cancer_expr[d['cancer']]['tumor'].append(d['expr'])
for sid in tcga_normal:
    d = all_tcga[sid]
    cancer_expr[d['cancer']]['normal'].append(d['expr'])

de_results = []
for cancer, data in sorted(cancer_expr.items()):
    tv, nv = data['tumor'], data['normal']
    if len(tv) >= 10 and len(nv) >= 3:  # min sample thresholds
        tm, nm = np.mean(tv), np.mean(nv)
        l2fc = tm - nm
        from scipy.stats import mannwhitneyu
        _, p = mannwhitneyu(tv, nv, alternative='two-sided')
        de_results.append({'cancer': cancer, 'n_tumor': len(tv), 'n_normal': len(nv),
                          'tumor_mean': round(tm, 3), 'normal_mean': round(nm, 3),
                          'log2FC': round(l2fc, 3), 'pval': p,
                          'fdr': min(p * 100, 1.0)})  # approximate

de_results.sort(key=lambda x: x['log2FC'], reverse=True)
print(f'{cancer:10s} {"N_t":>5s} {"N_n":>5s} {"T_mean":>7s} {"N_mean":>7s} {"Log2FC":>8s} {"P":>10s}')
print('-'*55)
for r in de_results[:20]:
    sig = '*' if r['pval'] < 0.05 else ''
    print(f'{r["cancer"]:10s} {r["n_tumor"]:5d} {r["n_normal"]:5d} '
          f'{r["tumor_mean"]:7.3f} {r["normal_mean"]:7.3f} '
          f'{r["log2FC"]:+8.3f} {r["pval"]:10.2e}{sig}')

outpath = f'{OUT}/outputs/pancancer_de.csv'
pd.DataFrame(de_results).to_csv(outpath, index=False)
print(f'Saved: {outpath}')

# 5. Survival analysis per cancer type
print('\n=== SURVIVAL ANALYSIS ===')
conn = sqlite3.connect(f'{CRC}/TCGA.survival.sqlite')
tables = pd.read_sql('SELECT name FROM sqlite_master WHERE type="table"', conn)
surv_results = []

for tbl_name in tables['name']:
    abbrev = tbl_name.split('.')[1]
    # Map cancer type to TCGA study abbreviation
    surv = pd.read_sql(f'SELECT * FROM "{tbl_name}"', conn)
    surv.columns = ['pid', 'OS', 'OS_event', 'RFS', 'RFS_event']
    # Already in months
    # Map sample IDs to patient IDs (match on TCGA-XX-XXXX)
    pid_set = set(surv['pid'].str[:12])
    matched = 0
    expr_list = []
    for sid, d in all_tcga.items():
        sp = sid[:12]
        if sp in pid_set:
            expr_list.append({'pid': sid[:12], 'cancer': d['cancer'], 'expr': d['expr']})
            matched += 1

    if matched < 30:
        continue

    expr_df = pd.DataFrame(expr_list)
    merged = surv.merge(expr_df, on='pid')
    merged = merged.dropna(subset=['OS', 'OS_event', 'expr']).copy()
    n = len(merged)
    if n < 30:
        continue

    med = merged['expr'].median()
    merged['grp'] = np.where(merged['expr'] >= med, 'High', 'Low')
    h = merged[merged['grp'] == 'High']
    l = merged[merged['grp'] == 'Low']

    kh = KaplanMeierFitter().fit(h['OS'], h['OS_event'])
    kl = KaplanMeierFitter().fit(l['OS'], l['OS_event'])
    lr = logrank_test(h['OS'], l['OS'],
                      event_observed_A=h['OS_event'], event_observed_B=l['OS_event'])

    cph = CoxPHFitter()
    cx = merged[['OS', 'OS_event', 'expr']].rename(columns={'OS': 't', 'OS_event': 'e', 'expr': 'x'})
    try:
        cph.fit(cx, duration_col='t', event_col='e')
        hr = float(np.exp(cph.params_.iloc[0]))
        hl = float(np.exp(cph.confidence_intervals_.iloc[0, 0]))
        hu = float(np.exp(cph.confidence_intervals_.iloc[0, 1]))
        cp = float(cph.summary.p.iloc[0])
    except:
        hr = hl = hu = cp = float('nan')

    sig = lr.p_value < 0.05
    surv_results.append({
        'cancer': abbrev, 'n': n, 'events': int(merged['OS_event'].sum()),
        'med_hi': kh.median_survival_time_, 'med_lo': kl.median_survival_time_,
        'logrank': lr.p_value, 'hr': hr, 'hr_lo': hl, 'hr_hi': hu, 'cox_p': cp, 'sig': sig
    })
    print(f'{abbrev:6s} n={n:4d} ev={int(merged["OS_event"].sum()):3d} '
          f'med={kh.median_survival_time_:5.1f}/{kl.median_survival_time_:5.1f} '
          f'P={lr.p_value:.4f} HR={hr:.4f}" + " *" if sig else ""')

    # KM plot
    fig, ax = plt.subplots(figsize=(7, 5))
    kh.plot(ax=ax, color='#d62728', lw=1.5, label=f'High (n={len(h)})')
    kl.plot(ax=ax, color='#1f77b4', lw=1.5, label=f'Low (n={len(l)})')
    ax.set_xlabel('Months'); ax.set_ylabel('OS')
    ax.set_title(f'{abbrev}: SLC5A6\nHR={hr:.3f}, P={lr.p_value:.4f}')
    ax.legend()
    fig.savefig(f'{OUT}/figures/KM_{abbrev}_SMVT_survival.pdf', dpi=300, bbox_inches='tight')
    plt.close()

pd.DataFrame(surv_results).to_csv(f'{OUT}/outputs/pancancer_survival.csv', index=False)

sig_cts = [r['cancer'] for r in surv_results if r['sig']]
print(f'\nTotal cancers with analysis: {len(surv_results)}')
print(f'Significant (P<0.05): {sig_cts}')
print(f'Results saved to {OUT}/outputs/')
