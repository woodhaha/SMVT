"""Recompute pan-cancer SLC5A6 DE with BH-FDR (unpaired Mann-Whitney U) — manuscript Table 1 source."""
import json
import numpy as np
from scipy import stats

T = json.load(open('D:/Researching/SMVT/03_Analysis/data/slc5a6_tn.json'))
tumors, normals = T['tumors'], T['normals']

rows = []
for ab in tumors:
    t = np.array(tumors[ab], dtype=float)
    n = np.array(normals.get(ab, []), dtype=float)
    if len(t) < 10 or len(n) < 3:
        continue
    fc = np.median(t) - np.median(n)
    p = stats.mannwhitneyu(t, n, alternative='two-sided').pvalue
    rows.append({'cancer': ab, 'n_t': len(t), 'n_n': len(n), 'log2FC': fc, 'pval': p})

# BH FDR correction (proper order)
ps = np.array([r['pval'] for r in rows])
m = len(ps)
order = np.argsort(ps)
ps_sorted = ps[order]
fdr_sorted = ps_sorted * m / np.arange(1, m + 1)
fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
fdr = np.empty_like(fdr_sorted)
fdr[order] = fdr_sorted
for r, f in zip(rows, fdr):
    r['FDR'] = f
    r['sig'] = '**' if f < 0.001 else '*' if f < 0.05 else ''

rows.sort(key=lambda r: -r['log2FC'])

print(f"{'cancer':<6}{'n_t':>5}{'n_n':>4}{'log2FC':>8}{'pval':>11}{'FDR':>11}  sig")
n_sig_up = 0
for r in rows:
    print(f"{r['cancer']:<6}{r['n_t']:>5}{r['n_n']:>4}{r['log2FC']:>+8.2f}{r['pval']:>11.2e}{r['FDR']:>11.2e}  {r['sig']}")
    if r['log2FC'] > 0 and r['FDR'] < 0.05:
        n_sig_up += 1
print(f"\nSignificant UP (FDR<0.05, log2FC>0): {n_sig_up} / {len(rows)}")

json.dump(rows, open('D:/Researching/SMVT/03_Analysis/outputs/de_table_fdr.json', 'w'), indent=1, default=float)
print('saved -> 03_Analysis/outputs/de_table_fdr.json')
