#!/usr/bin/env python3
"""Regenerate FigS16 (survival vs expression) from REAL data.
Scatter: per-cancer median SLC5A6 expression (33 types) vs Cox HR.
Box: median expression in FDR-significant (q<0.05, 6 types) vs non-significant (27).
"""
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

OUT = '04_Manuscript/figures'
expr = json.load(open('03_Analysis/data/slc5a6_by_cancer.json'))  # {cancer: [per-sample log2 expr]}
surv = list(csv.DictReader(open('03_Analysis/outputs/pancan_survival.csv')))

# per-cancer median expression + HR + FDR q
rows = []
for r in surv:
    c = r['cancer']
    if c not in expr:
        continue
    med = float(np.median(expr[c]))
    rows.append({'c': c, 'med_expr': med, 'hr': float(r['hr']), 'cox_p': float(r['cox_p'])})
rows.sort(key=lambda x: x['cox_p'])
n = len(rows)
ps = np.array([x['cox_p'] for x in rows])
order = np.argsort(ps); rank = np.empty(n); rank[order] = np.arange(1, n+1)
q = ps * n / rank; qs = q[order]
for i in range(n-2, -1, -1): qs[i] = min(qs[i], qs[i+1])
q[order] = qs
for i, x in enumerate(rows): x['q'] = q[i]

sig = [x for x in rows if x['q'] < 0.05]
nsig = [x for x in rows if x['q'] >= 0.05]
print(f'{len(sig)} significant, {len(nsig)} not (n={len(rows)})')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# (a) Scatter: HR (log) vs median expression
hrs = np.array([x['hr'] for x in rows]); mex = np.array([x['med_expr'] for x in rows])
sig_flag = np.array([x['q'] < 0.05 for x in rows])
ax1.scatter(mex[~sig_flag], hrs[~sig_flag], c='#9e9e9e', s=40, alpha=0.7, label='q≥0.05')
ax1.scatter(mex[sig_flag], hrs[sig_flag], c='#c62828', s=60, zorder=5, label='q<0.05')
for x in rows:
    ax1.annotate(x['c'], (x['med_expr'], x['hr']), fontsize=6, alpha=0.6,
                 textcoords='offset points', xytext=(3, 3))
rho, p = stats.spearmanr(mex, hrs)
ax1.set_yscale('log'); ax1.set_xlabel('Median SLC5A6 expression (log2)')
ax1.set_ylabel('Cox HR (log scale)')
ax1.set_title(f'(a) HR vs expression (ρ={rho:.2f}, P={p:.3g})', fontsize=10)
ax1.axhline(1, color='gray', ls='--', lw=0.8)
ax1.legend(fontsize=8); ax1.grid(alpha=0.2)

# (b) Box plot: significant vs not
sig_med = np.array([x['med_expr'] for x in sig]); nsig_med = np.array([x['med_expr'] for x in nsig])
bp = ax2.boxplot([nsig_med, sig_med], labels=[f'Not significant (n={len(nsig)})', f'q<0.05 (n={len(sig)})'],
                 patch_artist=True, widths=0.5)
for patch, col in zip(bp['boxes'], ['#e0e0e0', '#ef9a9a']):
    patch.set_facecolor(col)
for x, y in zip([1]*len(nsig_med), nsig_med):
    ax2.scatter(x + np.random.normal(0, 0.05, 1), y, c='#9e9e9e', s=30, alpha=0.6)
for x, y in zip([2]*len(sig_med), sig_med):
    ax2.scatter(x + np.random.normal(0, 0.05, 1), y, c='#c62828', s=30, alpha=0.6)
mw, mp = stats.mannwhitneyu(nsig_med, sig_med, alternative='two-sided')
ax2.set_ylabel('Median SLC5A6 expression (log2)')
ax2.set_title(f'(b) Expression by survival status (Mann-Whitney P={mp:.2f})', fontsize=10)
ax2.grid(alpha=0.2, axis='y')

fig.suptitle('Fig. S16 | SLC5A6 survival association is strongest in low-baseline-expression types', fontsize=11)
plt.tight_layout()
for fmt in ('png', 'pdf'):
    fig.savefig(f'{OUT}/FigS16_survival_vs_expression.{fmt}', dpi=300, bbox_inches='tight')
plt.close()
print('saved FigS16 (real data)')
