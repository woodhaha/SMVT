"""Compare ESM-2 35M vs 650M LLR for SLC5A6 missense. Robustness check of 'passenger' claim."""
import pandas as pd
import numpy as np
from scipy import stats

small = pd.read_csv('D:/Researching/SMVT/03_Analysis/outputs/esm_pathogenicity_mc3.csv')
big = pd.read_csv('D:/Researching/SMVT/03_Analysis/outputs/esm650_mc3.csv')
struct = pd.read_csv('D:/Researching/SMVT/03_Analysis/outputs/esm_mutations_structural.csv')

# merge on (pos, ref, alt)
m = pd.merge(small[['pos','ref','alt','hgvsp','esm_llr','esm_zscore','blosum62','composite','composite_percentile']],
             big[['pos','ref','alt','esm650_llr','esm650_percentile']],
             on=['pos','ref','alt'])
m = m.merge(struct[['pos','ref','alt','struct_dist']], on=['pos','ref','alt'], how='left')
m['domain'] = small.set_index(['pos','ref','alt'])['domain'].reindex(pd.MultiIndex.from_frame(m[['pos','ref','alt']])).values
print(f'merged: {len(m)} rows (expect 87 incl control)')

# 1. model agreement
r, p = stats.spearmanr(m['esm_llr'], m['esm650_llr'])
print(f'\n[1] 35M vs 650M LLR Spearman rho={r:.3f}, p={p:.2e}')

# 2. is the distribution still uniform (passenger) in 650M?
from scipy.stats import shapiro, kstest
sk, sp = shapiro(m['esm650_llr'])
ks, kp = kstest(m['esm650_llr'], 'norm')
print(f'[2] 650M LLR normality: Shapiro W={sk:.3f} p={sp:.3f}; KS p={kp:.3f}')

# 3. control rank in 650M (match by pos 123 / R>L)
ctl = m[(m['pos'] == 123) & (m['ref'] == 'R') & (m['alt'] == 'L')]
if len(ctl):
    print(f'[3] R123L control: 35M percentile={ctl.iloc[0]["composite_percentile"]:.0f}, 650M percentile={ctl.iloc[0]["esm650_percentile"]:.0f}')

# 4. structural avoidance: distance from biotin site vs pathogenicity
ms = m[m['struct_dist'].notna() & (m['hgvsp'] != 'p.R123L')]
if len(ms) > 5:
    rr, rp = stats.spearmanr(ms['esm650_llr'], ms['struct_dist'])
    print(f'[4] 650M: pathogenicity vs biotin-site distance rho={rr:.3f} p={rp:.3f} (n={len(ms)})')
    # fraction of most-destabilizing that are >10A
    thresh = ms['esm650_llr'].quantile(0.25)
    far = ms[ms['esm650_llr'] <= thresh]
    print(f'    bottom-25% destabilizing: {len(far)} variants, {len(far[far["struct_dist"]>10])}/{len(far)} >10A from biotin site')

# 5. directional consistency (top/bottom quartile overlap)
q_low_35 = set(m[m['composite_percentile'] <= 25]['pos'])
q_low_65 = set(m[m['esm650_percentile'] <= 25]['pos'])
inter = len(q_low_35 & q_low_65)
print(f'[5] bottom-quartile overlap: {inter}/{len(q_low_35)} (35M) ∩ {len(q_low_65)} (650M)')

print('\n=== Top 10 most destabilizing (650M consensus) ===')
print(m.sort_values('esm650_llr').head(10)[['hgvsp','pos','domain','esm_llr','esm650_llr','composite_percentile','esm650_percentile','struct_dist']].to_string(index=False))
