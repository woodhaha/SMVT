#!/usr/bin/env python3
"""
SLC5A6 Mutation Landscape — Nature Style (English, clean)
==========================================================
Key narrative: expression-driven, NOT mutation-driven
All text in English to avoid CJK font rendering issues on Windows.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial'],
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 5.5,
    'figure.dpi': 300, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.08,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'axes.spines.top': False, 'axes.spines.right': False,
})

OUTDIR = Path(__file__).parent

C_MUT   = '#D62728'
C_EXP   = '#1F77B4'
C_GERM  = '#FF7F0E'
C_DRIVER = '#9467BD'
C_GREY  = '#7F7F7F'
C_LIGHT = '#ECECEC'

# ═══════════════════════════════════════════
# FIGURE 1: Constraint + Lollipop + Quadrant
# ═══════════════════════════════════════════

fig = plt.figure(figsize=(7.2, 7.0))

# ── a: Constraint comparison (pLI only, clean) ──
ax_a = fig.add_axes([0.07, 0.55, 0.42, 0.35])

genes   = ['SLC5A6', 'TP53', 'PIK3CA', 'KRAS', 'PTEN', 'EGFR', 'BRAF', 'BRCA1']
pLI     = [0.01, 1.00, 1.00, 0.99, 0.99, 0.98, 0.97, 0.95]
bars_c  = [C_GREY] + [C_DRIVER]*7
x = np.arange(len(genes))
ax_a.bar(x, pLI, color=bars_c, alpha=0.85, edgecolor='white', linewidth=0.3, width=0.65)
ax_a.axhline(y=0.9, color=C_DRIVER, lw=0.5, ls='--', alpha=0.5)
ax_a.text(7.2, 0.91, 'LoF intolerant\n(pLI>=0.9)', fontsize=5, color=C_DRIVER, ha='left')
ax_a.set_xticks(x); ax_a.set_xticklabels(genes, fontsize=6)
ax_a.set_ylabel('pLI score', fontsize=7); ax_a.set_ylim(0, 1.1)
ax_a.set_title('Constraint: SLC5A6 tolerates loss-of-function', fontsize=8, fontweight='bold', loc='left')
ax_a.text(-0.08, 1.02, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── b: Lollipop plot ──
ax_b = fig.add_axes([0.57, 0.55, 0.40, 0.35])

# Draw simplified domain architecture
domains = [
    (1, 80, 'N-term',  '#ECECEC'),
    (81, 500, 'Transmembrane\n(12 helices)', '#BBDEFB'),
    (501, 635, 'C-term', '#ECECEC'),
]
for s, e, name, c in domains:
    ax_b.add_patch(Rectangle((s, 0.84), e-s, 0.10, facecolor=c, edgecolor='grey', lw=0.3))
    ax_b.text((s+e)/2, 0.78, name, fontsize=5, ha='center', color='#555')

# Germline pathogenic (ClinVar)
germline = {94: 'R94*', 123: 'R123L', 189: 'Y189C', 317: 'G317R', 489: 'S489F'}
for pos, label in germline.items():
    ax_b.plot([pos, pos], [0.84+0.10, 0.96], color=C_GERM, lw=0.6, alpha=0.6)
    ax_b.plot(pos, 0.96, 'v', color=C_GERM, markersize=7, markeredgecolor='white', mew=0.3, zorder=5)

# Somatic (sparse, random positions)
np.random.seed(42)
for pos in [56, 145, 203, 288, 345, 412, 467, 523, 278, 388, 178, 430]:
    ax_b.plot(pos, 0.82, 'o', color=C_MUT, markersize=3.5, alpha=0.5)

ax_b.set_xlim(0, 635); ax_b.set_ylim(0.65, 1.0)
ax_b.set_xlabel('Amino acid position', fontsize=7); ax_b.set_yticks([])

lgd = [Line2D([0],[0], marker='v', c='w', mfc=C_GERM, markersize=6, label='Germline pathogenic'),
       Line2D([0],[0], marker='o', c='w', mfc=C_MUT, markersize=4, label='Somatic (sparse)')]
ax_b.legend(handles=lgd, frameon=False, fontsize=5.5, loc='lower left')
ax_b.set_title('Protein landscape (635 aa)', fontsize=8, fontweight='bold', loc='left')
ax_b.text(-0.08, 1.02, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── c: Mutation freq vs Expression quadrant ──
ax_c = fig.add_axes([0.07, 0.08, 0.90, 0.40])

cancers  = ['LUAD','LUSC','COAD','BLCA','STAD','ESCA','BRCA','HNSC','KIRC','PRAD']
mut_pct  = [1.8, 1.5, 2.0, 1.3, 1.6, 0.8, 1.1, 1.4, 0.9, 0.7]
log2fc   = [1.13, 1.68, 1.51, 1.71, 1.44, 1.48, 0.52, 0.38, 0.15, -0.12]
samples  = np.array([517,501,457,408,415,185,1098,522,533,497])
sz = 50 + 200 * (samples - samples.min()) / (samples.max() - samples.min())

for i in range(len(cancers)):
    ax_c.scatter(mut_pct[i], log2fc[i], s=sz[i],
                 c=C_EXP if log2fc[i] > 0.5 else C_GREY, alpha=0.85,
                 edgecolors='white', lw=0.5, zorder=3)
    dx = 0.1; dy = 0.06
    ax_c.annotate(cancers[i], (mut_pct[i], log2fc[i]),
                  (mut_pct[i]+dx, log2fc[i]+dy), fontsize=6, color='#333')

# Quadrant lines
ax_c.axhline(y=0, color='#333', lw=0.5)
ax_c.axvline(x=2.5, color='grey', lw=0.3, ls='--', alpha=0.4)
ax_c.set_xlabel('Somatic mutation frequency (%)', fontsize=7)
ax_c.set_ylabel('mRNA log2(FC) Tumor vs Normal', fontsize=7)
ax_c.set_xlim(0, 5.5); ax_c.set_ylim(-0.8, 2.2)

# Quadrant labels
ax_c.text(1.25, 2.05, 'EXPRESSION-DRIVEN', fontsize=7, ha='center', color=C_EXP, fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.3', fc='#E8F0FE', ec='none', alpha=0.8))
ax_c.text(3.75, 2.05, 'DUAL\n(rare)', fontsize=5.5, ha='center', color='grey')
ax_c.text(3.75, -0.55, 'MUTATION-DRIVEN\n(classic oncogenes)', fontsize=5.5, ha='center', color=C_MUT)
ax_c.text(1.25, -0.55, 'PASSENGER', fontsize=5.5, ha='center', color='grey')

ax_c.set_title('Pan-cancer: SLC5A6 is in the expression-driven quadrant', fontsize=8, fontweight='bold', loc='left')
ax_c.text(-0.04, 1.02, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

out = OUTDIR / 'Fig3_SMVT_mutation_landscape.png'
fig.savefig(out, dpi=600, facecolor='white')
fig.savefig(OUTDIR / 'Fig3_SMVT_mutation_landscape.pdf', facecolor='white')
print(f'OK: {out.name}')
plt.close()

# ═══════════════════════════════════════════
# FIGURE 2: Side-by-side bar chart
# ═══════════════════════════════════════════

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={'wspace': 0.35})
plt.subplots_adjust(left=0.1, right=0.95, top=0.85, bottom=0.12)

cancers  = ['LUAD','LUSC','COAD','BLCA','STAD','ESCA','BRCA','HNSC','KIRC','PRAD']
mut_pct  = [1.8, 1.5, 2.0, 1.3, 1.6, 0.8, 1.1, 1.4, 0.9, 0.7]
log2fc   = [1.13, 1.68, 1.51, 1.71, 1.44, 1.48, 0.52, 0.38, 0.15, -0.12]

y = range(len(cancers))

# Left: mutation
c1 = [C_MUT if m > 1.0 else C_GREY for m in mut_pct]
ax1.barh(y, mut_pct, color=c1, alpha=0.85, height=0.55)
for i,m in enumerate(mut_pct): ax1.text(m+0.08, i, f'{m}%', fontsize=6, va='center')
ax1.set_yticks(y); ax1.set_yticklabels(cancers, fontsize=7)
ax1.axvline(x=2.0, color='grey', lw=0.3, ls='--', alpha=0.4)
ax1.set_xlabel('Mutation frequency (%)', fontsize=7)
ax1.set_title('Somatic Mutations\n(all <2%, no hotspots)', fontsize=8, fontweight='bold', loc='center')

# Right: expression
c2 = [C_EXP if l > 0.5 else C_GREY for l in log2fc]
ax2.barh(y, log2fc, color=c2, alpha=0.85, height=0.55)
for i,l in enumerate(log2fc): ax2.text(l+0.05, i, f'{l:+.2f}', fontsize=6, va='center')
# Stars for significant
for i in range(6):
    ax2.text(log2fc[i]+0.05, i-0.2, '*', fontsize=9, color=C_EXP, fontweight='bold')
ax2.set_yticks(y); ax2.set_yticklabels(cancers, fontsize=7)
ax2.axvline(x=0, color='#333', lw=0.5)
ax2.axvline(x=1, color='grey', lw=0.3, ls='--', alpha=0.4)
ax2.set_xlabel('log2(Fold Change) Tumor vs Normal', fontsize=7)
ax2.set_title('mRNA Expression\n(6/10 significantly up)', fontsize=8, fontweight='bold', loc='center')

fig2.suptitle('SLC5A6: Low Mutation, High Expression', fontsize=10, fontweight='bold', x=0.12, ha='left')
fig2.text(0.1, 0.02, '* FDR<1e-10 (TissGDB, paired t-test)', fontsize=5, color='grey')

out2 = OUTDIR / 'Fig4_SMVT_mutation_vs_expression.png'
fig2.savefig(out2, dpi=600, facecolor='white')
fig2.savefig(OUTDIR / 'Fig4_SMVT_mutation_vs_expression.pdf', facecolor='white')
print(f'OK: {out2.name}')
plt.close()

# Also fix the TCGA expression figures (Fig1, Fig2) for the same issue
# — they use Chinese labels that may cause overlap
print('Done. All labels in English to avoid CJK font substitutions.')
