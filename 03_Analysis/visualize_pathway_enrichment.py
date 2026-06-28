#!/usr/bin/env python3
"""
SMVT Pathway Enrichment — Nature Style Visualization
=====================================================
Reactome (17) + GO (52) + KEGG (13) enrichment results
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
    'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 6.5,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5, 'legend.fontsize': 5,
    'figure.dpi': 300, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.08,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'axes.spines.top': False, 'axes.spines.right': False,
})

OUTDIR = Path(__file__).parent

# ═══ Colors ═══
C_REACTOME = '#E41A1C'
C_GO       = '#377EB8'
C_KEGG     = '#4DAF4A'
C_SLC      = '#1F77B4'
C_METAB    = '#D62728'
C_LIPID    = '#FF7F0E'
C_OTHER    = '#7F7F7F'

# ═══ Reactome data (17 terms) ═══
reactome_data = [
    ('Biotin transport and metabolism',          2.76e-16, 6, C_METAB),
    ('Metabolism of water-soluble vitamins',     1.29e-15, 9, C_METAB),
    ('Metabolism of vitamins and cofactors',     8.83e-14, 9, C_METAB),
    ('SLC-mediated transmembrane transport',     1.45e-07, 6, C_SLC),
    ('Vitamin B5 (pantothenate) metabolism',     1.35e-06, 3, C_METAB),
    ('Defects in vitamin and cofactor metabolism',1.82e-06, 3, C_METAB),
    ('Activation of SREBP (SREBF)',              1.34e-05, 3, C_LIPID),
    ('Regulation of cholesterol by SREBP',       3.03e-05, 3, C_LIPID),
    ('Carnitine shuttle',                        1.11e-04, 2, C_LIPID),
    ('SLC transporter disorders',                1.65e-04, 3, C_SLC),
    ('Integration of energy metabolism',         2.27e-04, 3, C_LIPID),
    ('Transport of vitamins/nucleosides',        5.59e-04, 2, C_METAB),
    ('Metabolism of steroids',                   6.44e-04, 3, C_LIPID),
    ('Disorders of transmembrane transporters',  7.47e-04, 3, C_SLC),
    ('Fatty acyl-CoA biosynthesis',              7.97e-04, 2, C_LIPID),
    ('SLC-mediated organic anion transport',     8.41e-04, 2, C_SLC),
    ('Fatty acid metabolism',                    9.35e-04, 3, C_LIPID),
]

# ═══ GO data (top 18 from 52) ═══
go_data = [
    ('Sulfur compound metabolic process',       1.53e-14, 11, C_METAB),
    ('Sulfur compound biosynthetic process',    8.97e-09, 6,  C_METAB),
    ('Organic anion transport',                 2.79e-09, 9,  C_SLC),
    ('Organic acid transport',                  7.12e-09, 8,  C_SLC),
    ('Carboxylic acid transport',               7.12e-09, 8,  C_SLC),
    ('Acyl-CoA biosynthetic process',           9.81e-09, 5,  C_LIPID),
    ('Vascular transport',                      1.03e-07, 5,  C_SLC),
    ('Branched-chain amino acid metab. process', 2.99e-08, 4, C_METAB),
    ('Small molecule catabolic process',        2.93e-07, 7,  C_METAB),
    ('Monocarboxylic acid metabolic process',   9.98e-07, 7,  C_METAB),
    ('Biotin metabolic process',                1.78e-06, 3,  C_METAB),
    ('Fatty acid metabolic process',            2.14e-05, 5,  C_LIPID),
    ('Lipid biosynthetic process',              4.92e-05, 6,  C_LIPID),
    ('Coenzyme metabolic process',              6.08e-05, 6,  C_METAB),
    ('Anion transmembrane transport',           8.13e-05, 5,  C_SLC),
    ('Acetyl-CoA metabolic process',            1.11e-04, 2,  C_LIPID),
    ('Biotin carboxyl carrier activity',        1.32e-04, 2,  C_METAB),
    ('Solute:inorganic anion antiporter',       2.78e-04, 2,  C_SLC),
]

# ═══ KEGG data (all 13) ═══
kegg_data = [
    ('Vitamin digestion and absorption',        2.15e-07, 4, C_METAB),
    ('Propanoate metabolism',                   5.14e-07, 4, C_LIPID),
    ('Valine, leucine, isoleucine degradation', 2.72e-06, 4, C_LIPID),
    ('Fatty acid biosynthesis',                 6.57e-06, 3, C_LIPID),
    ('AMPK signaling pathway',                  1.11e-04, 4, C_LIPID),
    ('Pyruvate metabolism',                     1.51e-04, 3, C_LIPID),
    ('Alcoholic liver disease',                 2.11e-04, 4, C_OTHER),
    ('Fatty acid metabolism',                   2.23e-04, 3, C_LIPID),
    ('Non-alcoholic fatty liver disease',       4.37e-03, 3, C_OTHER),
    ('Biosynthesis of cofactors',               5.09e-03, 3, C_METAB),
    ('Biotin metabolism',                       1.24e-02, 1, C_METAB),
    ('Citrate cycle (TCA cycle)',               5.83e-02, 1, C_LIPID),
    ('Central carbon metabolism in cancer',     1.01e-01, 3, C_OTHER),
]

# ═══════════════════════════════════════════════════
# FIGURE: 3-panel pathway enrichment
# ═══════════════════════════════════════════════════

fig = plt.figure(figsize=(9.0, 10.5))

# ── Panel a: Reactome (top 17) ──
ax_a = fig.add_axes([0.06, 0.62, 0.42, 0.35])
r_data = list(reversed(reactome_data[:15]))
labels_a = [d[0][:55] for d in r_data]
pvals_a  = [-np.log10(d[1]) for d in r_data]
counts_a = [d[2] for d in r_data]
colors_a = [d[3] for d in r_data]

y = range(len(labels_a))
ax_a.barh(y, pvals_a, height=0.7, color=colors_a, alpha=0.85, edgecolor='white', lw=0.3)
for i, (p, cnt) in enumerate(zip(pvals_a, counts_a)):
    ax_a.text(p + 0.3, i, f'n={cnt}', fontsize=4.5, va='center', color='#666')

ax_a.set_yticks(y); ax_a.set_yticklabels(labels_a, fontsize=5)
ax_a.set_xlabel('-log10(P-value)', fontsize=6.5)
ax_a.axvline(x=-np.log10(0.05), color='grey', lw=0.3, ls='--', alpha=0.4)
ax_a.text(-np.log10(0.05)+0.1, len(labels_a)-1, 'P=0.05', fontsize=4.5, color='grey')

leg = [Line2D([0],[0],c=C_METAB,lw=4,label='Vitamin/Cofactor'),
       Line2D([0],[0],c=C_SLC,lw=4,label='Transport'),
       Line2D([0],[0],c=C_LIPID,lw=4,label='Lipid/Energy')]
ax_a.legend(handles=leg, frameon=False, fontsize=5, loc='lower right')
ax_a.set_title('Reactome (17 pathways)', fontsize=8, fontweight='bold', loc='left')
ax_a.text(-0.08, 1.03, 'a', transform=ax_a.transAxes, fontsize=12, fontweight='bold')

# ── Panel b: GO (top 18 from 52) ──
ax_b = fig.add_axes([0.55, 0.62, 0.42, 0.35])
g_data = list(reversed(go_data))
labels_b = [d[0][:55] for d in g_data]
pvals_b  = [-np.log10(d[1]) for d in g_data]
counts_b = [d[2] for d in g_data]
colors_b = [d[3] for d in g_data]

yb = range(len(labels_b))
ax_b.barh(yb, pvals_b, height=0.7, color=colors_b, alpha=0.85, edgecolor='white', lw=0.3)
for i, (p, cnt) in enumerate(zip(pvals_b, counts_b)):
    ax_b.text(p + 0.3, i, f'n={cnt}', fontsize=4.5, va='center', color='#666')

ax_b.set_yticks(yb); ax_b.set_yticklabels(labels_b, fontsize=5)
ax_b.set_xlabel('-log10(P-value)', fontsize=6.5)
ax_b.axvline(x=-np.log10(0.05), color='grey', lw=0.3, ls='--', alpha=0.4)

leg2 = [Line2D([0],[0],c=C_METAB,lw=4,label='Metabolism'),
        Line2D([0],[0],c=C_SLC,lw=4,label='Transport'),
        Line2D([0],[0],c=C_LIPID,lw=4,label='Lipid/Acyl-CoA')]
ax_b.legend(handles=leg2, frameon=False, fontsize=5, loc='lower right')
ax_b.set_title('GO Biological Process (top 18 of 52)', fontsize=8, fontweight='bold', loc='left')
ax_b.text(-0.08, 1.03, 'b', transform=ax_b.transAxes, fontsize=12, fontweight='bold')

# ── Panel c: KEGG (all 13) ──
ax_c = fig.add_axes([0.06, 0.38, 0.42, 0.20])
k_data = list(reversed(kegg_data))
labels_c = [d[0][:50] for d in k_data]
pvals_c  = [-np.log10(d[1]) for d in k_data]
counts_c = [d[2] for d in k_data]
colors_c = [d[3] for d in k_data]

yc = range(len(labels_c))
ax_c.barh(yc, pvals_c, height=0.65, color=colors_c, alpha=0.85, edgecolor='white', lw=0.3)
for i, (p, cnt) in enumerate(zip(pvals_c, counts_c)):
    ax_c.text(p + 0.15, i, f'n={cnt}', fontsize=4.5, va='center', color='#666')

ax_c.set_yticks(yc); ax_c.set_yticklabels(labels_c, fontsize=5)
ax_c.set_xlabel('-log10(P-value)', fontsize=6.5)
ax_c.axvline(x=-np.log10(0.05), color='grey', lw=0.3, ls='--', alpha=0.4)
ax_c.set_title('KEGG (13 pathways)', fontsize=8, fontweight='bold', loc='left')
ax_c.text(-0.08, 1.03, 'c', transform=ax_c.transAxes, fontsize=12, fontweight='bold')

# ── Panel d: Cross-database convergence diagram ──
ax_d = fig.add_axes([0.55, 0.38, 0.42, 0.20])
ax_d.axis('off')

# Draw convergence schematic
boxes = [
    (0.25, 0.80, 'Reactome\nBiotin transport\nSLC-mediated\ntransport', C_REACTOME, 'white'),
    (0.50, 0.80, 'GO\nOrganic anion\ntransport\nAcyl-CoA biosynth', C_GO, 'white'),
    (0.75, 0.80, 'KEGG\nVitamin digestion\nFatty acid biosynth\nAMPK signaling', C_KEGG, 'white'),
]
for cx, cy, text, bg, fg in boxes:
    rect = FancyBboxPatch((cx-0.20, cy-0.15), 0.40, 0.28, boxstyle="round,pad=0.04",
                          fc=bg, alpha=0.15, ec=bg, lw=1.2, transform=ax_d.transAxes)
    ax_d.add_patch(rect)
    ax_d.text(cx, cy+0.17, text.split('\n')[0], transform=ax_d.transAxes,
              ha='center', fontsize=6, fontweight='bold', color=bg)
    ax_d.text(cx, cy-0.05, '\n'.join(text.split('\n')[1:]), transform=ax_d.transAxes,
              ha='center', fontsize=4.5, color='#555')

# Convergence point
ax_d.text(0.50, 0.35, 'SLC5A6 (SMVT)\nVitamin Transport → Lipid Synthesis → Tumor',
          transform=ax_d.transAxes, ha='center', fontsize=7, fontweight='bold', color='#333',
          bbox=dict(boxstyle='round,pad=0.3', fc='#FFF5E8', ec=C_METAB, lw=1.5))

# Arrows converging
for cx in [0.25, 0.50, 0.75]:
    ax_d.annotate('', xy=(0.50, 0.38), xytext=(cx, 0.62),
                  transform=ax_d.transAxes, arrowprops=dict(arrowstyle='->', color='#999', lw=0.8))

ax_d.set_title('Cross-Database Convergence', fontsize=8, fontweight='bold', loc='left')
ax_d.text(-0.04, 1.03, 'd', transform=ax_d.transAxes, fontsize=12, fontweight='bold')

# ── Panel e: SREBP mechanism schematic ──
ax_e = fig.add_axes([0.06, 0.10, 0.42, 0.24])
ax_e.axis('off')

pathway = [
    (0.15, 0.82, 'SMVT\nSLC5A6', C_SLC, 'white'),
    (0.15, 0.58, 'Biotin\nUptake', C_METAB, 'white'),
    (0.15, 0.34, 'ACC\nACACA/B', C_LIPID, 'white'),
    (0.50, 0.58, 'Malonyl-CoA\n(FASN substrate)', C_LIPID, '#333'),
    (0.50, 0.34, 'SREBP1\nCleavage', C_LIPID, 'white'),
    (0.50, 0.10, 'FASN, SCD, ACLY\nTranscription', C_LIPID, '#333'),
]
for cx, cy, text, bg, fg in pathway:
    is_box = fg == 'white'
    if is_box:
        rect = FancyBboxPatch((cx-0.10, cy-0.09), 0.20, 0.17,
                              boxstyle="round,pad=0.03", fc=bg, alpha=0.88, ec='white', lw=0.5,
                              transform=ax_e.transAxes)
        ax_e.add_patch(rect)
    ax_e.text(cx, cy, text, transform=ax_e.transAxes, ha='center', va='center',
              fontsize=5.5 if is_box else 5, color=fg, fontweight='bold' if is_box else 'normal')

# Arrows
arrows = [(0.15,0.73,0.15,0.67),(0.15,0.49,0.15,0.43),(0.25,0.46,0.40,0.52),
          (0.50,0.49,0.50,0.43)]
for x1, y1, x2, y2 in arrows:
    ax_e.annotate('', xy=(x2, y2), xytext=(x1, y1),
                  transform=ax_e.transAxes, arrowprops=dict(arrowstyle='->', color='#555', lw=0.8))

# SREBP callout
ax_e.text(0.78, 0.76, 'KEY FINDING:', fontsize=6, fontweight='bold', color=C_LIPID, transform=ax_e.transAxes)
ax_e.text(0.78, 0.64, 'SREBP links SMVT\nvitamin transport to\nlipid transcription\nprogram (Reactome #7\np=1.3e-5)', fontsize=4.8,
          color='#555', transform=ax_e.transAxes, fontstyle='italic')

ax_e.set_title('SMVT→Biotin→ACC→SREBP→FASN Axis', fontsize=8, fontweight='bold', loc='left')
ax_e.text(-0.04, 1.03, 'e', transform=ax_e.transAxes, fontsize=12, fontweight='bold')

# ── Panel f: Top hits comparison ──
ax_f = fig.add_axes([0.55, 0.10, 0.42, 0.24])

# Top-5 from each DB
top_hits = [
    ('Biotin transport\nand metabolism', 15.56, C_REACTOME, 'Reactome'),
    ('Sulfur compound\nmetabolic process', 13.82, C_GO, 'GO'),
    ('Vitamin digestion\nand absorption', 6.67, C_KEGG, 'KEGG'),
    ('Metab. of water-\nsoluble vitamins', 14.89, C_REACTOME, 'Reactome'),
    ('Organic anion\ntransport', 8.55, C_GO, 'GO'),
    ('Propanoate\nmetabolism', 6.29, C_KEGG, 'KEGG'),
    ('SLC-mediated\ntransmembrane transp.', 6.84, C_REACTOME, 'Reactome'),
    ('Acyl-CoA\nbiosynthetic process', 8.01, C_GO, 'GO'),
    ('Fatty acid\nbiosynthesis', 5.18, C_KEGG, 'KEGG'),
]
labels_f = [t[0] for t in top_hits]
pvals_f  = [t[1] for t in top_hits]
colors_f = [t[2] for t in top_hits]
dbs_f    = [t[3] for t in top_hits]

yf = range(len(labels_f)-1, -1, -1)
ax_f.barh(yf, pvals_f, height=0.6, color=colors_f, alpha=0.85, edgecolor='white', lw=0.3)
ax_f.set_yticks(yf); ax_f.set_yticklabels(labels_f, fontsize=5)

for i, (p, db) in enumerate(zip(pvals_f, dbs_f)):
    ax_f.text(p + 0.2, yf[i], db, fontsize=4.5, va='center', color='#666', fontstyle='italic')

ax_f.set_xlabel('-log10(P-value)', fontsize=6.5)
ax_f.set_title('Top-9 Cross-Database Hits', fontsize=8, fontweight='bold', loc='left')
ax_f.text(-0.04, 1.03, 'f', transform=ax_f.transAxes, fontsize=12, fontweight='bold')

# ── Global Title ──
fig.suptitle('SLC5A6 Interactome — Pathway Enrichment (Reactome + GO + KEGG)',
             fontsize=10.5, fontweight='bold', x=0.06, y=0.998, ha='left')
fig.text(0.06, 0.005, 'pathlinkR v1.8.0 + clusterProfiler + ReactomePA  |  SLC5A6 interactome (25 genes)  |  P-value cutoff: 0.5 (Reactome), 0.3 (GO), 0.3 (KEGG)',
         fontsize=4.5, color='grey', style='italic')

out = OUTDIR / 'Fig_pathlinkR_master.png'
fig.savefig(out, dpi=600, facecolor='white')
fig.savefig(OUTDIR / 'Fig_pathlinkR_master.pdf', facecolor='white')
print(f'OK: {out.name}')
plt.close()
print('Pathway enrichment visualization complete.')
