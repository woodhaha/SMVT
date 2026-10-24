"""FigS4 (GO enrichment) + FigS5 (KEGG + Reactome) — unified publication style.

Splits the combined pathlinkR master into the two manuscript figures.
Data mirror the manuscript Results section (pathlinkR, SLC5A6 interactome).
"""
import sys
sys.path.insert(0, 'D:/Researching/SMVT/03_Analysis')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig_style import apply_style, PALETTE
apply_style()

OUT = r'D:\Researching\SMVT\04_Manuscript\figures'

# category colors (unified palette)
C_METAB = PALETTE['red']      # vitamin/cofactor metabolism
C_TRANS = PALETTE['blue']     # transport
C_LIPID = PALETTE['orange']   # lipid/energy
C_OTHER = PALETTE['gray']

go_data = [
    ('Sulfur compound metabolic process',       1.53e-14, 11, C_METAB),
    ('Sulfur compound biosynthetic process',    8.97e-09, 6,  C_METAB),
    ('Organic anion transport',                 2.79e-09, 9,  C_TRANS),
    ('Organic acid transport',                  7.12e-09, 8,  C_TRANS),
    ('Carboxylic acid transport',               7.12e-09, 8,  C_TRANS),
    ('Acyl-CoA biosynthetic process',           9.81e-09, 5,  C_LIPID),
    ('Vascular transport',                      1.03e-07, 5,  C_TRANS),
    ('Branched-chain amino acid metab. process', 2.99e-08, 4, C_METAB),
    ('Small molecule catabolic process',        2.93e-07, 7,  C_METAB),
    ('Monocarboxylic acid metabolic process',   9.98e-07, 7,  C_METAB),
    ('Biotin metabolic process',                1.78e-06, 3,  C_METAB),
    ('Fatty acid metabolic process',            2.14e-05, 5,  C_LIPID),
    ('Lipid biosynthetic process',              4.92e-05, 6,  C_LIPID),
    ('Coenzyme metabolic process',              6.08e-05, 6,  C_METAB),
    ('Anion transmembrane transport',           8.13e-05, 5,  C_TRANS),
    ('Acetyl-CoA metabolic process',            1.11e-04, 2,  C_LIPID),
    ('Biotin carboxyl carrier activity',        1.32e-04, 2,  C_METAB),
    ('Solute:inorganic anion antiporter',       2.78e-04, 2,  C_TRANS),
]

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

reactome_data = [
    ('Biotin transport and metabolism',          2.76e-16, 6, C_METAB),
    ('Metabolism of water-soluble vitamins',     1.29e-15, 9, C_METAB),
    ('Metabolism of vitamins and cofactors',     8.83e-14, 9, C_METAB),
    ('SLC-mediated transmembrane transport',     1.45e-07, 6, C_TRANS),
    ('Vitamin B5 (pantothenate) metabolism',     1.35e-06, 3, C_METAB),
    ('Defects in vitamin and cofactor metabolism',1.82e-06, 3, C_METAB),
    ('Activation of SREBP (SREBF)',              1.34e-05, 3, C_LIPID),
    ('Regulation of cholesterol by SREBP',       3.03e-05, 3, C_LIPID),
    ('Carnitine shuttle',                        1.11e-04, 2, C_LIPID),
    ('SLC transporter disorders',                1.65e-04, 3, C_TRANS),
    ('Integration of energy metabolism',         2.27e-04, 3, C_LIPID),
    ('Transport of vitamins/nucleosides',        5.59e-04, 2, C_METAB),
    ('Metabolism of steroids',                   6.44e-04, 3, C_LIPID),
    ('Disorders of transmembrane transporters',  7.47e-04, 3, C_TRANS),
    ('Fatty acyl-CoA biosynthesis',              7.97e-04, 2, C_LIPID),
    ('SLC-mediated organic anion transport',     8.41e-04, 2, C_TRANS),
    ('Fatty acid metabolism',                    9.35e-04, 3, C_LIPID),
]


def barh_terms(ax, data, title):
    d = list(reversed(data))
    labels = [t[0][:55] for t in d]
    pvals = [-np.log10(t[1]) for t in d]
    colors = [t[3] for t in d]
    y = range(len(labels))
    ax.barh(y, pvals, height=0.7, color=colors, alpha=0.85, edgecolor='white', lw=0.3)
    for i, t in enumerate(d):
        ax.text(-np.log10(t[1]) + 0.25, i, f'n={t[2]}', fontsize=6, va='center', color='#555')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('−log$_{10}$(P-value)')
    ax.axvline(x=-np.log10(0.05), color='grey', lw=0.6, ls='--', alpha=0.4)
    ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
    ax.tick_params(width=0.6)
    return ax


def cat_legend(ax):
    handles = [Line2D([0], [0], c=C_METAB, lw=4, label='Vitamin/Cofactor'),
               Line2D([0], [0], c=C_TRANS, lw=4, label='Transport'),
               Line2D([0], [0], c=C_LIPID, lw=4, label='Lipid/Energy')]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc='lower right')


# ── FigS4: GO enrichment (single panel) ──
fig, ax = plt.subplots(figsize=(6.5, 5.5))
barh_terms(ax, go_data, 'GO Biological Process (top 18 of 52)')
cat_legend(ax)
fig.tight_layout()
for fmt in ('png', 'pdf'):
    fig.savefig(f'{OUT}/FigS4_SMVT_GO_enrichment.{fmt}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('saved FigS4')

# ── FigS5: KEGG (left) + Reactome (right) ──
fig, (ax_k, ax_r) = plt.subplots(1, 2, figsize=(11, 6))
barh_terms(ax_k, kegg_data, 'KEGG pathways (13)')
barh_terms(ax_r, reactome_data[:15], 'Reactome pathways (top 15 of 17)')
cat_legend(ax_k)
cat_legend(ax_r)
fig.tight_layout()
for fmt in ('png', 'pdf'):
    fig.savefig(f'{OUT}/FigS5_SMVT_KEGG_Reactome.{fmt}', dpi=300, bbox_inches='tight')
plt.close(fig)
print('saved FigS5')
