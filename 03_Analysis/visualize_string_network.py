#!/usr/bin/env python3
"""
SLC5A6 (SMVT) STRING Network — Nature Style Visualization
===========================================================
Hub-and-spoke PPI network + evidence breakdown + functional modules
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Arc, FancyArrowPatch, Circle, Wedge
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
C_CENTER  = '#D62728'  # SLC5A6
C_ANCHOR  = '#E76F51'  # PDZD11
C_BIOTIN  = '#2CA02C'  # biotin module
C_SLC     = '#1F77B4'  # SLC cluster
C_OTHER   = '#7F7F7F'
C_BG      = '#FAFAFA'

# ═══ Data ═══
partners = [
    ('PDZD11',   0.969, 'Anchor',    'PDZ scaffold, apical membrane'),
    ('HLCS',     0.635, 'Biotin',    'Holocarboxylase synthetase'),
    ('SLC5A7',   0.655, 'SLC',       'Choline transporter'),
    ('SLC22A12', 0.631, 'SLC',       'URAT1, urate transporter'),
    ('SLC26A4',  0.578, 'SLC',       'Pendrin, anion exchanger'),
    ('BTD',      0.555, 'Biotin',    'Biotinidase, biotin recycling'),
    ('SLC5A3',   0.511, 'SLC',       'Myo-inositol transporter'),
    ('SLC23A1',  0.493, 'SLC',       'Vitamin C transporter'),
    ('DPH2',     0.471, 'Other',     'Diphthamide biosynthesis'),
    ('SLC19A2',  0.469, 'SLC',       'Thiamine transporter B1'),
]

# Evidence sub-scores (from STRING JSON)
evidence = {
    'PDZD11':   {'exp':0.292, 'db':0.900, 'txt':0.609, 'coe':0.000},
    'SLC5A7':   {'exp':0.000, 'db':0.000, 'txt':0.581, 'coe':0.000, 'phy':0.212},
    'HLCS':     {'exp':0.000, 'db':0.000, 'txt':0.625, 'coe':0.067, 'ngh':0.042},
    'SLC22A12': {'exp':0.045, 'db':0.000, 'txt':0.629, 'coe':0.042},
    'SLC26A4':  {'exp':0.000, 'db':0.000, 'txt':0.574, 'coe':0.050},
    'BTD':      {'exp':0.000, 'db':0.000, 'txt':0.551, 'coe':0.049},
    'SLC5A3':   {'exp':0.137, 'db':0.000, 'txt':0.316, 'coe':0.054, 'phy':0.228},
    'SLC23A1':  {'exp':0.000, 'db':0.000, 'txt':0.454, 'coe':0.110, 'fus':0.003},
    'DPH2':     {'exp':0.000, 'db':0.000, 'txt':0.000, 'coe':0.471},
    'SLC19A2':  {'exp':0.000, 'db':0.000, 'txt':0.459, 'coe':0.059},
}

module_color = {'Anchor': C_ANCHOR, 'Biotin': C_BIOTIN, 'SLC': C_SLC, 'Other': C_OTHER}

# ═══════════════════════════════════════════
# FIGURE 1: Multi-panel Network Composite
# ═══════════════════════════════════════════

fig = plt.figure(figsize=(9.0, 7.5))

# ── Panel a: Hub-and-Spoke Network Diagram ──
ax_a = fig.add_axes([0.02, 0.25, 0.52, 0.72])
ax_a.set_xlim(-3.5, 3.5); ax_a.set_ylim(-3.5, 3.5)
ax_a.axis('off')
ax_a.set_aspect('equal')

# Central node — SLC5A6
center = Circle((0, 0), 0.55, facecolor=C_CENTER, edgecolor='white', lw=2, zorder=10)
ax_a.add_patch(center)
ax_a.text(0, 0, 'SLC5A6\n(SMVT)', ha='center', va='center', fontsize=7, color='white', fontweight='bold')

# Position partners in a circle
n = len(partners)
angles = np.linspace(np.pi/2, np.pi/2 + 2*np.pi, n, endpoint=False)
radii = 2.8

for i, (name, score, module, desc) in enumerate(partners):
    x = radii * np.cos(angles[i])
    y = radii * np.sin(angles[i])

    # Node size proportional to score
    r = 0.22 + 0.28 * (score - 0.45)
    node = Circle((x, y), r, facecolor=module_color[module], edgecolor='white', lw=1.5, zorder=5, alpha=0.92)
    ax_a.add_patch(node)

    # Edge thickness proportional to score
    lw = 0.5 + 8 * (score - 0.45)
    ax_a.plot([0, x], [0, y], color='#999', lw=lw, alpha=0.5, zorder=1, solid_capstyle='round')

    # Label
    offset_x = 0.35 * np.cos(angles[i])
    offset_y = 0.35 * np.sin(angles[i])
    label_x = x + offset_x
    label_y = y + offset_y

    fs = 5.5 if score < 0.6 else 6.5
    fw = 'bold' if score > 0.9 else 'normal'
    ax_a.text(x, y, name, ha='center', va='center', fontsize=fs, color='white', fontweight=fw)

# Ring annotations for modules
for angle_start, angle_end, color, label, label_angle in [
    (angles[1]-0.15, angles[3]+0.15, C_SLC, 'SLC Cluster', angles[2]),
    (angles[4]-0.15, angles[5]+0.15, C_BIOTIN, '', 0),
    (angles[0]-0.15, angles[0]+0.15, C_ANCHOR, '', 0),
]:
    arc = Arc((0, 0), radii*2 + 0.6, radii*2 + 0.6, angle=0,
              theta1=np.degrees(angle_start), theta2=np.degrees(angle_end),
              color=color, lw=2, alpha=0.4, zorder=0)
    ax_a.add_patch(arc)

# Legend
legend_elements = [
    Line2D([0],[0], marker='o', c='w', mfc=C_CENTER, markersize=10, label='SLC5A6 (SMVT)'),
    Line2D([0],[0], marker='o', c='w', mfc=C_ANCHOR, markersize=7, label='Anchor (PDZD11, 0.969)'),
    Line2D([0],[0], marker='o', c='w', mfc=C_BIOTIN, markersize=7, label='Biotin module'),
    Line2D([0],[0], marker='o', c='w', mfc=C_SLC, markersize=7, label='SLC transporter cluster'),
]
ax_a.legend(handles=legend_elements, frameon=False, fontsize=5.5, loc='lower left',
            bbox_to_anchor=(-0.02, -0.08))

ax_a.set_title('STRING PPI Network (score >= 0.4)', fontsize=8, fontweight='bold', loc='left')
ax_a.text(0.02, 0.97, 'a', transform=ax_a.transAxes, fontsize=12, fontweight='bold')

# ── Panel b: Evidence Breakdown ──
ax_b = fig.add_axes([0.58, 0.55, 0.40, 0.42])

gene_names = [p[0] for p in partners]
scores = [p[1] for p in partners]
modules = [p[2] for p in partners]
colors_b = [module_color[m] for m in modules]

# Horizontal bar chart
y_pos = range(len(gene_names)-1, -1, -1)
bars = ax_b.barh(y_pos, scores, height=0.6, color=colors_b, alpha=0.85, edgecolor='white', lw=0.3)

# Add evidence sub-score annotations
for i, (name, score) in enumerate(zip(gene_names, scores)):
    if name in evidence:
        ev = evidence[name]
        # Draw stacked sub-bars as small dots
        sub_labels = []
        if ev.get('exp', 0) > 0: sub_labels.append('E')
        if ev.get('db', 0) > 0: sub_labels.append('D')
        if ev.get('txt', 0) > 0: sub_labels.append('T')
        if ev.get('coe', 0) > 0: sub_labels.append('C')
        if ev.get('phy', 0) > 0: sub_labels.append('P')
        sub_text = ' '.join(sub_labels)
        ax_b.text(score + 0.03, y_pos[i], sub_text, fontsize=4.5, va='center', color='#666')

ax_b.set_yticks(list(y_pos))
ax_b.set_yticklabels(gene_names, fontsize=6, fontweight='bold')
ax_b.set_xlabel('Combined STRING score', fontsize=7)
ax_b.set_xlim(0, 1.15)
ax_b.axvline(x=0.7, color='grey', lw=0.3, ls='--', alpha=0.4)
ax_b.text(0.71, 9.3, 'high\nconf.', fontsize=4.5, color='grey')

ax_b.set_title('Interaction Scores & Evidence', fontsize=8, fontweight='bold', loc='left')
ax_b.text(-0.05, 1.02, 'b', transform=ax_b.transAxes, fontsize=12, fontweight='bold')

# ── Panel c: Evidence Legend + Module Description ──
ax_c = fig.add_axes([0.58, 0.28, 0.40, 0.22])
ax_c.axis('off')

descriptions = [
    'E  Experimental (biochemical, Co-IP, etc.)',
    'D  Database (curated pathway/PPI databases)',
    'T  Text mining (co-mentioned in abstracts)',
    'C  Co-expression (RNA correlation across tissues)',
    'P  Phylogenetic co-occurrence',
]
for i, d in enumerate(descriptions):
    ax_c.text(0.02, 0.85 - i*0.16, d, fontsize=5, transform=ax_c.transAxes, color='#555')

# HLCS regulatory annotation
ax_c.text(0.02, 0.05, 'HLCS feedback: HLCS biotinylates H4K12 at SLC5A6 promoter\n'
         '--> silences SLC5A6 when biotin is abundant\n'
         '--> loss of this loop may explain cancer overexpression',
         fontsize=4.8, transform=ax_c.transAxes, color=C_BIOTIN, fontstyle='italic',
         bbox=dict(boxstyle='round,pad=0.3', fc='#F0FFF0', ec='none', alpha=0.7))

ax_c.set_title('Evidence Codes & Key Regulatory Circuit', fontsize=7.5, fontweight='bold', loc='left')
ax_c.text(-0.05, 1.02, 'c', transform=ax_c.transAxes, fontsize=12, fontweight='bold')

# ── Panel d: PDZ Interaction Detail ──
ax_d = fig.add_axes([0.02, 0.03, 0.52, 0.18])
ax_d.axis('off')

# Draw PDZ binding schematic
# SLC5A6 C-terminus motif
ax_d.text(0.08, 0.65, 'SLC5A6 C-terminus', fontsize=6.5, fontweight='bold', color=C_CENTER)
# Draw the motif
motif = '...S-E-R-T-L-COOH (Class I PDZ motif: x-S/T-x-phi)'
ax_d.text(0.08, 0.42, motif, fontsize=6, color='#333', fontfamily='monospace')

# PDZD11
ax_d.text(0.62, 0.65, 'PDZD11', fontsize=6.5, fontweight='bold', color=C_ANCHOR)
ax_d.text(0.62, 0.42, 'PDZ domain protein\nApical membrane scaffold\nBrush border localization', fontsize=5.5, color='#555')

# Arrow between them
ax_d.annotate('', xy=(0.58, 0.55), xytext=(0.48, 0.55),
              transform=ax_d.transAxes,
              arrowprops=dict(arrowstyle='<->', color='#333', lw=1.5))
ax_d.text(0.53, 0.62, 'PDZ-ligand\nbinding', fontsize=5, ha='center', color='#666', fontstyle='italic')

# Box around them
rect_l = FancyBboxPatch((0.04, 0.15), 0.40, 0.68,
                         boxstyle="round,pad=0.06", facecolor='#FFF5F5', edgecolor=C_CENTER, lw=0.8,
                         transform=ax_d.transAxes)
ax_d.add_patch(rect_l)
rect_r = FancyBboxPatch((0.56, 0.15), 0.40, 0.68,
                         boxstyle="round,pad=0.06", facecolor='#FFF0E8', edgecolor=C_ANCHOR, lw=0.8,
                         transform=ax_d.transAxes)
ax_d.add_patch(rect_r)

ax_d.text(0.04, 0.88, 'd', transform=ax_d.transAxes, fontsize=12, fontweight='bold')
ax_d.set_title('Strongest Interaction (score 0.969)', fontsize=7.5, fontweight='bold', loc='left')

# ── Global title ──
fig.suptitle('SLC5A6 (SMVT) Protein-Protein Interaction Network — STRING v12',
             fontsize=10, fontweight='bold', x=0.02, y=0.995, ha='left')

fig.text(0.02, 0.01, 'Data: STRING v12 (https://string-db.org/), Homo sapiens, combined score >= 0.4',
         fontsize=5, color='grey', style='italic')

out = OUTDIR / 'Fig5_SMVT_STRING_network.png'
fig.savefig(out, dpi=600, facecolor='white')
fig.savefig(OUTDIR / 'Fig5_SMVT_STRING_network.pdf', facecolor='white')
print(f'OK: {out.name}')
plt.close()

# ═══════════════════════════════════════════
# FIGURE 2: Simple evidence heatmap
# ═══════════════════════════════════════════

fig2, ax = plt.subplots(figsize=(6.0, 3.5))
plt.subplots_adjust(left=0.2, right=0.95, top=0.88, bottom=0.08)

# Build evidence matrix
ev_types = ['exp', 'db', 'txt', 'coe', 'phy', 'ngh', 'fus']
ev_labels = ['Experimental', 'Database', 'Text Mining', 'Co-expression', 'Phylo. co-occ.', 'Neighborhood', 'Gene Fusion']
ev_colors = ['#D62728', '#9467BD', '#1F77B4', '#2CA02C', '#FF7F0E', '#8C564B', '#E377C2']

matrix = np.zeros((len(partners), len(ev_types)))
for i, (name, score, module, desc) in enumerate(partners):
    if name in evidence:
        for j, evt in enumerate(ev_types):
            matrix[i, j] = evidence[name].get(evt, 0)

# Heatmap
im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)

# Labels
ax.set_yticks(range(len(partners)))
ax.set_yticklabels([p[0] for p in partners], fontsize=6.5, fontweight='bold')
ax.set_xticks(range(len(ev_types)))
ax.set_xticklabels(ev_labels, fontsize=6, rotation=25, ha='right')

# Annotate cells
for i in range(len(partners)):
    for j in range(len(ev_types)):
        if matrix[i, j] > 0:
            val = matrix[i, j]
            color = 'white' if val > 0.4 else '#333'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=4.5, color=color,
                    fontweight='bold' if val > 0.5 else 'normal')

# Add combined score column as text
for i, (name, score, module, desc) in enumerate(partners):
    color = module_color[module]
    ax.text(len(ev_types) + 0.5, i, f'{score:.3f}', fontsize=5.5, va='center', ha='center',
            fontweight='bold', color=color,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=color, lw=0.5, alpha=0.8))

# Additional column header
ax.text(len(ev_types) + 0.5, -1.2, 'Score', fontsize=5.5, ha='center', fontweight='bold', rotation=25)

fig2.suptitle('SLC5A6 STRING Interaction Evidence Matrix', fontsize=9, fontweight='bold', x=0.2, ha='left')

out2 = OUTDIR / 'Fig6_SMVT_STRING_evidence_heatmap.png'
fig2.savefig(out2, dpi=600, facecolor='white')
fig2.savefig(OUTDIR / 'Fig6_SMVT_STRING_evidence_heatmap.pdf', facecolor='white')
print(f'OK: {out2.name}')
plt.close()

print('\nDone: STRING network visualization complete.')
