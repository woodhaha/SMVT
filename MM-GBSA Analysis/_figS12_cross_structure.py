"""FigS12 — Cross-structure docking validation: AF2 vs cryo-EM (26va)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
                     'figure.dpi': 300, 'savefig.dpi': 300, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial']})

OUT = r'D:\Researching\SMVT\04_Manuscript\figures'
COMPOUNDS = ['Biotin', 'Naftazone', 'Phenobarbital', 'Hydromorphone',
             'Esketamine', 'Furosemide', 'Gabapentin\nEnacarbil', 'Riboflavin']

# Docking scores: AF2 vs cryo-EM (from manuscript data)
af2_scores  = [-6.82, -8.03, -7.31, -7.72, -7.36, -7.56, -6.97, -7.56]
cryo_scores = [-6.37, -7.58, -7.49, -3.74, -5.90, -6.91, -5.73, -6.47]

af2_ranks  = [4, 1, 2, 2, 5, 3, 7, 4][:8]
cryo_ranks = [4, 1, 2, 7, 6, 5, 8, 3]

colors = ['#3498db','#e67e22','#f39c12','#9b59b6','#1abc9c','#e74c3c','#2ecc71','#95a5a6']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

# Panel a: Paired bar chart
x = np.arange(len(COMPOUNDS))
w = 0.35
bars1 = ax1.bar(x - w/2, af2_scores, w, label='AlphaFold2 (AF-Q9Y289-F1)', color='#5D7DB3', edgecolor='white', linewidth=0.5)
bars2 = ax1.bar(x + w/2, cryo_scores, w, label='Cryo-EM (PDB: 26va)', color='#E8834A', edgecolor='white', linewidth=0.5)

for bar in bars1:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, h + 0.1, f'{h:.1f}', ha='center', va='bottom', fontsize=7, fontweight='bold')
for bar in bars2:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, h + 0.1, f'{h:.1f}', ha='center', va='bottom', fontsize=7, fontweight='bold')

ax1.set_xticks(x)
ax1.set_xticklabels(COMPOUNDS, rotation=30, ha='right', fontsize=8)
ax1.set_ylabel('Vina Score (kcal/mol)', fontsize=11)
ax1.set_title('Cross-Structure Docking Score Comparison', fontsize=12, fontweight='bold')
ax1.legend(fontsize=8, framealpha=0.9)
ax1.axhline(y=-6.76, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax1.text(-0.3, -6.6, 'Biotin baseline (−6.76)', fontsize=7, color='gray', fontstyle='italic')
ax1.invert_yaxis()
ax1.set_ylim(-9, -2)

# Inset: score difference
ax_inset = ax1.inset_axes([0.55, 0.55, 0.4, 0.35])
diffs = [c - a for a, c in zip(af2_scores, cryo_scores)]
ax_inset.barh(range(len(COMPOUNDS)), diffs, color=colors, height=0.6)
ax_inset.set_yticks(range(len(COMPOUNDS)))
ax_inset.set_yticklabels([c.replace('\n', ' ') for c in COMPOUNDS], fontsize=6)
ax_inset.set_xlabel('Δ Score (Cryo-EM − AF2)', fontsize=8)
ax_inset.axvline(x=0, color='black', linewidth=0.5)
ax_inset.tick_params(labelsize=7)

# Panel b: Pocket residue comparison
pocket_af2  = {'PHE79', 'PHE80', 'TYR84', 'GLU91', 'PHE98', 'LEU99', 'ASN132', 'ILE133',
               'LEU135', 'ALA138', 'ALA139', 'SER213', 'GLY214', 'ILE217', 'THR265',
               'MET267', 'TYR271', 'GLN276', 'GLN277', 'ALA299', 'ALA300', 'GLN301', 'ILE501'}
pocket_cryo = {'PHE79', 'PHE80', 'TYR84', 'GLU91', 'PHE98', 'LEU99', 'ASN132', 'LEU135',
               'ALA138', 'ALA139', 'GLY214', 'ILE217', 'THR265', 'MET267', 'TYR271',
               'GLN276', 'GLN277', 'ALA299', 'GLN300', 'GLN301', 'ILE501'}

shared = pocket_af2 & pocket_cryo
only_af2 = pocket_af2 - pocket_cryo
only_cryo = pocket_cryo - pocket_af2

precision = len(shared) / (len(shared) + len(only_af2))
recall = len(shared) / (len(shared) + len(only_cryo))
f1 = 2 * precision * recall / (precision + recall)

# Draw Venn-like representation
from matplotlib.patches import Circle, FancyBboxPatch
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 6)
ax2.set_aspect('equal')
ax2.axis('off')

c1 = Circle((3.3, 3), 1.8, facecolor='#5D7DB3', alpha=0.3, edgecolor='#5D7DB3', linewidth=2)
c2 = Circle((6.7, 3), 1.8, facecolor='#E8834A', alpha=0.3, edgecolor='#E8834A', linewidth=2)
ax2.add_patch(c1)
ax2.add_patch(c2)

ax2.text(3.3, 3, f'Shared\n{len(shared)} residues', ha='center', va='center', fontsize=10, fontweight='bold')
ax2.text(0.8, 4.5, f'AF2-only\n{len(only_af2)}', ha='center', fontsize=8, color='#5D7DB3')
ax2.text(9.2, 4.5, f'Cryo-EM-only\n{len(only_cryo)}', ha='center', fontsize=8, color='#E8834A')

# Metrics table
ax2.text(5, 0.4, f'Precision: {precision:.0%}   Recall: {recall:.0%}   F1: {f1:.0%}',
         ha='center', va='center', fontsize=10, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))

# Shared residues text
shared_str = ', '.join(sorted(shared))
ax2.text(5, 5.6, f'Shared residues: {shared_str}', ha='center', va='top',
         fontsize=6, fontstyle='italic', color='gray', wrap=True)

ax2.set_title('Binding Pocket Residue Overlap', fontsize=12, fontweight='bold')

fig.savefig(f'{OUT}/FigS12_cross_structure_docking.png', bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}/FigS12_cross_structure_docking.pdf', bbox_inches='tight', facecolor='white')
print(f'✅ FigS12 saved')
print(f'   Precision={precision:.0%}, Recall={recall:.0%}, F1={f1:.0%}')
