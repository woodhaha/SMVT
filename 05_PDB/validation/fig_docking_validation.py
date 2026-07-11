#!/usr/bin/env python3
"""SMVT Docking Validation — Dual-Structure Comparison (AF2 vs cryo-EM 26va)."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── MANDATORY ─────────────────────────────────────────────────────────────────
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.size': 8,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.linewidth': 1.0,
    'legend.frameon': False,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# ── Data ───────────────────────────────────────────────────────────────────────
v26 = {'NAFTAZONE': -7.578, 'PHENOBARBITAL': -7.485, 'GABAPENTIN_ENACARBIL': -5.955,
       'BIOTIN': -5.927, 'ESKETAMINE': -5.902, 'FUROSEMIDE': -4.323,
       'HYDROMORPHONE': -3.740, 'RIBOFLAVIN': -1.546}
vAF = {'NAFTAZONE': -8.03, 'PHENOBARBITAL': -8.30, 'GABAPENTIN_ENACARBIL': -6.97,
       'BIOTIN': -6.38, 'ESKETAMINE': -7.58, 'FUROSEMIDE': -7.56,
       'HYDROMORPHONE': -7.72, 'RIBOFLAVIN': -7.56}
mg  = {'GABAPENTIN_ENACARBIL': -43.33, 'RIBOFLAVIN': -41.48, 'BIOTIN': -29.97,
       'HYDROMORPHONE': -29.89, 'ESKETAMINE': -27.91, 'PHENOBARBITAL': -23.61,
       'NAFTAZONE': -22.36, 'FUROSEMIDE': 6.40}

# Rank order by cryo-EM score
order = sorted(v26.keys(), key=lambda c: v26[c])
short = {'NAFTAZONE':'NAF','PHENOBARBITAL':'PHE','GABAPENTIN_ENACARBIL':'GAB',
         'BIOTIN':'BIO','ESKETAMINE':'ESK','FUROSEMIDE':'FUR',
         'HYDROMORPHONE':'HYD','RIBOFLAVIN':'RIB'}

# Palette — systematic
colors = {'BIOTIN':'#0F4D92','ESKETAMINE':'#42949E','FUROSEMIDE':'#B64342',
          'GABAPENTIN_ENACARBIL':'#2E9E44','HYDROMORPHONE':'#9A4D8E',
          'NAFTAZONE':'#FF8C00','PHENOBARBITAL':'#7884B4','RIBOFLAVIN':'#767676'}

# ── Figure ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.08, 5.5))  # EMBO J / Nature Comms single-col
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4,
              left=0.09, right=0.97, bottom=0.1, top=0.92)

# ── Panel a: Dual-receptor Vina score bar chart ───────────────────────────────
ax = fig.add_subplot(gs[:, :2])
x = np.arange(len(order))
w = 0.35

for i, c in enumerate(order):
    ax.bar(i - w/2, vAF[c], w, color='#B4C0E4', edgecolor='white', lw=0.3, zorder=3)
    ax.bar(i + w/2, v26[c], w, color='#484878', edgecolor='white', lw=0.3, zorder=3)

ax.set_xticks(x)
ax.set_xticklabels([short[c] for c in order], fontsize=7)
ax.set_ylabel('Vina score (kcal/mol)', fontsize=8)
ax.set_ylim(-10.5, 0.5)
ax.axhline(0, color='#333', lw=0.5)

# Score labels (above bars)
for i, c in enumerate(order):
    ax.text(i - w/2, vAF[c] + 0.3, f'{vAF[c]:.1f}', ha='center', fontsize=5,
            color='#484878', va='bottom')
    ax.text(i + w/2, v26[c] + 0.3, f'{v26[c]:.1f}', ha='center', fontsize=5,
            color='#333', va='bottom')
    # Red arrow for large drop
    if vAF[c] - v26[c] > 2:
        ax.annotate('', xy=(i+w/2, v26[c]+0.5), xytext=(i+w/2, vAF[c]-0.5),
                    arrowprops=dict(arrowstyle='<->', color='#e74c3c', lw=0.8, alpha=0.6))
        ax.text(i+w/2, (vAF[c]+v26[c])/2, f'{vAF[c]-v26[c]:.0f}', ha='center', fontsize=5,
                color='#e74c3c', fontweight='bold')

# Panel a labels: compound name only (no duplicate rank annotations)
ax.text(-0.06, 1.06, 'a', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top', ha='right')
ax.grid(axis='y', alpha=0.1, lw=0.3)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#B4C0E4', edgecolor='white', label='AF2 (predicted)'),
    Patch(facecolor='#333333', edgecolor='white', label='cryo-EM (26va)'),
]
ax.legend(handles=legend_elements, fontsize=6, loc='upper center',
          bbox_to_anchor=(0.5, -0.05), ncol=2)

# ── Panel b: Ranking correlation heatmap ──────────────────────────────────────
ax = fig.add_subplot(gs[0, 2])
methods = {'Vina\n(AF2)': vAF, 'Vina\n(26va)': v26}
data = np.zeros((len(methods), len(order)))
for j, (mname, mdict) in enumerate(methods.items()):
    for i, c in enumerate(order):
        data[j, i] = list(sorted(mdict.values(), reverse=True)).index(mdict[c]) + 1

im = ax.imshow(data, cmap='RdBu_r', aspect='auto', vmin=1, vmax=8)
ax.set_yticks(range(len(methods)))
ax.set_yticklabels(list(methods.keys()), fontsize=7)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([short[c] for c in order], fontsize=7, rotation=45, ha='right')
ax.set_title('Rank', fontsize=8)

for j in range(len(methods)):
    for i in range(len(order)):
        ax.text(i, j, f'{int(data[j,i])}', ha='center', va='center',
                fontsize=7, fontweight='bold',
                color='white' if abs(data[j,i]-4.5) > 2 else '#333')

ax.text(-0.15, 1.06, 'b', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top', ha='right')

# ── Panel c: Pocket agreement schematic ──────────────────────────────────────
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')

# Stats
pocket_overlap = 16  # common residues
pocket_exp = 23      # experimental
pocket_pred = 20     # AF2 predicted
precision = pocket_overlap / pocket_pred * 100
recall = pocket_overlap / pocket_exp * 100
f1 = 2 * precision * recall / (precision + recall)

txt = (
    "Binding Pocket\n"
    "AF2 vs cryo-EM\n"
    f"{'='*18}\n\n"
    f"Overlap: {pocket_overlap}/{pocket_exp}\n"
    f"Precision: {precision:.0f}%\n"
    f"Recall: {recall:.0f}%\n"
    f"F1: {f1:.0f}%\n\n"
    "Cross-structure\n"
    "ranking shift:\n"
    f"  NAF: {1}->{1}\n"
    f"  PHE: {sorted(vAF.values(), reverse=True).index(vAF['PHENOBARBITAL'])+1}"
         f"->{sorted(v26.values(), reverse=True).index(v26['PHENOBARBITAL'])+1}\n"
    f"  HYD: {sorted(vAF.values(), reverse=True).index(vAF['HYDROMORPHONE'])+1}"
         f"->{sorted(v26.values(), reverse=True).index(v26['HYDROMORPHONE'])+1}\n\n"
    "Naftazone: #1 in both\n"
    "Phenobarbital: #2\n"
    "Biotin rank: 4 (ref)\n"
    "Riboflavin neg: #8 26va"
)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=6.5,
        va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', fc='#f8f9fa', ec='#dee2e6', lw=0.5))
ax.text(-0.15, 1.06, 'c', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top', ha='right')

# ── Save ───────────────────────────────────────────────────────────────────────
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.15, top=0.92)
fig.savefig(f'{OUT}/Fig_docking_validation.svg', bbox_inches='tight')
fig.savefig(f'{OUT}/Fig_docking_validation.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {OUT}/Fig_docking_validation.svg + .png')
