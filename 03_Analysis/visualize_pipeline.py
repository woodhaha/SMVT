"""
SMVT Virtual Screening Pipeline — Comprehensive Visualization
Shows: funnel, multi-round strategy, pharmacophore model, top hits
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, Polygon
import numpy as np
import os

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/figures", exist_ok=True)

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9,
                     'axes.titlesize': 12, 'figure.dpi': 200,
                     'savefig.dpi': 300, 'savefig.bbox': 'tight'})

fig, axes = plt.subplots(2, 2, figsize=(22, 16))
(ax_funnel, ax_family, ax_pharm, ax_hits) = axes.flatten()

# ═══════════════════════════════════════════════════════
# PANEL A: SCREENING FUNNEL
# ═══════════════════════════════════════════════════════
ax_funnel.set_xlim(0, 10); ax_funnel.set_ylim(0, 10)
ax_funnel.axis('off')
ax_funnel.set_title('A. Screening Funnel', fontweight='bold', loc='left', fontsize=14)

stages = [
    (5.0, 9.2, 8.0, 'ZINC Database\n~230M compounds', '#BDC3C7'),
    (5.0, 7.8, 7.0, 'FDA / ChEMBL Filter\n~3,300 approved drugs', '#95A5A6'),
    (5.0, 6.5, 6.0, 'Pharmacophore ML Pre-screen\nECFP4 RF (AUC=0.888)', '#7F8C8D'),
    (5.0, 5.2, 5.0, 'Diversity Selection\n264 Murcko scaffolds → 500', '#95A5A6'),
    (5.0, 3.9, 4.0, 'AutoDock Vina\n4 rounds, exh=8-32', '#3498DB'),
    (5.0, 2.6, 3.6, 'Z-score Normalization\nbatch effect correction', '#2980B9'),
    (5.0, 1.3, 2.8, '8 Elite Hits (ΔG < −8.0)\n+ Barbiturate class confirmation', '#E74C3C'),
]

for cx, cy, w, label, color in stages:
    rect = FancyBboxPatch((cx - w/2, cy - 0.45), w, 0.85,
                           boxstyle="round,pad=0.15", facecolor=color,
                           edgecolor='white', linewidth=2, alpha=0.85)
    ax_funnel.add_patch(rect)
    # shadow effect
    if 'Elite' in label:
        ax_funnel.text(cx, cy, label, ha='center', va='center', fontsize=9,
                       fontweight='bold', color='white')
    else:
        ax_funnel.text(cx, cy, label, ha='center', va='center', fontsize=8,
                       color='white')

# Connect stages with tapered lines
for i in range(len(stages)-1):
    x1, y1, w1, _, _ = stages[i]
    x2, y2, w2, _, _ = stages[i+1]
    # Draw trapezoid-like connector
    points = [(x1-w1/2+0.3, y1-0.5), (x1+w1/2-0.3, y1-0.5),
              (x2+w2/2-0.3, y2+0.5), (x2-w2/2+0.3, y2+0.5)]
    ax_funnel.fill(points, color='#ECF0F1', alpha=0.5, edgecolor='#BDC3C7', linewidth=0.5)

# Annotations
ax_funnel.annotate('230M → 3,300', xy=(9.3, 8.7), fontsize=7, color='#7F8C8D')
ax_funnel.annotate('ML filter', xy=(8.8, 6.2), fontsize=7, color='#7F8C8D')
ax_funnel.annotate('500 docked', xy=(8.5, 4.7), fontsize=7, color='#7F8C8D')
ax_funnel.annotate('702 total', xy=(8.0, 3.4), fontsize=7, color='#7F8C8D')
ax_funnel.annotate('ΔG < −8.0', xy=(8.3, 1.8), fontsize=7, color='#E74C3C', fontweight='bold')


# ═══════════════════════════════════════════════════════
# PANEL B: CHEMICAL FAMILY BREAKDOWN
# ═══════════════════════════════════════════════════════
ax_family.set_title('B. Chemical Family Hit Rates', fontweight='bold', loc='left', fontsize=14)

families = ['Barbiturate', 'NSAID', 'Opioid', 'Vitamin/Cofactor', 'Sulfonamide', 'Other']
n_tested = [8, 12, 5, 8, 6, 663]
n_hits = [8, 6, 2, 2, 1, 31]  # ΔG < -7.0
hit_rate = [100, 50, 40, 25, 17, 4.7]

colors_fam = ['#E74C3C', '#E67E22', '#9B59B6', '#2ECC71', '#3498DB', '#95A5A6']
x = np.arange(len(families))
width = 0.55

bars = ax_family.bar(x, hit_rate, width, color=colors_fam, edgecolor='white', linewidth=1.5, alpha=0.9)

# Add count labels
for i, (rate, tested, hits) in enumerate(zip(hit_rate, n_tested, n_hits)):
    ax_family.text(i, rate + 1.5, f'{hits}/{tested}\n({rate}%)',
                   ha='center', fontsize=9, fontweight='bold')

ax_family.set_xticks(x)
ax_family.set_xticklabels(families, rotation=25, ha='right', fontsize=10)
ax_family.set_ylabel('Hit Rate (%)', fontsize=11)
ax_family.set_ylim(0, 115)
ax_family.axhline(y=50, color='gray', linestyle=':', alpha=0.5)

# Highlight best
best_bar = bars[0]
best_bar.set_edgecolor('#C0392B'); best_bar.set_linewidth(3)

# Star marker for 100%
ax_family.annotate('★ 100%', xy=(0, 102), ha='center', fontsize=11,
                   color='#C0392B', fontweight='bold')


# ═══════════════════════════════════════════════════════
# PANEL C: PHARMACOPHORE MODEL
# ═══════════════════════════════════════════════════════
ax_pharm.set_xlim(0, 10); ax_pharm.set_ylim(0, 10)
ax_pharm.axis('off')
ax_pharm.set_title('C. Pharmacophore Model & SAR Rules', fontweight='bold', loc='left', fontsize=14)

# Biotin reference structure (simplified)
ax_pharm.text(5, 9.5, 'Biotin (Natural Substrate) — Reference', ha='center', fontsize=10,
              fontweight='bold', color='#2C3E50')

# Draw simplified biotin ureido ring + carboxyl + thiophene
# Ureido ring
circle1 = plt.Circle((4.2, 8.05), 0.35, facecolor='#3498DB', alpha=0.4, edgecolor='#2980B9', linewidth=2)
circle2 = plt.Circle((5.8, 8.05), 0.35, facecolor='#3498DB', alpha=0.4, edgecolor='#2980B9', linewidth=2)
ax_pharm.add_patch(circle1); ax_pharm.add_patch(circle2)
ax_pharm.text(5, 7.8, 'Ureido Ring\n(HBD+HBA)', ha='center', fontsize=7, color='#2980B9')

# Carboxyl
circle3 = plt.Circle((2.5, 7.3), 0.3, facecolor='#E74C3C', alpha=0.35, edgecolor='#C0392B', linewidth=2)
ax_pharm.add_patch(circle3)
ax_pharm.text(2.5, 6.9, 'COOH', ha='center', fontsize=7, color='#C0392B')

# Thiophene
circle4 = plt.Circle((3.3, 8.6), 0.25, facecolor='#F39C12', alpha=0.35, edgecolor='#E67E22', linewidth=2)
ax_pharm.add_patch(circle4)
ax_pharm.text(3.3, 8.95, 'Thioether', ha='center', fontsize=6, color='#E67E22')

# Connect
ax_pharm.plot([4.5, 5.5], [8.05, 8.05], 'k-', linewidth=2)
ax_pharm.plot([2.8, 3.85], [7.4, 7.8], 'k-', linewidth=1.5)
ax_pharm.plot([3.5, 3.85], [8.5, 8.3], 'k-', linewidth=1)

# ===== Barbiturate =====
ax_pharm.text(5, 6.2, 'Barbiturate (Best Hit Class) — 100% Hit Rate', ha='center', fontsize=10,
              fontweight='bold', color='#E74C3C')

# Barbituric acid core
circle5 = plt.Circle((4.2, 4.8), 0.32, facecolor='#E74C3C', alpha=0.4, edgecolor='#C0392B', linewidth=2.5)
circle6 = plt.Circle((5.8, 4.8), 0.32, facecolor='#E74C3C', alpha=0.4, edgecolor='#C0392B', linewidth=2.5)
ax_pharm.add_patch(circle5); ax_pharm.add_patch(circle6)
ax_pharm.text(5, 4.5, 'Malonylurea\n(Mimics Ureido)', ha='center', fontsize=7, color='#C0392B', fontweight='bold')

# Phenyl ring
circle7 = plt.Circle((2.8, 3.9), 0.4, facecolor='#9B59B6', alpha=0.35, edgecolor='#8E44AD', linewidth=2)
ax_pharm.add_patch(circle7)
ax_pharm.text(2.8, 3.5, 'Aromatic\nRing', ha='center', fontsize=7, color='#8E44AD')

# Carbonyls
for xo in [3.4, 6.6]:
    ax_pharm.scatter([xo], [5.3], s=80, marker='o', facecolor='#F39C12', edgecolor='#E67E22', linewidth=1.5)
ax_pharm.text(5, 5.6, 'C=O (HBA)', ha='center', fontsize=7, color='#E67E22')

# Connections
ax_pharm.plot([4.5, 5.5], [4.8, 4.8], 'k-', linewidth=2.5)
ax_pharm.plot([3.2, 3.9], [3.9, 4.55], 'k-', linewidth=1.5)

# ===== SAR Rules Box =====
rules_box = FancyBboxPatch((0.3, 0.3), 9.4, 2.8, boxstyle="round,pad=0.3",
                            facecolor='#2C3E50', edgecolor='#1A252F', linewidth=2, alpha=0.92)
ax_pharm.add_patch(rules_box)

rules_text = [
    ('Rule 1:', '#E74C3C', ' Cyclic ureide/carboxamide core is the key pharmacophore'),
    ('Rule 2:', '#F39C12', ' Carboxyl group is NOT essential (barbiturates lack it entirely)'),
    ('Rule 3:', '#3498DB', ' Aromatic ring enhances affinity (hydrophobic contact)'),
    ('Rule 4:', '#2ECC71', ' H-bond acceptors (C=O) ⩾2 for optimal binding'),
    ('Rule 5:', '#9B59B6', ' Halogen substituents (Cl) tolerated, add specificity'),
]
for i, (label, color, desc) in enumerate(rules_text):
    y = 2.8 - i * 0.48
    ax_pharm.text(0.6, y, label, fontsize=8.5, fontweight='bold', color=color)
    ax_pharm.text(2.1, y, desc, fontsize=8.5, color='white')

ax_pharm.text(0.6, 0.5, '→ Expands chemical space beyond carboxylic acid mimics',
              fontsize=8, color='#BDC3C7', fontstyle='italic')


# ═══════════════════════════════════════════════════════
# PANEL D: TOP HITS RANKING
# ═══════════════════════════════════════════════════════
ax_hits.set_title('D. Top 10 Hits & Controls', fontweight='bold', loc='left', fontsize=14)

top_hits = [
    ('Hydromorphone', -8.58, 'Opioid', '#9B59B6'),
    ('Furosemide', -8.36, 'Sulfonamide', '#3498DB'),
    ('Naftazone', -8.34, 'Naphthoquinone', '#1ABC9C'),
    ('Phenobarbital', -8.30, 'Barbiturate', '#E74C3C'),
    ('Pentobarbital', -8.18, 'Barbiturate', '#E74C3C'),
    ('Diclofenac', -8.07, 'NSAID', '#E67E22'),
    ('Carprofen', -8.04, 'NSAID', '#E67E22'),
    ('Butabarbital', -8.02, 'Barbiturate', '#E74C3C'),
    ('Toloxatone', -7.92, 'MAO Inhibitor', '#95A5A6'),
    ('Lenalidomide', -7.85, 'Imide', '#95A5A6'),
]

# Controls
controls = [
    ('Biotin (天然底物)', -6.76, '#2ECC71'),
    ('Gabapentin (FDA药)', -6.63, '#27AE60'),
    ('Riboflavin (阴性)', -0.01, '#BDC3C7'),
]

names = [h[0] for h in top_hits] + [c[0] for c in controls]
values = [-h[1] for h in top_hits] + [-c[1] for c in controls]
colors_hits = [h[3] for h in top_hits]
colors_ctrl = [c[2] for c in controls]
colors = colors_hits + colors_ctrl

y_pos = range(len(names))
bars = ax_hits.barh(y_pos, values, color=colors, edgecolor='white', linewidth=1.2, height=0.7, alpha=0.9)

# Hit threshold
ax_hits.axvline(x=7.0, color='#E74C3C', linestyle='--', linewidth=1.5, alpha=0.6)
ax_hits.text(7.05, len(names)-0.3, '← Hit threshold (ΔG < −7.0)', fontsize=8, color='#E74C3C', alpha=0.8)
ax_hits.axvline(x=8.0, color='#C0392B', linestyle=':', linewidth=1, alpha=0.4)
ax_hits.text(8.05, len(names)-0.9, 'Elite (ΔG < −8.0)', fontsize=7, color='#C0392B', alpha=0.7)

# Labels
ax_hits.set_yticks(y_pos)
ax_hits.set_yticklabels(names, fontsize=9, fontfamily='monospace')
ax_hits.invert_yaxis()
ax_hits.set_xlabel('Binding Affinity −ΔG (kcal/mol)', fontsize=11)

# Value labels
all_compounds = [(n, v, f, c) for n, v, f, c in top_hits] + [(c[0], c[1], '', c[2]) for c in controls]
for i, (name, val, family, color) in enumerate(all_compounds):
    label = f'{val:.2f}'
    ax_hits.text(-val + 0.1, i, label, va='center', fontsize=8, fontweight='bold')
    if family:
        ax_hits.text(-val + 0.1, i - 0.25, family, va='top', fontsize=6, color='#7F8C8D')

# Highlight MD targets
md_targets = ['Hydromorphone', 'Furosemide', 'Naftazone', 'Phenobarbital']
for i, name in enumerate(names):
    if name in md_targets:
        ax_hits.annotate('→MD', xy=(-values[i] + 0.3, i + 0.15),
                         fontsize=7, fontweight='bold', color='#C0392B', va='center')

# Separator
ax_hits.axhline(y=9.5, color='#BDC3C7', linewidth=1, linestyle='-', alpha=0.5)
ax_hits.text(0.3, 10.8, 'Controls:', fontsize=8, fontweight='bold', color='#27AE60')

# Legend
leg = [mpatches.Patch(color='#E74C3C', label='Barbiturate'),
       mpatches.Patch(color='#E67E22', label='NSAID'),
       mpatches.Patch(color='#9B59B6', label='Opioid'),
       mpatches.Patch(color='#3498DB', label='Sulfonamide'),
       mpatches.Patch(color='#2ECC71', label='Control')]
ax_hits.legend(handles=leg, loc='lower right', fontsize=7, ncol=2)

# ═══════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════
fig.suptitle('SMVT (SLC5A6) Virtual Screening Pipeline\nAutoDock Vina · Pharmacophore ML · 702 Compounds → 8 Elite Hits',
             fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig("03_Analysis/figures/Fig_Screening_Pipeline.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_Screening_Pipeline.pdf", facecolor='white')
plt.close()
print("✓ Fig_Screening_Pipeline saved")

# ═══════════════════════════════════════════════════════
# BONUS: Chemical space plot (logP vs MW colored by affinity)
# ═══════════════════════════════════════════════════════
try:
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    df = pd.read_csv("03_Analysis/docking/docking_expanded_results.csv")
    df = df[df['best_affinity'].notna()].copy()

    # Add SMILES from the pharmacophore script's dictionary
    smiles_map = {
        "Biotin": "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2",
        "Lipoic_Acid": "C1CSSC1CCCCC(=O)O",
        "Pantothenic_Acid": "CC(C)(CO)C(O)C(=O)NCCC(=O)O",
        "Gabapentin_enacarbil": "CC(C)(C(=O)O)OC(=O)NCC1(CCCCC1)CC(=O)O",
        "Desthiobiotin": "CC1C(NC(=O)N1)CCCCCC(=O)O",
        "Biotin_Methyl_Ester": "C1C2C(C(S1)CCCCC(=O)OC)NC(=O)N2",
        "Biotin_Sulfoxide": "C1C2C(C(S1(=O))CCCCC(=O)O)NC(=O)N2",
        "Biotin_Sulfone": "C1C2C(C(S1(=O)=O)CCCCC(=O)O)NC(=O)N2",
        "Norbiotin": "C1C2C(C(S1)CCCC(=O)O)NC(=O)N2",
        "Homobiotin": "C1C2C(C(S1)CCCCCC(=O)O)NC(=O)N2",
        "Indomethacin": "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=C(C=C2)OC)CC(=O)O",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Diclofenac": "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl",
        "Ketoprofen": "CC(C1=CC=CC(=C1)C(=O)C2=CC=CC=C2)C(=O)O",
        "Flurbiprofen": "CC1=CC(=CC=C1)C2=CC=C(C=C2)C(C)C(=O)O",
        "Phenylbutazone": "CCCCC1C(=O)N(N(C1=O)C2=CC=CC=C2)C3=CC=CC=C3",
        "Naproxen": "CC(C1=CC2=C(C=C1)C=CC(=C2)OC)C(=O)O",
        "Celecoxib": "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",
        "Piroxicam": "CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3",
        "Meloxicam": "CC1=CC=CS(=O)(=O)N1C(=O)NC2=C(C3=NC=C(S3)C)O2",
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Sulindac": "CC1=C(C2=C(C=CC(=C2)F)C(=C1)CC(=O)O)C3=CC=C(C=C3)S(=O)C",
        "Mefenamic_Acid": "CC1=CC=CC(=C1NC2=CC=CC=C2C(=O)O)C",
        "Phenobarbital": "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2",
        "Pentobarbital": "CCCC(C)C1(C(=O)NC(=O)NC1=O)CC",
        "Secobarbital": "CCCC(C)C1(C(=O)NC(=O)NC1=O)CC=C",
        "Butabarbital": "CCCC(C)C1(C(=O)NC(=O)NC1=O)CC",
        "Naftazone": "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N",
        "Furosemide": "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl",
        "Hydromorphone": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5",
        "Esketamine": "CNC1(C2=CC=CC=C2Cl)CCCCC1=O",
        "Toloxatone": "O=C1OCCN1c1ccc(OC)cc1C",
        "Lenalidomide": "O=C1CCC(N1C(=O)c1cccc2c1CNC2=O)N",
        "Carprofen": "Clc1cc2c(c(Cl)c1)nc1c2ccc(C(C)C(=O)O)c1",
        "Oxymorphone": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5",
    }

    mw_list, logp_list, aff_list, colors_list = [], [], [], []
    for _, row in df.iterrows():
        smi = smiles_map.get(row['name'])
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        aff = -row['best_affinity']
        mw_list.append(mw)
        logp_list.append(logp)
        aff_list.append(aff)
        if aff > 8:
            colors_list.append('#E74C3C')
        elif aff > 7:
            colors_list.append('#F39C12')
        elif aff > 6:
            colors_list.append('#3498DB')
        else:
            colors_list.append('#BDC3C7')

    fig2, ax2 = plt.subplots(figsize=(10, 7))
    scatter = ax2.scatter(mw_list, logp_list, c=colors_list, s=[max(40, a*20) for a in aff_list],
                          alpha=0.8, edgecolors='white', linewidth=0.5, zorder=5)

    # Highlight top 8
    for i in range(min(8, len(df))):
        row = df.iloc[i]
        smi = smiles_map.get(row['name'])
        if not smi: continue
        mol = Chem.MolFromSmiles(smi)
        if not mol: continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        ax2.annotate(row['name'].replace('_', '\n'), (mw, logp),
                     fontsize=6, fontweight='bold', ha='center',
                     color='#2C3E50')

    ax2.set_xlabel('Molecular Weight (Da)', fontsize=12)
    ax2.set_ylabel('logP', fontsize=12)
    ax2.set_title('E. Chemical Space: MW vs logP (colored by affinity)', fontweight='bold', loc='left', fontsize=13)

    leg = [mpatches.Patch(color='#E74C3C', label='Elite (ΔG < −8.0)'),
           mpatches.Patch(color='#F39C12', label='Hit (ΔG < −7.0)'),
           mpatches.Patch(color='#3498DB', label='Moderate (ΔG < −6.0)'),
           mpatches.Patch(color='#BDC3C7', label='Weak (ΔG > −6.0)')]
    ax2.legend(handles=leg, loc='upper right', fontsize=9)

    plt.tight_layout()
    fig2.savefig("03_Analysis/figures/Fig_Chemical_Space.png", facecolor='white')
    fig2.savefig("03_Analysis/figures/Fig_Chemical_Space.pdf", facecolor='white')
    plt.close()
    print("✓ Fig_Chemical_Space saved")
except Exception as e:
    print(f"Chemical space plot skipped: {e}")

print("\n✓ ALL DONE — 5 panel screening pipeline visualization")
