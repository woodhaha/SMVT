#!/usr/bin/env python3
"""SMVT paper figures — large, readable fonts, Nature-clean style."""
import os, warnings, json, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial'],
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 14,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'xtick.direction': 'out', 'ytick.direction': 'out',
})

DIR = r"D:\Researching\SMVT\SMVT-MD-Analysis"
OUT = r"D:\Researching\SMVT\Paper_Figures"
os.makedirs(OUT, exist_ok=True); os.chdir(DIR)

C8 = ['#2166AC','#D6604D','#4DAF4A','#FF7F00','#984EA3','#A65628','#F781BF','#999999']
C = ['NAFTAZONE','BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL','HYDROMORPHONE','PHENOBARBITAL','RIBOFLAVIN']
S = ['Naftazone','Biotin','Esketamine','Furosemide','Gab.Enacarbil','Hydromorphone','Phenobarbital','Riboflavin']

with open(r'D:\Researching\SMVT\MM-GBSA Analysis\energy_decomposition.json') as f:
    E = json.load(f)

rmsd = {}
for c in C:
    rmsd[c] = pd.read_csv(f'{c}_100ns_full.csv')
    rmsd[c]['t'] = np.arange(len(rmsd[c])) * 100.0 / 1000

pca = {}
for c in C:
    try: pca[c] = np.load(f'{c}_pca_full.npz')
    except: pass

def ax_style(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ══════════════════════════════════
# Figure 1 — Virtual Screening
# ══════════════════════════════════
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('Figure 1 | Virtual Screening Results', fontsize=16, fontweight='bold')
x = range(8); w = 0.6
V = [-8.34,-6.76,-7.58,-8.36,-6.63,-8.58,-8.30,-0.01]

b1 = a1.bar(x, V, w, color=C8, ec='white', lw=0.4)
for b,v in zip(b1,V):
    a1.text(b.get_x()+w/2, b.get_height()-0.2, f'{v:.1f}', ha='center', va='top',
            fontsize=9, color='white', fontweight='bold')
a1.axhline(-6.76, color='#888', ls='--', lw=1, dashes=(4,2))
a1.text(7.4, -6.3, 'Biotin ref', fontsize=10, color='#888', ha='right')
a1.set_xticks(list(x)); a1.set_xticklabels(S, fontsize=10, rotation=25, ha='right')
ax_style(a1); a1.set_ylabel('Vina Score (kcal/mol)'); a1.set_title('A | Docking Score')

mg = [E[c]['dG_kcal'] for c in C]
b2 = a2.bar(x, mg, w, color=C8, ec='white', lw=0.4)
for b,v in zip(b2,mg):
    off, cl = (-1.2, 'white') if v<0 else (0.5, '#555')
    a2.text(b.get_x()+w/2, b.get_height()+off, f'{v:.0f}', ha='center', fontsize=9, color=cl, fontweight='bold')
a2.axhline(-29.97, color='#888', ls='--', lw=1, dashes=(4,2))
a2.set_xticks(list(x)); a2.set_xticklabels(S, fontsize=10, rotation=25, ha='right')
ax_style(a2); a2.set_ylabel('MM-GBSA (kcal/mol)'); a2.set_title('B | Binding Free Energy')

plt.tight_layout(pad=1.5)
plt.savefig(os.path.join(OUT, 'Figure_1_screening.pdf'), bbox_inches='tight'); plt.close()
print('OK Figure_1')

# ══════════════════════════════════
# Figure 2 — MD RMSD (100 ns)
# ══════════════════════════════════
fig, axes = plt.subplots(2, 4, figsize=(18, 9))
fig.suptitle('Figure 2 | MD Simulation Stability (100 ns)', fontsize=16, fontweight='bold')
for i,c in enumerate(C):
    ax = axes.flat[i]; df = rmsd[c]
    raw = df['RMSD_A']; ma = raw.rolling(51, center=True).mean()
    ax.plot(df['t'], raw, color=C8[i], lw=0.4, alpha=0.2)
    ax.plot(df['t'], ma, color=C8[i], lw=2)
    ax.axhline(raw.mean(), color=C8[i], ls='--', lw=0.8, alpha=0.3)
    ax_style(ax); ax.set_xlim(0, 100); ax.set_ylim(0, 10)
    ax.set_xlabel('Time (ns)'); ax.set_ylabel('RMSD (A)')
    ax.set_title(S[i], fontsize=12, fontweight='bold', color=C8[i])
plt.tight_layout(pad=1.5)
plt.savefig(os.path.join(OUT, 'Figure_2_MD.pdf'), bbox_inches='tight'); plt.close()
print('OK Figure_2')

# ══════════════════════════════════
# Figure 3 — FEL (top 4)
# ══════════════════════════════════
KT = 0.596
TOP = ['HYDROMORPHONE','FUROSEMIDE','NAFTAZONE','PHENOBARBITAL']
TL = ['Hydromorphone (best hit)', 'Furosemide', 'Naftazone', 'Phenobarbital (barbiturate)']
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
fig.suptitle('Figure 3 | Free Energy Landscape (PC1-PC2)', fontsize=16, fontweight='bold')
for i,c in enumerate(TOP):
    ax = axes.flat[i]; d = pca[c]
    pc1,pc2 = d['pc1'][:5000], d['pc2'][:5000]
    r = np.sqrt(pc1**2+pc2**2); m = r < np.percentile(r,99); pc1,pc2 = pc1[m], pc2[m]
    H,xe,ye = np.histogram2d(pc1,pc2,bins=60,density=True)
    H = gaussian_filter(H,sigma=1.2); H[H<1e-12]=1e-12; FE = -KT*np.log(H); FE -= FE.min()
    X,Y = np.meshgrid((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,indexing='ij')
    cnt = ax.contourf(X,Y,FE.T,levels=np.linspace(0,7,15),cmap='viridis',extend='max')
    ax_style(ax); ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
    ax.set_title(TL[i], fontsize=12, fontweight='bold')
    if i==3:
        cb = plt.colorbar(cnt,ax=ax,shrink=0.6,pad=0.05)
        cb.set_label('Free Energy (kcal/mol)', fontsize=10)
plt.tight_layout(pad=1.2)
plt.savefig(os.path.join(OUT, 'Figure_3_FEL.pdf'), bbox_inches='tight'); plt.close()
print('OK Figure_3')

# ══════════════════════════════════
# Figure 4 — Binding Energetics
# ══════════════════════════════════
fig = plt.figure(figsize=(14, 5.5))
gs = fig.add_gridspec(1, 3, wspace=0.45, left=0.06, right=0.97, bottom=0.18, top=0.88)
fig.suptitle('Figure 4 | Binding Energetics & Mode Classification', fontsize=16, fontweight='bold')

x = range(8); w = 0.55
ax_a = fig.add_subplot(gs[0,0])
ma_abs = [abs(E[c]['dG_kcal']) for c in C]
mc = [E[c]['dE_coulomb_LJ'] for c in C]
mg = [E[c]['dG_GB_polar'] for c in C]
ax_a.bar(x, ma_abs, w, color='#B2B2B2', ec='white', lw=0.4, label='Total')
ax_a.bar(x, [abs(v/4) for v in mc], w*0.6, color='#2166AC', ec='none', alpha=0.7, label='vdW+Elec /4')
ax_a.bar(x, [abs(v/4) for v in mg], w*0.3, color='#D6604D', ec='none', alpha=0.7, label='GB Polar /4')
ax_a.set_xticks(list(x)); ax_a.set_xticklabels(S, fontsize=10, rotation=25, ha='right')
ax_style(ax_a); ax_a.set_ylabel('|Binding Energy| (kcal/mol)'); ax_a.set_title('A | MM-GBSA Decomposition')
ax_a.legend(fontsize=9)

ax_b = fig.add_subplot(gs[0,1])
rat = [abs(mc[i])/max(abs(mg[i]),1) for i in range(8)]
bb = ax_b.barh(range(8), rat, 0.6, color=C8, ec='white', lw=0.4)
for b,v in zip(bb,rat):
    ax_b.text(b.get_width()+0.1, b.get_y()+b.get_height()/2, f'{v:.1f}x', va='center', fontsize=10, fontweight='bold')
ax_b.set_yticks(range(8)); ax_b.set_yticklabels(S, fontsize=10); ax_b.invert_yaxis()
ax_style(ax_b); ax_b.set_xlabel('|vdW| / |GB| Ratio'); ax_b.set_title('B | Energy Balance')

ax_c = fig.add_subplot(gs[0,2]); ax_c.axis('off')
CL = ['Lock&Key','Induced Fit','Lock&Key','Induced Fit',
      'Lock&Key','Induced Fit','Lock&Key','Induced Fit']
CT = [75,85,64,56,75,66,81,80]; CK = [561,538,571,513,543,223,703,480]
td = [[f'{S[i]}', f'{V[i]:.1f}', f'{ma_abs[i]:.0f}', CL[i], f'{CT[i]}%'] for i in range(8)]
cols = ['Compound','Vina','MM-GBSA','Mode','Contact']
tbl = ax_c.table(cellText=td, colLabels=cols, loc='center', cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1,1.6)
for i in range(8):
    tbl[i+1,0].set_text_props(color=C8[i], fontweight='bold')
    tbl[i+1,3].set_color('#3498DB' if CL[i]=='Lock&Key' else '#E74C3C')
ax_c.set_title('C | Summary', fontsize=14, fontweight='bold', pad=15)

plt.savefig(os.path.join(OUT, 'Figure_4_energetics.pdf'), bbox_inches='tight'); plt.close()
print('OK Figure_4')

for f in ['Figure_1_screening.pdf','Figure_2_MD.pdf','Figure_3_FEL.pdf','Figure_4_energetics.pdf']:
    p = os.path.join(OUT, f)
    print(f'  {f}: {os.path.getsize(p)/1024:.0f} KB')
print('All done!')
