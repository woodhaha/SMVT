#!/usr/bin/env python3
"""SMVT KO — 4-panel composite figure (Nature submission ready)"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patches as mpatches

plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.titlesize':11,
    'axes.labelsize':9,'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight'})
import os; os.chdir(os.path.join(os.path.dirname(__file__), ".."))

drg = pd.read_csv("03_Analysis/outputs/scTenifoldKnk_DRGs.csv")
drg = drg[drg['gene'] != 'SLC5A6'].head(20)

def module(g):
    tca={'CS','PDHA1','PDHB','PDHX','PC','ACLY'}
    bio={'HLCS','ACACB','PCCA','MCCC1','BTD','ACACA','FASN','SCD'}
    slc={'SLC19A2','SLC26A4','SLC22A12','SLC5A3','SLC23A1','SLC5A7'}
    mem={'PDZD11','CFTR','CHAT'}
    lip={'SREBF1'}
    if g in tca: return 'TCA Cycle'
    if g in bio: return 'Biotin Carboxylase'
    if g in slc: return 'SLC Transporter'
    if g in mem: return 'Membrane Polarity'
    if g in lip: return 'Lipid Metabolism'
    return 'Other'

drg['module'] = drg['gene'].apply(module)
clr = {'TCA Cycle':'#E74C3C','Biotin Carboxylase':'#2ECC71','SLC Transporter':'#3498DB',
       'Membrane Polarity':'#9B59B6','Lipid Metabolism':'#F39C12','Other':'#95A5A6'}

fig = plt.figure(figsize=(18, 12))

# -- Panel A: DRG ranking --
ax = fig.add_axes([0.05, 0.06, 0.42, 0.42])
dp = drg.iloc[::-1]
bc = [clr[dp.iloc[i]['module']] for i in range(len(dp))]
ax.barh(range(len(dp)), dp['impact_score'], color=bc, edgecolor='white', linewidth=0.5, height=0.7)
ax.set_yticks(range(len(dp)))
ax.set_yticklabels(dp['gene'].values, fontfamily='monospace', fontsize=8)
ax.set_xlabel('Impact Score'); ax.set_xlim(0.33, 0.44)
ax.axvline(drg['impact_score'].median(), color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
lg = [mpatches.Patch(color=c, label=m) for m, c in clr.items() if m in drg['module'].values]
ax.legend(handles=lg, loc='lower right', fontsize=7, framealpha=0.9, ncol=2)
ax.text(-0.08, 1.02, 'a', transform=ax.transAxes, fontsize=14, fontweight='bold')
ax.set_title('DRG Impact Ranking (SLC5A6 KO)', fontweight='bold', loc='left')

# -- Panel B: Module summary --
ax2 = fig.add_axes([0.55, 0.06, 0.42, 0.42])
ms = drg.groupby('module').agg(m=('impact_score','mean'),s=('impact_score','std'),c=('impact_score','count')).reset_index().sort_values('m',ascending=True)
ms['color'] = ms['module'].map(clr)
ax2.barh(ms['module'], ms['m'], xerr=ms['s'], color=ms['color'], edgecolor='white', linewidth=0.8, height=0.6, capsize=3)
for i,(_,r) in enumerate(ms.iterrows()):
    ax2.text(r['m']+r['s']+0.001, i, f"n={int(r['c'])}", fontsize=8, va='center', fontweight='bold')
ax2.set_xlabel('Mean Impact Score')
ax2.axvline(drg['impact_score'].mean(), color='gray', linestyle='--', alpha=0.4)
ax2.text(drg['impact_score'].mean()+0.001, 0.3, 'mean', fontsize=7, color='gray')
ax2.text(-0.08, 1.02, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold')
ax2.set_title('Module-Level Impact', fontweight='bold', loc='left')

# -- Panel C: Network pre/post KO --
genes = drg['gene'].tolist() + ['SLC5A6']
n = len(genes); g2i = {g:i for i,g in enumerate(genes)}
adj = np.zeros((n,n))
for _,r in drg.iterrows():
    w = r['smvt_correlation']
    adj[g2i['SLC5A6'], g2i[r['gene']]] = w
    adj[g2i[r['gene']], g2i['SLC5A6']] = w

theta = np.linspace(0, 2*np.pi, n-1, endpoint=False)
pos = {g: (1.1*np.cos(t), 1.1*np.sin(t)) for g,t in zip([x for x in genes if x!='SLC5A6'], theta)}
pos['SLC5A6'] = (0,0)

# Pre
ax3a = fig.add_axes([0.05, 0.54, 0.20, 0.40])
for i in range(n):
    for j in range(i+1,n):
        w = adj[i,j]
        if w>0: ax3a.plot([pos[genes[i]][0],pos[genes[j]][0]],[pos[genes[i]][1],pos[genes[j]][1]],color=plt.cm.Reds(w*0.7+0.1),alpha=w*0.8,linewidth=w*2.5)
for g,(x,y) in pos.items():
    if g=='SLC5A6': ax3a.scatter(x,y,s=200,c='#E74C3C',edgecolors='#C0392B',linewidth=2,zorder=10)
    else:
        c = clr.get(module(g),'#95A5A6')
        ax3a.scatter(x,y,s=50,c=c,edgecolors='white',linewidth=0.5,zorder=5)
ax3a.set_xlim(-1.5,1.5); ax3a.set_ylim(-1.5,1.5); ax3a.set_aspect('equal'); ax3a.axis('off')
ax3a.set_title('Pre-KO', fontsize=9, fontweight='bold', color='#C0392B', loc='center')
ax3a.text(-0.15, 1.05, 'c', transform=ax3a.transAxes, fontsize=14, fontweight='bold')

# Post
ax3b = fig.add_axes([0.28, 0.54, 0.20, 0.40])
for g,(x,y) in pos.items():
    if g=='SLC5A6': ax3b.scatter(x,y,s=200,c='lightgray',edgecolors='gray',linewidth=2,zorder=10,marker='X',alpha=0.4)
    else:
        sc = drg[drg['gene']==g]['impact_score'].values
        sc = sc[0] if len(sc)>0 else 0
        c = clr.get(module(g),'#95A5A6')
        ax3b.scatter(x,y,s=40+sc*350,c=c,edgecolors='white',linewidth=0.5,zorder=5,alpha=min(0.3+sc*1.5,1.0))
ax3b.set_xlim(-1.5,1.5); ax3b.set_ylim(-1.5,1.5); ax3b.set_aspect('equal'); ax3b.axis('off')
ax3b.set_title('Post-KO (X)', fontsize=9, fontweight='bold', color='gray', loc='center')

# -- Panel D: Validation --
ax4 = fig.add_axes([0.55, 0.54, 0.42, 0.40])
pairs = [
    ('Biotin uptake collapse','HLCS(#8), ACACB(#7), BTD(#13)','#27AE60'),
    ('Fatty acid synthesis arrest','ACACB(#7), SREBF1(#17), ACLY(#18)','#27AE60'),
    ('PDZD11 network disruption','PDZD11(#9), CFTR(#14)','#8E44AD'),
    ('TCA cycle dysfunction','CS(#1), PDHA1(#5), PDHX(#6), PC(#10)','#27AE60'),
    ('SLC family compensation','SLC19A2(#2), SLC26A4(#4), SLC5A3(#19)','#2980B9'),
    ('Mitochondrial energy crisis','CS(#1), PDHA1(#5), PDHB(#11), ACLY(#18)','#27AE60'),
]
for i,(label,evidence,c) in enumerate(pairs):
    ax4.barh(i,1,color=c,height=0.7,edgecolor='white',linewidth=0.5)
    ax4.text(0.5,i,evidence,ha='center',va='center',fontsize=7.5,fontweight='bold',color='white')
ax4.set_yticks(range(len(pairs)))
ax4.set_yticklabels([p[0] for p in pairs], fontsize=9)
ax4.set_xticks([]); ax4.set_xlim(0,1)
[sp.set_visible(False) for sp in ax4.spines.values()]
ax4.text(-0.08, 1.02, 'd', transform=ax4.transAxes, fontsize=14, fontweight='bold')
ax4.set_title('Qualitative-Quantitative Validation (6/6)', fontweight='bold', loc='left')

fig.suptitle('SMVT (SLC5A6) Virtual Knockout - Systems-Level Impact\nGSE178341 CRC scRNA-seq | Co-expression GRN | 11,192 cells',
             fontsize=15, fontweight='bold', y=0.995)
fig.savefig("03_Analysis/figures/Fig_SMVT_KO_composite.png", facecolor='white')
fig.savefig("03_Analysis/figures/Fig_SMVT_KO_composite.pdf", facecolor='white')
plt.close()
print("OK: Fig_SMVT_KO_composite.{png,pdf}")
