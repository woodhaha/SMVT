#!/usr/bin/env python3
"""Summary visualization: centroid clusters + binding classes + key metrics."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

DATA = {
    "NAFTAZONE":         {"clusters": [561,122,301,15],         "k":4, "rmsd":5.52, "pc":96.6, "class":"Lock-and-Key", "color":"#3498db"},
    "BIOTIN":            {"clusters": [538,162,106,67,42,28,56], "k":7, "rmsd":4.90, "pc":84.7, "class":"Induced Fit", "color":"#e74c3c"},
    "ESKETAMINE":        {"clusters": [571,386,42],              "k":3, "rmsd":5.03, "pc":97.2, "class":"Lock-and-Key", "color":"#3498db"},
    "FUROSEMIDE":        {"clusters": [513,225,77,70,57,49,8],   "k":7, "rmsd":4.94, "pc":86.0, "class":"Induced Fit", "color":"#e74c3c"},
    "GABAPENTIN_ENACARBIL": {"clusters": [543,259,197],         "k":3, "rmsd":4.94, "pc":88.7, "class":"Lock-and-Key", "color":"#3498db"},
    "HYDROMORPHONE":     {"clusters": [223,23,192,269,49,98,116,29], "k":8, "rmsd":7.03, "pc":82.5, "class":"Induced Fit", "color":"#e74c3c"},
    "PHENOBARBITAL":     {"clusters": [703,215,57,24],           "k":4, "rmsd":6.02, "pc":94.6, "class":"Lock-and-Key", "color":"#3498db"},
    "RIBOFLAVIN":        {"clusters": [480,245,138,122,8,4,2],   "k":7, "rmsd":4.63, "pc":82.9, "class":"Induced Fit", "color":"#e74c3c"},
}

COMPOUNDS = list(DATA.keys())
N = len(COMPOUNDS)

# ═════════════════════════════════════════════════╗
#  FIGURE 1: Centroid cluster composition (stacked bar + key metrics)
# ═════════════════════════════════════════════════╝
fig, ax = plt.subplots(figsize=(14, 6))
fig.suptitle('SMVT 100ns MD — Optimal Conformational State Landscape', fontsize=14, fontweight='bold')

# Stacked bar: cluster size distribution normalized
colors = plt.cm.Set2(np.linspace(0.05, 0.95, 8))
bottom = np.zeros(N)

bars = []
for ci in range(8):
    vals = []
    for c in COMPOUNDS:
        cl = DATA[c]["clusters"]
        vals.append(cl[ci] if ci < len(cl) else 0)
    if ci == 0:
        b = ax.bar(range(N), vals, color=colors[ci], edgecolor='white', linewidth=0.5, label='Cluster 1 (largest)')
    else:
        b = ax.bar(range(N), vals, bottom=bottom, color=colors[ci], edgecolor='white', linewidth=0.5,
                   label=f'Cluster {ci+1}' if ci < 7 else '_nolegend_')
    bottom += vals

# Labels
short = {"NAFTAZONE":"NAF","BIOTIN":"BIO","ESKETAMINE":"ESK","FUROSEMIDE":"FUR",
         "GABAPENTIN_ENACARBIL":"GAB","HYDROMORPHONE":"HYD","PHENOBARBITAL":"PHE","RIBOFLAVIN":"RIB"}
ax.set_xticks(range(N))
ax.set_xticklabels([short[c] for c in COMPOUNDS], fontsize=9)
ax.set_ylabel('Number of frames')
ax.set_title('Cluster composition (k-means on PC1-PC2)', fontsize=12)
ax.legend(fontsize=6, ncol=4, loc='upper left')

# Annotate top cluster %
for i, c in enumerate(COMPOUNDS):
    pct = DATA[c]["clusters"][0] / sum(DATA[c]["clusters"]) * 100
    ax.text(i, DATA[c]["clusters"][0]/2, f'{pct:.0f}%', ha='center', va='center',
            fontsize=10, fontweight='bold', color='white')
    ax.text(i, sum(DATA[c]["clusters"]) + 15, f'k={DATA[c]["k"]}', ha='center', fontsize=8, color='#555')

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig('ALL_cluster_composition.png', dpi=150, bbox_inches='tight')
plt.close()
print('1. ALL_cluster_composition.png')

# ═════════════════════════════════════════════════╗
#  FIGURE 2: Binding class classification dashboard
# ═════════════════════════════════════════════════╝
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('SMVT 100ns MD — Binding Mode Classification & Key Metrics', fontsize=14, fontweight='bold')

# Panel 1: Clusters vs RMSD scatter
ax = axes[0,0]
for c in COMPOUNDS:
    d = DATA[c]
    cls = "Lock-and-Key" if d["class"] == "Lock-and-Key" else "Induced Fit"
    mk = 'o' if cls == "Lock-and-Key" else 'D'
    clr = '#3498db' if cls == "Lock-and-Key" else '#e74c3c'
    top_pct = d["clusters"][0] / sum(d["clusters"]) * 100
    ax.scatter(d["k"], d["rmsd"], s=top_pct*8, c=clr, marker=mk, alpha=0.7, edgecolors='k', linewidths=0.5, zorder=5)
    ax.annotate(short[c], (d["k"], d["rmsd"]), fontsize=7, ha='center', va='bottom',
                xytext=(0, 5), textcoords='offset points')

ax.set_xlabel('Number of clusters (k)'); ax.set_ylabel('Mean RMSD (Å)')
ax.set_title('Cluster count vs RMSD (size = top cluster %)')
ax.grid(alpha=0.2)

# Legend
lk = mpatches.Patch(color='#3498db', label='Lock-and-Key', alpha=0.7)
idf = mpatches.Patch(color='#e74c3c', label='Induced Fit', alpha=0.7)
ax.legend(handles=[lk, idf], loc='lower right', fontsize=9)

# Panel 2: PC variance bar
ax = axes[0,1]
x = range(N)
w = 0.35
lock_idx = [i for i,c in enumerate(COMPOUNDS) if DATA[c]["class"] == "Lock-and-Key"]
ind_idx = [i for i,c in enumerate(COMPOUNDS) if DATA[c]["class"] != "Lock-and-Key"]

ax.bar([i for i in lock_idx], [DATA[COMPOUNDS[i]]["pc"] for i in lock_idx],
       w, color='#3498db', alpha=0.8, edgecolor='k', linewidth=0.5, label='Lock-and-Key')
ax.bar([i for i in ind_idx], [DATA[COMPOUNDS[i]]["pc"] for i in ind_idx],
       w, color='#e74c3c', alpha=0.8, edgecolor='k', linewidth=0.5, label='Induced Fit')
ax.set_xticks(range(N))
ax.set_xticklabels([short[c] for c in COMPOUNDS], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('PC1+PC2 cumulative variance (%)')
ax.set_title('PCA Convergence Strength')
ax.legend(fontsize=8)
ax.grid(axis='y', alpha=0.3)

# Panel 3: Top cluster % bar
ax = axes[1,0]
top_pcts = [DATA[c]["clusters"][0]/sum(DATA[c]["clusters"])*100 for c in COMPOUNDS]
colors_bar = ['#3498db' if DATA[c]["class"]=="Lock-and-Key" else '#e74c3c' for c in COMPOUNDS]
bars = ax.bar(range(N), top_pcts, color=colors_bar, alpha=0.8, edgecolor='k', linewidth=0.5)
ax.set_xticks(range(N)); ax.set_xticklabels([short[c] for c in COMPOUNDS], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Top cluster (%)'); ax.set_title('Conformational Concentration')
ax.grid(axis='y', alpha=0.3)
for bar, pct in zip(bars, top_pcts):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1, f'{pct:.0f}%',
            ha='center', fontsize=8, fontweight='bold')

# Panel 4: Classification table as heatmap
ax = axes[1,1]
ax.axis('off')
table_data = []
for c in COMPOUNDS:
    d = DATA[c]
    top_pct = d["clusters"][0] / sum(d["clusters"]) * 100
    table_data.append([c, str(d["k"]), f'{top_pct:.0f}%',
                       f'{d["rmsd"]:.1f} A', f'{d["pc"]:.1f}%',
                       d["class"]])

col_labels = ['Compound', 'Clusters', 'Top%', 'RMSD', 'PC1+2', 'Binding Class']
cell_text = table_data
cell_colors = [['#f8f9fa']*6 for _ in range(N)]
for i, c in enumerate(COMPOUNDS):
    if DATA[c]["class"] == "Lock-and-Key":
        cell_colors[i] = ['#e8f4f8', '#e8f4f8', '#e8f4f8', '#e8f4f8', '#e8f4f8', '#d4eef7']
    else:
        cell_colors[i] = ['#fde8e8', '#fde8e8', '#fde8e8', '#fde8e8', '#fde8e8', '#fcc']

table = ax.table(cellText=cell_text, colLabels=col_labels, cellColours=cell_colors,
                 loc='center', cellLoc='center')
table.auto_set_font_size(False); table.set_fontsize(8)
table.scale(1, 1.4)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_fontsize(9); cell.set_text_props(weight='bold')
ax.set_title('Classification Summary', fontsize=12)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('ALL_classification_dashboard.png', dpi=150, bbox_inches='tight')
plt.close()
print('2. ALL_classification_dashboard.png')

# ═════════════════════════════════════════════════╗
#  FIGURE 3: Centroid frame timeline
# ═════════════════════════════════════════════════╝
fig, ax = plt.subplots(figsize=(12, 5))
fig.suptitle('Time Distribution of Best Conformations (Cluster Centroid Frames)', fontsize=14, fontweight='bold')

centroid_times = {
    "NAFTAZONE": 85.7, "BIOTIN": 43.4, "ESKETAMINE": 42.4,
    "FUROSEMIDE": 19.1, "GABAPENTIN_ENACARBIL": 41.1,
    "HYDROMORPHONE": 71.5, "PHENOBARBITAL": 30.8, "RIBOFLAVIN": 33.2
}

colors_t = ['#3498db','#e74c3c','#3498db','#e74c3c','#3498db','#e74c3c','#3498db','#e74c3c']

for i, c in enumerate(COMPOUNDS):
    t = centroid_times[c]
    marker = 'o' if DATA[c]["class"] == "Lock-and-Key" else 'D'
    ax.scatter(t, i, s=200, c=colors_t[i], marker=marker, edgecolors='k', linewidths=1, zorder=5)
    ax.annotate(f'  {c} ({t:.0f}ns)', (t, i), fontsize=9, va='center')

ax.set_yticks(range(N)); ax.set_yticklabels([])
ax.set_xlabel('Simulation Time (ns)')
ax.set_title('Each dot = the most representative frame for that compound')
# Shade the first 44ns region
ax.axvspan(0, 44, alpha=0.08, color='blue', label='First 44ns')
ax.axvspan(44, 88, alpha=0.08, color='orange', label='Last 44ns')
ax.legend(fontsize=8)
ax.set_xlim(0, 95)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig('ALL_centroid_timeline.png', dpi=150, bbox_inches='tight')
plt.close()
print('3. ALL_centroid_timeline.png')

print('\nALL SUMMARY FIGURES DONE')
print('Files: ALL_cluster_composition.png, ALL_classification_dashboard.png, ALL_centroid_timeline.png')
