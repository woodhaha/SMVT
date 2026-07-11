"""2D schematic: Vina docked binding modes for top compounds."""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Load docked pose coords for Hydromorphone and Biotin
compounds = ['HYDROMORPHONE', 'BIOTIN']
palette = {'HYDROMORPHONE': '#9b59b6', 'BIOTIN': '#3498db'}

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
fig.suptitle('Vina Docked Binding Modes (Top 2 Cross-Validated Candidates)',
             fontsize=13, fontweight='bold')

for idx, name in enumerate(compounds):
    ax = axes[idx]
    fn = os.path.join(ANALYSIS, f'{name}_docked.pdb')
    if not os.path.exists(fn): continue

    # Parse PDB coords (robust: split on whitespace)
    coords = []
    atom_names = []
    for l in open(fn):
        l = l.rstrip()
        if l.startswith(('ATOM', 'HETATM')):
            parts = l.split()
            # parts: ['HETATM', '1', 'C', 'LIG', 'A', '1', '3.226', '3.321', '-2.367', ...]
            # find xyz - they're always at positions -5, -4, -3 from end
            x = float(parts[-5]); y = float(parts[-4]); z = float(parts[-3])
            coords.append([x, y, z])
            atom_names.append(parts[2])

    coords = np.array(coords)
    center = coords.mean(axis=0)

    # 2D projection: view along Z axis (top-down)
    ax.scatter(coords[:,0], coords[:,1], s=30, c=palette[name],
               edgecolors='k', linewidths=0.5, zorder=5, alpha=0.8)

    # Draw bonds between nearby atoms (AutoDock Vina tree)
    # Use distance-based heuristic: connect atoms < 1.8A apart
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            d = np.linalg.norm(coords[i] - coords[j])
            if d < 1.8:
                ax.plot([coords[i,0], coords[j,0]],
                       [coords[i,1], coords[j,1]],
                       '-', color=palette[name], linewidth=1.5, alpha=0.6, zorder=3)

    # Annotate center
    ax.text(center[0], center[1], f'  {name}', fontsize=8, fontweight='bold',
            color=palette[name])

    # Look up Vina score
    import subprocess
    log_fn = os.path.join(ANALYSIS, f'{name}_vina_log.txt')
    score = None
    if os.path.exists(log_fn):
        for ln in open(log_fn):
            parts = ln.strip().split()
            if len(parts) >= 2 and parts[0] == '1':
                score = float(parts[1]); break
    if score:
        mg_fn = os.path.join(ANALYSIS, name, 'mmgbsa_results.json')
        mg = ''
        if os.path.exists(mg_fn):
            mg_v = json.load(open(mg_fn))['dG_kcal']
            mg = f'MM-GBSA={mg_v:.1f}'
        ax.set_title(f'{name}\nVina={score:.2f} kcal/mol {mg}', fontsize=9)

    ax.set_xlabel('X (A)')
    ax.set_ylabel('Y (A)')
    ax.grid(alpha=0.15)
    ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join(ANALYSIS, 'ALL_binding_poses.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join(ANALYSIS, 'ALL_binding_poses.pdf'), dpi=150, bbox_inches='tight')
plt.close()
print('ALL_binding_poses.png + .pdf done')
