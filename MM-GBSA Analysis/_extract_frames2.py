"""Extract frames from DCD using OpenMM's built-in reader."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm.app as app
import openmm.unit as unit

COMPOUNDS = ['BIOTIN']
TRAJ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

c = 'BIOTIN'
dcd = os.path.join(TRAJ_DIR, c, f'{c}_100ns.dcd')
top_pdb = os.path.join(TRAJ_DIR, c, f'{c}_final.pdb')
outdir = os.path.join(ANALYSIS_DIR, c, f'{c}_multiframes')
os.makedirs(outdir, exist_ok=True)

print(f'Loading DCD: {dcd}')
pdb = app.PDBFile(top_pdb)
dcd_file = app.DCDFile(dcd, pdb.topology, dt=0.1*unit.picoseconds)

n_frames = dcd_file.getNumFrames()
print(f'N frames: {n_frames}, N atoms: {pdb.topology.getNumAtoms()}')

step = max(1, n_frames // 10)
selected = list(range(0, n_frames, step))[:10]

for i, idx in enumerate(selected):
    dcd_file.seek(idx)
    positions = dcd_file.readFrame() * unit.nanometers  # DCD is in nm

    out_fn = os.path.join(outdir, f'frame_{i}.pdb')
    with open(out_fn, 'w') as f:
        app.PDBFile.writeFile(pdb.topology, positions, f)
    print(f'  frame {i}: idx {idx} -> {out_fn}')
    sys.stdout.flush()

dcd_file.close()
print('DONE')
