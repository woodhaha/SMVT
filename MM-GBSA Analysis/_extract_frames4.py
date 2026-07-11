"""Extract frames using mdtraj for load + OpenMM for PDB write."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import mdtraj as md
import openmm.app as app
import openmm.unit as unit
import numpy as np

TRAJ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
COMPOUNDS = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
             'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']

for c in COMPOUNDS:
    dcd = os.path.join(TRAJ_DIR, c, f'{c}_100ns.dcd')
    top = os.path.join(TRAJ_DIR, c, f'{c}_final.pdb')
    outdir = os.path.join(ANALYSIS_DIR, c, f'{c}_multiframes')
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(dcd):
        print(f'{c}: skip (no DCD)'); continue

    print(f'{c}: loading...', end=' ', flush=True)
    sys.stdout.flush()
    t = md.load_dcd(dcd, top=top)
    n = t.n_frames
    print(f'{n} frames')

    # Load OpenMM topology from the same PDB
    pdb = app.PDBFile(top)

    step = max(1, n // 10)
    for i in range(10):
        idx = min(i * step, n - 1)
        # Get positions in nanometers (OpenMM uses nm)
        positions = t.xyz[idx] * unit.nanometers
        fn = os.path.join(outdir, f'frame_{i}.pdb')
        with open(fn, 'w') as f:
            app.PDBFile.writeFile(pdb.topology, positions, f)
        sz = os.path.getsize(fn)
        print(f'  frame {i} (idx {idx}): {sz} bytes', flush=True)

print('ALL DONE')
