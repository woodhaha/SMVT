"""Extract 10 evenly spaced frames from each trajectory DCD for multi-frame MM-GBSA."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import mdtraj as md

COMPOUNDS = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
             'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']
TRAJ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
N_FRAMES = 10

for c in COMPOUNDS:
    dcd = os.path.join(TRAJ_DIR, c, f'{c}_100ns.dcd')
    top = os.path.join(TRAJ_DIR, c, f'{c}_final.pdb')
    outdir = os.path.join(ANALYSIS_DIR, c, f'{c}_multiframes')
    os.makedirs(outdir, exist_ok=True)

    if not os.path.exists(dcd):
        print(f'{c}: DCD not found, skipping')
        continue

    print(f'{c}: loading DCD...', end=' ', flush=True)
    t = md.load_dcd(dcd, top=top)
    print(f'{t.n_frames} frames, {t.n_atoms} atoms')

    step = max(1, t.n_frames // N_FRAMES)
    for i in range(N_FRAMES):
        idx = i * step if i * step < t.n_frames else t.n_frames - 1
        fn = os.path.join(outdir, f'frame_{i}.pdb')
        # Save as PDB
        t[idx].save(fn)
        print(f'  frame {i} (idx {idx})')

print('ALL DONE')
