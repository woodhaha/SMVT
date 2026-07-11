"""Extract frames from DCD using MDTraj with reduced topology."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import mdtraj as md

COMPOUNDS = ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
             'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']
TRAJ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

for c in COMPOUNDS:
    dcd = os.path.join(TRAJ_DIR, c, f'{c}_100ns.dcd')
    top = os.path.join(TRAJ_DIR, c, f'{c}_final.pdb')
    outdir = os.path.join(ANALYSIS_DIR, c, f'{c}_multiframes')
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(dcd):
        print(f'{c}: skip (no DCD)'); continue

    print(f'{c}: loading...', end=' ', flush=True)
    sys.stdout.flush()

    # Load with topology from prmtop if available, else use PDB
    try:
        t = md.load_dcd(dcd, top=top)
    except Exception as e:
        print(f'ERROR: {e}'); continue

    n = t.n_frames
    print(f'{n} frames', flush=True)
    step = max(1, n // 10)

    for i in range(10):
        idx = min(i * step, n - 1)
        fn = os.path.join(outdir, f'frame_{i}.pdb')
        t[idx].save(fn)
        sz = os.path.getsize(fn)
        print(f'  frame {i} (idx {idx}): {sz} bytes', flush=True)

print('ALL DONE')
