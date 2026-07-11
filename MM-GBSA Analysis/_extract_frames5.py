"""Extract frames from DCD using pure OpenMM."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import openmm.app as app
import openmm.unit as unit

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

    print(f'{c}: loading OpenMM topology...')
    pdb = app.PDBFile(top)
    natoms = pdb.topology.getNumAtoms()
    print(f'{c}: reading DCD ({os.path.getsize(dcd)//1024**2} MB)...')
    sys.stdout.flush()

    # Read DCD as binary, parse header manually for frame count
    with open(dcd, 'rb') as f:
        # CHARMM DCD header: first 4 bytes = "CORD", then various ints
        import struct
        f.read(4)  # "CORD" or marker
        # NSET (number of frames) at offset 8 (big-endian in CHARMM DCD)
        f.read(4)  # skip
        n_frames = struct.unpack('i', f.read(4))[0]
        f.read(4)  # skip
        n_atoms_dcd = struct.unpack('i', f.read(4))[0]
        print(f'  DCD header: {n_frames} frames, {n_atoms_dcd} atoms (topology: {natoms})')
        sys.stdout.flush()

    # Use OpenMM's DCD reader with file handle
    with open(dcd, 'rb') as f:
        dcd_file = app.DCDFile(f, pdb.topology, dt=0.1*unit.picoseconds)
        n = dcd_file.getNumFrames()
        print(f'  OpenMM reads: {n} frames')
        sys.stdout.flush()

        step = max(1, n // 10)
        for i in range(10):
            idx = min(i * step, n - 1)
            dcd_file.seek(idx)
            # DCD stores in angstroms, but OpenMM converts; positions may be in nm
            frame_positions = dcd_file.readFrame()
            # Convert from nm to nm (already correct for OpenMM)
            fn = os.path.join(outdir, f'frame_{i}.pdb')
            with open(fn, 'w') as fout:
                app.PDBFile.writeFile(pdb.topology, frame_positions, fout)
            print(f'  frame {i} (idx {idx}): {os.path.getsize(fn)} bytes')
            sys.stdout.flush()

    print(f'{c}: done')

print('ALL DONE')
