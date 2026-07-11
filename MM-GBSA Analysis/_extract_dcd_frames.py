"""Extract frames from DCD using raw binary parsing."""
import os, sys, struct, warnings
import numpy as np
warnings.filterwarnings("ignore")
import openmm.app as app
import openmm.unit as unit

TRAJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))
COMPOUNDS = ['BIOTIN']
N = 5

for c in COMPOUNDS:
    dcd = os.path.join(TRAJ, c, f'{c}_100ns.dcd')
    top = os.path.join(TRAJ, c, f'{c}_final.pdb')
    outdir = os.path.join(ANALYSIS, c, f'{c}_multiframes')
    os.makedirs(outdir, exist_ok=True)

    print(f'{c}:...', end=' ', flush=True)
    pdb = app.PDBFile(top)
    natoms = pdb.topology.getNumAtoms()

    with open(dcd, 'rb') as f:
        # Parse header
        reclen = struct.unpack('<i', f.read(4))[0]
        assert reclen == 84
        magic = f.read(4)  # CORD
        nset = struct.unpack('<i', f.read(4))[0]
        f.read(4*5)  # ISTART, NSAVC, (3 unused)
        box_flag = struct.unpack('<i', f.read(4))[0]
        f.read(reclen - 32)  # rest of header
        trailing = struct.unpack('<i', f.read(4))[0]

        print(f'{nset} frames, {natoms} atoms, box={box_flag}')

        step = max(1, nset // N)
        selected = list(range(0, nset, step))[:N]

        for fi, frame_idx in enumerate(selected):
            # Seek to frame start: header was 92 bytes + 4 trailing = 96
            f.seek(0)
            # Skip header
            reclen = struct.unpack('<i', f.read(4))[0]
            # Skip past header (reclen bytes + 4 trailing)
            f.seek(84 + 4, os.SEEK_CUR)

            # Skip frames before target
            for skip_idx in range(frame_idx):
                if box_flag:
                    rc = struct.unpack('<i', f.read(4))[0]
                    f.read(rc)
                    struct.unpack('<i', f.read(4))[0]
                for _ in range(3):
                    rc = struct.unpack('<i', f.read(4))[0]
                    f.read(rc)
                    struct.unpack('<i', f.read(4))[0]

            # Read target frame
            if box_flag:
                rc = struct.unpack('<i', f.read(4))[0]
                cell = struct.unpack('<6d', f.read(rc))
                struct.unpack('<i', f.read(4))[0]

            x = np.zeros(natoms, dtype=np.float32)
            y = np.zeros(natoms, dtype=np.float32)
            z = np.zeros(natoms, dtype=np.float32)

            for arr in [x, y, z]:
                rc = struct.unpack('<i', f.read(4))[0]
                f.readinto(arr)
                struct.unpack('<i', f.read(4))[0]

            # Positions in Angstroms → nm for OpenMM
            pos_nm = np.column_stack([x, y, z]) * 0.1 * unit.nanometers

            out_fn = os.path.join(outdir, f'frame_{fi}.pdb')
            with open(out_fn, 'w') as fout:
                app.PDBFile.writeFile(pdb.topology, pos_nm, fout)
            print(f'  f{fi}@{frame_idx} OK')
    print(f'  Done: {outdir}')

print('ALL DONE')
