"""Minimal DCD reader + frame extractor using OpenMM for PDB output."""
import os, sys, struct, warnings
import numpy as np
warnings.filterwarnings("ignore")
import openmm.app as app
import openmm.unit as unit

TRAJ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SMVT_MD', 'trajectories')
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
N_FRAMES = 5  # 5 frames per compound for averaging

def read_dcd_header(f):
    """Read DCD header, return (n_frames, natoms, box_flag)."""
    reclen = struct.unpack('<i', f.read(4))[0]
    magic = f.read(4)
    nset = struct.unpack('<i', f.read(4))[0]
    istart = struct.unpack('<i', f.read(4))[0]
    nsavc = struct.unpack('<i', f.read(4))[0]
    ncrd = struct.unpack('<i', f.read(4))[0]
    f.read(4)
    # Byte 24-27: box flag at offset 20 in the 84-byte header
    # Actually let me just scan for it
    f.read(4)  # skip
    box_flag = struct.unpack('<i', f.read(4))[0]
    # Skip rest of header
    f.seek(4, os.SEEK_CUR)  # trailing reclen
    return nset, box_flag

def read_dcd_frame_xyz(f, natoms, box_flag):
    """Read one frame from DCD, return (natoms, 3) positions in Angstroms."""
    if box_flag:
        # Unit cell: 6 doubles
        reclen = struct.unpack('<i', f.read(4))[0]
        cell = struct.unpack('<6d', f.read(reclen))
        reclen2 = struct.unpack('<i', f.read(4))[0]

    x = np.zeros(natoms, dtype=np.float32)
    y = np.zeros(natoms, dtype=np.float32)
    z = np.zeros(natoms, dtype=np.float32)

    for name, arr in [('X', x), ('Y', y), ('Z', z)]:
        reclen = struct.unpack('<i', f.read(4))[0]
        assert reclen == 4 * natoms, f'{name}: reclen={reclen} != {4*natoms}'
        f.readinto(arr)
        reclen2 = struct.unpack('<i', f.read(4))[0]

    return np.column_stack([x, y, z])

for c in ['BIOTIN','ESKETAMINE','FUROSEMIDE','GABAPENTIN_ENACARBIL',
          'HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','RIBOFLAVIN']:
    dcd = os.path.join(TRAJ_DIR, c, f'{c}_100ns.dcd')
    top = os.path.join(TRAJ_DIR, c, f'{c}_final.pdb')
    outdir = os.path.join(ANALYSIS_DIR, c, f'{c}_multiframes')
    os.makedirs(outdir, exist_ok=True)

    if not os.path.exists(dcd):
        print(f'{c}: skip (no DCD)')
        continue

    print(f'{c}:', end=' ', flush=True)
    pdb = app.PDBFile(top)
    natoms = pdb.topology.getNumAtoms()
    print(f'Topology: {natoms} atoms, DCD: {os.path.getsize(dcd)//1024//1024} MB', end=' ')

    with open(dcd, 'rb') as f:
        n_frames, box_flag = read_dcd_header(f)
        print(f'{n_frames} frames, box={box_flag}')

        step = max(1, n_frames // N_FRAMES)
        selected = list(range(0, n_frames, step))[:N_FRAMES]

        for fi, frame_idx in enumerate(selected):
            # Seek to frame position
            f.seek(0)
            read_dcd_header(f)  # re-read header
            offset = 0
            if box_flag:
                for _ in range(frame_idx):
                    reclen = struct.unpack('<i', f.read(4))[0]
                    f.read(reclen + 4)
                    for _ in range(3):  # X, Y, Z
                        reclen = struct.unpack('<i', f.read(4))[0]
                        f.seek(reclen + 4, os.SEEK_CUR)
                positions_A = read_dcd_frame_xyz(f, natoms, box_flag)
            else:
                for _ in range(frame_idx):
                    for _ in range(3):
                        reclen = struct.unpack('<i', f.read(4))[0]
                        f.seek(reclen + 4, os.SEEK_CUR)
                positions_A = read_dcd_frame_xyz(f, natoms, box_flag)

            # Convert Angstroms to nm for OpenMM
            pos_nm = positions_A * 0.1 * unit.nanometers

            fn = os.path.join(outdir, f'frame_{fi}.pdb')
            with open(fn, 'w') as fout:
                app.PDBFile.writeFile(pdb.topology, pos_nm, fout)
            sz = os.path.getsize(fn)
            print(f'  [{fi}@{frame_idx}] {sz} bytes', flush=True)

print('ALL DONE')
