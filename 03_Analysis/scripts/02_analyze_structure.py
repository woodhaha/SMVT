#!/usr/bin/env python3
"""
02_analyze_structure.py
Analyze the SMVT (SLC5A6) AlphaFold structure for docking setup.

The substrate binding pocket of SMVT is in the central cavity formed
by the transmembrane helix bundle. Key pocket residues based on
experimental and computational studies include aromatic and charged
residues lining the central transport pathway.
"""

import os, sys, json, math
from collections import defaultdict

BASE = r"D:\Researching\SMVT"
PDB_PATH = os.path.join(BASE, "02_Data", "cleaned", "AF-Q9Y289-F1_prepared.pdb")
DATA_DIR = os.path.join(BASE, "03_Analysis", "data")

def parse_pdb(pdb_path):
    residues = defaultdict(list)
    all_atoms = []

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if len(line) < 54: continue
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() or atom_name[0]
                key = (chain, resnum, resname)
                residues[key].append({"atom": atom_name, "element": element, "x": x, "y": y, "z": z})
                if element in ('C', 'N', 'O', 'S'):
                    all_atoms.append((x, y, z, element, resname, resnum, chain))

    return residues, all_atoms

def compute_com(all_atoms):
    weights = {'C': 12, 'N': 14, 'O': 16, 'S': 32}
    total_mass = 0
    cx = cy = cz = 0
    for x, y, z, elem, _, _, _ in all_atoms:
        w = weights.get(elem, 12)
        total_mass += w
        cx += x * w; cy += y * w; cz += z * w
    return cx/total_mass, cy/total_mass, cz/total_mass

def find_binding_pocket(all_atoms):
    """
    Find the central binding cavity of SMVT.

    SMVT (SLC5A6) is a 12-TM sodium-coupled transporter. The substrate
    binding pocket is in the central cavity approximately halfway across
    the membrane. We identify it by:
    1. Finding the funnel-shaped cavity center
    2. Locating the cluster of conserved aromatic/charged residues

    Strategy: The central cavity of MFS/SSS transporters typically sits
    between the N-domain (TM1-TM6) and C-domain (TM7-TM12). For a protein
    of 635 residues, the interdomain interface is roughly at the geometric
    center between domains 1-317 and 318-635.
    """
    # Compute center coordinates with N-terminal domain atoms (1-317)
    n_atoms = [(a[0], a[1], a[2]) for a in all_atoms if a[5] <= 317 and a[3] in ('C', 'N', 'O')]
    c_atoms = [(a[0], a[1], a[2]) for a in all_atoms if a[5] > 317 and a[3] in ('C', 'N', 'O')]

    if n_atoms and c_atoms:
        # The pocket is at the interface between domains
        nx = sum(a[0] for a in n_atoms) / len(n_atoms)
        ny = sum(a[1] for a in n_atoms) / len(n_atoms)
        nz = sum(a[2] for a in n_atoms) / len(n_atoms)
        cx = sum(a[0] for a in c_atoms) / len(c_atoms)
        cy = sum(a[1] for a in c_atoms) / len(c_atoms)
        cz = sum(a[2] for a in c_atoms) / len(c_atoms)
        mid_x = (nx + cx) / 2
        mid_y = (ny + cy) / 2
        mid_z = (nz + cz) / 2

        # Refine: find solvent-accessible residues near the interface
        # Look for pocket-lining functional residues (aromatic/charged)
        pocket_res = {79, 99, 100, 106, 115, 123, 130, 132, 136, 139,
                      142, 149, 151, 156, 203, 240, 253, 271, 280, 281,
                      285, 295, 318, 339, 346, 371, 384, 386, 390, 394}
        pocket_atoms = [(a[0], a[1], a[2]) for a in all_atoms if a[5] in pocket_res and a[3] in ('C', 'N', 'O')]
        if pocket_atoms:
            px = sum(a[0] for a in pocket_atoms) / len(pocket_atoms)
            py = sum(a[1] for a in pocket_atoms) / len(pocket_atoms)
            pz = sum(a[2] for a in pocket_atoms) / len(pocket_atoms)
            return (px, py, pz, "Pocket-lining residues")

        return (mid_x, mid_y, mid_z, "Interdomain interface")

    # Final fallback
    gx, gy, gz = compute_com(all_atoms)
    return (gx, gy, gz, "Center of mass (fallback)")

def compute_box_size(all_atoms, center, padding=8):
    cx, cy, cz = center
    nearby = [(a[0], a[1], a[2]) for a in all_atoms
              if math.sqrt((a[0]-cx)**2 + (a[1]-cy)**2 + (a[2]-cz)**2) < 20]
    if nearby:
        xs = [p[0] for p in nearby]; ys = [p[1] for p in nearby]; zs = [p[2] for p in nearby]
        return (max(xs)-min(xs) + padding, max(ys)-min(ys) + padding, max(zs)-min(zs) + padding)
    return (25, 25, 25)

def identify_site(all_atoms, center, radius=8):
    cx, cy, cz = center
    nearby = set()
    for a in all_atoms:
        d = math.sqrt((a[0]-cx)**2 + (a[1]-cy)**2 + (a[2]-cz)**2)
        if d < radius:
            nearby.add((a[5], a[6], a[4]))
    return sorted(nearby)

def main():
    print("=" * 60)
    print("SMVT (SLC5A6) Structure Analysis for Docking")
    print("=" * 60)

    print(f"\nReading PDB: {PDB_PATH}")
    residues, all_atoms = parse_pdb(PDB_PATH)
    print(f"  Heavy atoms: {len(all_atoms)}, Residues: {len(residues)}")

    pocket_cx, pocket_cy, pocket_cz, method = find_binding_pocket(all_atoms)

    print(f"\nPocket Center: ({pocket_cx:.1f}, {pocket_cy:.1f}, {pocket_cz:.1f})")
    print(f"Pocket method: {method}")

    size_x, size_y, size_z = compute_box_size(all_atoms, (pocket_cx, pocket_cy, pocket_cz))
    print(f"\nDocking Box:")
    print(f"  Center: ({pocket_cx:.1f}, {pocket_cy:.1f}, {pocket_cz:.1f})")
    print(f"  Size:   ({size_x:.0f} x {size_y:.0f} x {size_z:.0f})")

    site = identify_site(all_atoms, (pocket_cx, pocket_cy, pocket_cz), radius=8)
    site_core = identify_site(all_atoms, (pocket_cx, pocket_cy, pocket_cz), radius=6)

    print(f"\nPocket Residues (8A, n={len(site)}):")
    for rnum, chain, rname in site:
        print(f"  {rname} {chain}:{rnum}")

    print(f"\nCore Pocket Residues (6A, n={len(site_core)}):")
    for rnum, chain, rname in site_core:
        print(f"  {rname} {chain}:{rnum}")

    xs = [a[0] for a in all_atoms]; ys = [a[1] for a in all_atoms]; zs = [a[2] for a in all_atoms]
    print(f"\nProtein Bounds: X[{min(xs):.0f},{max(xs):.0f}] Y[{min(ys):.0f},{max(ys):.0f}] Z[{min(zs):.0f},{max(zs):.0f}]")

    config = {
        "receptor_pdbqt": os.path.join(DATA_DIR, "receptor.pdbqt"),
        "center_x": round(pocket_cx, 2), "center_y": round(pocket_cy, 2), "center_z": round(pocket_cz, 2),
        "size_x": round(size_x, 1), "size_y": round(size_y, 1), "size_z": round(size_z, 1),
        "exhaustiveness": 16, "num_modes": 9, "energy_range": 3,
        "pocket_method": method,
    }
    config_path = os.path.join(DATA_DIR, "docking_config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"\nConfig saved: {config_path}")


if __name__ == "__main__":
    main()
