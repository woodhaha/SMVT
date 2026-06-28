"""
Analyze Vina docking pose for steric clashes between Esketamine and SMVT protein.
Clashes at 150K+ could explain why NVT fails at exactly 200K.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Load protein PDB and extract atom positions
protein_atoms = []
with open(r"D:\Researching\SMVT\SMVT_Colab_MD_Files\AF-Q9Y289-F1.pdb") as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip() if len(line) > 78 else line[12:16].strip()
                protein_atoms.append((elem, np.array([x, y, z])))
            except:
                continue
print(f"Protein: {len(protein_atoms)} atoms")

# Parse Vina PDBQT for ligand pose (first model)
ligand_pos = []
with open(r"D:\Researching\SMVT\SMVT_Colab_MD_Files\ESKETAMINE_docked.pdbqt") as f:
    in_model = False
    for line in f:
        if line.startswith("MODEL") and int(line.split()[1]) == 1:
            in_model = True
        elif line.startswith("ENDMDL") and in_model:
            break
        elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip() if len(line) > 78 else line[12:16].strip()
                ligand_pos.append((elem, np.array([x, y, z])))
            except:
                continue
print(f"Ligand (Vina pose): {len(ligand_pos)} atoms")

# Also get RDKit-generated ligand atoms for comparison
smiles = "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer()

print(f"\n{'='*60}")
print("Steric Clash Analysis: Vina Pose vs SMVT Protein")
print(f"{'='*60}")

# vdW radii (Å) — common values
vdw_radii = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'P': 1.80, 'S': 1.80, 'CL': 1.75, 'Cl': 1.75,
}

def get_vdw(elem):
    return vdw_radii.get(elem.upper(), 1.70)

# Find clashes: any ligand atom < 0.8 * (vdW_i + vdW_j) from any protein atom
# (0.8 = generous threshold; <1.0 = overlap, severe if <0.7)
clashes = []
close_contacts = []
for i, (l_elem, l_pos) in enumerate(ligand_pos):
    l_vdw = get_vdw(l_elem)
    for j, (p_elem, p_pos) in enumerate(protein_atoms):
        p_vdw = get_vdw(p_elem)
        dist = np.linalg.norm(l_pos - p_pos)
        vdw_sum = l_vdw + p_vdw
        ratio = dist / vdw_sum if vdw_sum > 0 else 999
        if ratio < 0.8:  # severe clash
            clashes.append((i, l_elem, j, p_elem, dist, ratio))
        elif ratio < 1.0:  # close contact
            close_contacts.append((i, l_elem, j, p_elem, dist, ratio))

print(f"\nSevere clashes (dist < 0.8×vdW sum): {len(clashes)}")
print(f"Close contacts (0.8×vdW ≤ dist < 1.0×vdW): {len(close_contacts)}")

if clashes:
    print("\n⚠️  SEVERE CLASHES:")
    clashes.sort(key=lambda x: x[4])  # sort by distance
    for i, l_elem, j, p_elem, dist, ratio in clashes[:20]:
        print(f"  Ligand[{i}] {l_elem} — Protein[{j}] {p_elem}: {dist:.2f}Å (ratio={ratio:.2f})")

if close_contacts and not clashes:
    print("\nClose contacts (not severe):")
    close_contacts.sort(key=lambda x: x[4])
    for i, l_elem, j, p_elem, dist, ratio in close_contacts[:10]:
        print(f"  Ligand[{i}] {l_elem} — Protein[{j}] {p_elem}: {dist:.2f}Å (ratio={ratio:.2f})")

if not clashes and not close_contacts:
    print("✓ No steric issues detected.")

# Also check ligand internal geometry
print(f"\nLigand internal distances (min-max):")
internals = []
for i in range(len(ligand_pos)):
    for j in range(i+1, len(ligand_pos)):
        d = np.linalg.norm(ligand_pos[i][1] - ligand_pos[j][1])
        internals.append(d)
print(f"  Min: {min(internals):.3f}Å  Max: {max(internals):.3f}Å  Median: {np.median(internals):.3f}Å")

# Check for zero or near-zero distances (overlapping atoms in ligand itself)
bad = [d for d in internals if d < 0.5]
if bad:
    print(f"  ⚠️ {len(bad)} near-zero distances (<0.5Å) in ligand! Overlapping atoms in Vina pose?")
else:
    print("  ✓ No self-overlaps in ligand")
