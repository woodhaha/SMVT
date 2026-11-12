"""RQ4 — MM-GBSA (fixed atom mapping). meeko SMILES-IDX -> RDKit order, re-add H.

Key fix: meeko reorders atoms + strips non-polar H; gen_gaff_v3 template uses RDKit
AddHs order (all H). We map the docked heavy-atom coordinates back to RDKit order via
the meeko PDBQT "REMARK SMILES IDX" remark, then re-embed H with heavy atoms fixed.
"""
import os, sys, re, subprocess, tempfile

SMVT = r"D:\Researching\SMVT"
sys.path.insert(0, os.path.join(SMVT, "SMVT_MD", "scripts"))
from gen_gaff_v3 import generate_template
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry

VINA = r"D:\Researching\SMVT\03_Analysis\tools\vina.exe"
MK_PREP_LIG = r"C:\anaconda3\Scripts\mk_prepare_ligand.exe"
RECEPTOR = os.path.join(SMVT, "MM-GBSA Analysis", "receptor_clean.pdbqt")
CENTER = (3.9, 2.0, -7.6)
BOX = 22
ANALYSIS = os.path.join(SMVT, "MM-GBSA Analysis")
LIGDIR = os.path.join(SMVT, "SMVT_MD", "ligands")
PROT_SRC = os.path.join(ANALYSIS, "BIOTIN_complex.pdb")

COMPOUNDS = [
    ("BISCF3", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("CF3",    "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("AC",     "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)C)C2"),
    ("BIOTIN", "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"),
]


def prep_and_dock(smiles, workdir):
    """meeko prep -> vina dock. Returns (meeko_pdbqt, vina_out) paths."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) != 0:
        return None, None
    AllChem.MMFFOptimizeMolecule(mol)
    sdf = os.path.join(workdir, "l.sdf")
    Chem.MolToMolFile(mol, sdf)
    pdbqt = os.path.join(workdir, "l.pdbqt")
    r = subprocess.run([MK_PREP_LIG, "-i", sdf, "-o", pdbqt], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(pdbqt):
        return None, None
    out = os.path.join(workdir, "docked.pdbqt")
    cmd = [VINA, "--receptor", RECEPTOR, "--ligand", pdbqt,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
           "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
           "--out", out, "--exhaustiveness", "16", "--num_modes", "9", "--seed", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return (pdbqt, out) if r.returncode == 0 and os.path.exists(out) else (None, None)


def build_ligand_lines(smiles, meeko_pdbqt, vina_out, n_prot=0):
    """Map docked pose to RDKit AddHs order, re-embed H, return (PDB lines, bonds)."""
    # 1. parse meeko "SMILES IDX": smiles_atom(1-based) -> meeko_pdbqt_atom(1-based)
    smi2mek = {}
    for ln in open(meeko_pdbqt):
        if "SMILES IDX" in ln:
            nums = [int(x) for x in ln.split()[3:]]
            for i in range(0, len(nums), 2):
                smi2mek[nums[i]] = nums[i + 1]
    # 2. docked MODEL 1 coords: meeko_atom(1-based) -> (x,y,z)
    mek_coord = {}
    in_m = False
    for ln in open(vina_out):
        if ln.startswith("MODEL") and ln.split()[1] == "1":
            in_m = True
            continue
        if ln.startswith("ENDMDL"):
            break
        if in_m and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            mek_coord[int(ln[6:11])] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    # 3. build heavy-only mol, set docked heavy coords, then AddHs (addCoords) to place H
    heavy = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(heavy.GetNumAtoms())
    n_set = 0
    for smi_idx, mek_idx in smi2mek.items():
        if mek_idx in mek_coord:
            x, y, z = mek_coord[mek_idx]
            conf.SetAtomPosition(smi_idx - 1, Geometry.Point3D(x, y, z))
            n_set += 1
    if n_set < 5:
        return None, None
    heavy.AddConformer(conf)
    mol = Chem.AddHs(heavy, addCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    c = mol.GetConformer()
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        name = f"{atom.GetSymbol()}{i + 1}"
        p = c.GetAtomPosition(i)
        lines.append(f"HETATM{n_prot + i + 1:5d} {name:<4s} LIG B   1    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {atom.GetSymbol():>2s}")
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    return lines, bonds


def main():
    # templates
    for name, smi in COMPOUNDS:
        xml, n, tc = generate_template(smi, "LIG")
        with open(os.path.join(LIGDIR, f"{name}_template.xml"), "w") as f:
            f.write(xml)
        print(f"template {name}: {n} atoms")

    # dock + build complex (correct atom mapping + CONECT bonds)
    prot_lines = [l for l in open(PROT_SRC) if l.startswith("ATOM")]
    n_prot = len(prot_lines)
    for name, smi in COMPOUNDS:
        wd = tempfile.mkdtemp(prefix="mmgbsa_")
        mek, out = prep_and_dock(smi, wd)
        if not out:
            print(f"{name}: DOCK FAILED")
            continue
        lig, bonds = build_ligand_lines(smi, mek, out, n_prot)
        if not lig:
            print(f"{name}: MAPPING FAILED")
            continue
        # build CONECT for ligand (absolute serials = n_prot + lig_idx + 1)
        adj = {}
        for i, j in bonds:
            adj.setdefault(i, []).append(j)
            adj.setdefault(j, []).append(i)
        conect = []
        for i in sorted(adj):
            conect.append("CONECT" + "".join(f"{n_prot + x + 1:5d}" for x in [i] + sorted(adj[i])))
        with open(os.path.join(ANALYSIS, f"{name}_complex.pdb"), "w") as f:
            f.writelines(prot_lines)
            for l in lig:
                f.write(l + "\n")
            for cl in conect:
                f.write(cl + "\n")
            f.write("END\n")
        print(f"{name}: complex built ({len(lig)} lig atoms, {len(conect)} CONECT)")

    # MM-GBSA
    print("\n=== MM-GBSA ===")
    sys.path.insert(0, ANALYSIS)
    import _run_mmgbsa as mg
    for name, _ in COMPOUNDS:
        try:
            mg.run(name)
        except Exception as e:
            print(f"{name}: MM-GBSA ERROR {str(e)[:150]}")


if __name__ == "__main__":
    main()
