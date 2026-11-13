#!/usr/bin/env python3
"""
Co-crystal pose RMSD validation: redocked biotin (26va) vs experimental biotin pose (26va).
Same coordinate frame (receptor = 26va, cocrystal BTN = 26va) -> direct heavy-atom RMSD, no alignment.
Uses the already-redocked pose from rq2 (biotin -6.54) and re-docks once cleanly to be safe.
"""
import os, subprocess, re, tempfile, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

BASE = r"D:\Researching\SMVT"
VINA = r"D:\Researching\SMVT\03_Analysis\tools\vina.exe"
OBABEL = r"C:\anaconda3\Scripts\obabel.exe"
MK_PREP_LIG = r"C:\anaconda3\Scripts\mk_prepare_ligand.exe"
RECEPTOR = r"D:\Researching\SMVT\autoresearch\data\26va_receptor_clean.pdbqt"
PDB_26VA = r"D:\Researching\SMVT\05_PDB\26va.pdb"
CENTER = (109.09, 109.19, 107.32)   # biotin centroid in 26va frame (from rq2)
SIZE = (25.0, 25.0, 25.0)
EXH = 16
NMODES = 9
BIOTIN_SMILES = "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"

workdir = os.path.join(BASE, "03_Analysis", "cocrystal_rmsd_out")
os.makedirs(workdir, exist_ok=True)

# ── Step 1: extract co-crystal biotin (BTN) from 26va.pdb ──
cocrystal_pdb = os.path.join(workdir, "biotin_cocrystal.pdb")
with open(PDB_26VA) as f, open(cocrystal_pdb, "w") as g:
    n = 0
    for line in f:
        if line.startswith("HETATM") and "BTN" in line[17:20]:
            g.write(line); n += 1
print(f"[1] cocrystal BTN atoms: {n}")

cocrystal_sdf = os.path.join(workdir, "biotin_cocrystal.sdf")
r = subprocess.run([OBABEL, cocrystal_pdb, "-O", cocrystal_sdf], capture_output=True, text=True)
if r.returncode != 0:
    print("obabel cocrystal fail:", r.stderr[:200]); sys.exit(1)
ref = Chem.MolFromMolFile(cocrystal_sdf, removeHs=False)
ref = Chem.RemoveHs(ref)
print(f"[1] cocrystal mol: {ref.GetNumAtoms()} heavy atoms, formula {Chem.rdMolDescriptors.CalcMolFormula(ref)}")

# ── Step 2: prepare biotin ligand (meeko, same as rq2) + re-dock ──
lig_sdf = os.path.join(workdir, "biotin.sdf")
smiles_mol = Chem.MolFromSmiles(BIOTIN_SMILES)
smiles_mol = Chem.AddHs(smiles_mol)
Chem.AllChem.EmbedMolecule(smiles_mol, randomSeed=42)
Chem.AllChem.MMFFOptimizeMolecule(smiles_mol)
w = Chem.SDWriter(lig_sdf); w.write(smiles_mol); w.close()
lig_pdbqt = os.path.join(workdir, "biotin.pdbqt")
r = subprocess.run([MK_PREP_LIG, "-i", lig_sdf, "-o", lig_pdbqt], capture_output=True, text=True)
if r.returncode != 0:
    print("mk_prepare_ligand fail:", r.stderr[:300]); sys.exit(1)
print("[2] biotin ligand prepared")

out_pdbqt = os.path.join(workdir, "biotin_docked.pdbqt")
cmd = [VINA, "--receptor", RECEPTOR, "--ligand", lig_pdbqt,
       "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
       "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
       "--out", out_pdbqt, "--exhaustiveness", str(EXH), "--num_modes", str(NMODES), "--seed", "42"]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
print("[3] vina stderr tail:", r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "(none)")
scores = [float(m.group(2)) for m in re.finditer(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", r.stdout, re.M)]
print(f"[3] docked {len(scores)} modes, best {min(scores):.2f}")

# ── Step 3: split modes, compute heavy-atom RMSD vs cocrystal ──
def read_models(pdbqt_path):
    models = []
    with open(pdbqt_path) as f:
        cur = []
        for line in f:
            if line.startswith("MODEL"):
                cur = []
            elif line.startswith("ENDMDL"):
                if cur: models.append("".join(cur))
            elif line.startswith(("ATOM", "HETATM")):
                cur.append(line)
    if not models:  # single-model output
        with open(pdbqt_path) as f:
            models = ["".join(l for l in f if l.startswith(("ATOM", "HETATM")))]
    return models

def mol_from_pdbqt_block(block, idx):
    p = os.path.join(workdir, f"mode{idx}.pdbqt")
    with open(p, "w") as f:
        f.write("MODEL 1\n" + block + "ENDMDL\n")
    sdf = os.path.join(workdir, f"mode{idx}.sdf")
    subprocess.run([OBABEL, p, "-O", sdf], capture_output=True, text=True)
    m = Chem.MolFromMolFile(sdf, removeHs=False)
    return Chem.RemoveHs(m) if m else None

def heavy_rmsd(prb, ref):
    match = prb.GetSubstructMatch(ref)  # prb atom idx matching each ref atom
    if not match:
        try:
            return rdMolAlign.GetBestRMS(prb, ref)
        except Exception as e:
            return float("nan")
    rc = ref.GetConformer(); pc = prb.GetConformer()
    d2 = 0.0; n = 0
    for i in range(ref.GetNumAtoms()):
        if ref.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue
        j = match[i]
        d2 += (rc.GetAtomPosition(i) - pc.GetAtomPosition(j)).LengthSq()
        n += 1
    return (d2 / n) ** 0.5 if n else float("nan")

models = read_models(out_pdbqt)
print(f"[4] parsed {len(models)} modes")
rows = []
for i, blk in enumerate(models):
    m = mol_from_pdbqt_block(blk, i)
    if m is None:
        rows.append((i+1, scores[i] if i < len(scores) else float("nan"), float("nan")))
        continue
    rmsd = heavy_rmsd(m, ref)
    sc = scores[i] if i < len(scores) else float("nan")
    rows.append((i+1, sc, rmsd))
    print(f"    mode {i+1}: score {sc:7.2f}  RMSD {rmsd:6.2f} A")

top1 = rows[0][2] if rows else float("nan")
best = min(r[2] for r in rows) if rows else float("nan")
best_mode = min(rows, key=lambda r: r[2])[0] if rows else 0
print("\n=== RESULT ===")
print(f"Top-1 pose RMSD : {top1:.2f} A")
print(f"Best pose RMSD  : {best:.2f} A (mode {best_mode})")
print(f"Success (<2.0 A): top1={'YES' if top1<2.0 else 'NO'}, best={'YES' if best<2.0 else 'NO'}")
