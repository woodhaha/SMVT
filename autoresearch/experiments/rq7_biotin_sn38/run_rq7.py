"""RQ7 — biotin-SN38 vs COOH-bisCF3-SN38: does esterifying biotin's carboxyl (consuming the
free-carboxyl anchor) change rejection? Dock biotin-SN38 on AF2 + 26va + 26vd.

Biotin-SN38 built by esterification (biotin acid + SN-38 C10 phenol) via RDKit reaction.
Compare to RQ6 (bisCF3-SN38 retains free carboxyl; rejected by size only).
"""
import sys, os, csv, re, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import _to_pdbqt, VINA
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import numpy as np

BIOTIN = "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"
SN38 = "CCC1=C2C=C(O)C=CC2=NC2=C1CN1C(=O)C3=C(C=C21)C(CC)(O)COC(=O)3"

RECEPTORS = [
    ("AF2", r"D:\Researching\SMVT\03_Analysis\docking\SMVT_receptor.pdbqt", (5.03, 1.5, -2.44)),
    ("26va", r"D:\Researching\SMVT\autoresearch\data\26va_receptor_clean.pdbqt", (109.09, 109.19, 107.32)),
    ("26vd", r"D:\Researching\SMVT\autoresearch\data\26vd_receptor_clean.pdbqt", (104.99, 103.44, 105.33)),
]

# esterification: biotin acid + SN38 C10 phenol -> ester
rxn = AllChem.ReactionFromSmarts("[C:1](=[O:2])[O:3].[O:4][c:5]>>[C:1](=[O:2])[O:4][c:5]")
prods = rxn.RunReactants((Chem.MolFromSmiles(BIOTIN), Chem.MolFromSmiles(SN38)))
m = prods[0][0]
Chem.SanitizeMol(m)
biotin_sn38 = Chem.MolToSmiles(m)
formula = rdMolDescriptors.CalcMolFormula(m)
mw = rdMolDescriptors.CalcExactMolWt(m)
print(f"biotin-SN38 formula: {formula}  MW {mw:.1f}")
print(f"biotin-SN38 SMILES:  {biotin_sn38}")
assert formula == "C32H34N4O7S", f"formula mismatch: {formula}"
print("formula OK (expect C32H34N4O7S)\n")


def max_dimension(smiles):
    mm = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mm, randomSeed=42, useRandomCoords=True, maxAttempts=200)
    AllChem.MMFFOptimizeMolecule(mm)
    c = mm.GetConformer()
    pts = np.array([list(c.GetAtomPosition(i)) for i in range(mm.GetNumAtoms())])
    return float(np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1).max())


def dock(smiles, name, receptor, center, workdir, seed=42):
    ligand = _to_pdbqt(smiles, workdir)
    if ligand is None:
        return {"name": name, "best_affinity": None, "ok": False, "error": "prep failed"}
    out = os.path.join(workdir, "o.pdbqt")
    cmd = [VINA, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", "25", "--size_y", "25", "--size_z", "25",
           "--out", out, "--exhaustiveness", "16", "--num_modes", "9", "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    scores = [float(mm.group(2)) for mm in re.finditer(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", r.stdout, re.M)]
    if not scores:
        return {"name": name, "best_affinity": None, "ok": False, "error": r.stderr[:120]}
    return {"name": name, "best_affinity": min(scores), "ok": True, "error": ""}


print(f"biotin-SN38 max dimension: {max_dimension(biotin_sn38):.1f} A\n")
base = os.path.dirname(os.path.abspath(__file__))
rows = []
for rname, rec, center in RECEPTORS:
    r = dock(biotin_sn38, "biotin-SN38", rec, center, os.path.join(base, "results", rname))
    rows.append({"receptor": rname, "best_affinity": r["best_affinity"], "ok": r["ok"], "error": r["error"]})
    print(f"{rname:5s}  biotin-SN38   {str(r['best_affinity']):>10s}  {'OK' if r['ok'] else 'FAIL:'+r['error']}")

with open(os.path.join(base, "rq7_results.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "best_affinity", "ok", "error"])
    w.writeheader()
    w.writerows(rows)
print(f"\nsaved: {os.path.join(base, 'rq7_results.csv')}")
