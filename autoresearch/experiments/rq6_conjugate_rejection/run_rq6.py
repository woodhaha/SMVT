"""RQ6 — cross-structure refutation: dock the full COOH-bisCF3-SN38 conjugate (763.7 Da)
against AF2 + cryo-EM 26va (occluded) + 26vd (outward-open). Expected: no productive pose
(conjugate ~22 A extends beyond the ~300 Da pocket). Completes the 3-state refutation table.

Also re-docks the HEAD (COOH-bisCF3, MW 345) as a within-run positive control.
"""
import sys, os, csv, re, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import _to_pdbqt, VINA
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

CONJUGATE = "CCC1=C2C=C(OC(=O)NCC34CC5(C(F)(F)F)CC(CC(C5)(C3)C(F)(F)F)C4C(=O)O)C=CC2=NC2=C1CN1C(=O)C3=C(C=C21)C(CC)(O)COC(=O)3"
HEAD = "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"

RECEPTORS = [
    ("AF2", r"D:\Researching\SMVT\03_Analysis\docking\SMVT_receptor.pdbqt", (5.03, 1.5, -2.44)),
    ("26va", r"D:\Researching\SMVT\autoresearch\data\26va_receptor_clean.pdbqt", (109.09, 109.19, 107.32)),
    ("26vd", r"D:\Researching\SMVT\autoresearch\data\26vd_receptor_clean.pdbqt", (104.99, 103.44, 105.33)),
]


def max_dimension(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=42, useRandomCoords=True, maxAttempts=200)
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    pts = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    return float(d.max())


def dock(smiles, name, receptor, center, workdir, seed=42):
    ligand = _to_pdbqt(smiles, workdir)
    if ligand is None:
        return {"name": name, "best_affinity": None, "ok": False, "error": "ligand prep failed"}
    out = os.path.join(workdir, "o.pdbqt")
    cmd = [VINA, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", "25", "--size_y", "25", "--size_z", "25",
           "--out", out, "--exhaustiveness", "16", "--num_modes", "9", "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    scores = [float(m.group(2)) for m in re.finditer(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", r.stdout, re.M)]
    if not scores:
        return {"name": name, "best_affinity": None, "ok": False, "error": r.stderr[:120]}
    return {"name": name, "best_affinity": min(scores), "ok": True, "error": ""}


base = os.path.dirname(os.path.abspath(__file__))
print(f"conjugate max dimension: {max_dimension(CONJUGATE):.1f} A")
print(f"head max dimension:      {max_dimension(HEAD):.1f} A\n")

rows = []
for rname, rec, center in RECEPTORS:
    for ligand, lname in [(CONJUGATE, "COOH-bisCF3-SN38"), (HEAD, "COOH-bisCF3 (head)")]:
        wd = os.path.join(base, "results", rname, lname)
        r = dock(ligand, lname, rec, center, wd)
        rows.append({"receptor": rname, "name": lname, "best_affinity": r["best_affinity"],
                     "ok": r["ok"], "error": r["error"]})
        print(f"{rname:5s}  {lname:22s}  {str(r['best_affinity']):>10s}  {'OK' if r['ok'] else 'FAIL:'+r['error']}")

out_csv = os.path.join(base, "rq6_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "name", "best_affinity", "ok", "error"])
    w.writeheader()
    w.writerows(rows)
print(f"\nsaved: {out_csv}")
