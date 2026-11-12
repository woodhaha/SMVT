"""RQ8 — negative control: low-affinity SMVT head (1-carboxy-3-aminomethyladamantane,
i.e. COOH-bisCF3 minus the two CF3 groups) conjugated to SN-38, vs COOH-bisCF3-SN38 and
biotin-SN38. Isolates the head-affinity variable: same cage, same free carboxyl, same
carbamate handle, only the CF3 affinity enhancers removed.

Low-affinity head + conjugate derived by deleting the two CF3 groups from the bisCF3 series.
"""
import sys, os, csv, re, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import _to_pdbqt, VINA
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
import numpy as np

# low-affinity head = 1-carboxy-3-aminomethyladamantane (free amine)
HEAD = "OC(=O)C12CC3CC(CC(C3)(C1))C2(CN)"
# low-affinity conjugate = SN38 carbamate to the aminomethyl (bisCF3-SN38 minus 2 CF3)
CONJ = "CCC1=C2C=C(OC(=O)NCC34CC5CC(CC(C5)(C3))C4C(=O)O)C=CC2=NC2=C1CN1C(=O)C3=C(C=C21)C(CC)(O)COC(=O)3"

RECEPTORS = [
    ("AF2", r"D:\Researching\SMVT\03_Analysis\docking\SMVT_receptor.pdbqt", (5.03, 1.5, -2.44)),
    ("26va", r"D:\Researching\SMVT\autoresearch\data\26va_receptor_clean.pdbqt", (109.09, 109.19, 107.32)),
    ("26vd", r"D:\Researching\SMVT\autoresearch\data\26vd_receptor_clean.pdbqt", (104.99, 103.44, 105.33)),
]


def check(name, smi):
    m = Chem.MolFromSmiles(smi)
    print(f"{name}: {rdMolDescriptors.CalcMolFormula(m)}  MW {rdMolDescriptors.CalcExactMolWt(m):.1f}")
    return Chem.MolToSmiles(m)


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


head_smi = check("low-aff head (1-COOH-3-CH2NH2-adam)", HEAD)
conj_smi = check("low-aff conjugate", CONJ)
print(f"head max dim: {max_dimension(head_smi):.1f} A   conjugate max dim: {max_dimension(conj_smi):.1f} A\n")

base = os.path.dirname(os.path.abspath(__file__))
rows = []
for rname, rec, center in RECEPTORS:
    for smi, nm in [(head_smi, "low-aff head"), (conj_smi, "low-aff-SN38")]:
        r = dock(smi, nm, rec, center, os.path.join(base, "results", rname, nm))
        rows.append({"receptor": rname, "name": nm, "best_affinity": r["best_affinity"], "ok": r["ok"]})
        print(f"{rname:5s}  {nm:16s}  {str(r['best_affinity']):>10s}  {'OK' if r['ok'] else 'FAIL:'+r['error']}")

with open(os.path.join(base, "rq8_results.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "name", "best_affinity", "ok"])
    w.writeheader()
    w.writerows(rows)
print(f"\nsaved: {os.path.join(base, 'rq8_results.csv')}")
