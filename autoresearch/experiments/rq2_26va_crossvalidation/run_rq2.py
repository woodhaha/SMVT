"""RQ2 — cross-structure validation: re-dock the CF3 hits against cryo-EM 26va (occluded).

Question: does the CF3 enhancement (COOH-bisCF3 -9.68 on AF2) survive on the experimental
cryo-EM structure, or is it AF2-overfitting?
Calibration: biotin re-docks -6.54 (manuscript 26va biotin -6.49).

Receptor: autoresearch/data/26va_receptor_clean.pdbqt  (obabel -xr, cleaned)
Box: center (109.09, 109.19, 107.32) = biotin centroid, size 25^3.
"""
import sys, os, csv, re, subprocess, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import _to_pdbqt, VINA, MK_PREP_LIG

RECEPTOR = r"D:\Researching\SMVT\autoresearch\data\26va_receptor_clean.pdbqt"
CENTER = (109.09, 109.19, 107.32)
SIZE = (25.0, 25.0, 25.0)
EXH = 16
NMODES = 9

# CF3 hits + references (name, SMILES)
LIGANDS = [
    ("biotin (calib -6.49)",      "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"),
    ("COOH-bisCF3",               "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("COOH-CF3",                  "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("diCOOH-CF3",                "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(=O)O)C2"),
    ("COOH-OCF3",                 "OC(=O)C12CC3CC(CC(C3)(C1)OC(F)(F)F)C2"),
    ("COOH-CH2CF3",               "OC(=O)C12CC3CC(CC(C3)(C1)CC(F)(F)F)C2"),
    ("COOH-Ac (prior -8.84)",     "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)C)C2"),
    ("adamantane-1,3-diCOOH",     "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)O)C2"),
    ("1-adamantanecarboxylic",    "OC(=O)C12CC3CC(CC(C3)C1)C2"),
]


def dock26(smiles, name, workdir, seed=42):
    ligand = _to_pdbqt(smiles, workdir)
    if ligand is None:
        return {"name": name, "best_affinity": 99.0, "ok": False, "error": "ligand prep failed"}
    out = os.path.join(workdir, "o.pdbqt")
    cmd = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
           "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
           "--out", out, "--exhaustiveness", str(EXH), "--num_modes", str(NMODES), "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    scores = [float(m.group(2)) for m in re.finditer(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", r.stdout, re.M)]
    if not scores:
        return {"name": name, "best_affinity": 99.0, "ok": False, "error": r.stderr[:100]}
    return {"name": name, "best_affinity": min(scores), "ok": True, "error": ""}


workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)
results = [dock26(smi, name, os.path.join(workdir, name)) for name, smi in LIGANDS]
results.sort(key=lambda r: r["best_affinity"])

out_csv = os.path.join(os.path.dirname(__file__), "rq2_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], r["ok"], r["error"]])

print("\n=== RQ2 26va cross-validation (ranked) ===")
for r in results:
    print(f"{r['name']:28s} {r['best_affinity']:8.2f}  {'' if r['ok'] else 'FAIL:'+r['error']}")
print(f"\nsaved: {out_csv}")
