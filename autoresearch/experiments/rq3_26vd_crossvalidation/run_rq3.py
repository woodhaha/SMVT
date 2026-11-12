"""RQ3 — cross-structure validation on cryo-EM 26vd (outward-open state).

Completes the 3-state check (AF2 + 26va occluded + 26vd outward-open).
Question: does the CF3 ranking hold on the outward-open state too?
Receptor: autoresearch/data/26vd_receptor_clean.pdbqt (obabel -xr, cleaned).
Box: pocket-residue centroid (104.99, 103.44, 105.33), size 25^3.
"""
import sys, os, csv, re, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import _to_pdbqt, VINA

RECEPTOR = r"D:\Researching\SMVT\autoresearch\data\26vd_receptor_clean.pdbqt"
CENTER = (104.99, 103.44, 105.33)
SIZE = (25.0, 25.0, 25.0)

LIGANDS = [
    ("biotin",                  "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"),
    ("COOH-bisCF3",             "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("COOH-CF3",                "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("diCOOH-CF3",              "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(=O)O)C2"),
    ("COOH-OCF3",               "OC(=O)C12CC3CC(CC(C3)(C1)OC(F)(F)F)C2"),
    ("COOH-CH2CF3",             "OC(=O)C12CC3CC(CC(C3)(C1)CC(F)(F)F)C2"),
    ("COOH-Ac",                 "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)C)C2"),
    ("adamantane-1,3-diCOOH",   "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)O)C2"),
]


def dock(smiles, name, workdir, seed=42):
    ligand = _to_pdbqt(smiles, workdir)
    if ligand is None:
        return {"name": name, "best_affinity": 99.0, "ok": False, "error": "ligand prep failed"}
    out = os.path.join(workdir, "o.pdbqt")
    cmd = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
           "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
           "--out", out, "--exhaustiveness", "16", "--num_modes", "9", "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    scores = [float(m.group(2)) for m in re.finditer(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", r.stdout, re.M)]
    if not scores:
        return {"name": name, "best_affinity": 99.0, "ok": False, "error": r.stderr[:100]}
    return {"name": name, "best_affinity": min(scores), "ok": True, "error": ""}


workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)
results = [dock(smi, name, os.path.join(workdir, name)) for name, smi in LIGANDS]
results.sort(key=lambda r: r["best_affinity"])

with open(os.path.join(os.path.dirname(__file__), "rq3_results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], r["ok"], r["error"]])

print("\n=== RQ3 26vd (outward-open) cross-validation ===")
for r in results:
    print(f"{r['name']:26s} {r['best_affinity']:8.2f}  {'' if r['ok'] else 'FAIL:'+r['error']}")
