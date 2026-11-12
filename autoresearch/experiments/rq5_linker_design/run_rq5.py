"""RQ5 — nanoparticle targeting head: can the 4th bridgehead (7-position) of COOH-bisCF3
carry a linker (for NP conjugation) without losing SMVT binding?

Template: OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(L)  — L at bridgehead 7.
Target: find linkers that keep score < -8.5 (preserving the -9.68 bisCF3 affinity).
"""
import sys, os, csv
from rdkit import Chem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_batch

BASE = "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"

LIGANDS = [
    ("base bisCF3 (ref -9.68)", BASE),
    # smallest control
    ("7-CH3",        "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(C)"),
    # polar linkers (nanoparticle conjugation handles)
    ("7-CH2OH",      "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CO)"),
    ("7-CH2NH2",     "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CN)"),
    ("7-CH2N3",      "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CN=[N+]=[N-])"),
    ("7-CH2SH",      "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CS)"),
    ("7-CC#C",       "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CC#C)"),
    # extended linkers
    ("7-CH2CH2OH",   "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CCO)"),
    ("7-CH2CH2NH2",  "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CCN)"),
    ("7-CH2CH2CH2NH2","OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CCCN)"),
    # PEG + second carboxyl comparison
    ("7-CH2OCH2CH2NH2","OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(COCCN)"),
    ("7-COOH (diCOOH)","OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(C(=O)O)"),
]

print("=== SMILES validation ===")
valid = []
for name, smi in LIGANDS:
    m = Chem.MolFromSmiles(smi)
    print(f"  {name:26s} {'OK' if m else 'FAIL'}")
    if m:
        valid.append((name, smi))

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)
results = dock_batch(valid, workdir)

with open(os.path.join(os.path.dirname(__file__), "rq5_results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "all_affinities", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]])

print("\n=== RQ5 linker design (ranked) ===")
for r in results:
    keep = " ***" if r["best_affinity"] < -8.5 else ""
    print(f"{r['name']:26s} {r['best_affinity']:8.2f}{keep}  {'' if r['ok'] else 'FAIL:'+r['error']}")

print(f"\nlinkers preserving < -8.5: {sum(1 for r in results if r['best_affinity'] < -8.5)}/{len(results)}")
