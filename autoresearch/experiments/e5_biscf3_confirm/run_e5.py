"""E5 — confirm COOH-bisCF3 -9.69 (3 seeds) + push further (tri-CF3, diCOOH-bisCF3)."""
import sys, os, csv
from rdkit import Chem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_smiles, dock_batch

BIS_CF3 = "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"
TRI_CF3 = "OC(=O)C12CC3(C(F)(F)F)CC(C(F)(F)F)(CC(C3)(C1)C(F)(F)F)C2"
DI_BIS = "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(=O)O)C2"

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)

print("=== COOH-bisCF3 confirmation (3 seeds) ===")
conf = []
for seed in [1, 42, 2024]:
    r = dock_smiles(BIS_CF3, f"bisCF3_seed{seed}", workdir, seed=seed)
    conf.append(r)
    print(f"  seed {seed}: {r['best_affinity']:.2f}")

print("\n=== push further ===")
LIGANDS = [
    ("bisCF3 (ref)", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("triCF3",        "OC(=O)C12CC3(C(F)(F)F)CC(C(F)(F)F)(CC(C3)(C1)C(F)(F)F)C2"),
    ("diCOOH-bisCF3", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(=O)O)C2"),
]
valid = []
for name, smi in LIGANDS:
    m = Chem.MolFromSmiles(smi)
    print(f"  {name:20s} {'OK' if m else 'FAIL'}")
    if m:
        valid.append((name, smi))

results = dock_batch(valid, workdir)

rows = [("name", "best_affinity", "all_affinities", "ok", "error")]
for r in results:
    rows.append((r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]))
for r in conf:
    rows.append((r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]))
with open(os.path.join(os.path.dirname(__file__), "e5_results.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows)

print("\n=== E5 RESULTS ===")
for r in sorted(results, key=lambda x: x["best_affinity"]):
    mark = " ***" if r["best_affinity"] < -9.0 else ""
    print(f"{r['name']:20s} {r['best_affinity']:8.2f}{mark}")
conf_str = [f"{r['best_affinity']:.2f}" for r in conf]
print(f"\nbisCF3 confirmation (3 seeds): {conf_str}")
