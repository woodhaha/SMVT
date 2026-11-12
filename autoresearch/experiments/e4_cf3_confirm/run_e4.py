"""E4 — confirm COOH-CF3 -9.02 (3 seeds) + expand the CF3 SAR.

COOH-CF3 (-9.02) beat the -8.84 ceiling in E3. Confirm it's reproducible across Vina
seeds, then explore fluorinated variants to map the CF3 effect:
diCOOH-CF3 (2 COOH + CF3), CHF2, CH2CF3, OCF3, and bis-CF3.
"""
import sys, os, csv
from rdkit import Chem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_smiles, dock_batch

COOH_CF3 = "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)

# ---- 1. confirmation: COOH-CF3 at 3 fixed seeds ----
print("=== COOH-CF3 confirmation (3 seeds) ===")
conf = []
for seed in [1, 42, 2024]:
    r = dock_smiles(COOH_CF3, f"COOH-CF3_seed{seed}", workdir, seed=seed)
    conf.append(r)
    print(f"  seed {seed}: {r['best_affinity']:.2f}")

# ---- 2. fluorinated variants ----
LIGANDS = [
    ("COOH-CF3 (ref)",       "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("diCOOH-CF3",           "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)O)C2"),  # placeholder replaced below
    ("COOH-CHF2",            "OC(=O)C12CC3CC(CC(C3)(C1)C(F)F)C2"),
    ("COOH-CH2CF3",          "OC(=O)C12CC3CC(CC(C3)(C1)CC(F)(F)F)C2"),
    ("COOH-OCF3",            "OC(=O)C12CC3CC(CC(C3)(C1)OC(F)(F)F)C2"),
    ("COOH-bisCF3",          "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),  # placeholder
]

# build diCOOH-CF3 (1,3-diCOOH-5-CF3) and bis-CF3 carefully
LIGANDS[1] = ("diCOOH-CF3", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(=O)O)C2")
LIGANDS[5] = ("COOH-bisCF3", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2")

print("\n=== SMILES validation ===")
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
# also append the 3 confirmation runs
for r in conf:
    rows.append((r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]))

with open(os.path.join(os.path.dirname(__file__), "e4_results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerows(rows)

print("\n=== E4 RESULTS (ranked) ===")
for r in sorted(results, key=lambda x: x["best_affinity"]):
    mark = " ***" if r["best_affinity"] < -8.84 else ""
    print(f"{r['name']:20s} {r['best_affinity']:8.2f}{mark}")

conf_str = [f"{r['best_affinity']:.2f}" for r in conf]
print(f"\nCF3 confirmation (3 seeds): {conf_str}")
print(f"hits < -8.84: {sum(1 for r in results if r['best_affinity'] < -8.84)}/{len(results)}")
