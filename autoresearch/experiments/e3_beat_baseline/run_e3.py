"""E3 — final push to beat adamantane-COOH-Ac -8.84.

The prior focused SAR found COOH-Ac (-8.84) and diCOOH-F/CN/Me (-8.63 to -8.67).
E3 tests the untested combos: carboxyl + acetyl + small group, and carboxyl + iPr/cPr
(not in the prior Me/Et/F/CN/Ac/OH partner set). Target: < -8.84.
"""
import sys, os, csv
from rdkit import Chem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_batch

LIGANDS = [
    # ---- reproduce the -8.84 baseline (1-COOH-3-Ac-adamantane) ----
    ("COOH-Ac (baseline -8.84)", "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)C)C2"),
    ("adamantane-1,3-diCOOH (ref)", "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)O)C2"),
    # ---- untested small partners at the C3 bridgehead ----
    ("COOH-iPr",   "OC(=O)C12CC3CC(CC(C3)(C1)C(C)C)C2"),
    ("COOH-cPr",   "OC(=O)C12CC3CC(CC(C3)(C1)C1CC1)C2"),
    ("COOH-CF3",   "OC(=O)C12CC3CC(CC(C3)(C1)C(F)(F)F)C2"),
    ("COOH-CN",    "OC(=O)C12CC3CC(CC(C3)(C1)C#N)C2"),
    ("COOH-CH2Ac", "OC(=O)C12CC3CC(CC(C3)(C1)CC(=O)C)C2"),
]

print("=== SMILES validation ===")
valid = []
for name, smi in LIGANDS:
    m = Chem.MolFromSmiles(smi)
    print(f"  {name:30s} {'OK' if m else 'FAIL'}")
    if m:
        valid.append((name, smi))

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)
results = dock_batch(valid, workdir)

out_csv = os.path.join(os.path.dirname(__file__), "e3_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "all_affinities", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]])

print("\n=== E3 RESULTS (ranked) ===")
for r in results:
    mark = " ***" if r["best_affinity"] < -8.84 else ""
    print(f"{r['name']:30s} {r['best_affinity']:8.2f}{mark}  {'FAIL:'+r['error'] if not r['ok'] else ''}")

print(f"\nhits < -8.84: {sum(1 for r in results if r['best_affinity'] < -8.84)}/{len(results)}")
