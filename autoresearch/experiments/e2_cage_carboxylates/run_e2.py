"""E2 — diamondoid/bicyclic cage carboxylates (H4).

Aromatic chemotype capped at ~-7.7 (E1). The adamantane cage is where -8.84 lives.
Test NEW cage scaffolds (noradamantane, bicyclo[2.2.2]octane, cubane, homoadamantane)
+ small-EWG adamantane variants, targeting < -8.84. Internal calibration: adamantane
-1,3,5-tricarboxylic acid (manuscript -8.51).
"""
import sys, os, csv
from rdkit import Chem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_batch

LIGANDS = [
    # ---- internal calibration (known baseline) ----
    ("adamantane-1,3,5-triCOOH", "OC(=O)C12CC3CC(C(=O)O)(CC(C(=O)O)(C3)C1)C2"),
    ("1-adamantanecarboxylic",   "OC(=O)C12CC3CC(CC(C3)C1)C2"),
    ("adamantane-1,3-diCOOH",    "OC(=O)C12CC3CC(CC(C3)(C1)C(=O)O)C2"),
    # ---- small-EWG adamantane (push below -8.84) ----
    ("3-methyl-adamantane-1-COOH",  "OC(=O)C12CC3CC(CC(C)(C3)C1)C2"),
    ("3,5-dimethyl-adam-1-COOH",    "OC(=O)C12CC3CC(C)(CC(C)(C3)C1)C2"),
    ("1-adamantaneacetic-acid",     "OC(=O)CC12CC3CC(CC(C3)C1)C2"),
    ("3,5,7-trimethyl-adam-COOH",   "OC(=O)C12CC3CC(C)(CC(C)(C3)(C1)C)C2"),
    # ---- new cages (H4) ----
    ("bicyclo[2.2.2]octane-1-COOH",      "OC(=O)C12CCC(CC1)CC2"),
    ("bicyclo[2.2.2]octane-1,4-diCOOH",  "OC(=O)C12CCC(CC1)(CC2)C(=O)O"),
    ("cubane-1-COOH",                    "OC(=O)C12C3C4C1C5C2C3C45"),
    ("cubane-1,4-diCOOH",                "OC(=O)C12C3C4C1C5C2C3C45C(=O)O"),
    ("noradamantane-1-COOH",             "OC(=O)C12CC3CC(C1)C(C2)C3"),
    ("homoadamantane-1-COOH",            "OC(=O)C12CC3CC(CC(C3)CC1)C2"),
    # ---- naphthyl optimization (push E1 winner -7.74) ----
    ("1-naphthylalanine",      "NC(Cc1cccc2ccccc12)C(=O)O"),
    ("4-F-1-naphthylalanine",  "NC(Cc1ccc(F)c2ccccc12)C(=O)O"),
    ("1-naphthylglycine",      "NC(C(=O)O)c1cccc2ccccc12"),
    ("1-naphthaleneacetic",    "OC(=O)Cc1cccc2ccccc12"),
    ("2-naphthaleneacetic",    "OC(=O)Cc1ccc2ccccc2c1"),
]

# pre-validate SMILES
print("=== SMILES validation ===")
valid = []
for name, smi in LIGANDS:
    m = Chem.MolFromSmiles(smi)
    tag = "OK" if m else "FAIL"
    if m:
        valid.append((name, smi))
    print(f"  {name:28s} {tag}")

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)
results = dock_batch(valid, workdir)

out_csv = os.path.join(os.path.dirname(__file__), "e2_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "all_affinities", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]])

print("\n=== E2 RESULTS (ranked, more negative = better) ===")
for r in results:
    mark = " ***" if r["best_affinity"] < -8.0 else ""
    print(f"{r['name']:28s} {r['best_affinity']:8.2f}{mark}  {'FAIL:'+r['error'] if not r['ok'] else ''}")

print(f"\nsaved: {out_csv}")
print(f"hits < -8.0: {sum(1 for r in results if r['best_affinity'] < -8.0)}/{len(results)}")
print(f"hits < -8.84: {sum(1 for r in results if r['best_affinity'] < -8.84)}/{len(results)}")
