"""E1 — aromatic carboxyl-bearing analogs (H3: ring extension; partial H1: EWG).

Dock ~28 aromatic amino-acid / indole-carboxyl analogs against the SMVT AF2 pocket.
Prediction: F/CF3/CN-substituted and naphthalene/benzofuran-extended analogs beat
5-HTP (-7.66), some approaching -8.0.
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from dock import dock_batch

LIGANDS = [
    # ---- baselines / controls (known scores) ----
    ("5-HTP",            "NC(Cc1c[nH]c2ccc(O)cc12)C(=O)O"),
    ("Tryptophan",       "NC(Cc1c[nH]c2ccccc12)C(=O)O"),
    ("Phenylalanine",    "NC(Cc1ccccc1)C(=O)O"),
    ("Tyrosine",         "NC(Cc1ccc(O)cc1)C(=O)O"),
    ("Biotin",           "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"),
    # ---- H1: EWG-substituted tryptophans ----
    ("5-F-tryptophan",   "NC(Cc1c[nH]c2ccc(F)cc12)C(=O)O"),
    ("5-Cl-tryptophan",  "NC(Cc1c[nH]c2ccc(Cl)cc12)C(=O)O"),
    ("5-Br-tryptophan",  "NC(Cc1c[nH]c2ccc(Br)cc12)C(=O)O"),
    ("5-CF3-tryptophan", "NC(Cc1c[nH]c2ccc(C(F)(F)F)cc12)C(=O)O"),
    ("5-CN-tryptophan",  "NC(Cc1c[nH]c2ccc(C#N)cc12)C(=O)O"),
    ("6-F-tryptophan",   "NC(Cc1c[nH]c2cc(F)ccc12)C(=O)O"),
    ("7-F-tryptophan",   "NC(Cc1c[nH]c2cccc(F)c12)C(=O)O"),
    # ---- H3: ring-extended aromatic amino acids ----
    ("1-naphthylalanine",  "NC(Cc1cccc2ccccc12)C(=O)O"),
    ("2-naphthylalanine",  "NC(Cc1ccc2ccccc2c1)C(=O)O"),
    ("benzofuran-alanine", "NC(Cc1coc2ccccc12)C(=O)O"),
    ("benzothiophene-ala", "NC(Cc1csc2ccccc12)C(=O)O"),
    ("4-F-phenylalanine",  "NC(Cc1ccc(F)cc1)C(=O)O"),
    ("4-CF3-phenylalanine","NC(Cc1ccc(C(F)(F)F)cc1)C(=O)O"),
    ("3,5-diF-phenylala",  "NC(Cc1cc(F)cc(F)c1)C(=O)O"),
    ("3-F-tyrosine",       "NC(Cc1ccc(O)c(F)c1)C(=O)O"),
    # ---- indole-carboxyl (auxin-like, no amine) ----
    ("indole-3-acetic",      "OC(=O)Cc1c[nH]c2ccccc12"),
    ("indole-3-propionic",   "OC(=O)CCc1c[nH]c2ccccc12"),
    ("indole-3-butyric",     "OC(=O)CCCc1c[nH]c2ccccc12"),
    ("5-OH-indoleacetic",    "OC(=O)Cc1c[nH]c2ccc(O)cc12"),
    ("5-F-indoleacetic",     "OC(=O)Cc1c[nH]c2ccc(F)cc12"),
    ("6-Cl-indoleacetic",    "OC(=O)Cc1c[nH]c2cc(Cl)ccc12"),
    # ---- multi-carboxyl (H2 probe) ----
    ("xanthurenic_acid",     "OC(=O)c1cc(O)c2[nH]c(C(=O)O)cc2c1"),
    ("kynurenic_acid",       "OC(=O)c1c[nH]c2ccccc2c1=O"),
]

workdir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(workdir, exist_ok=True)

results = dock_batch(LIGANDS, workdir)

out_csv = os.path.join(os.path.dirname(__file__), "e1_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "best_affinity", "all_affinities", "ok", "error"])
    for r in results:
        w.writerow([r["name"], r["best_affinity"], str(r["all_affinities"]), r["ok"], r["error"]])

print("\n=== E1 RESULTS (ranked, more negative = better) ===")
for r in results:
    mark = " ***" if r["best_affinity"] < -8.0 else ""
    print(f"{r['name']:24s} {r['best_affinity']:8.2f}{mark}  {'FAIL:'+r['error'] if not r['ok'] else ''}")

print(f"\nsaved: {out_csv}")
print(f"hits < -8.0: {sum(1 for r in results if r['best_affinity'] < -8.0)}/{len(results)}")
