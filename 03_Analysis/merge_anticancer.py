#!/usr/bin/env python3
"""Merge anticancer library results:
1. Extract anticancer compounds from the 10569-combined run (overlap).
2. Add the dedicated anticancer docking run (new compounds).
3. Apply perfluoro false-positive filter (nF>=6 & F/heavy>=0.4).
4. Rank and save.
"""
import json, os
import pandas as pd
from rdkit import Chem

BASE = "/root/docking"
LIB = f"{BASE}/anticancer_library_full.csv"
CK_MAIN = f"{BASE}/results.csv.checkpoint.json"
CK_AC = f"{BASE}/anticancer_results.csv.checkpoint.json"
OUT = f"{BASE}/anticancer_combined_results.csv"

lib = pd.read_csv(LIB, dtype=str)
lib_map = dict(zip(lib["inchikey"], lib["smiles"]))
print(f"anticancer library: {len(lib_map)}")

rows = []
seen = set()

# 1. extract overlap from 10569 main run
with open(CK_MAIN) as f:
    ck = json.load(f)
n_overlap = 0
for r in ck.get("results", []):
    k = r["key"]
    if k in lib_map and k not in seen:
        rows.append({"name": r["name"], "smiles": lib_map[k],
                     "best_affinity": r["best_affinity"], "inchikey": k, "source": "10569_overlap"})
        seen.add(k); n_overlap += 1
print(f"overlap extracted from 10569: {n_overlap}")

# 2. add dedicated anticancer run
if os.path.exists(CK_AC):
    with open(CK_AC) as f:
        ck2 = json.load(f)
    for r in ck2.get("results", []):
        k = r["key"]
        if k in seen:
            continue
        rows.append({"name": r["name"], "smiles": lib_map.get(k, ""),
                     "best_affinity": r["best_affinity"], "inchikey": k, "source": "anticancer_new"})
        seen.add(k)
    print(f"anticancer new results: {len(ck2.get('results', []))}")
else:
    print("anticancer checkpoint not found yet")

# 3. perfluoro filter
def nF_nHeavy(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return 0, 0
    nF = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "F")
    nH = sum(1 for a in m.GetAtoms() if a.GetSymbol() != "H")
    return nF, nH

df = pd.DataFrame(rows)
df["nF"], df["nHeavy"] = zip(*df["smiles"].map(nF_nHeavy))
df["perfluoro"] = (df["nF"] >= 6) & (df["nF"] / df["nHeavy"].clip(lower=1) >= 0.4)
print(f"perfluoro-flagged (excluded): {len(df[df['perfluoro']])}")
df_clean = df[~df["perfluoro"]].sort_values("best_affinity")
df_clean.to_csv(OUT, index=False)
print(f"final anticancer results (clean): {len(df_clean)} -> {OUT}")
for _, r in df_clean.head(10).iterrows():
    print(f"  {r['best_affinity']:>7}  {str(r['name'])[:35]:35s}  {r['source']}")
