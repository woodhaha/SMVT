"""
Download FDA-approved drugs from ChEMBL API with retry + resume.
"""
import os, sys, ssl, json, time
from urllib.request import Request, urlopen

ssl_ctx = ssl._create_unverified_context()
DATA = r"D:\Researching\SMVT\03_Analysis\data"
SAVE_PATH = os.path.join(DATA, "chembl_fda_approved.csv")
os.makedirs(DATA, exist_ok=True)

# Resume: load already-downloaded
existing = set()
if os.path.exists(SAVE_PATH):
    with open(SAVE_PATH, "r") as f:
        next(f)  # skip header
        for line in f:
            cid = line.split(",")[0]
            if cid:
                existing.add(cid)
print(f"Already have: {len(existing)} drugs")

BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PARAMS = "max_phase=4&molecule_type=Small%20molecule&limit=100"
HEADERS = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}

all_drugs = list(existing)
total = 3475
offset = 0
max_retries = 5

f = open(SAVE_PATH, "a", encoding="utf-8")
if not existing:
    f.write("chembl_id,name,smiles,mw,alogp,ro5_violations,first_approval\n")

while offset < total:
    url = f"{BASE}?{PARAMS}&offset={offset}"

    for attempt in range(max_retries):
        try:
            req = Request(url, headers=HEADERS)
            resp = urlopen(req, timeout=90, context=ssl_ctx)
            data = json.loads(resp.read())
            break
        except Exception as e:
            wait = 2 ** attempt
            print(f"\n  Retry {attempt+1}/{max_retries} offset={offset}: {e} (wait {wait}s)")
            time.sleep(wait)
    else:
        print(f"\n  SKIP offset {offset} after {max_retries} retries")
        offset += 100
        continue

    molecules = data.get("molecules", [])
    new_count = 0
    for mol in molecules:
        cid = mol.get("molecule_chembl_id", "")
        if cid in all_drugs:
            continue
        structures = mol.get("molecule_structures") or {}
        smiles = structures.get("canonical_smiles", "")
        if not smiles:
            continue
        name = mol.get("pref_name", "").replace(",", ";")
        props = mol.get("molecule_properties") or {}
        mw = props.get("full_mwt", "")
        alogp = props.get("alogp", "")
        ro5 = props.get("num_ro5_violations", "")
        fa = mol.get("first_approval", "")
        f.write(f'{cid},{name},{smiles},{mw},{alogp},{ro5},{fa}\n')
        all_drugs.append(cid)
        new_count += 1

    pct = min(100, int(len(all_drugs) / total * 100))
    print(f"  [{len(all_drugs)}/{total}] +{new_count} | page {offset//100+1}/35 ({pct}%)", flush=True)

    if len(molecules) < 100:
        break
    offset += 100
    time.sleep(0.2)

f.close()
sz = os.path.getsize(SAVE_PATH)
print(f"\nDone: {len(all_drugs)} drugs | File: {SAVE_PATH} ({sz/1024:.0f} KB)")
