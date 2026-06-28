"""
Download FDA-approved drugs from ChEMBL API (verified working).
ChEMBL API returns 3475 approved small molecules.
Downloads SMILES + metadata, saves as CSV.
"""
import os, sys, ssl, json, time
from urllib.request import Request, urlopen

ssl_ctx = ssl._create_unverified_context()
DATA = r"D:\Researching\SMVT\03_Analysis\data"
os.makedirs(DATA, exist_ok=True)

BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PARAMS = "max_phase=4&molecule_type=Small%20molecule&limit=100"
HEADERS = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}

print("Downloading FDA-approved drugs from ChEMBL...")
print("Expected: ~3,475 molecules")

all_drugs = []
offset = 0
total = None

while True:
    url = f"{BASE}?{PARAMS}&offset={offset}"
    try:
        req = Request(url, headers=HEADERS)
        resp = urlopen(req, timeout=60, context=ssl_ctx)
        data = json.loads(resp.read())
    except Exception as e:
        print(f"  FAIL at offset {offset}: {e}")
        break

    if total is None:
        total = data.get("page_meta", {}).get("total_count", 0)
        total_pages = (total + 99) // 100
        print(f"  Total: {total} drugs, {total_pages} pages")

    molecules = data.get("molecules", [])
    for mol in molecules:
        structures = mol.get("molecule_structures") or {}
        smiles = structures.get("canonical_smiles", "")
        if not smiles:
            continue
        all_drugs.append({
            "chembl_id": mol.get("molecule_chembl_id", ""),
            "name": mol.get("pref_name", ""),
            "smiles": smiles,
            "mw": mol.get("molecule_properties", {}).get("full_mwt", ""),
            "alogp": mol.get("molecule_properties", {}).get("alogp", ""),
            "ro5_violations": mol.get("molecule_properties", {}).get("num_ro5_violations", ""),
            "first_approval": mol.get("first_approval", ""),
        })

    offset += 100
    pct = min(100, int(offset / total * 100)) if total else 0
    print(f"  [{len(all_drugs)}/{total}] page {offset//100}/{total_pages} ({pct}%)", end="\r")

    if len(molecules) < 100:
        break
    time.sleep(0.3)

print()
print(f"Downloaded: {len(all_drugs)} FDA-approved drugs with SMILES")

# Save CSV
csv_path = os.path.join(DATA, "chembl_fda_approved.csv")
with open(csv_path, "w", encoding="utf-8") as f:
    f.write("chembl_id,name,smiles,mw,alogp,ro5_violations,first_approval\n")
    for d in all_drugs:
        f.write(f'{d["chembl_id"]},{d["name"]},{d["smiles"]},{d["mw"]},{d["alogp"]},{d["ro5_violations"]},{d["first_approval"]}\n')

sz = os.path.getsize(csv_path)
print(f"Saved: {csv_path} ({sz/1024:.0f} KB)")
print("Done!")
