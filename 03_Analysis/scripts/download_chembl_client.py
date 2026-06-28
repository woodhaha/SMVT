"""Download FDA-approved drugs using official ChEMBL client (better caching/retry)."""
import os, sys, time, csv

DATA = r"D:\Researching\SMVT\03_Analysis\data"
SAVE_PATH = os.path.join(DATA, "chembl_fda_approved.csv")

from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule

# Get total count
print("Counting approved drugs...")
approved = molecule.filter(max_phase=4, molecule_type="Small molecule")
total = len(approved)
print(f"Total: {total} approved small molecules")

# Download all
all_drugs = []
batch_size = 500
offset = 0

print(f"Downloading in batches of {batch_size}...")
while offset < total:
    for attempt in range(5):
        try:
            batch = molecule.filter(
                max_phase=4,
                molecule_type="Small molecule",
                offset=offset,
                limit=batch_size,
            )
            break
        except Exception as e:
            wait = 2 ** attempt
            print(f"  Retry {attempt+1}: {e} (wait {wait}s)")
            time.sleep(wait)
    else:
        print(f"  SKIP offset {offset}")
        offset += batch_size
        continue

    for mol in batch:
        smiles = None
        if mol.get("molecule_structures"):
            smiles = mol["molecule_structures"].get("canonical_smiles")
        if not smiles:
            continue

        name = (mol.get("pref_name") or "").replace(",", ";")
        props = mol.get("molecule_properties") or {}
        all_drugs.append({
            "chembl_id": mol["molecule_chembl_id"],
            "name": name,
            "smiles": smiles,
            "mw": props.get("full_mwt", ""),
            "alogp": props.get("alogp", ""),
            "ro5": props.get("num_ro5_violations", ""),
            "first_approval": mol.get("first_approval", ""),
        })

    offset += batch_size
    pct = min(100, len(all_drugs) * 100 // total)
    print(f"  [{len(all_drugs)}/{total}] ({pct}%)", flush=True)
    time.sleep(0.5)

# Save
with open(SAVE_PATH, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["chembl_id","name","smiles","mw","alogp","ro5","first_approval"])
    w.writeheader()
    w.writerows(all_drugs)

sz = os.path.getsize(SAVE_PATH)
print(f"\nDone: {len(all_drugs)} drugs | {SAVE_PATH} ({sz/1024:.0f} KB)")
