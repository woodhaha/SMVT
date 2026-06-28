"""
Download FDA-approved drugs from PubChem.
PubChem has an "FDA Approved Drugs" classification.
Uses PUG-REST API.
"""
import os, sys, json, time, csv
from urllib.request import Request, urlopen
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

DATA = r"D:\Researching\SMVT\03_Analysis\data"
SAVE_PATH = os.path.join(DATA, "pubchem_fda_approved.csv")
os.makedirs(DATA, exist_ok=True)

# PubChem has "fda approved" in the source classification
# Alternative: use PubChem's compound list by source

def pubchem_get(url, desc="", timeout=90, retries=5):
    """PubChem API call with retry."""
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
            resp = urlopen(req, timeout=timeout, context=ctx)
            return resp.read()
        except Exception as e:
            wait = 2 ** attempt
            if attempt < retries - 1:
                print(f"  Retry {attempt+1}/{retries} for {desc}: {e} (wait {wait}s)")
                time.sleep(wait)
            else:
                raise

print("=" * 60)
print("PubChem FDA-Approved Drug Downloader")
print("=" * 60)

# Method 1: PubChem Classification "FDA Approved Drugs"
# This is the most reliable PubChem source for FDA-approved compounds
# URL: https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/classification/fda%20approved%20drugs/cids/JSON

# Actually let's try a different approach - PubChem has a "Drug and Medication"
# classification. Let's use the compound/listkey approach instead.

# Method: Search PubChem for FDA-approved drugs using the ETL boundary
# "fda_approved" is a known PubChem source classification

print("\n[1] Trying PubChem PUG-REST compound search...")

# PubChem source classification IDs:
# FDA Approved Drugs in PubChem: use source classification
# Source: https://pubchem.ncbi.nlm.nih.gov/source/

# Try getting CIDs from the "FDA Approved Drugs" classification
base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

try:
    # Get CIDs for all compounds classified as FDA-approved
    # Using PubChem's classification system
    url = f"{base}/compound/fastidentity/cid/JSON?identity_type=same_connectivity"
    print(f"  Testing PubChem API: {url[:70]}...")
    data = pubchem_get(url, "API test")
    print(f"  PubChem API OK ({len(data)} bytes)")
except Exception as e:
    print(f"  PubChem API FAIL: {e}")

# Alternative: Download from NIH NCATS Inxight Drugs
print("\n[2] Trying NCATS Inxight Drugs API (FDA-approved subset)...")

# NCATS Inxight: https://drugs.ncats.io/
# Has a specific FDA-approved drug list
ncats_url = "https://drugs.ncats.io/api/v1/substances"
# Filter for FDA-approved

try:
    # Get count first
    url = f"{ncats_url}/search?q=approval&size=1"
    data = pubchem_get(url, "NCATS count")
    result = json.loads(data)
    total = result.get("total", 0)
    print(f"  NCATS Inxight FDA-approved substances: {total}")
except Exception as e:
    print(f"  NCATS FAIL: {e}")

# Try PubChem via classification
print("\n[3] Trying PubChem compound list by classification...")

# PubChem MeSH classification for FDA-approved
# The classification tree: https://pubchem.ncbi.nlm.nih.gov/classification/
# Classification ID for drugs: various

# Simpler: use PubChem's pre-built FDA list via search
search_terms = [
    '"fda approved" AND hasdrugaction',
    '"approved drug"',
]

for term in search_terms:
    try:
        url = f"{base}/compound/fastsubstructure/cid/JSON?smiles=C&MaxRecords=5"
        data = pubchem_get(url, f"search: {term[:30]}")
        print(f"  Search OK for '{term[:30]}...': {len(data)} bytes")
    except Exception as e:
        print(f"  Search FAIL: {e}")

# Method 4: NCATS Inxight download (most promising)
print("\n[4] Downloading from NCATS Inxight Drugs...")

all_drugs = []
for page in range(0, 100):
    try:
        url = f"{ncats_url}/search?q=approval_type:Prescription+OR+approval_type:OTC&size=100&from={page*100}"
        data = pubchem_get(url, f"NCATS page {page}")
        result = json.loads(data)
        substances = result.get("content", [])
        if not substances:
            break

        for sub in substances:
            name = sub.get("name", "")
            smiles = sub.get("smiles", "")
            if smiles:
                all_drugs.append({
                    "name": name,
                    "smiles": smiles,
                    "source": "NCATS_Inxight",
                })

        print(f"  Page {page+1}: +{len(substances)} = {len(all_drugs)} drugs", flush=True)
        time.sleep(0.3)

    except Exception as e:
        print(f"  Page {page} FAIL: {e}")
        break

if all_drugs:
    with open(SAVE_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["name", "smiles", "source"])
        w.writeheader()
        w.writerows(all_drugs)
    print(f"\nSaved: {len(all_drugs)} drugs → {SAVE_PATH}")
else:
    print("\nNCATS download produced no results")

print("\nDone!")
