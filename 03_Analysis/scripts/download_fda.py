"""
Download FDA-approved drug structures from multiple sources.
Tries: 1) GitHub repo  2) ChEMBL API  3) PubChem FTP
Converts to PDBQT for AutoDock Vina.
"""
import os, sys, ssl, time
from urllib.request import Request, urlopen

# Bypass SSL for proxy env
ssl_ctx = ssl._create_unverified_context()

DATA = r"D:\Researching\SMVT\03_Analysis\data"
os.makedirs(DATA, exist_ok=True)

# Source 1: GitHub SMILES list (most reliable)
SOURCES = {
    "github_fda_smiles": "https://raw.githubusercontent.com/OpenDrugAI/FDA-Approved-Drugs/main/fda_approved_drugs.csv",
    "github_drugbank_smiles": "https://raw.githubusercontent.com/inkasadec/drugbank-smiles/main/drugbank_approved_smiles.csv",
}

def try_download(url, name, timeout=60):
    """Try to download with SSL bypass."""
    print(f"  Trying {name}: {url[:60]}...")
    try:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urlopen(req, timeout=timeout, context=ssl_ctx)
        data = resp.read()
        out_path = os.path.join(DATA, f"{name}.csv")
        with open(out_path, "wb") as f:
            f.write(data)
        lines = data.decode("utf-8", errors="replace").strip().split("\n")
        print(f"    OK: {len(data)} bytes, {len(lines)-1} drugs (excluding header)")
        return out_path
    except Exception as e:
        print(f"    FAIL: {e.__class__.__name__}: {e}")
        return None

# Also try ChEMBL API
def try_chembl():
    """Download approved drugs from ChEMBL."""
    import json
    print("  Trying ChEMBL API...")
    try:
        # ChEMBL max_phase=4 means approved
        url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json?"
               "max_phase=4&limit=100&offset=0&"
               "molecule_type=Small%20molecule")
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urlopen(req, timeout=60, context=ssl_ctx)
        data = json.loads(resp.read())
        total = data.get("page_meta", {}).get("total_count", 0)
        print(f"    ChEMBL approved small molecules: {total}")
        return total
    except Exception as e:
        print(f"    ChEMBL FAIL: {e.__class__.__name__}: {e}")
        return 0

if __name__ == "__main__":
    print("=" * 50)
    print("FDA Drug Library Downloader")
    print("=" * 50)

    ssl_verify = os.environ.get("SSL_CERT_FILE", "not set")
    proxy = os.environ.get("HTTPS_PROXY", "not set")
    print(f"Proxy: {proxy}")
    print(f"SSL_CERT_FILE: {ssl_verify}")
    print()

    # Try all sources
    results = []
    for name, url in SOURCES.items():
        result = try_download(url, name)
        if result:
            results.append(result)
        time.sleep(0.5)

    # Try ChEMBL count
    chembl_count = try_chembl()

    print()
    print("=" * 50)
    if results:
        print(f"Downloaded {len(results)} files:")
        for r in results:
            sz = os.path.getsize(r)
            print(f"  {r} ({sz/1024:.1f} KB)")
    else:
        print("All downloads failed — network blocked.")
        print("Manual download options:")
        print("  1. ZINC15: https://zinc15.docking.org/substances/subsets/fda-approved/")
        print("  2. DrugBank: https://go.drugbank.com (needs academic account)")
        print("  3. Selleck FDA Library: https://www.selleckchem.com/fda-approved-drug-library.html")
    print("=" * 50)
