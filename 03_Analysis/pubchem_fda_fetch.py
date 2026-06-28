#!/usr/bin/env python3
"""
Quick PubChem FDA-approved drug fetcher → ML screen → docking candidates.
PubChem PUG REST is the most reliable FDA drug source (no 403 blocks like ZINC).
"""
import os, sys, time, json, pickle, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator

os.chdir("D:/Researching/SMVT")
os.makedirs("03_Analysis/outputs", exist_ok=True)
os.makedirs("03_Analysis/models", exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("03_Analysis/models/pubchem_fetch.log", mode="w"),
              logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

PUBCHEM_CACHE = "03_Analysis/models/pubchem_fda_raw.json"
MODEL_PATH = "03_Analysis/models/smvt_rf_regressor.pkl"
OUT_CSV = "03_Analysis/outputs/pubchem_fda_top500.csv"
TOP_N = 500

# ═══ Step 1: Get FDA-approved drug CIDs from PubChem ═══
if os.path.exists(PUBCHEM_CACHE):
    with open(PUBCHEM_CACHE) as f:
        drugs = json.load(f)
    log.info(f"Loaded {len(drugs)} drugs from cache")
else:
    log.info("Querying PubChem for FDA-approved drug CIDs...")

    # PubChem classification: 'fda_approved_drugs'
    # Using PUG REST to search by classification
    CID_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cids/JSON"
               "?classification=fda_approved_drugs&list_return=list")

    try:
        r = requests.get(CID_URL, timeout=60)
        r.raise_for_status()
        data = r.json()
        cid_list = data.get("IdentifierList", {}).get("CID", [])
        log.info(f"PubChem returned {len(cid_list)} FDA-approved CIDs")
    except Exception as e:
        log.error(f"CID fetch failed: {e}")
        # Fallback: use the more basic classification query
        alt_url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cids/JSON"
                   "?list_return=list&compound_type=approved")
        r = requests.get(alt_url, timeout=60)
        data = r.json()
        cid_list = data.get("IdentifierList", {}).get("CID", [])
        log.info(f"Fallback returned {len(cid_list)} CIDs")

    if not cid_list:
        log.error("No CIDs found. Exiting.")
        sys.exit(1)

    # Step 2: Batch-fetch properties (SMILES, MW, names) — PubChem allows 100 CIDs per request
    drugs = []
    batch_size = 100

    for i in range(0, len(cid_list), batch_size):
        batch = cid_list[i:i+batch_size]
        cid_str = ",".join(str(c) for c in batch)

        prop_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/"
            f"property/CanonicalSMILES,MolecularWeight,XLogP,Charge,HBondDonorCount,"
            f"HBondAcceptorCount,RotatableBondCount,IUPACName/JSON"
        )

        try:
            r = requests.get(prop_url, timeout=30)
            r.raise_for_status()
            props = r.json().get("PropertyTable", {}).get("Properties", [])
            for p in props:
                smi = p.get("CanonicalSMILES")
                if smi and len(smi) > 5:
                    drugs.append({
                        "cid": int(p.get("CID", 0)),
                        "smiles": smi,
                        "name": p.get("IUPACName", f"CID{p.get('CID','?')}"),
                        "mw": float(p.get("MolecularWeight", 0)),
                        "logp": float(p.get("XLogP", 0)),
                    })

            log.info(f"  Batch {i//batch_size + 1}: {len(props)} properties "
                     f"(total drugs: {len(drugs)})")
            time.sleep(0.3)  # Be nice to PubChem

        except Exception as e:
            log.warning(f"  Batch {i//batch_size + 1} failed: {e}")
            time.sleep(2)
            continue

    # Save cache
    with open(PUBCHEM_CACHE, 'w') as f:
        json.dump(drugs, f, indent=2)
    log.info(f"Cached {len(drugs)} drugs to {PUBCHEM_CACHE}")

# ═══ Step 2: Filter drug-like + deduplicate ═══
log.info(f"\nFiltering and deduplicating {len(drugs)} compounds...")

# Remove invalid / too small / too large
valid = []
seen_smiles = set()
for d in drugs:
    smi = d.get("smiles", "")
    if not smi or smi in seen_smiles:
        continue
    mw = d.get("mw", 0)
    if mw < 100 or mw > 800:  # typical drug range
        continue
    # Exclude metals, salts (contains . in SMILES)
    if '.' in smi:
        continue
    seen_smiles.add(smi)
    valid.append(d)

log.info(f"  After drug-like filter: {len(valid)} compounds")

# Deduplicate against already-screened ChEMBL set
if os.path.exists("03_Analysis/outputs/screening_master_results.csv"):
    screened = pd.read_csv("03_Analysis/outputs/screening_master_results.csv")
    screened_smiles = set(screened.get("smiles", [])) if "smiles" in screened.columns else set()
    # Also check by name
    screened_names = set(screened.get("compound", [])) | set(screened.get("name", []))

    new_drugs = []
    for d in valid:
        if d["smiles"] not in screened_smiles:
            new_drugs.append(d)

    log.info(f"  After dedup vs screened: {len(new_drugs)} NEW compounds")
    valid = new_drugs
else:
    log.info("  No screening history found — all compounds are new")

# ═══ Step 3: ML pre-screen with existing RF model ═══
log.info(f"\n[Phase 2] ML pre-screening with pharmacophore RF model...")

if not os.path.exists(MODEL_PATH):
    log.error(f"Model not found: {MODEL_PATH}. Run pharmacophore_ml_screen.py first.")
    sys.exit(1)

with open(MODEL_PATH, 'rb') as f:
    model = pickle.load(f)

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

scored = []
for d in valid:
    mol = Chem.MolFromSmiles(d["smiles"])
    if mol is None:
        continue
    try:
        fp = fp_gen.GetFingerprint(mol)
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        score = float(model.predict([arr])[0])
        d["ml_score"] = score
        scored.append(d)
    except:
        continue

scored.sort(key=lambda x: x["ml_score"])
log.info(f"  Scored {len(scored)} compounds")
log.info(f"  Score range: {scored[0]['ml_score']:.2f} to {scored[-1]['ml_score']:.2f}")

# ═══ Step 4: Show top candidates ═══
top = scored[:TOP_N]
log.info(f"\n{'='*65}")
log.info(f"Top {TOP_N} ML-predicted SMVT ligands:")
log.info(f"{'Rank':<6} {'CID':<10} {'Name':<28} {'Score':>8} {'MW':>6}")
log.info(f"{'-'*65}")
for i, d in enumerate(top[:25]):
    name = d["name"][:26] if d["name"] else f"CID{d['cid']}"
    log.info(f"{i+1:<6} {d['cid']:<10} {name:<28} {d['ml_score']:>8.2f} {d['mw']:>6.0f}")

# ═══ Step 5: Save for docking ═══
df_top = pd.DataFrame(top)
# Ensure columns match what docking_batch_screen.py expects
df_top["name"] = df_top.apply(
    lambda r: f"{r['name'][:20] if r['name'] else 'CID'+str(r['cid'])}_CID{r['cid']}", axis=1
)
df_top = df_top.rename(columns={"smiles": "smiles"})

# Select columns for docking
dock_cols = ["name", "smiles", "ml_score", "mw", "logp"]
df_dock = df_top[[c for c in dock_cols if c in df_top.columns]]
df_dock.to_csv(OUT_CSV, index=False)
log.info(f"\nSaved top {len(top)} for docking → {OUT_CSV}")
log.info(f"Next: python docking_batch_screen.py --input {OUT_CSV}")

# Also save full scored library
full_csv = "03_Analysis/outputs/pubchem_fda_scored_full.csv"
df_full = pd.DataFrame(scored)
df_full["name"] = df_full.apply(
    lambda r: f"PubChem_CID{r['cid']}", axis=1
)
df_full.to_csv(full_csv, index=False)
log.info(f"Full library ({len(scored)} compounds) → {full_csv}")
