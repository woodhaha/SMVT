#!/usr/bin/env python3
"""
Path B — ZINC FDA Library → ML Pre-Screen → Top 500 for Docking
================================================================
Phase 1: Download ZINC FDA-approved drug SMILES (~1,600 compounds)
Phase 2: Score with existing RF pharmacophore model (AUC=0.888)
Phase 3: Output top 500 CSV → ready for docking_batch_screen.py

Usage:
  python zinc_fda_pipeline.py                  # full pipeline
  python zinc_fda_pipeline.py --fetch-only     # only download (Phase 1)
  python zinc_fda_pipeline.py --top 500        # custom top-N
"""
import os, sys, time, json, pickle, logging, argparse, gzip, io, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/outputs", exist_ok=True)
os.makedirs("03_Analysis/models", exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("03_Analysis/models/zinc_pipeline.log", mode="w"),
              logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument("--fetch-only", action="store_true")
parser.add_argument("--top", type=int, default=500)
args = parser.parse_args()

# ═══════════════════════════════════════════════════════════════
# PHASE 1: Download ZINC FDA Library
# ═══════════════════════════════════════════════════════════════
ZINC_CACHE = "03_Analysis/models/zinc_fda_raw.json"

def fetch_zinc_fda():
    """Fetch FDA-approved drug SMILES from ZINC20 API.
    Uses ZINC's tranche browser — FDA-approved = drugs that have
    'fda' flag in ZINC's annotation. Falls back to ZINC InMan
    (in-stock + FDA annotation) if tranche API is slow.

    Strategy: Use ZINC20 substance API with 'fda' property filter.
    Page size 500, parallel page fetches for speed.
    """
    if os.path.exists(ZINC_CACHE):
        with open(ZINC_CACHE) as f:
            data = json.load(f)
        log.info(f"Loaded {len(data)} ZINC compounds from cache: {ZINC_CACHE}")
        return data

    log.info("Fetching ZINC FDA-approved library...")

    # Approach: ZINC20 REST API — /substances with fda=approved filter
    # Base URL for ZINC20 substances
    BASE = "https://zinc20.docking.org/substances.json"

    all_drugs = []
    page = 0

    # First, get count
    try:
        r = requests.get(BASE, params={"fda": "approved", "count": "true"}, timeout=30)
        if r.status_code == 200:
            total = r.json().get("count", 1600)
            log.info(f"ZINC FDA-approved count: {total}")
        else:
            total = 1600
            log.warning(f"Could not get count (HTTP {r.status_code}), assuming ~1600")
    except Exception as e:
        total = 1600
        log.warning(f"Count endpoint failed: {e}, assuming ~1600")

    while len(all_drugs) < total:
        params = {
            "fda": "approved",
            "page": page,
            "output_fields": "zinc_id smiles pref_name mol_weight logp"
        }
        try:
            r = requests.get(BASE, params=params, timeout=60)
            if r.status_code == 200:
                data = r.json()
                substances = data.get("substances", [])
                if not substances:
                    log.info(f"Page {page}: empty — ZINC API exhausted")
                    break
                all_drugs.extend(substances)
                log.info(f"Page {page}: +{len(substances)} compounds (total: {len(all_drugs)})")
                page += 1
            elif r.status_code == 404:
                log.info(f"Page {page}: 404 — ZINC API exhausted")
                break
            else:
                log.warning(f"Page {page}: HTTP {r.status_code}, retrying...")
                time.sleep(2)
                continue
        except Exception as e:
            log.warning(f"Page {page}: {e}, retrying in 5s...")
            time.sleep(5)
            continue

    # Save cache
    with open(ZINC_CACHE, 'w') as f:
        json.dump(all_drugs, f, indent=2)
    log.info(f"Saved {len(all_drugs)} compounds to {ZINC_CACHE}")
    return all_drugs


def fetch_zinc_fallback():
    """Fallback: Use ZINC's ready-to-dock FDA subset if API fails.
    This downloads the pre-computed ZINC FDA library in SMILES format.
    Source: https://zinc.docking.org/substances/subsets/fda.csv
    """
    if os.path.exists(ZINC_CACHE):
        with open(ZINC_CACHE) as f:
            data = json.load(f)
        log.info(f"Loaded {len(data)} from cache")
        return data

    log.info("API approach failed. Using ZINC bulk download fallback...")

    # ZINC FDA-approved subset (direct CSV download)
    urls = [
        "https://zinc20.docking.org/substances/subsets/fda-approved.csv",
        "https://files.docking.org/zinc20/fda/fda_approved.tsv",
    ]

    all_drugs = []

    for url in urls:
        try:
            log.info(f"Trying: {url}")
            r = requests.get(url, timeout=120, stream=True)
            if r.status_code == 200:
                content = r.text
                lines = content.strip().split('\n')
                header = lines[0].strip()
                log.info(f"Downloaded {len(lines)-1} lines from {url}")

                # Auto-detect format
                if '\t' in header:
                    sep = '\t'
                elif ',' in header:
                    sep = ','
                else:
                    sep = '\t'

                cols = header.lower().split(sep)
                zinc_idx = next((i for i, c in enumerate(cols) if 'zinc' in c or 'id' in c), 0)
                smiles_idx = next((i for i, c in enumerate(cols) if 'smiles' in c), 1)
                name_idx = next((i for i, c in enumerate(cols) if 'name' in c or 'pref' in c),
                               2 if len(cols) > 2 else None)

                for line in lines[1:]:
                    if not line.strip():
                        continue
                    parts = line.strip().split(sep)
                    if len(parts) > max(zinc_idx, smiles_idx):
                        entry = {
                            "zinc_id": parts[zinc_idx] if zinc_idx < len(parts) else "",
                            "smiles": parts[smiles_idx] if smiles_idx < len(parts) else "",
                        }
                        if name_idx and name_idx < len(parts):
                            entry["pref_name"] = parts[name_idx]
                        else:
                            entry["pref_name"] = entry["zinc_id"]
                        all_drugs.append(entry)

                if len(all_drugs) > 500:
                    log.info(f"Got {len(all_drugs)} from {url} — sufficient, stopping")
                    break
            else:
                log.warning(f"HTTP {r.status_code} from {url}")
        except Exception as e:
            log.warning(f"Failed to fetch {url}: {e}")
            continue

    if not all_drugs:
        log.error("All ZINC download methods failed!")
        return []

    with open(ZINC_CACHE, 'w') as f:
        json.dump(all_drugs, f, indent=2)
    log.info(f"Saved {len(all_drugs)} compounds to cache")
    return all_drugs


# ═══════════════════════════════════════════════════════════════
# PHASE 2: ML Pre-Screen with existing RF model
# ═══════════════════════════════════════════════════════════════
def ml_screen(compounds, top_n=500):
    """Score all compounds with existing pharmacophore RF model."""
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, rdFingerprintGenerator

    MODEL_PATH = "03_Analysis/models/smvt_rf_regressor.pkl"
    if not os.path.exists(MODEL_PATH):
        log.error(f"Model not found: {MODEL_PATH}")
        log.info("Run pharmacophore_ml_screen.py first to train the model")
        return compounds[:top_n]

    with open(MODEL_PATH, 'rb') as f:
        model = pickle.load(f)
    log.info(f"Loaded RF model from {MODEL_PATH}")

    # Generate ECFP4 fingerprints
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    valid = []
    for c in compounds:
        smi = c.get("smiles", "")
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            fp = fp_gen.GetFingerprint(mol)
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            score = model.predict([arr])[0]
            c["ml_score"] = float(score)
            c["mw"] = Chem.Descriptors.MolWt(mol)
            c["logp"] = Chem.Descriptors.MolLogP(mol)
            valid.append(c)
        except Exception as e:
            continue
        if len(valid) % 200 == 0:
            log.info(f"  Scored {len(valid)} compounds...")

    log.info(f"Scored {len(valid)}/{len(compounds)} valid compounds")

    # Sort by predicted binding affinity (more negative = better)
    valid.sort(key=lambda x: x.get("ml_score", 0))

    # Save full scored library
    df = pd.DataFrame(valid)
    df.to_csv("03_Analysis/outputs/zinc_fda_scored.csv", index=False)
    log.info(f"Saved full scored library: 03_Analysis/outputs/zinc_fda_scored.csv")

    # Show top candidates
    top = valid[:top_n]
    log.info(f"\n{'='*60}")
    log.info(f"Top {top_n} ML-predicted candidates:")
    log.info(f"{'Rank':<6} {'Name':<30} {'Score':>8} {'MW':>6} {'logP':>5}")
    log.info(f"{'-'*60}")
    for i, c in enumerate(top[:20]):
        name = (c.get("pref_name") or c.get("zinc_id", "?"))[:28]
        log.info(f"{i+1:<6} {name:<30} {c['ml_score']:>8.2f} {c.get('mw',0):>6.0f} {c.get('logp',0):>5.1f}")
    log.info(f"{'='*60}")

    return top


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    log.info("=" * 60)
    log.info("Path B — ZINC FDA → ML Screen → Top 500 for Docking")
    log.info("=" * 60)

    t0 = time.time()

    # Phase 1: Download
    log.info("\n[Phase 1] Downloading ZINC FDA library...")
    compounds = fetch_zinc_fda()
    if not compounds:
        log.warning("Primary API failed, trying fallback...")
        compounds = fetch_zinc_fallback()

    if not compounds:
        log.error("No compounds downloaded. Exiting.")
        sys.exit(1)

    log.info(f"Phase 1 complete: {len(compounds)} compounds ({time.time()-t0:.0f}s)")

    if args.fetch_only:
        log.info("--fetch-only flag set. Exiting after download.")
        sys.exit(0)

    # Phase 2: ML Screen
    log.info(f"\n[Phase 2] ML pre-screening with pharmacophore RF model...")
    top = ml_screen(compounds, top_n=args.top)

    # Write top-N CSV for docking
    df_top = pd.DataFrame(top)
    docking_csv = f"03_Analysis/outputs/zinc_fda_top{args.top}_for_docking.csv"
    df_top.to_csv(docking_csv, index=False)
    log.info(f"Phase 2 complete: top {args.top} saved to {docking_csv} ({time.time()-t0:.0f}s)")

    log.info(f"\n{'='*60}")
    log.info(f"Pipeline complete in {time.time()-t0:.0f}s")
    log.info(f"Next: python docking_batch_screen.py --input {docking_csv}")
    log.info(f"{'='*60}")
