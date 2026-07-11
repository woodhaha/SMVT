#!/usr/bin/env python3
"""
Phase A Step 2 — Fetch Drug Candidates from ChEMBL + ML Pre-Screen (ROBUST)
============================================================================
- Paginated fetch with retry logic and checkpointing
- Resumes from saved checkpoint if connection drops
- ML pre-screening with SMVT pharmacophore model
"""

import os, sys, pickle, logging, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/outputs", exist_ok=True)
os.makedirs("03_Analysis/models", exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("03_Analysis/models/drugbank_fetch.log", mode="w"),
              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ═══ Step 1: Robust paginated fetch from ChEMBL ═══
CHECKPOINT = "03_Analysis/models/chembl_checkpoint.json"
PAGE_SIZE = 100
MAX_RETRIES = 5

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

session = requests.Session()
retry_strategy = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retry_strategy))

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PARAMS_TEMPLATE = {
    "max_phase": 4,
    "molecule_type": "Small molecule",
    "limit": PAGE_SIZE,
}

# Resume from checkpoint
if os.path.exists(CHECKPOINT):
    with open(CHECKPOINT) as f:
        ckpt = json.load(f)
    all_drugs = ckpt.get("drugs", [])
    offset = ckpt.get("offset", 0)
    log.info(f"Resuming from checkpoint: {len(all_drugs)} drugs, offset={offset}")
else:
    all_drugs = []
    offset = 0

# Get total count first
try:
    count_params = {**PARAMS_TEMPLATE, "limit": 1, "offset": 0}
    r = session.get(BASE_URL, params=count_params, timeout=30)
    total_count = r.json().get("page_meta", {}).get("total_count", 3000)
    log.info(f"ChEMBL total approved small molecules: {total_count}")
except:
    total_count = 3000
    log.warning(f"Could not get total count, assuming {total_count}")

# Fetch remaining pages
consecutive_failures = 0
while offset < total_count:
    if consecutive_failures >= MAX_RETRIES:
        log.error(f"Too many failures, stopping at offset {offset}")
        break

    params = {**PARAMS_TEMPLATE, "offset": offset}
    try:
        r = session.get(BASE_URL, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        molecules = data.get("molecules", [])

        if not molecules:
            log.info(f"No more molecules at offset {offset}, done.")
            break

        for mol in molecules:
            smiles = None
            structures = mol.get("molecule_structures")
            if structures:
                smiles = structures.get("canonical_smiles")
            if smiles:
                all_drugs.append({
                    "chembl_id": mol.get("molecule_chembl_id", ""),
                    "name": mol.get("pref_name") or mol.get("molecule_chembl_id", "UNK"),
                    "smiles": smiles,
                })

        offset += PAGE_SIZE
        consecutive_failures = 0

        # Save checkpoint every 500
        if offset % 500 == 0:
            with open(CHECKPOINT, "w") as f:
                json.dump({"drugs": all_drugs, "offset": offset}, f)
            log.info(f"  ... {len(all_drugs)} drugs fetched (offset={offset}/{total_count})")

        time.sleep(0.3)  # Be polite to ChEMBL servers

    except Exception as e:
        consecutive_failures += 1
        log.warning(f"Error at offset {offset} (fail {consecutive_failures}/{MAX_RETRIES}): {e}")
        time.sleep(5 * consecutive_failures)

# Final save
with open(CHECKPOINT, "w") as f:
    json.dump({"drugs": all_drugs, "offset": offset}, f)

log.info(f"Fetched {len(all_drugs)} drugs with SMILES")

# ═══ Step 2: Filter + deduplicate ═══
df_drugs = pd.DataFrame(all_drugs)
df_drugs = df_drugs.drop_duplicates(subset="smiles")
log.info(f"After dedup: {len(df_drugs)}")

def is_drug_like(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    if mw < 120 or mw > 800:
        return False
    logp = Descriptors.MolLogP(mol)
    if logp < -3 or logp > 7:
        return False
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    if hbd > 8 or hba > 15:
        return False
    return True

df_drugs["valid"] = df_drugs["smiles"].apply(is_drug_like)
df_filtered = df_drugs[df_drugs["valid"]].copy()
log.info(f"Drug-like filter passed: {len(df_filtered)}/{len(df_drugs)}")

# ═══ Step 3: Generate features + ML screen ═══
log.info("Loading SMVT ML model...")
with open("03_Analysis/models/smvt_ml_screen.pkl", "rb") as f:
    model = pickle.load(f)

rf_reg = model["rf_regressor"]
rf_clf = model["rf_classifier"]
scaler = model["scaler"]

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

fps, desc_list, valid_smiles = [], [], []
for smi in df_filtered["smiles"]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    fps.append(list(mfpgen.GetFingerprint(mol)))
    desc_list.append({
        "MW": Descriptors.MolWt(mol), "LogP": Descriptors.MolLogP(mol),
        "HBA": Descriptors.NumHAcceptors(mol), "HBD": Descriptors.NumHDonors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol), "TPSA": Descriptors.TPSA(mol),
        "RingCount": Descriptors.RingCount(mol), "AromaticRing": Descriptors.NumAromaticRings(mol),
        "Carboxyl": smi.count("C(=O)O"), "FractionCsp3": Descriptors.FractionCSP3(mol),
        "HeavyAtom": Descriptors.HeavyAtomCount(mol),
    })
    valid_smiles.append(smi)

fp_array = np.array(fps)
df_desc = pd.DataFrame(desc_list)
desc_scaled = scaler.transform(df_desc[model["descriptor_names"]])
X_screen = np.hstack([fp_array, desc_scaled])

df_ml = df_filtered[df_filtered["smiles"].isin(valid_smiles)].copy()
df_ml["predicted_affinity"] = rf_reg.predict(X_screen)
df_ml["predicted_hit_prob"] = rf_clf.predict_proba(X_screen)[:, 1]

log.info(f"ML scored: {len(df_ml)} compounds")
log.info(f"Predicted hits (ΔG < -6.5): {(df_ml['predicted_affinity'] < -6.5).sum()}")

# ═══ Step 4: Diversity selection ═══
def get_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except: return None

TOP_N = 500
df_ranked = df_ml.sort_values("predicted_affinity")
selected, scaffold_counts = [], {}

for _, row in df_ranked.iterrows():
    if len(selected) >= TOP_N: break
    scaff = get_scaffold(row["smiles"])
    cnt = scaffold_counts.get(scaff, 0)
    limit = 5 if len(selected) < 200 else 10
    if cnt >= limit: continue
    selected.append({"name": row["name"], "chembl_id": row["chembl_id"], "smiles": row["smiles"],
                     "predicted_affinity": row["predicted_affinity"],
                     "predicted_hit_prob": row["predicted_hit_prob"], "scaffold": scaff})
    scaffold_counts[scaff] = cnt + 1

df_selected = pd.DataFrame(selected).sort_values("predicted_affinity")
log.info(f"Selected {len(df_selected)} compounds ({len(scaffold_counts)} unique scaffolds)")

# ═══ Save outputs ═══
cols = ["name", "chembl_id", "smiles", "predicted_affinity", "predicted_hit_prob"]
df_ml.sort_values("predicted_affinity")[cols].to_csv(
    "03_Analysis/outputs/drugbank_candidates_scored.csv", index=False)
df_selected[cols + ["scaffold"]].to_csv(
    "03_Analysis/outputs/drugbank_top500_for_docking.csv", index=False)

log.info(f"Saved all scored → drugbank_candidates_scored.csv")
log.info(f"Saved top {len(df_selected)} → drugbank_top500_for_docking.csv")

# Top 20 preview
log.info(f"\nTOP 20 PREDICTED SMVT BINDERS:")
for i, (_, row) in enumerate(df_selected.head(20).iterrows()):
    log.info(f"  {i+1:2d}. {row['name'][:35]:35s} | ΔG={row['predicted_affinity']:.2f} | p={row['predicted_hit_prob']:.2f}")

n_pred_hits = (df_ml['predicted_affinity'] < -6.5).sum()
log.info(f"\nDone. {n_pred_hits} predicted hits from {len(df_ml)} screened. Ready for docking.")
