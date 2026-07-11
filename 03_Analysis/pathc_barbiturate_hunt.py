#!/usr/bin/env python3
"""
Path C — Targeted barbiturate/ureido substructure search + direct docking.
Searches PubChem for structural analogs of known barbiturate hits,
then docks directly (no ML pre-filter needed — we know the scaffold works).
"""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator
from openbabel import pybel

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/outputs", exist_ok=True)
os.makedirs("03_Analysis/docking", exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

VINA_BIN = "os.environ.get("VINA_BIN", "vina")"
RECEPTOR = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCK_DIR = "03_Analysis/docking"
CKPT = "03_Analysis/docking/pathc_checkpoint.json"
RESULTS = "03_Analysis/outputs/docking_pathc_results.csv"
CACHE = "03_Analysis/outputs/pathc_candidates.json"

CENTER = [-2.5, 1.0, -1.0]
BOX = [22, 22, 22]
EXHAUSTIVENESS = 16

# ═══ Known barbiturate SMILES (seed compounds) ═══
SEEDS = {
    "Phenobarbital": "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2",
    "Hexobarbital":  "CC1(C(=O)NC(=O)N(C1=O)C)C2=CC=CCC2",
    "Barbital":      "CCC1(CC)C(=O)NC(=O)NC1=O",
    "Primidone":     "C1CC(=O)NC(=O)C1(C2=CC=CC=C2)CC",
    "Cyclobarbital": "C1CC2(CC1)C(=O)NC(=O)NC2=O",
}

# Barbituric acid core SMARTS
BARB_CORE = "O=C1CC(=O)NC(=O)N1"
HYDANTOIN_CORE = "O=C1CNC(=O)N1"
UREIDO = "NC(=O)N"
CYCLIC_IMIDE = "[#6]1C(=O)NC(=O)[#6]1"

# ═══ Phase 1: Fetch analogs from PubChem ═══
def fetch_pubchem_analogs(smiles, threshold=85, max_results=200):
    """Similarity search from PubChem using a seed SMILES."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON"
    params = {"Threshold": threshold, "MaxRecords": max_results}
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code == 200:
            cids = r.json().get("IdentifierList", {}).get("CID", [])
            return cids
    except:
        pass
    return []

def fetch_properties(cids, batch_size=100):
    """Batch-fetch SMILES, MW, names from PubChem."""
    all_props = []
    for i in range(0, len(cids), batch_size):
        batch = cids[i:i+batch_size]
        cid_str = ",".join(str(c) for c in batch)
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/"
               f"property/CanonicalSMILES,MolecularWeight,XLogP,IUPACName/JSON")
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                props = r.json().get("PropertyTable", {}).get("Properties", [])
                all_props.extend(props)
        except:
            pass
        time.sleep(0.3)
    return all_props

log.info("=" * 60)
log.info("Path C — Targeted Barbiturate/Ureido Analog Search")
log.info("=" * 60)

if os.path.exists(CACHE):
    with open(CACHE) as f:
        candidates = json.load(f)
    log.info(f"Loaded {len(candidates)} candidates from cache")
else:
    # Search from each seed compound
    all_cids = set()
    for name, smi in SEEDS.items():
        cids = fetch_pubchem_analogs(smi, threshold=85, max_results=200)
        log.info(f"  {name}: {len(cids)} analogs (threshold 85%)")
        all_cids.update(cids)
        time.sleep(0.5)

    log.info(f"Total unique CIDs: {len(all_cids)}")

    # Fetch properties
    log.info("Fetching properties...")
    cid_list = list(all_cids)
    props = fetch_properties(cid_list)
    log.info(f"Got properties for {len(props)} compounds")

    # Filter
    candidates = []
    for p in props:
        smi = p.get("CanonicalSMILES", "")
        mw = p.get("MolecularWeight", 0)
        if not smi or mw < 120 or mw > 650:
            continue
        if "." in smi:  # Skip salts
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        # Check for at least ONE of the target substructures
        has_barb = mol.HasSubstructMatch(Chem.MolFromSmarts(BARB_CORE)) if Chem.MolFromSmarts(BARB_CORE) else False
        has_hyd = mol.HasSubstructMatch(Chem.MolFromSmarts(HYDANTOIN_CORE)) if Chem.MolFromSmarts(HYDANTOIN_CORE) else False
        has_ureido = mol.HasSubstructMatch(Chem.MolFromSmarts(UREIDO))
        has_imide = mol.HasSubstructMatch(Chem.MolFromSmarts(CYCLIC_IMIDE)) if Chem.MolFromSmarts(CYCLIC_IMIDE) else False

        if has_barb or has_hyd or has_ureido or has_imide:
            candidates.append({
                "cid": int(p.get("CID", 0)),
                "name": p.get("IUPACName", f"CID{p.get('CID','?')}")[:40],
                "smiles": smi,
                "mw": mw,
                "logp": p.get("XLogP", 0),
                "has_barb": has_barb,
                "has_hydantoin": has_hyd,
                "has_ureido": has_ureido,
                "has_imide": has_imide,
            })

    log.info(f"Filtered to {len(candidates)} substructure-matched compounds")

    # Deduplicate vs already screened
    seen_names = set()
    seen_smiles = set()

    # Load all existing results
    for f in ["03_Analysis/outputs/screening_master_results.csv",
              "03_Analysis/outputs/docking_round5_results.csv",
              "03_Analysis/outputs/docking_round5b_results.csv",
              "03_Analysis/outputs/docking_batch_results.csv",
              "03_Analysis/outputs/docking_results.csv",
              "03_Analysis/outputs/docking_expanded_results.csv"]:
        if os.path.exists(f):
            df = pd.read_csv(f)
            if "name" in df.columns:
                seen_names.update(df["name"].values)
            if "smiles" in df.columns:
                seen_smiles.update(s for s in df["smiles"].values if isinstance(s, str))

    new = []
    for c in candidates:
        if c["name"] not in seen_names and c["smiles"] not in seen_smiles:
            new.append(c)

    log.info(f"After dedup vs {len(seen_names)} known: {len(new)} NEW compounds")
    candidates = new

    with open(CACHE, 'w') as f:
        json.dump(candidates, f, indent=2)

# ═══ Phase 2: Direct Docking ═══
def prep_and_dock(smiles, name):
    """OpenBabel PDBQT prep → Vina docking."""
    safe = name.replace("/","_").replace("\\","_").replace(" ","_").replace("*","_")[:40]
    pdbqt = os.path.join(DOCK_DIR, f"{safe}_pathc.pdbqt")

    try:
        # Prepare
        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            mol = pybel.readstring("smi", smiles)
            mol.addh()
            mol.make3D(forcefield="mmff94", steps=500)
            mol.localopt(forcefield="mmff94", steps=200)
            sdf = os.path.join(DOCK_DIR, f"{safe}_pathc.sdf")
            mol.write("sdf", sdf, overwrite=True)
            subprocess.run(["obabel", sdf, "-O", pdbqt, "--gen3d"],
                         capture_output=True, timeout=30)

        # Dock
        out = os.path.join(DOCK_DIR, f"{safe}_pathc_docked.pdbqt")
        result = subprocess.run([
            VINA_BIN, "--receptor", RECEPTOR, "--ligand", pdbqt,
            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
            "--center_z", str(CENTER[2]),
            "--size_x", str(BOX[0]), "--size_y", str(BOX[1]), "--size_z", str(BOX[2]),
            "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", "5", "--out", out,
        ], capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            return {"name": name, "status": "vina_error"}

        affs = []
        for line in result.stdout.split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    v = float(parts[1])
                    if abs(v) < 100: affs.append(round(v, 3))
                except: pass

        if not affs:
            with open(out) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        for p in line.split():
                            try:
                                v = float(p)
                                if abs(v) < 100: affs.append(round(v, 3))
                            except: pass

        if not affs:
            return {"name": name, "status": "parse_error"}

        return {"name": name, "best_affinity": min(affs), "all_affinities": sorted(set(affs)), "status": "success"}
    except Exception as e:
        return {"name": name, "status": "error", "error": str(e)[:200]}

# ═══ Run ═══
log.info(f"\n{'='*60}")
log.info(f"Docking {len(candidates)} targeted analogs (Vina ex={EXHAUSTIVENESS})")
log.info(f"{'='*60}")

if os.path.exists(CKPT):
    ckpt = json.load(open(CKPT))
else:
    ckpt = {"completed": [], "failed": [], "results": []}

done = set(ckpt.get("completed", [])) | set(ckpt.get("failed", []))
pending = [c for c in candidates if c["name"] not in done]
log.info(f"Pending: {len(pending)} | Done: {len(done)}")

start = time.time()
for i, c in enumerate(pending):
    r = prep_and_dock(c["smiles"], c["name"])
    if r["status"] == "success":
        ckpt["completed"].append(c["name"])
        ckpt["results"].append(r)
        flag = " 🏆 HIT!" if r["best_affinity"] < -7.0 else (" ⭐" if r["best_affinity"] < -6.5 else "")
        log.info(f"  [{i+1}/{len(pending)}] {c['name'][:30]:<30} ΔG={r['best_affinity']:.2f}{flag}")
    else:
        ckpt["failed"].append(c["name"])
        log.warning(f"  [{i+1}/{len(pending)}] {c['name'][:30]:<30} FAILED: {r.get('status')}")

    if (i+1) % 20 == 0:
        json.dump(ckpt, open(CKPT, "w"), indent=2)

json.dump(ckpt, open(CKPT, "w"), indent=2)

# Summary
elapsed = (time.time() - start) / 60
n_ok, n_fail = len(ckpt["completed"]), len(ckpt["failed"])
log.info(f"\nDONE: {n_ok} docked, {n_fail} failed in {elapsed:.1f}min")

if ckpt["results"]:
    df_r = pd.DataFrame(ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv(RESULTS, index=False)
    best = df_r.iloc[0]
    log.info(f"Best: {best['name']} = {best['best_affinity']:.2f}")
    hits = sum(1 for r in ckpt["results"] if r["best_affinity"] < -7.0)
    log.info(f"Hits (<-7.0): {hits}")
    log.info(f"\nTop 10:")
    for i, (_, r) in enumerate(df_r.head(10).iterrows()):
        log.info(f"  {i+1}. {r['name'][:35]:<35} ΔG={r['best_affinity']:.2f}")
