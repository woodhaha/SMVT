#!/usr/bin/env python3
"""Path C v2 — PubChem similarity search → dedup → direct docking. No SMARTS filter."""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import pandas as pd
import requests
from rdkit import Chem
from openbabel import pybel

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

VINA = "os.environ.get("VINA_BIN", "vina")"
RECEPTOR = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCK_DIR = "03_Analysis/docking"
os.makedirs(DOCK_DIR, exist_ok=True)
os.makedirs("03_Analysis/outputs", exist_ok=True)

CKPT = "03_Analysis/docking/pathc_v2_checkpoint.json"
RESULTS = "03_Analysis/outputs/docking_pathc_v2_results.csv"
CACHE = "03_Analysis/outputs/pathc_v2_candidates.json"

CENTER = [-2.5, 1.0, -1.0]
BOX = [22, 22, 22]

SEEDS = {
    "Phenobarbital": "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2",
    "Hexobarbital": "CC1(C(=O)NC(=O)N(C1=O)C)C2=CC=CCC2",
    "Barbital": "CCC1(CC)C(=O)NC(=O)NC1=O",
    "Primidone": "C1CC(=O)NC(=O)C1(C2=CC=CC=C2)CC",
}

# ═══ Phase 1: Fetch candidates ═══
if not os.path.exists(CACHE):
    all_cids = set()
    for name, smi in SEEDS.items():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smi}/cids/JSON"
        try:
            r = requests.get(url, params={"Threshold": 85, "MaxRecords": 200}, timeout=30)
            if r.status_code == 200:
                cids = r.json().get("IdentifierList", {}).get("CID", [])
                log.info(f"  {name}: {len(cids)} analogs")
                all_cids.update(cids)
        except Exception as e:
            log.warning(f"  {name}: {e}")
        time.sleep(0.5)

    log.info(f"Total unique CIDs: {len(all_cids)}")
    cid_list = list(all_cids)

    candidates = []
    for i in range(0, len(cid_list), 100):
        batch = cid_list[i:i+100]
        cid_str = ",".join(str(c) for c in batch)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/CanonicalSMILES,MolecularWeight,XLogP,IUPACName/JSON"
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                for p in r.json().get("PropertyTable", {}).get("Properties", []):
                    smi = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES", "")
                    mw = float(p.get("MolecularWeight", 0))
                    if smi and 120 < mw < 650 and "." not in smi:
                        candidates.append({
                            "cid": int(p.get("CID", 0)),
                            "name": str(p.get("IUPACName", f"CID{p.get('CID','?')}"))[:45],
                            "smiles": smi,
                            "mw": mw,
                            "logp": p.get("XLogP", 0),
                        })
        except Exception as e:
            pass
        if (i // 100) % 3 == 0:
            log.info(f"  Properties: {len(candidates)} valid...")
        time.sleep(0.3)

    log.info(f"Valid candidates: {len(candidates)}")

    # Dedup against all existing screening results
    seen_smiles = set()
    result_files = [
        "03_Analysis/outputs/screening_master_results.csv",
        "03_Analysis/outputs/docking_round5_results.csv",
        "03_Analysis/outputs/docking_round5b_results.csv",
        "03_Analysis/outputs/docking_batch_results.csv",
        "03_Analysis/outputs/docking_results.csv",
        "03_Analysis/outputs/docking_expanded_results.csv",
    ]
    for f in result_files:
        if os.path.exists(f):
            df = pd.read_csv(f)
            if "smiles" in df.columns:
                for s in df["smiles"]:
                    if isinstance(s, str):
                        seen_smiles.add(s)

    new = [c for c in candidates if c["smiles"] not in seen_smiles]
    log.info(f"After dedup: {len(new)} NEW compounds")

    with open(CACHE, "w") as f:
        json.dump(new, f, indent=2)
else:
    new = json.load(open(CACHE))
    log.info(f"Loaded {len(new)} from cache")

if not new:
    log.info("No new candidates to dock. Exiting.")
    sys.exit(0)

# ═══ Phase 2: Dock ═══
log.info(f"\nDocking {len(new)} PubChem barbiturate analogs (Vina ex=16)...")

if os.path.exists(CKPT):
    ckpt = json.load(open(CKPT))
else:
    ckpt = {"completed": [], "failed": [], "results": []}

done = set(ckpt.get("completed", [])) | set(ckpt.get("failed", []))
pending = [c for c in new if c["name"] not in done]
log.info(f"Pending: {len(pending)} | Done: {len(done)}")

start = time.time()
for i, c in enumerate(pending):
    safe = c["name"][:30].replace("/", "_").replace("\\", "_").replace(" ", "_").replace("*", "_")
    pdbqt = os.path.join(DOCK_DIR, f"{safe}_pcv2.pdbqt")
    out_pdbqt = os.path.join(DOCK_DIR, f"{safe}_pcv2_docked.pdbqt")

    try:
        # Prepare with OpenBabel
        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            mol = pybel.readstring("smi", c["smiles"])
            mol.addh()
            mol.make3D(forcefield="mmff94", steps=500)
            mol.localopt(forcefield="mmff94", steps=200)
            sdf = os.path.join(DOCK_DIR, f"{safe}_pcv2.sdf")
            mol.write("sdf", sdf, overwrite=True)
            subprocess.run(["obabel", sdf, "-O", pdbqt, "--gen3d"], capture_output=True, timeout=30)

        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            ckpt["failed"].append(c["name"])
            continue

        # Dock with Vina
        r = subprocess.run([
            VINA, "--receptor", RECEPTOR, "--ligand", pdbqt,
            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
            "--size_x", str(BOX[0]), "--size_y", str(BOX[1]), "--size_z", str(BOX[2]),
            "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt,
        ], capture_output=True, text=True, timeout=300)

        if r.returncode != 0:
            ckpt["failed"].append(c["name"])
            continue

        # Parse affinities
        affs = []
        for line in r.stdout.split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    v = float(parts[1])
                    if abs(v) < 100:
                        affs.append(round(v, 3))
                except ValueError:
                    pass
        if not affs:
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        for p in line.split():
                            try:
                                v = float(p)
                                if abs(v) < 100:
                                    affs.append(round(v, 3))
                            except ValueError:
                                pass

        if affs:
            best = min(affs)
            ckpt["completed"].append(c["name"])
            ckpt["results"].append({"name": c["name"], "best_affinity": best, "all_affinities": sorted(set(affs))})
            flag = " 🏆" if best < -7.0 else (" ⭐" if best < -6.5 else "")
            log.info(f"  [{i+1}/{len(pending)}] {c['name'][:30]:<30} ΔG={best:.2f}{flag}")
        else:
            ckpt["failed"].append(c["name"])

    except Exception as e:
        ckpt["failed"].append(c["name"])

    if (i + 1) % 30 == 0:
        json.dump(ckpt, open(CKPT, "w"), indent=2)

# Final save
json.dump(ckpt, open(CKPT, "w"), indent=2)

elapsed = (time.time() - start) / 60
n_ok = len(ckpt["completed"])
n_fail = len(ckpt["failed"])
log.info(f"\nDONE: {n_ok} docked, {n_fail} failed in {elapsed:.0f}min")

if ckpt["results"]:
    df_r = pd.DataFrame(ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv(RESULTS, index=False)
    best = df_r.iloc[0]
    log.info(f"Best: {best['name']} = {best['best_affinity']:.2f} kcal/mol")
    hits = sum(1 for r in ckpt["results"] if r["best_affinity"] < -7.0)
    log.info(f"Hits (<-7.0): {hits}")

    log.info(f"\nTop 15:")
    for i, (_, r) in enumerate(df_r.head(15).iterrows()):
        log.info(f"  {i+1}. {r['name'][:35]:<35} ΔG={r['best_affinity']:.2f}")
