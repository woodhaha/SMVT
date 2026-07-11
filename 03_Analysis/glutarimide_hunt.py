#!/usr/bin/env python3
"""Targeted glutarimide/piperidine-2,6-dione scaffold expansion."""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import pandas as pd
import requests
from openbabel import pybel

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

VINA = "os.environ.get("VINA_BIN", "vina")"
RECEPTOR = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCK_DIR = "03_Analysis/docking"
os.makedirs(DOCK_DIR, exist_ok=True)
os.makedirs("03_Analysis/outputs", exist_ok=True)

CKPT = "03_Analysis/docking/glutarimide_checkpoint.json"
RESULTS = "03_Analysis/outputs/docking_glutarimide_results.csv"
CACHE = "03_Analysis/outputs/glutarimide_candidates.json"
CENTER = [-2.5, 1.0, -1.0]
BOX = [22, 22, 22]

# Seeds: top glutarimide hits from Path C + thalidomide family
SEEDS = {
    "glutarimide_hit1": "O=C1NC(=O)CCC1(Cc2ccccc2)Cl",    # -9.36 champion
    "glutarimide_hit2": "O=C1NC(=O)CCC1(c2ccccc2)C(F)(F)F", # -8.79
    "thalidomide": "O=C1NC(=O)CCC1N2C(=O)c3ccccc3C2=O",
    "lenalidomide": "NC(=O)C1=C(C)CN(C1=O)C2CCC(=O)NC2=O",
    "pomalidomide": "NC1=CC(=O)c2ccccc2C1=O",
}

# ═══ Phase 1: Fetch ═══
if not os.path.exists(CACHE):
    all_cids = set()
    for name, smi in SEEDS.items():
        for threshold in [85, 80]:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smi}/cids/JSON"
            try:
                r = requests.get(url, params={"Threshold": threshold, "MaxRecords": 250}, timeout=30)
                if r.status_code == 200:
                    cids = r.json().get("IdentifierList", {}).get("CID", [])
                    log.info(f"  {name} (t={threshold}%): {len(cids)} analogs")
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
                            "name": str(p.get("IUPACName", f"CID{p.get('CID','?')}"))[:50],
                            "smiles": smi, "mw": mw,
                            "logp": p.get("XLogP", 0),
                        })
        except: pass
        if (i // 100) % 3 == 0:
            log.info(f"  Props: {len(candidates)} valid...")
        time.sleep(0.3)

    log.info(f"Valid: {len(candidates)}")

    # Dedup
    seen = set()
    for f in ["03_Analysis/outputs/all_time_top_hits.csv",
              "03_Analysis/outputs/pathc_v2_candidates.json",
              "03_Analysis/outputs/drugbank_candidates_scored.csv"]:
        if not os.path.exists(f): continue
        if f.endswith(".json"):
            for c in json.load(open(f)):
                seen.add(c.get("smiles", ""))
        else:
            df = pd.read_csv(f)
            if "smiles" in df.columns:
                seen.update(s for s in df["smiles"] if isinstance(s, str))
    new = [c for c in candidates if c["smiles"] not in seen]
    log.info(f"After dedup: {len(new)} NEW")
    with open(CACHE, "w") as f: json.dump(new, f, indent=2)
else:
    new = json.load(open(CACHE))
    log.info(f"Loaded {len(new)} from cache")

if not new:
    log.info("No new candidates. Exiting."); sys.exit(0)

# ═══ Phase 2: Dock ═══
log.info(f"\nDocking {len(new)} glutarimide analogs...")
ckpt = json.load(open(CKPT)) if os.path.exists(CKPT) else {"completed": [], "failed": [], "results": []}
done = set(ckpt.get("completed", [])) | set(ckpt.get("failed", []))
pending = [c for c in new if c["name"] not in done]
log.info(f"Pending: {len(pending)} | Done: {len(done)}")

start = time.time()
for i, c in enumerate(pending):
    safe = c["name"][:35].replace("/","_").replace("\\","_").replace(" ","_").replace("*","_")
    pdbqt = os.path.join(DOCK_DIR, f"{safe}_glu.pdbqt")
    out = os.path.join(DOCK_DIR, f"{safe}_glu_docked.pdbqt")
    try:
        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            mol = pybel.readstring("smi", c["smiles"]); mol.addh()
            mol.make3D(forcefield="mmff94", steps=500)
            mol.localopt(forcefield="mmff94", steps=200)
            sdf = os.path.join(DOCK_DIR, f"{safe}_glu.sdf")
            mol.write("sdf", sdf, overwrite=True)
            subprocess.run(["obabel", sdf, "-O", pdbqt, "--gen3d"], capture_output=True, timeout=30)
        r = subprocess.run([
            VINA, "--receptor", RECEPTOR, "--ligand", pdbqt,
            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
            "--size_x", str(BOX[0]), "--size_y", str(BOX[1]), "--size_z", str(BOX[2]),
            "--exhaustiveness", "16", "--num_modes", "5", "--out", out,
        ], capture_output=True, text=True, timeout=300)
        affs = []
        for line in r.stdout.split("\n"):
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
        if affs:
            best = min(affs)
            ckpt["completed"].append(c["name"]); ckpt["results"].append({"name": c["name"], "best_affinity": best, "all_affinities": sorted(set(affs))})
            flag = " 🏆🏆" if best < -8.5 else (" 🏆" if best < -7.0 else (" ⭐" if best < -6.5 else ""))
            log.info(f"  [{i+1}/{len(pending)}] {c['name'][:30]:<30} DG={best:.2f}{flag}")
        else:
            ckpt["failed"].append(c["name"])
    except: ckpt["failed"].append(c["name"])
    if (i+1) % 50 == 0:
        json.dump(ckpt, open(CKPT, "w"), indent=2)

json.dump(ckpt, open(CKPT, "w"), indent=2)
elapsed = (time.time()-start)/60
log.info(f"DONE: {len(ckpt['completed'])} ok, {len(ckpt['failed'])} fail in {elapsed:.0f}min")
if ckpt["results"]:
    df_r = pd.DataFrame(ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv(RESULTS, index=False)
    b = df_r.iloc[0]
    log.info(f"Best: {b['name']} = {b['best_affinity']:.2f}")
    for i, (_, r) in enumerate(df_r.head(20).iterrows()):
        log.info(f"  {i+1}. {r['name'][:35]:<35} DG={r['best_affinity']:.2f}")
