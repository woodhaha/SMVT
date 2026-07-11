#!/usr/bin/env python3
"""Dock leftover FDA drugs - v2: RDKit 3D + obabel PDBQT only. Robust version."""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

VINA = "os.environ.get("VINA_BIN", "vina")"
RECEPTOR = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCK_DIR = "03_Analysis/docking"
os.makedirs(DOCK_DIR, exist_ok=True)
CKPT = "03_Analysis/docking/fda_leftover_checkpoint.json"
RESULTS = "03_Analysis/outputs/docking_fda_leftover_results.csv"
CENTER = [-2.5, 1.0, -1.0]
BOX = [22, 22, 22]

df = pd.read_csv("03_Analysis/outputs/undocked_fda_filtered.csv")
candidates = [{"name": str(r["name"]), "smiles": str(r["smiles"])} for _, r in df.iterrows()]
log.info(f"Loaded {len(candidates)}")

ckpt = json.load(open(CKPT)) if os.path.exists(CKPT) else {"completed": [], "failed": [], "results": []}
done = set(ckpt.get("completed", [])) | set(ckpt.get("failed", []))
pending = [c for c in candidates if c["name"] not in done]
log.info(f"Pending: {len(pending)} | Done: {len(done)}")

start = time.time()
for i, c in enumerate(pending):
    safe = c["name"][:35].replace("/","_").replace("\\","_").replace(" ","_").replace("*","_")
    pdbqt = os.path.join(DOCK_DIR, f"{safe}_fda.pdbqt")
    out_pdbqt = os.path.join(DOCK_DIR, f"{safe}_fda_docked.pdbqt")
    try:
        # Step 1: RDKit 3D conformer
        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            mol = Chem.MolFromSmiles(c["smiles"])
            if mol is None:
                ckpt["failed"].append(c["name"]); continue
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.numThreads = 1
            status = AllChem.EmbedMolecule(mol, params)
            if status != 0:
                params.useRandomCoords = True
                params.maxIterations = 5000
                status = AllChem.EmbedMolecule(mol, params)
                if status != 0:
                    ckpt["failed"].append(c["name"]); continue
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            sdf = os.path.join(DOCK_DIR, f"{safe}_fda.sdf")
            w = Chem.SDWriter(sdf)
            w.write(mol)
            w.close()
            # Step 2: obabel SDF → PDBQT only (fast, no 3D gen needed)
            subprocess.run(["obabel", sdf, "-O", pdbqt, "-p", "7.4"],
                         capture_output=True, timeout=30)
            if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 50:
                ckpt["failed"].append(c["name"]); continue

        # Step 3: Dock
        r = subprocess.run([
            VINA, "--receptor", RECEPTOR, "--ligand", pdbqt,
            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
            "--size_x", str(BOX[0]), "--size_y", str(BOX[1]), "--size_z", str(BOX[2]),
            "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt,
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
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        for p in line.split():
                            try:
                                v = float(p)
                                if abs(v) < 100: affs.append(round(v, 3))
                            except: pass
        if affs:
            best = min(affs)
            ckpt["completed"].append(c["name"])
            ckpt["results"].append({"name": c["name"], "best_affinity": best, "all_affinities": sorted(set(affs))})
            flag = " ELITE" if best < -8.0 else (" HIT" if best < -7.0 else (" TOP" if best < -6.5 else ""))
            log.info(f"  [{i+1}/{len(pending)}] {c['name'][:30]:<30} DG={best:.2f}{flag}")
        else:
            ckpt["failed"].append(c["name"])
    except Exception:
        ckpt["failed"].append(c["name"])
    if (i + 1) % 5 == 0:
        json.dump(ckpt, open(CKPT, "w"), indent=2)

json.dump(ckpt, open(CKPT, "w"), indent=2)
elapsed = (time.time() - start) / 60
log.info(f"DONE: {len(ckpt['completed'])} ok, {len(ckpt['failed'])} fail in {elapsed:.0f}min")
if ckpt["results"]:
    df_r = pd.DataFrame(ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv(RESULTS, index=False)
    b = df_r.iloc[0]
    log.info(f"Best: {b['name']} = {b['best_affinity']:.2f}")
    elite = sum(1 for r in ckpt["results"] if r["best_affinity"] < -8.0)
    hits = sum(1 for r in ckpt["results"] if r["best_affinity"] < -7.0)
    log.info(f"Elite (<-8.0): {elite} | Hits (<-7.0): {hits}")
    for i, (_, r) in enumerate(df_r.head(15).iterrows()):
        log.info(f"  {i+1}. {r['name'][:35]:<35} DG={r['best_affinity']:.2f}")
