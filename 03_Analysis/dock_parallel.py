#!/usr/bin/env python3
"""Parallel FDA leftover docking — 12 workers, ~6 min."""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

os.chdir("D:/Researching/SMVT")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

VINA = "C:/Users/woodh/bin/vina"
RECEPTOR = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCK_DIR = "03_Analysis/docking"
CKPT_FILE = "03_Analysis/docking/fda_leftover_checkpoint.json"
CENTER = [-2.5, 1.0, -1.0]
BOX = [22, 22, 22]
NWORKERS = 10  # ProcessPool: 10 Vina processes

os.makedirs(DOCK_DIR, exist_ok=True)

df = pd.read_csv("03_Analysis/outputs/undocked_fda_filtered.csv")
ckpt = json.load(open(CKPT_FILE))
done = set(ckpt.get("completed", [])) | set(ckpt.get("failed", []))
pending = [
    {"name": str(r["name"]), "smiles": str(r["smiles"])}
    for _, r in df.iterrows() if str(r["name"]) not in done
]
log.info(f"Pending: {len(pending)} | Workers: {NWORKERS} (ProcessPool)")

def dock_one(c):
    safe = c["name"][:35].replace("/", "_").replace("\\", "_").replace(" ", "_").replace("*", "_")
    pdbqt = os.path.join(DOCK_DIR, f"{safe}_fda.pdbqt")
    out = os.path.join(DOCK_DIR, f"{safe}_fda_docked.pdbqt")
    try:
        if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
            mol = Chem.MolFromSmiles(c["smiles"])
            if mol is None:
                return {"name": c["name"], "status": "bad_smiles"}
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.numThreads = 1
            if AllChem.EmbedMolecule(mol, params) != 0:
                params.useRandomCoords = True
                if AllChem.EmbedMolecule(mol, params) != 0:
                    return {"name": c["name"], "status": "no_conformer"}
            AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
            sdf = os.path.join(DOCK_DIR, f"{safe}_fda.sdf")
            w = Chem.SDWriter(sdf)
            w.write(mol)
            w.close()
            subprocess.run(
                ["obabel", sdf, "-O", pdbqt, "-p", "7.4"],
                capture_output=True, timeout=30
            )
            if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 50:
                return {"name": c["name"], "status": "bad_pdbqt"}

        r = subprocess.run(
            [
                VINA, "--receptor", RECEPTOR, "--ligand", pdbqt,
                "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
                "--center_z", str(CENTER[2]),
                "--size_x", str(BOX[0]), "--size_y", str(BOX[1]),
                "--size_z", str(BOX[2]),
                "--exhaustiveness", "16", "--num_modes", "5", "--out", out,
                "--cpu", "1",
            ],
            capture_output=True, text=True, timeout=300,
        )

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
            with open(out) as f:
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
            return {
                "name": c["name"], "status": "success",
                "best_affinity": min(affs),
                "all_affinities": sorted(set(affs)),
            }
        return {"name": c["name"], "status": "no_score"}
    except Exception:
        return {"name": c["name"], "status": "error"}

if __name__ == "__main__":
    start = time.time()
    with ProcessPoolExecutor(max_workers=NWORKERS) as ex:
        futures = {ex.submit(dock_one, c): c for c in pending}
        for i, f in enumerate(as_completed(futures)):
            r = f.result()
            if r["status"] == "success":
                ckpt["completed"].append(r["name"])
                ckpt["results"].append(r)
                flag = (
                    " ELITE" if r["best_affinity"] < -8.0
                    else (" HIT" if r["best_affinity"] < -7.0
                    else (" TOP" if r["best_affinity"] < -6.5 else ""))
                )
                log.info(f"  [{i+1}/{len(pending)}] {r['name'][:30]:<30} DG={r['best_affinity']:.2f}{flag}")
            else:
                ckpt["failed"].append(r["name"])
            if (i + 1) % 10 == 0:
                json.dump(ckpt, open(CKPT_FILE, "w"), indent=2)

    json.dump(ckpt, open(CKPT_FILE, "w"), indent=2)
    elapsed = (time.time() - start) / 60
    log.info(f"DONE: {len(ckpt['completed'])} ok in {elapsed:.0f}min")
    results = sorted(ckpt["results"], key=lambda x: x["best_affinity"])
    elite = sum(1 for r in results if r["best_affinity"] < -8.0)
    hits = sum(1 for r in results if r["best_affinity"] < -7.0)
    log.info(f"Elite: {elite} | Hits: {hits}")
    for r in results[:10]:
        log.info(f"  {r['name'][:35]:<35} DG={r['best_affinity']:.2f}")
    df_r = pd.DataFrame(ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv("03_Analysis/outputs/docking_fda_leftover_results.csv", index=False)
