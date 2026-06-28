#!/usr/bin/env python3
"""
Phase A Step 3 — Batch Docking Pipeline for SMVT Virtual Screening
===================================================================
- Loads top-N ML-predicted candidates from drugbank_top500_for_docking.csv
- Prepares ligands: RDKit 3D conformer → meeko PDBQT
- Docks with AutoDock Vina (exhaustiveness=16)
- Checkpoints after every compound — resumable
- Uses ThreadPoolExecutor (works on Windows; Vina is the CPU bottleneck, not Python)

Usage:
  python docking_batch_screen.py              # dock all pending
  python docking_batch_screen.py --resume     # resume from checkpoint
  python docking_batch_screen.py --test 5     # dock first 5 only (quick test)
"""

import os, sys, subprocess, time, json, logging, argparse, warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ═══ Configuration ═══
VINA_BIN = "C:/Users/woodh/bin/vina"
MK_PREP_LIG = "C:/anaconda3/Scripts/mk_prepare_ligand.exe"

RECEPTOR_PDBQT = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCKING_DIR = "03_Analysis/docking"
OUTPUT_DIR = "03_Analysis/outputs"
CHECKPOINT = "03_Analysis/docking/batch_checkpoint.json"
RESULTS_CSV = f"{OUTPUT_DIR}/docking_batch_results.csv"

CENTER = [-2.5, 1.0, -1.0]
BOX_SIZE = [22, 22, 22]
EXHAUSTIVENESS = 16
NUM_MODES = 5
NUM_WORKERS = max(1, mp.cpu_count() - 2)  # Thread-based, Vina subprocess is the bottleneck

# ═══ Prepare ═══
os.makedirs(DOCKING_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"{DOCKING_DIR}/batch_docking.log", mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)


def prepare_ligand(smiles, name):
    """SMILES → 3D SDF → PDBQT. Returns path to PDBQT or raises."""
    safe_name = name.replace("/", "_").replace("\\", "_").replace(" ", "_").replace("*", "_")
    sdf_path = os.path.join(DOCKING_DIR, f"{safe_name}.sdf")
    pdbqt_path = os.path.join(DOCKING_DIR, f"{safe_name}.pdbqt")

    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
        return pdbqt_path

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise ValueError(f"RDKit cannot embed 3D conformer")

    AllChem.MMFFOptimizeMolecule(mol)
    # Keep explicit Hs — meeko needs them
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    result = subprocess.run(
        [MK_PREP_LIG, "-i", sdf_path, "-o", pdbqt_path],
        capture_output=True, text=True, timeout=60
    )
    if result.returncode != 0 or not os.path.exists(pdbqt_path):
        raise RuntimeError(f"meeko prep failed: {result.stderr[:200]}")

    return pdbqt_path


def dock_one(row):
    """Dock a single compound. Returns dict with results."""
    name, smiles = str(row["name"]), str(row["smiles"])
    safe_name = name.replace("/", "_").replace("\\", "_").replace(" ", "_").replace("*", "_")

    try:
        lig_pdbqt = prepare_ligand(smiles, name)
        out_pdbqt = os.path.join(DOCKING_DIR, f"{safe_name}_docked.pdbqt")

        cmd = [
            VINA_BIN,
            "--receptor", RECEPTOR_PDBQT,
            "--ligand", lig_pdbqt,
            "--center_x", str(CENTER[0]),
            "--center_y", str(CENTER[1]),
            "--center_z", str(CENTER[2]),
            "--size_x", str(BOX_SIZE[0]),
            "--size_y", str(BOX_SIZE[1]),
            "--size_z", str(BOX_SIZE[2]),
            "--exhaustiveness", str(EXHAUSTIVENESS),
            "--num_modes", str(NUM_MODES),
            "--out", out_pdbqt,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            return {"name": name, "status": "vina_error", "error": result.stderr[:200]}

        # Parse Vina output
        affinities = []
        for line in result.stdout.split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    val = float(parts[1])  # Vina format: "1   -7.5   0.0   0.0"
                    if abs(val) < 100:
                        affinities.append(round(val, 3))
                except (ValueError, IndexError):
                    continue

        if not affinities:
            # Try parsing output PDBQT
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        parts = line.split()
                        for p in parts:
                            try:
                                affinities.append(round(float(p), 3))
                            except ValueError:
                                pass

        if not affinities:
            return {"name": name, "status": "parse_error", "error": "No affinities parsed"}

        best = min(affinities)
        affinities = sorted(set(affinities))
        return {"name": name, "best_affinity": best, "all_affinities": affinities, "status": "success"}

    except subprocess.TimeoutExpired:
        return {"name": name, "status": "timeout"}
    except Exception as e:
        return {"name": name, "status": "error", "error": str(e)[:200]}


def save_checkpoint(checkpoint):
    with open(CHECKPOINT, "w") as f:
        json.dump(checkpoint, f, indent=2)


def main():
    """Main entry point — guarded for Windows multiprocessing compatibility."""
    os.chdir("D:/Researching/SMVT")

    parser = argparse.ArgumentParser()
    parser.add_argument("--test", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--input", type=str, default=None,
                       help="Custom input CSV (default: drugbank_top500_for_docking.csv)")
    parser.add_argument("--output", type=str, default=None,
                       help="Custom output CSV")
    parser.add_argument("--checkpoint", type=str, default=None,
                       help="Custom checkpoint file")
    args = parser.parse_args()

    # Load candidates
    if args.input:
        candidates_csv = args.input
    else:
        candidates_csv = f"{OUTPUT_DIR}/drugbank_top500_for_docking.csv"

    # Custom output
    if args.output:
        global RESULTS_CSV
        RESULTS_CSV = args.output

    # Custom checkpoint
    if args.checkpoint:
        global CHECKPOINT
        CHECKPOINT = args.checkpoint
    if not os.path.exists(candidates_csv):
        log.error(f"Candidate file not found: {candidates_csv}")
        sys.exit(1)

    df_candidates = pd.read_csv(candidates_csv)
    log.info(f"Loaded {len(df_candidates)} candidates")

    if args.test:
        df_candidates = df_candidates.head(args.test)
        log.info(f"TEST MODE: docking first {args.test}")

    # Load checkpoint
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            checkpoint = json.load(f)
        log.info(f"Loaded checkpoint: {len(checkpoint.get('completed', []))} completed")
    else:
        checkpoint = {"completed": [], "failed": [], "results": []}

    completed_set = set(checkpoint.get("completed", []))
    failed_set = set(checkpoint.get("failed", []))

    df_pending = df_candidates[~df_candidates["name"].isin(completed_set)]
    df_pending = df_pending[~df_pending["name"].isin(failed_set)]
    log.info(f"Pending: {len(df_pending)} | Done: {len(completed_set)} | Failed: {len(failed_set)}")

    if len(df_pending) == 0:
        log.info("All done!")
        sys.exit(0)

    # Docking loop
    log.info(f"Starting batch docking with {NUM_WORKERS} workers (exhaustiveness={EXHAUSTIVENESS})")

    rows = list(df_pending.iterrows())
    n_total = len(rows)
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(dock_one, row): row for _, row in rows}

        for future in as_completed(futures):
            result = future.result()
            name = result["name"]

            if result.get("status") == "success":
                checkpoint["completed"].append(name)
                checkpoint["results"].append({
                    "name": name,
                    "best_affinity": result["best_affinity"],
                    "all_affinities": result["all_affinities"],
                })
            else:
                checkpoint["failed"].append(name)
                log.warning(f"  FAILED: {name} ({result.get('status')}: {result.get('error', '')[:100]})")

            n_proc = len(checkpoint["completed"]) + len(checkpoint["failed"])
            elapsed = max(0.01, (time.time() - start_time) / 60)
            rate = n_proc / elapsed

            if n_proc % 10 == 0 or n_proc <= 3:
                eta = (n_total - n_proc) / max(0.01, rate)
                aff = result.get("best_affinity", "N/A")
                log.info(f"  [{n_proc}/{n_total}] {name[:30]}: {aff} | {rate:.1f} cpd/min | ETA {eta:.0f} min")

            save_checkpoint(checkpoint)

    # Final summary
    elapsed = (time.time() - start_time) / 60
    n_success = len(checkpoint["completed"])
    n_failed = len(checkpoint["failed"])
    log.info(f"\n{'='*60}\nDOCKING COMPLETE\n{'='*60}")
    log.info(f"  Successful: {n_success} | Failed: {n_failed} | Time: {elapsed:.1f} min")

    if checkpoint["results"]:
        affs = [r["best_affinity"] for r in checkpoint["results"] if r["best_affinity"] is not None]
        if affs:
            log.info(f"  Best ΔG: {min(affs):.2f} | Mean: {np.mean(affs):.2f}")
            log.info(f"  Hits (<-7.0): {sum(1 for a in affs if a < -7.0)} | (<-6.5): {sum(1 for a in affs if a < -6.5)}")
            df_results = pd.DataFrame(checkpoint["results"]).sort_values("best_affinity")
            df_results.to_csv(RESULTS_CSV, index=False)
            log.info(f"  Results → {RESULTS_CSV}")
            log.info(f"\nTOP 20:")
            for i, (_, r) in enumerate(df_results.head(20).iterrows()):
                log.info(f"  {i+1:2d}. {r['name'][:35]:35s} | {r['best_affinity']:.2f} kcal/mol")

    log.info("\nDone.")


if __name__ == "__main__":
    main()
