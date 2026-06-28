#!/usr/bin/env python3
"""
Fix and redock the 53 failed compounds using OpenBabel instead of meeko.
OB has better 3D conformer generation and more permissive SMILES handling.
"""
import os, sys, time, json, subprocess, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import openbabel as ob
from openbabel import pybel

os.chdir("D:/Researching/SMVT")
os.makedirs("03_Analysis/docking", exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

VINA_BIN = "C:/Users/woodh/bin/vina"
RECEPTOR_PDBQT = "03_Analysis/docking/SMVT_receptor.pdbqt"
DOCKING_DIR = "03_Analysis/docking"
CHECKPOINT = "03_Analysis/docking/round5b_checkpoint.json"
RESULTS_CSV = "03_Analysis/outputs/docking_round5b_results.csv"

CENTER = [-2.5, 1.0, -1.0]
BOX_SIZE = [22, 22, 22]
EXHAUSTIVENESS = 16

# ═══ Load failed compounds from round 5 ═══
ckpt = json.load(open("03_Analysis/docking/round5_checkpoint.json"))
failed_names = set(ckpt.get("failed", []))
log.info(f"Loading {len(failed_names)} failed compounds from round 5")

# Get SMILES from original retry CSV
df = pd.read_csv("03_Analysis/outputs/fixed_retry.csv")
df_failed = df[df["name"].isin(failed_names)].copy()
log.info(f"Matched {len(df_failed)} compounds with SMILES")

# ═══ Prepare using OpenBabel ═══
def prepare_pdbqt_obabel(smiles, name):
    """SMILES → 3D SDF → PDBQT using OpenBabel."""
    safe = name.replace("/","_").replace("\\","_").replace(" ","_").replace("*","_")
    pdbqt = os.path.join(DOCKING_DIR, f"{safe}_ob.pdbqt")

    if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 100:
        return pdbqt

    # Step 1: SMILES → 3D SDF with OpenBabel
    sdf = os.path.join(DOCKING_DIR, f"{safe}_ob.sdf")

    try:
        # Use pybel for robust 3D generation
        mol = pybel.readstring("smi", smiles)
        mol.addh()  # Add hydrogens
        mol.make3D(forcefield="mmff94", steps=500)
        mol.localopt(forcefield="mmff94", steps=200)
        mol.write("sdf", sdf, overwrite=True)
    except Exception as e:
        # Fallback: RDKit 3D → SDF → obabel to PDBQT
        log.warning(f"  pybel 3D failed for {name}: {e}, trying RDKit fallback")
        mol_r = Chem.MolFromSmiles(smiles)
        if mol_r is None:
            raise ValueError(f"Cannot parse SMILES: {smiles[:50]}")
        mol_r = Chem.AddHs(mol_r)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol_r, params)
        if status != 0:
            params.useRandomCoords = True
            AllChem.EmbedMolecule(mol_r, params)
        AllChem.MMFFOptimizeMolecule(mol_r)
        w = Chem.SDWriter(sdf)
        w.write(mol_r)
        w.close()

    if not os.path.exists(sdf) or os.path.getsize(sdf) < 50:
        raise RuntimeError(f"SDF generation failed")

    # Step 2: SDF → PDBQT with obabel (Gasteiger charges + AD4 atom types)
    result = subprocess.run(
        ["obabel", sdf, "-O", pdbqt, "--gen3d", "--partialcharge", "gasteiger"],
        capture_output=True, text=True, timeout=30
    )

    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
        # Try without --partialcharge
        result = subprocess.run(
            ["obabel", sdf, "-O", pdbqt],
            capture_output=True, text=True, timeout=30
        )

    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) < 100:
        raise RuntimeError(f"obabel PDBQT conversion failed")

    return pdbqt


def dock_one(name, smiles):
    """Dock single compound with OpenBabel-prepared PDBQT."""
    safe = name.replace("/","_").replace("\\","_").replace(" ","_").replace("*","_")

    try:
        lig_pdbqt = prepare_pdbqt_obabel(smiles, name)
        out_pdbqt = os.path.join(DOCKING_DIR, f"{safe}_ob_docked.pdbqt")

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
            "--num_modes", "5",
            "--out", out_pdbqt,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            return {"name": name, "status": "vina_error", "error": result.stderr[:200]}

        # Parse affinities
        affinities = []
        for line in result.stdout.split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    v = float(parts[1])
                    if abs(v) < 100:
                        affinities.append(round(v, 3))
                except:
                    continue

        if not affinities:
            # Try parsing from output file
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        for p in line.split():
                            try:
                                v = float(p)
                                if abs(v) < 100:
                                    affinities.append(round(v, 3))
                            except:
                                pass

        if not affinities:
            return {"name": name, "status": "parse_error"}

        best = min(affinities)
        return {"name": name, "best_affinity": best, "all_affinities": sorted(set(affinities)), "status": "success"}

    except subprocess.TimeoutExpired:
        return {"name": name, "status": "timeout"}
    except Exception as e:
        return {"name": name, "status": "error", "error": str(e)[:200]}


# ═══ Main ═══
log.info(f"\n{'='*60}")
log.info(f"Redocking {len(df_failed)} compounds with OpenBabel (Vina ex={EXHAUSTIVENESS})")
log.info(f"{'='*60}")

# Load/resume checkpoint
if os.path.exists(CHECKPOINT):
    with open(CHECKPOINT) as f:
        ob_ckpt = json.load(f)
else:
    ob_ckpt = {"completed": [], "failed": [], "results": []}

done_set = set(ob_ckpt.get("completed", [])) | set(ob_ckpt.get("failed", []))

start = time.time()
for i, (_, row) in enumerate(df_failed.iterrows()):
    name, smi = row["name"], row["smiles"]
    if name in done_set:
        continue

    result = dock_one(name, smi)

    if result.get("status") == "success":
        aff = result["best_affinity"]
        ob_ckpt["completed"].append(name)
        ob_ckpt["results"].append(result)
        hit_flag = " 🏆 HIT!" if aff < -7.0 else (" ⭐ TOP" if aff < -6.5 else "")
        log.info(f"  [{i+1}/{len(df_failed)}] {name[:30]:<30} ΔG={aff:.2f}{hit_flag}")
    else:
        ob_ckpt["failed"].append(name)
        log.warning(f"  [{i+1}/{len(df_failed)}] {name[:30]:<30} FAILED: {result.get('status')}")

    # Save checkpoint periodically
    if (i + 1) % 10 == 0:
        with open(CHECKPOINT, "w") as f:
            json.dump(ob_ckpt, f, indent=2)

# Final save
with open(CHECKPOINT, "w") as f:
    json.dump(ob_ckpt, f, indent=2)

# Summary
elapsed = (time.time() - start) / 60
n_ok = len(ob_ckpt["completed"])
n_fail = len(ob_ckpt["failed"])
log.info(f"\n{'='*60}")
log.info(f"COMPLETE: {n_ok} docked, {n_fail} failed in {elapsed:.1f}min")

if ob_ckpt["results"]:
    affs = [r["best_affinity"] for r in ob_ckpt["results"]]
    df_r = pd.DataFrame(ob_ckpt["results"]).sort_values("best_affinity")
    df_r.to_csv(RESULTS_CSV, index=False)

    best = df_r.iloc[0]
    log.info(f"Best: {best['name']} = {best['best_affinity']:.2f} kcal/mol")
    log.info(f"Hits (<-7.0): {sum(1 for a in affs if a < -7.0)}")
    log.info(f"Top (<-6.5): {sum(1 for a in affs if a < -6.5)}")

    # Show top 10
    log.info(f"\nTop 10:")
    for i, (_, r) in enumerate(df_r.head(10).iterrows()):
        log.info(f"  {i+1}. {r['name'][:35]:<35} ΔG={r['best_affinity']:.2f}")
