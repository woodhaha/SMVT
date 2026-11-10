#!/usr/bin/env python3
"""Stage A — dock Drug Repurposing Hub candidates against SMVT.

Uses docking_config.json values (center 5.03,1.5,-2.44 / box 25^3 /
exhaustiveness 16 / num_modes 9). Ligand prep = RDKit ETKDG + MMFF -> meeko PDBQT.
Resumable via checkpoint JSON keyed by InChIKey.

Usage:
  python docking_rephub.py --input candidates.csv --output results.csv [--test N]
"""
import os, sys, subprocess, time, json, argparse, warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp

warnings.filterwarnings("ignore")
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ── config ──
BASE = os.path.dirname(os.path.abspath(__file__))
VINA_BIN = os.path.join(BASE, "tools", "vina.exe")
MK_PREP_LIG = "C:/anaconda3/Scripts/mk_prepare_ligand.exe"
RECEPTOR_PDBQT = os.path.join(BASE, "docking", "SMVT_receptor.pdbqt")
DOCK_DIR = os.path.join(BASE, "docking_rephub")
CENTER = [-2.5, 1.0, -1.0]
BOX = [22.0, 22.0, 22.0]
EXHAUSTIVENESS = 16
NUM_MODES = 5
NUM_WORKERS = max(1, mp.cpu_count() - 2)

os.makedirs(DOCK_DIR, exist_ok=True)


def prepare_ligand(smiles, key):
    """SMILES -> 3D SDF -> meeko PDBQT. Returns PDBQT path or raises."""
    sdf = os.path.join(DOCK_DIR, f"{key}.sdf")
    pdbqt = os.path.join(DOCK_DIR, f"{key}.pdbqt")
    if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 100:
        return pdbqt

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit cannot parse SMILES")
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    st = AllChem.EmbedMolecule(mol, p)
    if st != 0:
        p.useRandomCoords = True
        st = AllChem.EmbedMolecule(mol, p)
        if st != 0:
            raise ValueError("cannot embed 3D")
    AllChem.MMFFOptimizeMolecule(mol)
    w = Chem.SDWriter(sdf); w.write(mol); w.close()

    r = subprocess.run([MK_PREP_LIG, "-i", sdf, "-o", pdbqt],
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0 or not os.path.exists(pdbqt):
        raise RuntimeError(f"meeko prep failed: {r.stderr[:150]}")
    return pdbqt


def dock_one(key, name, smiles):
    try:
        lig = prepare_ligand(smiles, key)
        out = os.path.join(DOCK_DIR, f"{key}_docked.pdbqt")
        cmd = [VINA_BIN, "--receptor", RECEPTOR_PDBQT, "--ligand", lig,
               "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
               "--center_z", str(CENTER[2]),
               "--size_x", str(BOX[0]), "--size_y", str(BOX[1]), "--size_z", str(BOX[2]),
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
               "--out", out]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            return {"key": key, "name": name, "status": "vina_error", "error": r.stderr[:150]}

        affs = []
        for line in r.stdout.split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    v = float(parts[1])
                    if abs(v) < 100:
                        affs.append(round(v, 3))
                except (ValueError, IndexError):
                    continue
        if not affs:
            try:
                with open(out) as f:
                    for line in f:
                        if "VINA RESULT" in line:
                            for tok in line.split():
                                try:
                                    affs.append(round(float(tok), 3))
                                except ValueError:
                                    pass
            except Exception:
                pass
        if not affs:
            return {"key": key, "name": name, "status": "parse_error", "error": "no affinity"}

        return {"key": key, "name": name, "best_affinity": min(affs),
                "all_affinities": sorted(set(affs)), "status": "success"}
    except subprocess.TimeoutExpired:
        return {"key": key, "name": name, "status": "timeout"}
    except Exception as e:
        return {"key": key, "name": name, "status": "error", "error": str(e)[:150]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--checkpoint", default=None)
    ap.add_argument("--test", type=int, default=0)
    ap.add_argument("--workers", type=int, default=NUM_WORKERS)
    args = ap.parse_args()

    cp = args.checkpoint or (args.output + ".checkpoint.json")
    df = pd.read_csv(args.input, dtype=str)
    if args.test:
        df = df.head(args.test)

    if os.path.exists(cp):
        ck = json.load(open(cp))
    else:
        ck = {"completed": [], "failed": [], "results": []}
    done = set(ck["completed"]) | set(ck["failed"])
    df = df[~df["inchikey"].isin(done)]
    print(f"loaded {len(df)} pending (done {len(done)})", flush=True)

    if len(df) == 0:
        print("all done"); return

    t0 = time.time()
    rows = [(r["inchikey"], r["name"], r["smiles"]) for _, r in df.iterrows()]
    n = len(rows)

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, k, nm, sm): (k, nm) for k, nm, sm in rows}
        for fut in as_completed(futs):
            res = fut.result()
            if res["status"] == "success":
                ck["completed"].append(res["key"])
                ck["results"].append({"key": res["key"], "name": res["name"],
                                      "best_affinity": res["best_affinity"],
                                      "all_affinities": res["all_affinities"]})
            else:
                ck["failed"].append(res["key"])
            np_ = len(ck["completed"]) + len(ck["failed"])
            el = max(0.01, (time.time() - t0) / 60)
            if np_ % 20 == 0 or np_ <= 5:
                aff = res.get("best_affinity", "N/A")
                print(f"  [{np_}/{n}] {res['name'][:28]:28s} {aff} | {np_/el:.1f} cpd/min | ETA {(n-np_)/max(0.01,np_/el):.0f} min", flush=True)
            json.dump(ck, open(cp, "w"))

    el = (time.time() - t0) / 60
    affs = [r["best_affinity"] for r in ck["results"] if r["best_affinity"] is not None]
    print(f"\nDONE: {len(ck['completed'])} ok / {len(ck['failed'])} fail / {el:.1f} min")
    if affs:
        print(f"best {min(affs):.2f} | mean {sum(affs)/len(affs):.2f} | hits<-8 {sum(1 for a in affs if a < -8)} | <-7 {sum(1 for a in affs if a < -7)}")
        pd.DataFrame(ck["results"]).sort_values("best_affinity").to_csv(args.output, index=False)
        print(f"results -> {args.output}")


if __name__ == "__main__":
    main()
