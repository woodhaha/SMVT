#!/usr/bin/env python3
"""
Phase B — MD Binding Stability Assessment for Top 10 SMVT Hits
================================================================
Requires: conda activate smvt-md  (Python 3.11 + openff-toolkit + openmmforcefields)
Protocol: Minim → NVT (100ps) → NPT (200ps) → Production (10ns)
Top 3 (Naftazone, Phenobarbital, Esketamine): 50ns extended
"""

import os, sys, json, time, logging, argparse, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule as OFFMolecule
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

# ═══ Config ═══
PROJECT_DIR = "D:/Researching/SMVT"
MD_DIR = f"{PROJECT_DIR}/03_Analysis/md"
RECEPTOR_PDB = f"{PROJECT_DIR}/02_Data/raw/AF-Q9Y289-F1.pdb"
DOCKING_DIR = f"{PROJECT_DIR}/03_Analysis/docking"
CHECKPOINT = f"{MD_DIR}/md_checkpoint.json"

# MD parameters
TEMP = 300 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
FRICTION = 1.0 / unit.picosecond
DT = 2.0 * unit.femtoseconds
CUTOFF = 1.0 * unit.nanometer
PADDING = 0.8 * unit.nanometers
IONIC = 0.15 * unit.molar
EQ_NVT_PS, EQ_NPT_PS = 100, 100
SAVE_PS = 50  # Save trajectory every 50 ps

# Top 10 hits (updated 2026-06-27 with FDA leftover results)
TOP10 = [
    ("HYDROMORPHONE",       "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5", -8.58),
    ("FUROSEMIDE",          "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl", -8.36),
    ("NAFTAZONE",           "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N", -8.34),
    ("PHENOBARBITAL",       "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2", -8.30),
    ("LENALIDOMIDE",        "Nc1cccc2c1CN(C1CCC(=O)NC1=O)C2=O", -8.25),
    ("BUFEXAMAC",           "CCCCOc1ccc(CC(=O)NO)cc1", -8.06),
    ("OXYMORPHONE",         "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@]3(O)[C@H]1C5", -8.04),
    ("TOLOXATONE",          "Cc1cccc(N2CC(CO)OC2=O)c1", -8.02),
    ("AVIBACTAM",           "NC(=O)[C@@H]1CC[C@@H]2CN1C(=O)N2OS(=O)(=O)O", -7.95),
    ("CYCLOBARBITAL",       "C1CC2=C(C1)C(=O)NC(=O)N2", -7.83),
]
TOP3_SET = {"HYDROMORPHONE", "NAFTAZONE", "PHENOBARBITAL"}  # Extended 50ns MD targets

os.makedirs(MD_DIR, exist_ok=True)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(f"{MD_DIR}/md_pipeline.log", mode="w"),
              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)


def extract_vina_pose(pdbqt_path):
    """Extract MODEL 1 atom positions from Vina PDBQT. Returns list of [x,y,z] in Angstrom."""
    pos = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if int(line.split()[1]) == 1:
                    in_model = True; pos = []
                elif in_model:
                    break
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except (ValueError, IndexError):
                    continue
    return pos


def prepare_protein(pdb_path, out_path):
    """PDBFixer: add missing atoms/hydrogens, save cleaned PDB."""
    if os.path.exists(out_path):
        return out_path
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return out_path


def run_md(name, smiles, prod_ns):
    """Full MD pipeline for one protein-ligand complex."""
    tag = name
    out_dir = f"{MD_DIR}/trajectories/{tag}"
    os.makedirs(out_dir, exist_ok=True)

    chk_file = f"{out_dir}/{tag}_prod.chk"
    if os.path.exists(chk_file):
        log.info(f"[{name}] Already done. Skipping.")
        return {"name": name, "status": "already_done"}

    log.info(f"{'='*50}\n[{name}] {prod_ns}ns MD\n{'='*50}")

    try:
        # ── 1. Prepare protein ──
        protein_clean = f"{out_dir}/protein_prepared.pdb"
        prepare_protein(RECEPTOR_PDB, protein_clean)
        log.info(f"  Protein prepared")

        # ── 2. Extract Vina pose ──
        pdbqt = f"{DOCKING_DIR}/{name}_docked.pdbqt"
        vina_pos = extract_vina_pose(pdbqt)
        if not vina_pos:
            return {"name": name, "status": "failed", "reason": "empty_pose"}
        log.info(f"  Vina pose: {len(vina_pos)} atoms")

        # ── 3. Parameterize ligand with OpenFF ──
        # Create OpenFF Molecule from SMILES, pre-assign Gasteiger charges via RDKit
        off_mol = OFFMolecule.from_smiles(smiles, allow_undefined_stereo=True)
        off_mol.generate_conformers(n_conformers=1)
        # Set Vina coordinates
        off_conf = off_mol.conformers[0]
        from openff.units import unit as off_unit
        for i in range(min(off_mol.n_atoms, len(vina_pos))):
            off_conf[i] = np.array(vina_pos[i]) * off_unit.angstrom
        # Pre-assign charges so SMIRNOFF doesn't try am1bcc
        off_mol.assign_partial_charges(partial_charge_method='gasteiger')

        # Save ligand as PDB (not SDF) for OpenMM Modeller
        lig_pdb_path = f"{out_dir}/ligand.pdb"
        rd_mol = Chem.MolFromSmiles(smiles)
        rd_mol = Chem.AddHs(rd_mol)
        AllChem.EmbedMolecule(rd_mol, randomSeed=42)
        conf = rd_mol.GetConformer()
        # Overwrite with Vina coordinates
        for i, pos_angstrom in enumerate(vina_pos):
            if i < rd_mol.GetNumAtoms():
                conf.SetAtomPosition(i, (pos_angstrom[0]/10.0, pos_angstrom[1]/10.0, pos_angstrom[2]/10.0))
        Chem.MolToPDBFile(rd_mol, lig_pdb_path)

        # SystemGenerator with OpenFF Sage + RDKit Gasteiger charges
        log.info(f"  Parameterizing ligand with OpenFF Sage (Gasteiger charges)...")
        system_generator = SystemGenerator(
            forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
            small_molecule_forcefield='openff-2.1.0',
            molecules=[off_mol],
            cache=f"{out_dir}/openff_cache.json"
        )

        # ── 4. Build solvated system ──
        log.info(f"  Building system...")
        prot_pdb = app.PDBFile(protein_clean)
        lig_pdb = app.PDBFile(lig_pdb_path)

        modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
        modeller.add(lig_pdb.topology, lig_pdb.positions)
        modeller.addSolvent(system_generator.forcefield, model='tip3p',
                            padding=PADDING, ionicStrength=IONIC, neutralize=True)

        system = system_generator.create_system(modeller.topology, molecules=[off_mol])
        n_particles = system.getNumParticles()
        log.info(f"  System: {n_particles} particles")

        # ── 5. Integrator + Simulation ──
        integrator = mm.LangevinMiddleIntegrator(TEMP, FRICTION, DT)
        # Use fastest available platform
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
            props = {'Precision': 'mixed'}
        except:
            platform = mm.Platform.getPlatformByName('CPU')
            props = {}
        simulation = app.Simulation(modeller.topology, system, integrator, platform, props)
        simulation.context.setPositions(modeller.positions)

        # ── 6. Minimization ──
        log.info(f"  Minimizing (5000 steps)...")
        simulation.minimizeEnergy(maxIterations=5000)

        # ── 7. NVT with slow heating (50K→300K, compact) ──
        log.info(f"  NVT equil (50K→300K heating)...")
        simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
        temps = [50, 150, 300]
        batch_steps = int(10 * 1000 / 1.0)  # 10ps per batch at 1fs
        for t in temps:
            simulation.integrator.setTemperature(t * unit.kelvin)
            simulation.step(batch_steps)
        simulation.integrator.setTemperature(TEMP)
        simulation.integrator.setStepSize(DT)
        simulation.step(int(70 * 1000 / 2.0))  # remaining 70ps at 2fs

        # ── 8. NPT ──
        log.info(f"  NPT equil ({EQ_NPT_PS}ps)...")
        system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
        simulation.context.reinitialize(preserveState=True)
        npt_steps = int(EQ_NPT_PS * 1000 / (DT / unit.picosecond))
        simulation.step(npt_steps)

        # ── 9. Production ──
        prod_steps = int(prod_ns * 1_000_000 / (DT / unit.picosecond))
        save_freq = int(SAVE_PS * 1000 / (DT / unit.picosecond))

        log.info(f"  Production: {prod_ns}ns ({prod_steps} steps)...")
        simulation.reporters.append(app.DCDReporter(f"{out_dir}/{tag}_prod.dcd", save_freq))
        simulation.reporters.append(app.StateDataReporter(
            f"{out_dir}/{tag}_prod.csv", save_freq,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            temperature=True, volume=True, density=True
        ))
        simulation.reporters.append(app.CheckpointReporter(chk_file, save_freq * 10))

        t0 = time.time()
        simulation.step(prod_steps)
        elapsed = (time.time() - t0) / 60

        # Save final frame
        final = simulation.context.getState(getPositions=True)
        with open(f"{out_dir}/{tag}_final.pdb", "w") as f:
            app.PDBFile.writeFile(simulation.topology, final.getPositions(), f)
        simulation.saveCheckpoint(chk_file)

        rate = prod_ns / elapsed if elapsed > 0 else 0
        log.info(f"  [{name}] DONE in {elapsed:.1f}min ({rate:.1f} ns/min)")
        return {"name": name, "status": "completed", "time_min": round(elapsed, 1),
                "ns_per_min": round(rate, 2)}

    except Exception as e:
        log.error(f"  [{name}] FAILED: {e}", exc_info=True)
        return {"name": name, "status": "failed", "reason": str(e)[:200]}


# ═══════════════════════════════════════════════════════════════
def main():
    os.chdir(PROJECT_DIR)
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", type=int, default=0, help="Test first N at 1ns")
    parser.add_argument("--prod-ns", type=int, default=0, help="Override production ns")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    if os.path.exists(CHECKPOINT):
        ckpt = json.load(open(CHECKPOINT))
    else:
        ckpt = {"completed": [], "failed": [], "results": []}

    completed = set(ckpt.get("completed", []))

    if args.test:
        compounds = [(n, s, a) for n, s, a in TOP10[:args.test]]
        default_ns = 1
    else:
        compounds = [(n, s, a) for n, s, a in TOP10 if n not in completed]
        default_ns = args.prod_ns or 10

    log.info(f"Phase B: {len(compounds)} compounds | Platform: OpenCL | Default: {default_ns}ns")
    for name, smiles, score in compounds:
        if args.test:
            prod_ns = args.prod_ns or 1
        elif args.prod_ns > 0:
            prod_ns = args.prod_ns
        else:
            prod_ns = 2  # 2ns for all on CPU
        result = run_md(name, smiles, prod_ns)

        if result["status"] in ("completed", "already_done"):
            ckpt["completed"].append(name)
            ckpt.setdefault("results", []).append({
                "name": name, "docking_ΔG": score, "md_ns": prod_ns,
                "time_min": result.get("time_min"),
            })
        elif result["status"] == "failed":
            ckpt["failed"].append(name)

        json.dump(ckpt, open(CHECKPOINT, "w"), indent=2)

    n_ok = len(ckpt.get("completed", []))
    n_fail = len(ckpt.get("failed", []))
    log.info(f"\nPhase B done. {n_ok} completed, {n_fail} failed.")
    log.info(f"Checkpoint: {CHECKPOINT}")


if __name__ == "__main__":
    main()
