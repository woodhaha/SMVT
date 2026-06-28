#!/usr/bin/env python3
"""
SMVT MD — Top 3 + Controls, 50ns each, GPU
============================================
Uses OpenFF + SMIRNOFF (amber/ff14SB + openff-2.1.0 + tip3p)
Based on tested md_binding_stability.py from SMVT project.
"""

import os, sys, json, time, argparse, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR_RAW = os.path.join(PACKAGE_DIR, "receptor", "AF-Q9Y289-F1.pdb")
RECEPTOR_PREPARED = os.path.join(PACKAGE_DIR, "receptor", "SMVT_prepared.pdb")
LIGAND_DIR = os.path.join(PACKAGE_DIR, "ligands")
TRAJ_DIR = os.path.join(PACKAGE_DIR, "trajectories")

COMPOUNDS = json.load(open(os.path.join(PACKAGE_DIR, "compounds.json")))["compounds"]
os.makedirs(TRAJ_DIR, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def extract_vina_pose(pdbqt_path):
    """Extract MODEL 1 atom positions from Vina PDBQT."""
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
    """Add missing atoms/hydrogens via PDBFixer."""
    from pdbfixer import PDBFixer
    import openmm.app as app
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


def run_md(compound_id: str, production_ns: float = 50):
    """Run full MD for one compound."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmmforcefields.generators import SystemGenerator
    from rdkit import Chem
    from rdkit.Chem import AllChem

    comp = next(c for c in COMPOUNDS if c["id"] == compound_id)
    out_dir = os.path.join(TRAJ_DIR, compound_id)
    os.makedirs(out_dir, exist_ok=True)

    chk_file = os.path.join(out_dir, f"{compound_id}_prod.chk")
    if os.path.exists(chk_file):
        log.info(f"[{compound_id}] Already done. Skip.")
        return {"name": compound_id, "status": "already_done"}

    log.info(f"{'='*50}\n[{compound_id}] {comp['name']} — {comp['dG']:.2f} kcal/mol — {production_ns}ns\n{'='*50}")

    try:
        # ── 1. Prepare protein ──
        protein_clean = prepare_protein(RECEPTOR_PREPARED if os.path.exists(RECEPTOR_PREPARED) else RECEPTOR_RAW,
                                         os.path.join(out_dir, "protein_prepared.pdb"))
        log.info("  Protein prepared")

        # ── 2. Generate ligand 3D conformer (RDKit ETKDGv3) ──
        rd_mol = Chem.MolFromSmiles(comp["smiles"])
        if rd_mol is None:
            return {"name": compound_id, "status": "failed", "reason": "bad_smiles"}
        rd_mol = Chem.AddHs(rd_mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(rd_mol, params)
        if status != 0:
            params.useRandomCoords = True
            status = AllChem.EmbedMolecule(rd_mol, params)
        if status != 0:
            return {"name": compound_id, "status": "failed", "reason": "no_conformer"}
        AllChem.MMFFOptimizeMolecule(rd_mol, maxIters=500)

        # Get conformer positions
        conf = rd_mol.GetConformer()
        n_atoms_rd = rd_mol.GetNumAtoms()
        rd_positions = np.array([list(conf.GetAtomPosition(i)) for i in range(n_atoms_rd)])
        log.info(f"  RDKit conformer: {n_atoms_rd} atoms")

        # Save as SDF
        lig_sdf = os.path.join(out_dir, "ligand.sdf")
        w = Chem.SDWriter(lig_sdf); w.write(rd_mol); w.close()

        # Convert to PDB via obabel
        lig_pdb = os.path.join(out_dir, "ligand.pdb")
        import subprocess
        subprocess.run(["obabel", lig_sdf, "-O", lig_pdb], capture_output=True, timeout=30)
        if not os.path.exists(lig_pdb) or os.path.getsize(lig_pdb) < 50:
            return {"name": compound_id, "status": "failed", "reason": "bad_pdb"}

        # ── 3. Parameterize with GAFF (amber/ff14SB + gaff-2.11) ──
        log.info(f"  Ligand prepared: {n_atoms_rd} atoms")

        # ── 4. Build system ──
        log.info("  Building system (solvation + FF)...")
        forcefield_kwargs = {
            'constraints': app.HBonds,
            'rigidWater': True,
            'hydrogenMass': 1.5 * unit.amu,
        }

        system_generator = SystemGenerator(
            forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
            small_molecule_forcefield='gaff-2.11',
            forcefield_kwargs=forcefield_kwargs,
        )

        # Modeller: protein + ligand
        protein_pdb = app.PDBFile(protein_clean)
        ligand_pdb = app.PDBFile(lig_pdb)
        modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
        modeller.add(ligand_pdb.topology, ligand_pdb.positions)
        modeller.addSolvent(
            system_generator.forcefield,
            padding=12.0 * unit.angstrom,
            ionicStrength=0.15 * unit.molar,
            model='tip3p',
        )

        system = system_generator.create_system(modeller.topology, molecules=ligand_pdb.topology)
        n_atoms = modeller.topology.getNumAtoms()
        log.info(f"  System: {n_atoms} atoms")

        # ── 5. Platform ──
        platform = mm.Platform.getPlatformByName("CUDA")
        log.info("  Platform: CUDA")

        # ── 6. Minimize ──
        log.info("  Minimizing...")
        integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=5000)

        state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        log.info(f"  Minimized: {energy:.0f} kJ/mol")

        # Save minimized
        with open(os.path.join(out_dir, "minimized.pdb"), "w") as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

        # ── 7. NVT 100ps ──
        log.info("  NVT equilibration (100ps)...")
        pos0 = state.getPositions()
        # Restrain protein backbone
        restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", 5.0*unit.kilocalories_per_mole/unit.angstrom**2)
        restraint.addPerParticleParameter("x0"); restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        protein_atoms = [a.index for a in modeller.topology.atoms()
                        if a.residue.name not in ("HOH","NA","CL") and a.element.symbol != "H"]
        for idx in protein_atoms:
            if idx < len(pos0):
                restraint.addParticle(idx, pos0[idx].value_in_unit(unit.nanometer))
        system.addForce(restraint)
        simulation.context.reinitialize(preserveState=True)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(50000)  # 100ps
        log.info("  NVT done")

        # ── 8. NPT 200ps ──
        log.info("  NPT equilibration (200ps)...")
        system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300*unit.kelvin, 25))
        simulation.context.setParameter("k", 2.0*unit.kilocalories_per_mole/unit.angstrom**2)
        simulation.context.reinitialize(preserveState=True)
        simulation.step(100000)  # 200ps
        log.info("  NPT done")

        # ── 9. Production ──
        log.info(f"  Production ({production_ns}ns)...")
        system.removeForce(len(system.getForces())-2)  # remove old barostat and restraint
        system.removeForce(len(system.getForces())-1)
        system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300*unit.kelvin, 25))
        simulation.context.reinitialize(preserveState=True)

        prod_steps = int(production_ns * 1e6 / 2)
        save_interval = 25000  # 50ps

        simulation.reporters.clear()
        simulation.reporters.append(app.DCDReporter(
            os.path.join(out_dir, "production.dcd"), save_interval))
        simulation.reporters.append(app.StateDataReporter(
            os.path.join(out_dir, "production_log.txt"), 5000,
            step=True, temperature=True, potentialEnergy=True, density=True, speed=True))
        simulation.reporters.append(app.CheckpointReporter(
            os.path.join(out_dir, f"{compound_id}_prod.chk"), 5000000))

        t0 = time.time()
        simulation.step(prod_steps)
        elapsed = time.time() - t0

        # Save final
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        with open(os.path.join(out_dir, "final.pdb"), "w") as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

        ns_per_day = production_ns/(elapsed/86400)
        log.info(f"  Done in {elapsed/3600:.1f}h ({ns_per_day:.0f} ns/day)")

        return {"name": compound_id, "status": "success",
                "elapsed_h": elapsed/3600, "ns_per_day": ns_per_day}

    except Exception as e:
        log.error(f"  FAILED: {e}")
        import traceback; traceback.print_exc()
        return {"name": compound_id, "status": "failed", "reason": str(e)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compound", type=str, help="Compound ID")
    parser.add_argument("--ns", type=float, default=50, help="Production ns")
    args = parser.parse_args()

    if args.compound:
        result = run_md(args.compound, args.ns)
        print(json.dumps(result, indent=2))
    else:
        parser.print_help()
