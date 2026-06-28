"""
SMVT Top-3 MD Simulation — Colab T4 GPU — GAFF 2.11 (SMIRNOFF types)

Uses pre-generated GAFF templates (ligand_params_v2/) with SMIRNOFF-derived
atom types + MMFF94 charges.  This replaces the heuristic RDKit-based typing
that caused NaN during NVT heating.

Protocol: Minim 5000 → NVT 50→300K → NPT 100ps → Production 50ns
"""

import os, time, json
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

# ── GPU detection ──
plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        break
    except:
        continue

print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

# ── Paths ──
DATA_DIR = "/content"
OUT_DIR  = "/content/trajectories"
os.makedirs(OUT_DIR, exist_ok=True)
RECEPTOR_PDB = os.path.join(DATA_DIR, "AF-Q9Y289-F1.pdb")

# ── Top 3 ligands ──
TOP3 = [
    ("NAFTAZONE",     "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N",
     "NAFTAZONE_template.xml"),
    ("PHENOBARBITAL", "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2",
     "PHENOBARBITAL_template.xml"),
    ("ESKETAMINE",    "CNC1(C2=CC=CC=C2Cl)CCCCC1=O",
     "ESKETAMINE_template.xml"),
]

# ── MD parameters ──
PROD_NS   = 50
TEMP      = 300.0 * unit.kelvin
PRESSURE  = 1.0 * unit.atmosphere
DT        = 2.0 * unit.femtoseconds
PADDING   = 1.0 * unit.nanometers

print(f"\n{'='*60}")
print(f"SMVT MD: {len(TOP3)} compounds × {PROD_NS}ns on {plat.getName()}")
print(f"{'='*60}")


# ═══════════════════════════════════════════════════
#  Utilities
# ═══════════════════════════════════════════════════

def extract_vina_pose(pdbqt_path):
    """Extract coordinates from Vina docking pose (first MODEL only)."""
    pos = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if int(line.split()[1]) == 1:
                    in_model = True
                    pos = []
                elif in_model:
                    break
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    pos.append([float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])])
                except:
                    continue
    return pos


def prepare_protein(pdb_path, out_path):
    """PDBFixer: missing atoms + hydrogens at pH 7.4."""
    if os.path.exists(out_path):
        return out_path
    print("  Preparing protein (PDBFixer)...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return out_path


def prepare_ligand_pdb(smiles, vina_pos, out_path, res_name="LIG"):
    """Generate ligand PDB from SMILES + Vina pose coordinates."""
    if os.path.exists(out_path):
        return out_path
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    n = min(mol.GetNumAtoms(), len(vina_pos))
    for i in range(n):
        conf.SetAtomPosition(i, (vina_pos[i][0], vina_pos[i][1], vina_pos[i][2]))
    pdb_str = Chem.MolToPDBBlock(mol)
    # Rename residue
    fixed = []
    for line in pdb_str.split("\n"):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            line = line[:17] + res_name + line[20:]
        fixed.append(line)
    with open(out_path, "w") as f:
        f.write("\n".join(fixed))
    return out_path


# ═══════════════════════════════════════════════════
#  Main MD loop
# ═══════════════════════════════════════════════════

def run_gpu_md(name, smiles, template_file):
    tag = name
    out_dir = os.path.join(OUT_DIR, tag)
    os.makedirs(out_dir, exist_ok=True)

    chk_file = os.path.join(out_dir, f"{tag}_{PROD_NS}ns.chk")
    if os.path.exists(chk_file):
        print(f"\n[{name}] Already completed. Skip.")
        return

    sep = "=" * 50
    print(f"\n{sep}\n[{name}] {PROD_NS}ns MD (GAFF2 + SMIRNOFF types)\n{sep}")
    t_start = time.time()

    # 1. Prepare protein
    prot_clean = os.path.join(out_dir, "protein.pdb")
    prepare_protein(RECEPTOR_PDB, prot_clean)

    # 2. Extract Vina pose
    pdbqt = os.path.join(DATA_DIR, f"{name}_docked.pdbqt")
    vina_pos = extract_vina_pose(pdbqt)
    print(f"  Vina pose: {len(vina_pos)} atoms")

    # 3. Generate ligand PDB with Vina pose
    lig_pdb = os.path.join(out_dir, "ligand.pdb")
    prepare_ligand_pdb(smiles, vina_pos, lig_pdb, "LIG")

    # 4. Load force field: AMBER ff14SB + GAFF 2.11 + custom template + TIP3P
    template_xml = os.path.join(DATA_DIR, template_file)
    gaff_xml = (
        "/usr/local/lib/python3.12/dist-packages/openmmforcefields/"
        "ffxml/amber/gaff/ffxml/gaff-2.11.xml"
    )

    print(f"  ForceField: amber14 + gaff-2.11 + {template_file} + tip3p")
    forcefield = app.ForceField(
        "amber14-all.xml",
        gaff_xml,
        template_xml,
        "tip3p.xml",
    )

    # 5. Build system
    print("  Building system...")
    prot_pdb = app.PDBFile(prot_clean)
    lig = app.PDBFile(lig_pdb)
    modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
    modeller.add(lig.topology, lig.positions)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=PADDING,
        ionicStrength=0.0 * unit.molar,
        neutralize=False,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        ignoreExternalBonds=True,
    )
    n_particles = system.getNumParticles()
    print(f"  System: {n_particles} particles")

    # 6. Backbone position restraints (equilibration only)
    k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    restraint_force = mm.CustomExternalForce(
        "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    restraint_force.addGlobalParameter("k", k_rest)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    backbone = {"CA", "C", "N"}
    for atom in modeller.topology.atoms():
        if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
            if atom.name in backbone:
                p = modeller.positions[atom.index]
                restraint_force.addParticle(atom.index, [p.x, p.y, p.z])
    restraint_idx = system.addForce(restraint_force)

    # 7. Integrator + Simulation
    integrator = mm.LangevinMiddleIntegrator(
        TEMP, 1.0 / unit.picosecond, DT
    )
    simulation = app.Simulation(
        modeller.topology, system, integrator, plat,
        {"Precision": "mixed"},
    )
    simulation.context.setPositions(modeller.positions)

    # ── Phase 1: Minimize (restrained) ──
    print("  Minimizing (restrained, 5000 steps)...")
    simulation.minimizeEnergy(maxIterations=5000)

    # ── Phase 2: NVT heating 50→300K (1fs, restrained) ──
    print("  NVT heating (50→300K, 1fs, restrained)...")
    simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
    for t in [50, 100, 150, 200, 250, 300]:
        simulation.integrator.setTemperature(t * unit.kelvin)
        simulation.step(10_000)  # 10ps per stage

    simulation.integrator.setTemperature(TEMP)
    simulation.step(50_000)  # 50ps more at 300K

    # ── Phase 3: NPT equil (reducing restraints) ──
    print("  NPT equil (reducing restraints)...")
    simulation.integrator.setStepSize(DT)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
    simulation.context.reinitialize(preserveState=True)

    for k_val in [2.5, 0.5]:
        new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
        simulation.context.setParameter("k", new_k)
        simulation.step(25_000)  # 50ps each @ 2fs

    # ── Phase 4: Remove restraints, final equil ──
    print("  Removing restraints, final equil...")
    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(50_000)  # 100ps unrestrained equil

    # ── Phase 5: Production ──
    prod_steps = int(PROD_NS * 1_000_000 / 2.0)
    save_freq  = 50_000  # every 100ps
    print(f"  Production: {PROD_NS}ns ({prod_steps} steps)...")

    simulation.reporters.append(
        app.DCDReporter(
            os.path.join(out_dir, f"{tag}_{PROD_NS}ns.dcd"),
            save_freq,
        )
    )
    simulation.reporters.append(
        app.StateDataReporter(
            os.path.join(out_dir, f"{tag}_{PROD_NS}ns.csv"),
            save_freq,
            step=True, time=True, potentialEnergy=True,
            temperature=True, volume=True, density=True,
        )
    )
    simulation.reporters.append(
        app.CheckpointReporter(chk_file, save_freq * 10)
    )

    simulation.step(prod_steps)

    elapsed = (time.time() - t_start) / 60.0
    print(f"  [{name}] DONE in {elapsed:.1f} min")

    # Save final frame
    state = simulation.context.getState(getPositions=True)
    with open(os.path.join(out_dir, f"{tag}_final.pdb"), "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Copy to results
    import subprocess as sp
    os.makedirs("/content/results", exist_ok=True)
    sp.run(["cp", "-r", out_dir, "/content/results/"], check=False)


# ═══════════════════════════════════════════════════
#  Run all 3
# ═══════════════════════════════════════════════════

for name, smiles, tfile in TOP3:
    try:
        run_gpu_md(name, smiles, tfile)
    except Exception as e:
        print(f"\n[{name}] FAILED: {str(e)[:500]}")
        import traceback
        traceback.print_exc()

print(f"\n{'='*60}")
print("ALL DONE!")
print(f"Results: /content/results/")
print(f"{'='*60}")
