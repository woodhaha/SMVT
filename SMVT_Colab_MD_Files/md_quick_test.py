"""
SMVT MD Quick Test — Verify NVT heating doesn't NaN with v3 GAFF templates.

Focus: Esketamine (has Cl, was problematic) × 2ns production only.
If NVT passes, proceed to full 50ns runs.
"""

import os, sys, time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

# ══ GPU ══
plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        break
    except:
        continue
print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

# ══ Config ══
DATA_DIR = "/content"
os.makedirs("/content/trajectories", exist_ok=True)

COMPOUNDS = [
    # (name, smiles, template, prod_ns)
    ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O",
     "ESKETAMINE_template.xml", 5),  # 5ns test
]

TEMP     = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT       = 2.0 * unit.femtoseconds
PADDING  = 1.0 * unit.nanometers


def extract_vina_pose(pdbqt_path):
    pos = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL") and int(line.split()[1]) == 1:
                in_model = True
                pos = []
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except:
                    continue
    return pos


def prepare_protein(pdb_path, out_path):
    if os.path.exists(out_path):
        return out_path
    print("  PDBFixer: missing atoms + H @ pH 7.4...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return out_path


def prepare_ligand_pdb(smiles, vina_pos, out_path, res_name="LIG"):
    if os.path.exists(out_path):
        return out_path
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for i in range(min(mol.GetNumAtoms(), len(vina_pos))):
        conf.SetAtomPosition(i, tuple(vina_pos[i]))
    pdb_str = Chem.MolToPDBBlock(mol)
    fixed = []
    for line in pdb_str.split("\n"):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            line = line[:17] + res_name + line[20:]
        fixed.append(line)
    with open(out_path, "w") as f:
        f.write("\n".join(fixed))
    return out_path


def run_md(name, smiles, template_file, prod_ns):
    tag = name
    out_dir = f"/content/trajectories/{tag}"
    os.makedirs(out_dir, exist_ok=True)

    chk = f"{out_dir}/{tag}_{prod_ns}ns.chk"
    if os.path.exists(chk):
        print(f"\n[{name}] Already done. Skip.")
        return True

    print(f"\n{'='*55}")
    print(f"[{name}] {prod_ns}ns MD test — v3 GAFF template")
    print(f"{'='*55}")
    t0 = time.time()

    # 1. Protein
    prot_clean = prepare_protein(f"{DATA_DIR}/AF-Q9Y289-F1.pdb",
                                 f"{out_dir}/protein.pdb")

    # 2. Extract Vina pose
    vina_pos = extract_vina_pose(f"{DATA_DIR}/{name}_docked.pdbqt")
    print(f"  Vina atoms: {len(vina_pos)}")

    # 3. Ligand PDB
    lig_pdb = prepare_ligand_pdb(smiles, vina_pos, f"{out_dir}/ligand.pdb")

    # 4. ForceField: AMBER ff14SB + GAFF 2.11 + v3 template + TIP3P
    template_xml = f"{DATA_DIR}/{template_file}"

    # Find gaff-2.11.xml
    import glob
    gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
    gaff_xml = gaff_candidates[0] if gaff_candidates else "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"

    print(f"  ForceField: amber14 + gaff-2.11 + {template_file} + tip3p")
    forcefield = app.ForceField("amber14-all.xml", gaff_xml, template_xml, "tip3p.xml")

    # 5. Build system
    print("  Building system...")
    prot_pdb = app.PDBFile(prot_clean)
    lig = app.PDBFile(lig_pdb)
    modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
    modeller.add(lig.topology, lig.positions)
    modeller.addSolvent(forcefield, model="tip3p", padding=PADDING,
                        ionicStrength=0.0 * unit.molar, neutralize=False)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds, ignoreExternalBonds=True,
    )
    n_particles = system.getNumParticles()
    print(f"  System: {n_particles} particles")

    # 6. Backbone restraints
    k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    restraint_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint_force.addGlobalParameter("k", k_rest)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    backbone = {"CA", "C", "N"}
    n_restrained = 0
    for atom in modeller.topology.atoms():
        if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
            if atom.name in backbone:
                p = modeller.positions[atom.index]
                restraint_force.addParticle(atom.index, [p.x, p.y, p.z])
                n_restrained += 1
    restraint_idx = system.addForce(restraint_force)
    print(f"  Restrained {n_restrained} backbone atoms")

    # 7. Integrator + Simulation
    integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0 / unit.picosecond, DT)
    simulation = app.Simulation(modeller.topology, system, integrator, plat,
                                {"Precision": "mixed"})
    simulation.context.setPositions(modeller.positions)

    # ── Helper: check for NaN ──
    def check_nan(label):
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        if np.isnan(pe) or np.isinf(pe):
            print(f"  ❌ NaN/Inf detected at {label}! PE={pe}")
            return False
        # print(f"  ✓ {label}: PE={pe:.1f} kcal/mol")
        return True

    # ── Phase 1: Minimize ──
    print("  Phase 1: Minimize (5000 steps)...")
    try:
        simulation.minimizeEnergy(maxIterations=5000)
    except Exception as e:
        print(f"  ❌ Minimize failed: {e}")
        return False
    if not check_nan("minimize"):
        return False
    state = simulation.context.getState(getEnergy=True)
    print(f"  ✓ Minimize: PE={state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole):.0f} kcal/mol")

    # ── Phase 2: NVT heating (1fs) ──
    print("  Phase 2: NVT heating 50→300K...")
    simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
    for t in [50, 100, 150, 200, 250, 300]:
        simulation.integrator.setTemperature(t * unit.kelvin)
        simulation.step(10_000)
        if not check_nan(f"NVT {t}K"):
            return False
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"    {t}K: PE={pe:.0f} kcal/mol")

    simulation.integrator.setTemperature(TEMP)
    simulation.step(50_000)
    if not check_nan("NVT 300K eq"):
        return False
    print("  ✓ NVT heating: PASS (no NaN)")

    # ── Phase 3: NPT equil ──
    print("  Phase 3: NPT equil...")
    simulation.integrator.setStepSize(DT)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
    simulation.context.reinitialize(preserveState=True)

    for k_val in [2.5, 0.5]:
        new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
        simulation.context.setParameter("k", new_k)
        simulation.step(25_000)
        if not check_nan(f"NPT k={k_val}"):
            return False

    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(50_000)
    if not check_nan("NPT unrestrained"):
        return False
    print("  ✓ NPT equil: PASS")

    # ── Phase 4: Production ──
    prod_steps = int(prod_ns * 1_000_000 / 2.0)
    save_freq  = 25_000  # every 50ps
    print(f"  Phase 4: Production {prod_ns}ns ({prod_steps} steps)...")

    simulation.reporters.append(app.DCDReporter(f"{out_dir}/{tag}_{prod_ns}ns.dcd", save_freq))
    simulation.reporters.append(app.StateDataReporter(
        f"{out_dir}/{tag}_{prod_ns}ns.csv", save_freq,
        step=True, time=True, potentialEnergy=True,
        temperature=True, volume=True, density=True,
    ))
    simulation.reporters.append(app.CheckpointReporter(chk, save_freq * 10))

    simulation.step(prod_steps)

    if not check_nan("production end"):
        return False

    elapsed = (time.time() - t0) / 60.0
    # Final frame
    state = simulation.context.getState(getPositions=True)
    with open(f"{out_dir}/{tag}_final.pdb", "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    print(f"  ✓ [{name}] DONE in {elapsed:.1f} min")
    return True


# ══ Run ══
print(f"\n{'='*55}")
print("SMVT MD Quick Test — v3 GAFF Templates")
print(f"{'='*55}")

results = {}
for name, smiles, tfile, ns in COMPOUNDS:
    ok = run_md(name, smiles, tfile, ns)
    results[name] = "PASS" if ok else "FAIL"

print(f"\n{'='*55}")
print("RESULTS:")
for k, v in results.items():
    print(f"  {k}: {v}")
print(f"{'='*55}")
