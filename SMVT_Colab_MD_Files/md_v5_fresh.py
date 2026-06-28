"""
SMVT MD v5 — Fresh start. Gasteiger charges + v3 protocol (which got to 150K).
Deletes cached files, rebuilds from scratch.
"""
import os, time, shutil, glob
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
    except: continue
print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

DATA_DIR = "/content"
OUT_ROOT = "/content/trajectories"

# Clean cached files from previous runs
for d in os.listdir(OUT_ROOT) if os.path.exists(OUT_ROOT) else []:
    path = os.path.join(OUT_ROOT, d)
    if os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)
os.makedirs(OUT_ROOT, exist_ok=True)

# Use Gasteiger template (conservative charges, correct types)
TEMPLATE = "ESKETAMINE_template_gasteiger.xml"
SMILES = "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"
NAME = "ESKETAMINE"
PROD_NS = 5

TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT = 2.0 * unit.femtoseconds
PADDING = 1.0 * unit.nanometers


def extract_vina_pose(pdbqt_path):
    pos = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL") and int(line.split()[1]) == 1:
                in_model = True; pos = []
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: continue
    return pos


# ══ Main ══
print(f"\n{'='*55}")
print(f"SMVT MD v5 — Fresh + Gasteiger + v3 protocol")
print(f"{'='*55}")

t0 = time.time()
out_dir = os.path.join(OUT_ROOT, NAME)
os.makedirs(out_dir, exist_ok=True)

# 1. Protein (fresh)
print("Preparing protein (PDBFixer)...")
fixer = PDBFixer(filename=f"{DATA_DIR}/AF-Q9Y289-F1.pdb")
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
prot_pdb_path = f"{out_dir}/protein.pdb"
with open(prot_pdb_path, "w") as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

# 2. Extract Vina pose
vina_pos = extract_vina_pose(f"{DATA_DIR}/{NAME}_docked.pdbqt")
print(f"Vina pose: {len(vina_pos)} atoms")

# 3. Ligand PDB (fresh)
print("Preparing ligand...")
mol = Chem.MolFromSmiles(SMILES)
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
        line = line[:17] + "LIG" + line[20:]
    fixed.append(line)
lig_pdb_path = f"{out_dir}/ligand.pdb"
with open(lig_pdb_path, "w") as f:
    f.write("\n".join(fixed))

# 4. ForceField (Gasteiger template)
gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
gaff_xml = gaff_candidates[0] if gaff_candidates else \
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
template_xml = f"{DATA_DIR}/{TEMPLATE}"
print(f"ForceField: amber14 + gaff-2.11 + {TEMPLATE} + tip3p")
forcefield = app.ForceField("amber14-all.xml", gaff_xml, template_xml, "tip3p.xml")

# 5. Build system
print("Building system...")
prot_pdb = app.PDBFile(prot_pdb_path)
lig = app.PDBFile(lig_pdb_path)
modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
modeller.add(lig.topology, lig.positions)
modeller.addSolvent(forcefield, model="tip3p", padding=PADDING,
                    ionicStrength=0.0 * unit.molar, neutralize=False)
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds, ignoreExternalBonds=True)
n_particles = system.getNumParticles()
print(f"System: {n_particles} particles")

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
print(f"Restrained {n_restrained} backbone atoms")

# 7. Integrator + Simulation
integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0 / unit.picosecond, DT)
simulation = app.Simulation(modeller.topology, system, integrator, plat,
                            {"Precision": "mixed"})
simulation.context.setPositions(modeller.positions)

def check_nan(label):
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    if np.isnan(pe) or np.isinf(pe):
        print(f"  ❌ NaN/Inf at {label}!")
        return False
    return True

# Phase 1: Minimize (5K steps — v3 protocol)
print("Phase 1: Minimize (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)
if not check_nan("minimize"):
    print("FAILED"); exit(1)
state = simulation.context.getState(getEnergy=True)
pe = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
print(f"  ✓ Minimize: PE={pe:.0f} kcal/mol")

# Phase 2: NVT heating (v3 protocol: 50K × 10ps)
print("Phase 2: NVT heating 50→300K...")
simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
for t in [50, 100, 150, 200, 250, 300]:
    simulation.integrator.setTemperature(t * unit.kelvin)
    simulation.step(10_000)
    if not check_nan(f"NVT {t}K"):
        print(f"  ❌ NaN at {t}K!"); exit(1)
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"  {t}K: PE={pe:.0f} kcal/mol")

simulation.integrator.setTemperature(TEMP)
simulation.step(50_000)
if not check_nan("NVT 300K eq"):
    print("FAILED"); exit(1)
print("  ✓ NVT heating: PASS (no NaN!)")

# Phase 3: NPT equil
print("Phase 3: NPT equil...")
simulation.integrator.setStepSize(DT)
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
simulation.context.reinitialize(preserveState=True)

for k_val in [2.5, 0.5]:
    new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    simulation.context.setParameter("k", new_k)
    simulation.step(25_000)
    if not check_nan(f"NPT k={k_val}"):
        print("FAILED"); exit(1)

system.removeForce(restraint_idx)
simulation.context.reinitialize(preserveState=True)
simulation.step(50_000)
if not check_nan("NPT unrestrained"):
    print("FAILED"); exit(1)
print("  ✓ NPT equil: PASS")

# Phase 4: Production
prod_steps = int(PROD_NS * 1_000_000 / 2.0)
save_freq = 25_000
print(f"Phase 4: Production {PROD_NS}ns ({prod_steps} steps)...")

simulation.reporters.append(app.DCDReporter(f"{out_dir}/{NAME}_{PROD_NS}ns.dcd", save_freq))
simulation.reporters.append(app.StateDataReporter(
    f"{out_dir}/{NAME}_{PROD_NS}ns.csv", save_freq,
    step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, density=True))

simulation.step(prod_steps)

if not check_nan("production end"):
    print("FAILED"); exit(1)

elapsed = (time.time() - t0) / 60.0
state = simulation.context.getState(getPositions=True)
with open(f"{out_dir}/{NAME}_final.pdb", "w") as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

print(f"\n{'='*55}")
print(f"✓ {NAME}: COMPLETE in {elapsed:.1f} min!")
print(f"  NVT: PASS | NPT: PASS | Production: {PROD_NS}ns")
print(f"  Results: {out_dir}/")
print(f"{'='*55}")
