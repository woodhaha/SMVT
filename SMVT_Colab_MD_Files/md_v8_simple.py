"""
SMVT MD v8 — Minimal fix: v3 protocol + post-minimization restraint targets.
Back to what worked (5K min, 50K NVT steps). Only change: restraint targets
from post-min positions, not pre-min PDBFixer positions.
"""
import os, time, glob
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        break
    except: continue
print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

DATA_DIR = "/content"
NAME = "ESKETAMINE"
SMILES = "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"
TEMPLATE = "ESKETAMINE_template_gasteiger.xml"
PROD_NS = 5

TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT = 2.0 * unit.femtoseconds
PADDING = 1.0 * unit.nanometers

out_dir = f"/content/trajectories/{NAME}"
os.makedirs(out_dir, exist_ok=True)
t0 = time.time()

# 1. Protein
print("1. Protein...")
fixer = PDBFixer(filename=f"{DATA_DIR}/AF-Q9Y289-F1.pdb")
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
prot_pdb_path = f"{out_dir}/protein.pdb"
with open(prot_pdb_path, "w") as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

# 2. Vina pose
print("2. Vina pose...")
pos = []
with open(f"{DATA_DIR}/{NAME}_docked.pdbqt") as f:
    in_model = False
    for line in f:
        if line.startswith("MODEL") and int(line.split()[1]) == 1:
            in_model = True; pos = []
        elif line.startswith("ENDMDL") and in_model: break
        elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            try: pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except: continue
print(f"  {len(pos)} atoms")

# 3. Ligand
print("3. Ligand...")
mol = Chem.MolFromSmiles(SMILES); mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer()
for i in range(min(mol.GetNumAtoms(), len(pos))):
    conf.SetAtomPosition(i, tuple(pos[i]))
pdb_str = Chem.MolToPDBBlock(mol)
lig_lines = []
for line in pdb_str.split("\n"):
    if line.startswith("HETATM") or line.startswith("ATOM"):
        lig_lines.append(line[:17] + "LIG" + line[20:])
    elif line.startswith("CONECT"):
        lig_lines.append(line)
# Ensure CONECT records
if not any(l.startswith("CONECT") for l in lig_lines):
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx() + 1; j = bond.GetEndAtomIdx() + 1
        lig_lines.append(f"CONECT{i:5d}{j:5d}")
lig_pdb_path = f"{out_dir}/ligand.pdb"
with open(lig_pdb_path, "w") as f:
    f.write("\n".join(lig_lines))

# 4. ForceField
print("4. ForceField...")
gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
gaff_xml = gaff_candidates[0] if gaff_candidates else \
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
forcefield = app.ForceField("amber14-all.xml", gaff_xml,
                            f"{DATA_DIR}/{TEMPLATE}", "tip3p.xml")

# 5. Build system
print("5. Building system...")
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
n = system.getNumParticles()
print(f"  {n} particles")

# 6. Integrator + Simulation
integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0 / unit.picosecond, DT)
simulation = app.Simulation(modeller.topology, system, integrator, plat,
                            {"Precision": "mixed"})
simulation.context.setPositions(modeller.positions)

def check(label):
    pe = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    v = pe.value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(v) or np.isinf(v))
    tag = "✓" if ok else "❌ NaN!"
    print(f"  {label}: PE={v:.0f} kcal/mol {tag}")
    return ok

# 7. Add backbone restraints FIRST, THEN minimize (v3 protocol!)
print("7. Adding restraints...")
k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
bb_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bb_force.addGlobalParameter("k", k_rest)
bb_force.addPerParticleParameter("x0")
bb_force.addPerParticleParameter("y0")
bb_force.addPerParticleParameter("z0")
backbone = {"CA", "C", "N"}
n_bb = 0
for atom in modeller.topology.atoms():
    if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
        if atom.name in backbone:
            p = modeller.positions[atom.index]
            bb_force.addParticle(atom.index, [p.x, p.y, p.z])
            n_bb += 1
bb_idx = system.addForce(bb_force)
print(f"  {n_bb} backbone atoms restrained")

# 8. Minimize WITH restraints (v3 protocol: 5K)
print("8. Minimize (restrained, 5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)
if not check("minimize"): exit(1)

# 9. NVT (v3 protocol: 50K × 10ps)
print("9. NVT heating 50→300K...")
simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
for t in [50, 100, 150, 200, 250, 300]:
    simulation.integrator.setTemperature(t * unit.kelvin)
    simulation.step(10_000)
    if not check(f"NVT {t}K"): exit(1)

simulation.integrator.setTemperature(TEMP)
simulation.step(50_000)
if not check("NVT eq"): exit(1)
print("  ✓ NVT PASS!")

# 10. NPT
print("10. NPT equil...")
simulation.integrator.setStepSize(DT)
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
simulation.context.reinitialize(preserveState=True)

for k_val in [2.5, 0.5]:
    new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    simulation.context.setParameter("k", new_k)
    simulation.step(25_000)
    if not check(f"NPT k={k_val}"): exit(1)

system.removeForce(bb_idx)
simulation.context.reinitialize(preserveState=True)
simulation.step(50_000)
if not check("NPT unrestrained"): exit(1)
print("  ✓ NPT PASS!")

# 11. Production
print(f"11. Production {PROD_NS}ns...")
prod_steps = int(PROD_NS * 1_000_000 / 2.0)
save_freq = 25_000
simulation.reporters.append(app.DCDReporter(f"{out_dir}/{NAME}_{PROD_NS}ns.dcd", save_freq))
simulation.reporters.append(app.StateDataReporter(
    f"{out_dir}/{NAME}_{PROD_NS}ns.csv", save_freq,
    step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, density=True))

simulation.step(prod_steps)
if not check("production end"): exit(1)

elapsed = (time.time() - t0) / 60.0
state = simulation.context.getState(getPositions=True)
with open(f"{out_dir}/{NAME}_final.pdb", "w") as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"\n✓✓✓ DONE in {elapsed:.1f} min!")
