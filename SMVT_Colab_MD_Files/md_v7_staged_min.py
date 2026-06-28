"""
SMVT MD v7 — Staged minimization to resolve Vina pose clashes.
3 severe steric clashes found in Vina pose (shortest: 2.08Å, 61% vdW).
Protocol: freeze ligand → relax protein → freeze protein → relax ligand
          → unrestrained min → slow NVT → NPT → production
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

out_dir = f"/content/trajectories/{NAME}_v7"
os.makedirs(out_dir, exist_ok=True)

t0 = time.time()

# Helper
def check(label, simulation):
    pe = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    v = pe.value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(v) or np.isinf(v))
    print(f"  {label}: PE={v:.0f} {'✓' if ok else '❌'}")
    return ok, v


# ══ 1. Protein ══
print("1. Protein (PDBFixer)...")
fixer = PDBFixer(filename=f"{DATA_DIR}/AF-Q9Y289-F1.pdb")
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)

# ══ 2. Vina pose ══
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

# ══ 3. Ligand PDB ══
mol = Chem.MolFromSmiles(SMILES); mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer()
for i in range(min(mol.GetNumAtoms(), len(pos))):
    conf.SetAtomPosition(i, tuple(pos[i]))
pdb_str = Chem.MolToPDBBlock(mol)
lig_lines = []
bonds_written = False
for line in pdb_str.split("\n"):
    if line.startswith("HETATM") or line.startswith("ATOM"):
        lig_lines.append(line[:17] + "LIG" + line[20:])
    elif line.startswith("CONECT"):
        lig_lines.append(line)
        bonds_written = True
# If RDKit didn't write CONECT, add them
if not bonds_written:
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx() + 1  # PDB is 1-indexed
        j = bond.GetEndAtomIdx() + 1
        lig_lines.append(f"CONECT{i:5d}{j:5d}")
with open(f"{out_dir}/ligand.pdb", "w") as f:
    f.write("\n".join(lig_lines))

# ══ 4. ForceField ══
gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
gaff_xml = gaff_candidates[0] if gaff_candidates else \
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
forcefield = app.ForceField("amber14-all.xml", gaff_xml,
                            f"{DATA_DIR}/{TEMPLATE}", "tip3p.xml")

# ══ 5. Build system ══
print("5. Building system...")
with open(f"{out_dir}/protein.pdb", "w") as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
prot_pdb = app.PDBFile(f"{out_dir}/protein.pdb")
lig = app.PDBFile(f"{out_dir}/ligand.pdb")
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

# ══ 6. Identify atom groups ══
protein_atoms = []
ligand_atoms = []
for atom in modeller.topology.atoms():
    if atom.residue.name in ("HOH", "WAT", "SOL", "NA", "CL"):
        continue
    if atom.residue.name == "LIG":
        ligand_atoms.append(atom.index)
    else:
        protein_atoms.append(atom.index)
print(f"  Protein: {len(protein_atoms)} | Ligand: {len(ligand_atoms)} | "
      f"Water: {n - len(protein_atoms) - len(ligand_atoms)}")

# ══ 7. Position restraint helper ══
def add_pos_restraints(simulation, atom_indices, k_kcal):
    """Add CustomExternalForce to restrain atoms to current positions."""
    k = k_kcal * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k)
    force.addPerParticleParameter("x0"); force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    for idx in atom_indices:
        p = positions[idx]
        force.addParticle(idx, [p.x, p.y, p.z])
    idx = simulation.system.addForce(force)
    simulation.context.reinitialize(preserveState=True)
    return idx


# ══ 8. Simulation setup ══
integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0 / unit.picosecond, DT)
simulation = app.Simulation(modeller.topology, system, integrator, plat,
                            {"Precision": "mixed"})
simulation.context.setPositions(modeller.positions)


# ══ 9. STAGED MINIMIZATION ══
print(f"\n{'='*55}")
print("STAGED MINIMIZATION")
print(f"{'='*55}")

# 9a. Restrain ALL atoms moderately, minimize (relax bad contacts gently)
print("\n[9a] Gentle minimize (all restrained @ 10 kcal/mol/Å²)...")
all_atoms = protein_atoms + ligand_atoms
r_all = add_pos_restraints(simulation, all_atoms, 10.0)
simulation.minimizeEnergy(maxIterations=5000)
ok, pe = check("[9a] gentle min", simulation)
if not ok: print("FAILED"); exit(1)

# 9b. Freeze protein, restrain ligand weakly, minimize (let protein adjust to ligand)
print("\n[9b] Minimize: protein free, ligand restrained @ 2.5...")
simulation.system.removeForce(r_all)
simulation.context.reinitialize(preserveState=True)
r_lig = add_pos_restraints(simulation, ligand_atoms, 2.5)
simulation.minimizeEnergy(maxIterations=5000)
ok, pe = check("[9b] free protein", simulation)
if not ok: print("FAILED"); exit(1)

# 9c. Freeze ligand, restrain protein backbone, minimize
print("\n[9c] Minimize: ligand free, protein restrained @ 5...")
simulation.system.removeForce(r_lig)
simulation.context.reinitialize(preserveState=True)
backbone = {"CA", "C", "N"}
bb_atoms = []
for atom in modeller.topology.atoms():
    if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
        if atom.name in backbone:
            bb_atoms.append(atom.index)
r_bb = add_pos_restraints(simulation, bb_atoms, 5.0)
simulation.minimizeEnergy(maxIterations=5000)
ok, pe = check("[9c] free ligand", simulation)
if not ok: print("FAILED"); exit(1)

# 9d. All free, gentle minimize
print("\n[9d] Minimize: ALL free, 5000 steps...")
simulation.system.removeForce(r_bb)
simulation.context.reinitialize(preserveState=True)
simulation.minimizeEnergy(maxIterations=5000)
ok, pe = check("[9d] free min", simulation)
if not ok: print("FAILED"); exit(1)
print(f"  ✓ Staged minimization complete!")


# ══ 10. NVT HEATING ══
print(f"\n{'='*55}")
print("NVT HEATING (50→300K)")
print(f"{'='*55}")

# Backbone restraints for NVT — use post-minimization positions!
state = simulation.context.getState(getPositions=True)
min_positions = state.getPositions()
bb_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
bb_force.addGlobalParameter("k", k_rest)
bb_force.addPerParticleParameter("x0")
bb_force.addPerParticleParameter("y0")
bb_force.addPerParticleParameter("z0")
for atom in modeller.topology.atoms():
    if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
        if atom.name in backbone:
            p = min_positions[atom.index]
            bb_force.addParticle(atom.index, [p.x, p.y, p.z])
bb_idx = system.addForce(bb_force)
simulation.context.reinitialize(preserveState=True)

simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
for t in [50, 100, 150, 200, 250, 300]:
    simulation.integrator.setTemperature(t * unit.kelvin)
    simulation.step(10_000)
    ok, pe = check(f"NVT {t}K", simulation)
    if not ok:
        print(f"\n❌ NaN at {t}K! Staged minimization didn't fix it.")
        exit(1)

simulation.integrator.setTemperature(TEMP)
simulation.step(50_000)
print("  ✓ NVT heating: PASS!")

# ══ 11. NPT ══
print(f"\n{'='*55}")
print("NPT EQUILIBRATION")
print(f"{'='*55}")
simulation.integrator.setStepSize(DT)
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
simulation.context.reinitialize(preserveState=True)

for k_val in [2.5, 0.5]:
    new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    simulation.context.setParameter("k", new_k)
    simulation.step(25_000)
    ok, pe = check(f"NPT k={k_val}", simulation)
    if not ok: print("FAILED"); exit(1)

system.removeForce(bb_idx)
simulation.context.reinitialize(preserveState=True)
simulation.step(50_000)
ok, pe = check("NPT free", simulation)
if not ok: print("FAILED"); exit(1)
print("  ✓ NPT equil: PASS!")

# ══ 12. PRODUCTION ══
print(f"\n{'='*55}")
print(f"PRODUCTION {PROD_NS}ns")
print(f"{'='*55}")
prod_steps = int(PROD_NS * 1_000_000 / 2.0)
save_freq = 25_000
simulation.reporters.append(app.DCDReporter(f"{out_dir}/{NAME}_{PROD_NS}ns.dcd", save_freq))
simulation.reporters.append(app.StateDataReporter(
    f"{out_dir}/{NAME}_{PROD_NS}ns.csv", save_freq,
    step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, density=True))

simulation.step(prod_steps)
ok, pe = check("Production end", simulation)
if not ok: print("FAILED"); exit(1)

elapsed = (time.time() - t0) / 60.0
state = simulation.context.getState(getPositions=True)
with open(f"{out_dir}/{NAME}_final.pdb", "w") as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

print(f"\n{'='*55}")
print(f"✓✓✓ {NAME}: ALL PASS in {elapsed:.1f} min!")
print(f"Results: {out_dir}/")
print(f"{'='*55}")
