"""
SMVT MD — Phenobarbital self-contained run.
Downloads AF structure, embeds Vina pose + GAFF template inline.
No file uploads needed. Colab run --gpu T4 this script.
"""
import subprocess, sys
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
    "openmm", "pdbfixer", "rdkit", "openmmforcefields"])
import os, time, urllib.request, glob
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

# ══ Config ══
NAME = "PHENOBARBITAL"
SMILES = "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2"
PROD_NS = 5
TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT = 2.0 * unit.femtoseconds
PADDING = 1.0 * unit.nanometers

# ══ Phenobarbital Vina pose coordinates (19 atoms) ══
VINA_POSE = [
    [-2.400,  0.168,  0.632], [-1.671, -1.104,  0.239], [ 0.319, -2.123, -0.136],
    [-0.363, -3.389, -0.387], [-1.241, -1.129, -0.852], [ 0.770, -0.182,  1.630],
    [ 1.551,  0.855,  1.004], [ 2.302, -0.068, -0.073], [ 1.456, -1.250, -0.519],
    [ 0.521, -0.271, -0.272], [-3.196, -0.001,  1.486], [-3.071,  0.548, -0.289],
    [-2.308,  1.187,  1.083], [ 0.897, -0.287,  2.692], [ 2.489,  1.494,  1.213],
    [ 1.217,  1.600,  0.340], [ 2.419,  0.807,  1.994], [ 0.180,  0.358,  1.869],
    [-2.642,  2.686,  1.567],
]

# ══ GAFF template (Gasteiger, 29 atoms) ══
TEMPLATE_XML = '''<ForceField>
 <Residues>
  <Residue name="LIG">
   <Atom name="C1" type="c3" charge="-0.060"/>
   <Atom name="C2" type="c3" charge="0.179"/>
   <Atom name="C3" type="c" charge="0.462"/>
   <Atom name="N4" type="n" charge="-0.366"/>
   <Atom name="C5" type="c" charge="0.478"/>
   <Atom name="N6" type="n" charge="-0.366"/>
   <Atom name="C7" type="c" charge="0.462"/>
   <Atom name="O8" type="o" charge="-0.464"/>
   <Atom name="O9" type="o" charge="-0.465"/>
   <Atom name="O10" type="o" charge="-0.464"/>
   <Atom name="C11" type="ca" charge="0.039"/>
   <Atom name="C12" type="ca" charge="-0.107"/>
   <Atom name="C13" type="ca" charge="-0.103"/>
   <Atom name="C14" type="ca" charge="-0.107"/>
   <Atom name="C15" type="ca" charge="-0.104"/>
   <Atom name="C16" type="ca" charge="-0.103"/>
   <Atom name="H17" type="hc" charge="0.019"/>
   <Atom name="H18" type="hc" charge="0.019"/>
   <Atom name="H19" type="hc" charge="0.019"/>
   <Atom name="H20" type="hc" charge="0.019"/>
   <Atom name="H21" type="hc" charge="0.019"/>
   <Atom name="H22" type="hn" charge="0.113"/>
   <Atom name="H23" type="hn" charge="0.113"/>
   <Atom name="H24" type="ha" charge="0.133"/>
   <Atom name="H25" type="ha" charge="0.133"/>
   <Atom name="H26" type="ha" charge="0.133"/>
   <Atom name="H27" type="ha" charge="0.133"/>
   <Atom name="H28" type="ha" charge="0.999"/>
   <Atom name="H29" type="hc" charge="0.010"/>
   <Bond from="0" to="1"/><Bond from="0" to="10"/><Bond from="0" to="16"/><Bond from="0" to="17"/>
   <Bond from="1" to="2"/><Bond from="1" to="18"/><Bond from="1" to="19"/>
   <Bond from="2" to="3"/><Bond from="2" to="7"/>
   <Bond from="3" to="4"/><Bond from="3" to="21"/>
   <Bond from="4" to="5"/><Bond from="4" to="8"/>
   <Bond from="5" to="6"/><Bond from="5" to="22"/>
   <Bond from="6" to="1"/><Bond from="6" to="9"/>
   <Bond from="10" to="11"/><Bond from="10" to="28"/>
   <Bond from="11" to="12"/><Bond from="11" to="23"/>
   <Bond from="12" to="13"/><Bond from="12" to="24"/>
   <Bond from="13" to="14"/><Bond from="13" to="25"/>
   <Bond from="14" to="15"/><Bond from="14" to="26"/>
   <Bond from="15" to="10"/><Bond from="15" to="27"/>
  </Residue>
 </Residues>
</ForceField>'''

# ══ GPU ══
plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        break
    except: continue
print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

# ══ Download AF structure ══
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-Q9Y289-F1-model_v4.pdb"
AF_PATH = "/content/AF-Q9Y289-F1.pdb"
if not os.path.exists(AF_PATH):
    print(f"Downloading AF structure...")
    urllib.request.urlretrieve(AF_URL, AF_PATH)

# ══ Write template ══
TEMPLATE_PATH = "/content/template.xml"
with open(TEMPLATE_PATH, "w") as f:
    f.write(TEMPLATE_XML)

# ══ Prepare protein ══
print("Preparing protein...")
fixer = PDBFixer(filename=AF_PATH)
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
prot_path = "/content/protein.pdb"
with open(prot_path, "w") as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

# ══ Ligand from SMILES + Vina pose ══
print("Preparing ligand...")
mol = Chem.MolFromSmiles(SMILES); mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer()
for i in range(min(mol.GetNumAtoms(), len(VINA_POSE))):
    conf.SetAtomPosition(i, tuple(VINA_POSE[i]))
pdb_str = Chem.MolToPDBBlock(mol)
lig_lines = ["REMARK Phenobarbital docking pose"]
for line in pdb_str.split("\n"):
    if line.startswith("HETATM") or line.startswith("ATOM"):
        lig_lines.append(line[:17] + "LIG" + line[20:])
    elif line.startswith("CONECT"):
        lig_lines.append(line)
for bond in mol.GetBonds():
    i = bond.GetBeginAtomIdx() + 1; j = bond.GetEndAtomIdx() + 1
    lig_lines.append(f"CONECT{i:5d}{j:5d}")
lig_path = "/content/ligand.pdb"
with open(lig_path, "w") as f:
    f.write("\n".join(lig_lines))
n_lig = mol.GetNumAtoms()
print(f"  Ligand: {n_lig} atoms")

# ══ ForceField ══
gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
gaff_xml = gaff_candidates[0] if gaff_candidates else \
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
forcefield = app.ForceField("amber14-all.xml", gaff_xml, TEMPLATE_PATH, "tip3p.xml")

# ══ Build system ══
print("Building system...")
prot_pdb = app.PDBFile(prot_path)
lig_pdb = app.PDBFile(lig_path)
modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
modeller.add(lig_pdb.topology, lig_pdb.positions)
modeller.addSolvent(forcefield, model="tip3p", padding=PADDING,
                    ionicStrength=0.0 * unit.molar, neutralize=False)
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds, ignoreExternalBonds=True)
n = system.getNumParticles()
print(f"  {n} particles")

# ══ Backbone restraints + Integrator ══
k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
bb_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
bb_force.addGlobalParameter("k", k_rest)
bb_force.addPerParticleParameter("x0"); bb_force.addPerParticleParameter("y0")
bb_force.addPerParticleParameter("z0")
for atom in modeller.topology.atoms():
    if atom.residue.name not in ("HOH", "WAT", "LIG", "SOL", "NA", "CL"):
        if atom.name in ("CA", "C", "N"):
            p = modeller.positions[atom.index]
            bb_force.addParticle(atom.index, [p.x, p.y, p.z])
bb_idx = system.addForce(bb_force)

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

# ══ 1. Minimize ══
print("\n1. Minimize (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)
if not check("minimize"): exit(1)

# ══ 2. NVT ══
print("2. NVT heating 50→300K...")
simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
for t in [50, 100, 150, 200, 250, 300]:
    simulation.integrator.setTemperature(t * unit.kelvin)
    simulation.step(10_000)
    if not check(f"NVT {t}K"): exit(1)

simulation.integrator.setTemperature(TEMP)
simulation.step(50_000)
if not check("NVT eq"): exit(1)
print("  ✓ NVT PASS!")

# ══ 3. NPT ══
print("3. NPT equil...")
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

# ══ 4. Production ══
print(f"4. Production {PROD_NS}ns...")
prod_steps = int(PROD_NS * 1_000_000 / 2.0)
save_freq = 25_000
out_dir = f"/content/trajectories/{NAME}"
os.makedirs(out_dir, exist_ok=True)
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
print(f"\n✓✓✓ {NAME}: ALL PASS in {elapsed:.1f} min!")
print(f"  {out_dir}/")
