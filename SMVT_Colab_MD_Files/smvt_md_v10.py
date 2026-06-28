import os, time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try: plat = mm.Platform.getPlatformByName(pname); break
    except: continue
print("OpenMM", mm.__version__)
print("Platform:", plat.getName())

DATA_DIR = "/content"
OUT_DIR = "/content/trajectories"
os.makedirs(OUT_DIR, exist_ok=True)
RECEPTOR_PDB = DATA_DIR + "/AF-Q9Y289-F1.pdb"

TOP3 = [
    ("NAFTAZONE", "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N", "NAFTAZONE_template.xml"),
    ("PHENOBARBITAL", "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2", "PHENOBARBITAL_template.xml"),
    ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O", "ESKETAMINE_template.xml"),
]

PROD_NS = 50
TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT = 2.0 * unit.femtoseconds

print("GPU MD:", len(TOP3), "compounds x", PROD_NS, "ns each (GAFF2 + restraints)")

def extract_vina_pose(pdbqt_path):
    pos = []; in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if int(line.split()[1]) == 1: in_model = True; pos = []
                elif in_model: break
            elif line.startswith("ENDMDL") and in_model: break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try: pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: continue
    return pos

def prepare_protein(pdb_path, out_path):
    if os.path.exists(out_path): return out_path
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return out_path

def prepare_ligand_pdb(smiles, vina_pos, out_path, res_name="LIG"):
    if os.path.exists(out_path): return out_path
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for i in range(min(mol.GetNumAtoms(), len(vina_pos))):
        conf.SetAtomPosition(i, (vina_pos[i][0], vina_pos[i][1], vina_pos[i][2]))
    pdb_str = Chem.MolToPDBBlock(mol)
    fixed_lines = []
    for line in pdb_str.split("\n"):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            line = line[:17] + res_name + line[20:]
        fixed_lines.append(line)
    with open(out_path, "w") as f:
        f.write("\n".join(fixed_lines))
    return out_path

def run_gpu_md(name, smiles, template_file):
    tag = name; out_dir = OUT_DIR + "/" + tag; os.makedirs(out_dir, exist_ok=True)
    chk_file = out_dir + "/" + tag + "_" + str(PROD_NS) + "ns.chk"
    if os.path.exists(chk_file):
        print("[" + name + "] Done. Skip."); return

    sep = "=" * 50
    print("\n" + sep + "\n[" + name + "] " + str(PROD_NS) + "ns MD (GAFF2)\n" + sep)
    t_start = time.time()

    prot_clean = out_dir + "/protein.pdb"
    prepare_protein(RECEPTOR_PDB, prot_clean)

    pdbqt = DATA_DIR + "/" + name + "_docked.pdbqt"
    vina_pos = extract_vina_pose(pdbqt)
    print("  Vina pose:", len(vina_pos), "atoms")
    lig_pdb = out_dir + "/ligand.pdb"
    prepare_ligand_pdb(smiles, vina_pos, lig_pdb, "LIG")

    template_xml = DATA_DIR + "/" + template_file
    gaff_xml = "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
    forcefield = app.ForceField("amber14-all.xml", gaff_xml, template_xml, "tip3p.xml")

    print("  Building system...")
    prot_pdb = app.PDBFile(prot_clean)
    lig = app.PDBFile(lig_pdb)
    modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
    modeller.add(lig.topology, lig.positions)
    modeller.addSolvent(forcefield, model="tip3p",
                        padding=0.8*unit.nanometers, ionicStrength=0.0*unit.molar,
                        neutralize=False)
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometers,
                                     constraints=app.HBonds,
                                     ignoreExternalBonds=True)
    print("  System:", system.getNumParticles(), "particles")

    # Backbone position restraints during equilibration
    k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms * unit.angstroms)
    restraint_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint_force.addGlobalParameter("k", k_rest)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    backbone_atoms = set(["CA", "C", "N"])
    for atom in modeller.topology.atoms():
        if atom.residue.name not in ["HOH", "WAT", "LIG", "SOL", "NA", "CL"]:
            if atom.name in backbone_atoms:
                p = modeller.positions[atom.index]
                restraint_force.addParticle(atom.index, [p.x, p.y, p.z])
    restraint_idx = system.addForce(restraint_force)

    integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0/unit.picosecond, DT)
    simulation = app.Simulation(modeller.topology, system, integrator, plat,
                                {"Precision": "mixed"})
    simulation.context.setPositions(modeller.positions)

    # Phase 1: Minimize with restraints
    print("  Minimizing (restrained)...")
    simulation.minimizeEnergy(maxIterations=5000)

    # Phase 2: NVT heating with restraints (1fs timestep)
    print("  NVT heating (50->300K, 1fs, restrained)...")
    simulation.integrator.setStepSize(1.0*unit.femtoseconds)
    for t in [50, 100, 150, 200, 250, 300]:
        simulation.integrator.setTemperature(t*unit.kelvin)
        simulation.step(int(10*1000/1.0))  # 10ps per stage
    simulation.integrator.setTemperature(TEMP)
    simulation.step(int(50*1000/1.0))  # 50ps more at 300K

    # Phase 3: NPT with decreasing restraints
    print("  NPT equil (reducing restraints)...")
    simulation.integrator.setStepSize(DT)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
    simulation.context.reinitialize(preserveState=True)

    for k_val in [2.5, 0.5]:
        new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms * unit.angstroms)
        simulation.context.setParameter("k", new_k)
        simulation.step(int(50*1000/2.0))

    # Phase 4: Remove restraints, short unrestrained equil
    print("  Removing restraints, final equil...")
    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(int(100*1000/2.0))

    # Phase 5: Production (unrestrained)
    prod_steps = int(PROD_NS * 1_000_000 / 2.0)
    save_freq = int(100 * 1000 / 2.0)
    print("  Production:", PROD_NS, "ns (", prod_steps, "steps)...")
    simulation.reporters.append(app.DCDReporter(out_dir + "/" + tag + "_" + str(PROD_NS) + "ns.dcd", save_freq))
    simulation.reporters.append(app.StateDataReporter(
        out_dir + "/" + tag + "_" + str(PROD_NS) + "ns.csv", save_freq,
        step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True
    ))
    simulation.reporters.append(app.CheckpointReporter(chk_file, save_freq*10))

    simulation.step(prod_steps)
    elapsed = (time.time()-t_start)/60
    print("  [" + name + "] DONE in", round(elapsed, 1), "min")

    state = simulation.context.getState(getPositions=True)
    with open(out_dir + "/" + tag + "_final.pdb", "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    simulation.saveCheckpoint(chk_file)

    import subprocess as sp
    os.makedirs("/content/results", exist_ok=True)
    sp.run(["cp", "-r", out_dir, "/content/results/"], check=False)

for name, smiles, tfile in TOP3:
    try:
        run_gpu_md(name, smiles, tfile)
    except Exception as e:
        print("[" + name + "] FAILED:", str(e)[:400])
        import traceback; traceback.print_exc()

print("\nALL DONE!")
