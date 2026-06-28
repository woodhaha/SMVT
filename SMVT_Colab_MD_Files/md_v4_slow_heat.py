"""
SMVT MD v4 — Slower heating + charge fallback.
Key changes from v3:
- 25K temp increments (not 50K) for gentler heating
- 20ps per stage (not 10ps)
- Longer minimize (10K steps)
- Try MMFF94 first, fall back to Gasteiger if NaN
"""
import os, time
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

DATA_DIR = "/content"
os.makedirs("/content/trajectories", exist_ok=True)

COMPOUNDS = [
    ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O",
     "ESKETAMINE_template.xml",            # MMFF94
     "ESKETAMINE_template_gasteiger.xml",   # Gasteiger fallback
     5),  # prod_ns
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
                in_model = True; pos = []
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: continue
    return pos


def prepare_protein(pdb_path, out_path):
    if os.path.exists(out_path):
        return out_path
    print("  PDBFixer...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
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


def run_md_v4(name, smiles, mmff_template, gasteiger_template, prod_ns=5):
    """Try MMFF94 template first. Fall back to Gasteiger if NaN."""
    tag = name
    out_dir = f"/content/trajectories/{tag}"
    os.makedirs(out_dir, exist_ok=True)

    # Prep protein + ligand (reuse across attempts)
    prot_clean = prepare_protein(f"{DATA_DIR}/AF-Q9Y289-F1.pdb", f"{out_dir}/protein.pdb")
    vina_pos = extract_vina_pose(f"{DATA_DIR}/{name}_docked.pdbqt")
    lig_pdb = prepare_ligand_pdb(smiles, vina_pos, f"{out_dir}/ligand.pdb")

    # Try MMFF94 first, then Gasteiger
    for charge_model, template_file in [
        ("MMFF94", mmff_template),
        ("Gasteiger", gasteiger_template),
    ]:
        print(f"\n  ── Attempt: {charge_model} charges ──")
        ok = _run_once(name, tag, out_dir, prot_clean, lig_pdb,
                       template_file, charge_model, prod_ns)
        if ok:
            print(f"  ✓ SUCCESS with {charge_model} charges!")
            return True
        print(f"  ✗ {charge_model} failed. Trying next...")
    return False


def _run_once(name, tag, out_dir, prot_clean, lig_pdb,
              template_file, charge_model, prod_ns):
    import glob
    t0 = time.time()

    template_xml = f"{DATA_DIR}/{template_file}"
    gaff_candidates = glob.glob("/usr/local/lib/**/gaff-2.11.xml", recursive=True)
    gaff_xml = gaff_candidates[0] if gaff_candidates else \
        "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"

    print(f"  ForceField: amber14 + gaff-2.11 + {template_file}")
    forcefield = app.ForceField("amber14-all.xml", gaff_xml, template_xml, "tip3p.xml")

    # Build system
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

    # Backbone restraints
    k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
    restraint_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
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

    # Integrator
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

    # Phase 1: Minimize (10K steps)
    print("  Minimize (10000 steps)...")
    try:
        simulation.minimizeEnergy(maxIterations=10000)
    except Exception as e:
        print(f"  ❌ Minimize failed: {e}")
        return False
    if not check_nan("minimize"): return False
    state = simulation.context.getState(getEnergy=True)
    pe0 = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"  ✓ Minimize: PE={pe0:.0f} kcal/mol")

    # Phase 2: SLOW NVT heating (25K × 20ps)
    print("  NVT heating 50→300K (25K × 20ps)...")
    simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
    for t in range(50, 301, 25):  # 50, 75, 100, 125, ..., 300
        simulation.integrator.setTemperature(t * unit.kelvin)
        simulation.step(20_000)  # 20ps per stage
        if not check_nan(f"NVT {t}K"):
            return False
        pe = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print(f"    {t:>3}K: PE={pe.value_in_unit(unit.kilocalories_per_mole):.0f}")

    simulation.integrator.setTemperature(TEMP)
    simulation.step(50_000)  # 50ps eq at 300K
    if not check_nan("NVT 300K eq"): return False
    print("  ✓ NVT heating: PASS")

    # Phase 3: NPT equil
    print("  NPT equil...")
    simulation.integrator.setStepSize(DT)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
    simulation.context.reinitialize(preserveState=True)

    for k_val in [2.5, 0.5]:
        new_k = k_val * unit.kilocalories_per_mole / (unit.angstroms ** 2)
        simulation.context.setParameter("k", new_k)
        simulation.step(25_000)
        if not check_nan(f"NPT k={k_val}"): return False

    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(50_000)
    if not check_nan("NPT unrestrained"): return False
    print("  ✓ NPT equil: PASS")

    # Phase 4: Production
    prod_steps = int(prod_ns * 1_000_000 / 2.0)
    save_freq = 25_000
    print(f"  Production {prod_ns}ns ({prod_steps} steps)...")
    simulation.reporters.append(app.DCDReporter(f"{out_dir}/{tag}_5ns.dcd", save_freq))
    simulation.reporters.append(app.StateDataReporter(
        f"{out_dir}/{tag}_5ns.csv", save_freq,
        step=True, time=True, potentialEnergy=True,
        temperature=True, volume=True, density=True))
    simulation.step(prod_steps)
    if not check_nan("production"): return False

    elapsed = (time.time() - t0) / 60.0
    state = simulation.context.getState(getPositions=True)
    with open(f"{out_dir}/{tag}_final.pdb", "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print(f"  ✓ DONE in {elapsed:.1f} min")
    return True


# ══ Run ══
print(f"\n{'='*55}")
print("SMVT MD v4 — Slow heat + charge fallback")
print(f"{'='*55}")

for name, smiles, mmff_t, gasteiger_t, ns in COMPOUNDS:
    ok = run_md_v4(name, smiles, mmff_t, gasteiger_t, ns)
    print(f"\n{name}: {'PASS' if ok else 'FAIL'}")

print("DONE")
