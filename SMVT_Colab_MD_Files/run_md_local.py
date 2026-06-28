#!/usr/bin/env python3
"""
SMVT-Phenobarbital 100ns MD — Local GPU Execution
GPU: RTX 2000 Ada (8GB) | OpenMM OpenCL | Amber14SB + GAFF2 + TIP3P
"""
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem
import time, json, os, sys, argparse
from pathlib import Path

# ═══ Config ═══
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "02_Data" / "cleaned"
DOCKING_DIR = PROJECT_ROOT / "03_Analysis" / "docking"
OUT_BASE = PROJECT_ROOT / "03_Analysis" / "md_output"
os.makedirs(OUT_BASE, exist_ok=True)

RECEPTOR_PDB = str(DATA_DIR / "AF-Q9Y289-F1_prepared.pdb")  # OpenMM minimized
COMPOUND = ('PHENOBARBITAL', 'CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2')

PROD_NS = 100
TEMP = 310 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
DT = 2.0 * unit.femtoseconds
CUTOFF = 1.0 * unit.nanometer
SAVE_PS = 100  # frame every 100ps → 1000 frames for 100ns

# ═══ Platform selection ═══
def get_platform():
    """Try CUDA first, then OpenCL, then CPU."""
    for name in ['CUDA', 'OpenCL', 'CPU']:
        try:
            p = mm.Platform.getPlatformByName(name)
            print(f'Using platform: {p.getName()}')
            return p
        except Exception:
            continue
    raise RuntimeError('No OpenMM platform available')

PLATFORM = get_platform()
OUT_DIR = str(OUT_BASE / f"phenobarbital_{PROD_NS}ns")
os.makedirs(OUT_DIR, exist_ok=True)

print(f'GPU MD: {COMPOUND[0]} x {PROD_NS}ns')
print(f'Platform: {PLATFORM.getName()}, Temp: {TEMP}, Pressure: {PRESSURE}')
print(f'Output: {OUT_DIR}')

# ═══ Helper functions ═══
def extract_vina_pose(pdbqt_path):
    """Extract first-model coordinates from Vina PDBQT."""
    pos = []; in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith('MODEL'):
                if int(line.split()[1]) == 1: in_model = True; pos = []
                elif in_model: break
            elif line.startswith('ENDMDL') and in_model: break
            elif in_model and (line.startswith('ATOM') or line.startswith('HETATM')):
                try: pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: continue
    return pos

def prepare_ligand_pdb(smiles, vina_pos, out_path):
    """Generate ligand PDB from SMILES + Vina docking pose."""
    if os.path.exists(out_path): return out_path
    print('  Building ligand from SMILES + Vina pose...')
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for i in range(min(mol.GetNumAtoms(), len(vina_pos))):
        conf.SetAtomPosition(i, (vina_pos[i][0]/10.0, vina_pos[i][1]/10.0, vina_pos[i][2]/10.0))
    Chem.MolToPDBFile(mol, out_path)
    return out_path

# ═══ Main MD ═══
name, smiles = COMPOUND
tag = name
chk_file = f'{OUT_DIR}/{tag}_{PROD_NS}ns.chk'
t_start = time.time()

sep = '=' * 60
print(f'\n{sep}')
print(f'[{name}] {PROD_NS}ns GPU MD (GAFF2, TIP3P, 310K)')
print(f'{sep}')

# 1. Prepare protein
prot_clean = f'{OUT_DIR}/protein.pdb'
if not os.path.exists(prot_clean):
    print('  Preparing protein...')
    fixer = PDBFixer(filename=RECEPTOR_PDB)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(prot_clean, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
else:
    print('  Protein already prepared.')

# 2. Extract Vina pose → ligand PDB
pdbqt = str(DOCKING_DIR / f'{name}_docked.pdbqt')
assert os.path.exists(pdbqt), f'PDBQT missing: {pdbqt}'
vina_pos = extract_vina_pose(pdbqt)
print(f'  Vina pose: {len(vina_pos)} atoms')
lig_pdb = f'{OUT_DIR}/ligand.pdb'
prepare_ligand_pdb(smiles, vina_pos, lig_pdb)

# 3. SystemGenerator (GAFF2)
print('  Parameterizing with GAFF2...')
system_gen = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='gaff-2.11',
    cache=f'{OUT_DIR}/cache.json'
)

# 4. Build solvated system
print('  Building solvated system (TIP3P, 0.15M NaCl, 8A pad)...')
prot_pdb = app.PDBFile(prot_clean)
lig = app.PDBFile(lig_pdb)
modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
modeller.add(lig.topology, lig.positions)
modeller.addSolvent(system_gen.forcefield, model='tip3p',
                    padding=0.8*unit.nanometers, ionicStrength=0.15*unit.molar,
                    neutralize=True)
system = system_gen.create_system(modeller.topology)
print(f'  System: {system.getNumParticles()} particles')

# 5. Simulation setup
sim_params = {'Precision': 'mixed'} if PLATFORM.getName() == 'CUDA' else {}
integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0/unit.picosecond, DT)
simulation = app.Simulation(modeller.topology, system, integrator, PLATFORM, sim_params)
simulation.context.setPositions(modeller.positions)

# 6. Energy minimization
print('  Minimizing...')
simulation.minimizeEnergy(maxIterations=5000)
energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f'  Minimized energy: {energy}')

# 7. NVT heating (50→150→310K, 100ps total with 1fs timestep)
print('  NVT heating (50→150→310K)...')
simulation.integrator.setStepSize(1.0*unit.femtoseconds)
for t in [50, 150, 310]:
    simulation.integrator.setTemperature(t*unit.kelvin)
    simulation.step(int(33*1000/1.0))
simulation.integrator.setTemperature(TEMP)
simulation.integrator.setStepSize(DT)

# 8. NPT equilibration (200ps)
print('  NPT equil (200ps)...')
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
simulation.context.reinitialize(preserveState=True)
simulation.step(int(200*1000/2.0))

# 9. Production (100ns)
prod_steps = int(PROD_NS * 1_000_000 / 2.0)
save_freq = int(SAVE_PS * 1000 / 2.0)
print(f'  Production: {PROD_NS}ns ({prod_steps:,} steps, save every {SAVE_PS}ps)...')
print(f'  T4 estimated: ~12h | RTX 2000 Ada estimated: ~8-10h')

simulation.reporters.append(app.DCDReporter(f'{OUT_DIR}/{tag}_{PROD_NS}ns.dcd', save_freq))
simulation.reporters.append(app.StateDataReporter(
    f'{OUT_DIR}/{tag}_{PROD_NS}ns.csv', save_freq,
    step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True
))
simulation.reporters.append(app.CheckpointReporter(chk_file, save_freq*10))

# Progress callback
last_report = [0]
def progress_callback(step, total):
    pct = step / total * 100
    if step - last_report[0] >= total // 20:  # report every 5%
        elapsed = (time.time() - t_start) / 60
        eta = elapsed / (step/total) - elapsed if step > 0 else 0
        print(f'  [{step/total*100:.0f}%] {step:,}/{total:,} steps | elapsed: {elapsed:.0f}min | ETA: {eta:.0f}min')
        last_report[0] = step

try:
    # Run with progress tracking
    chunk_size = int(prod_steps / 20)
    for i in range(20):
        simulation.step(chunk_size)
        elapsed = (time.time() - t_start) / 60
        eta = elapsed / (i+1) * (20-i-1)
        print(f'  [{(i+1)*5}%] {(i+1)*chunk_size:,}/{prod_steps:,} | elapsed: {elapsed:.0f}min | ETA: {eta:.0f}min')
except KeyboardInterrupt:
    print('\n  Interrupted — saving checkpoint...')
    simulation.saveCheckpoint(chk_file)
    print(f'  Checkpoint saved: {chk_file}')
    sys.exit(0)

elapsed = (time.time() - t_start) / 60
ns_per_day = PROD_NS / (elapsed / 60 / 24) if elapsed > 0 else 0
print(f'  [{name}] DONE in {elapsed:.1f}min (~{ns_per_day:.0f} ns/day)')

# Save final frame
state = simulation.context.getState(getPositions=True)
with open(f'{OUT_DIR}/{tag}_final.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
simulation.saveCheckpoint(chk_file)

# ═══ Save metadata ═══
metadata = {
    'compound': name,
    'smiles': smiles,
    'production_ns': PROD_NS,
    'temperature_k': 310,
    'force_field': 'Amber14SB + GAFF2 + TIP3P',
    'platform': PLATFORM.getName(),
    'elapsed_min': round(elapsed, 1),
    'ns_per_day': round(ns_per_day),
    'docking_score': -8.30,
    'receptor': RECEPTOR_PDB,
    'output_dir': OUT_DIR,
}
with open(f'{OUT_DIR}/metadata.json', 'w') as f:
    json.dump(metadata, f, indent=2)

print(f'\n{sep}')
print(f'Trajectory:  {OUT_DIR}/{tag}_{PROD_NS}ns.dcd')
print(f'State data:  {OUT_DIR}/{tag}_{PROD_NS}ns.csv')
print(f'Checkpoint:  {chk_file}')
print(f'Final PDB:   {OUT_DIR}/{tag}_final.pdb')
print(f'Metadata:    {OUT_DIR}/metadata.json')
print(f'{sep}')
print(f'Completed in {elapsed:.1f} minutes (~{ns_per_day:.0f} ns/day)')
print(f'Now run: python SMVT_Colab_MD_Files/analyze_md.py --dir {OUT_DIR}')
