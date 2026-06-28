#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Run HYDROMORPHONE 100ns MD on AutoDL RTX 4090
# ═══════════════════════════════════════════════════════════════
set -e

source activate smvt-md

echo "========================================="
echo " SMVT MD — HYDROMORPHONE 100ns"
echo " GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo " Start: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================="

cd "$(dirname "$0")"

python -c "
import os, sys, time, logging, warnings
warnings.filterwarnings('ignore')
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule as OFFMolecule
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.units import unit as off_unit

# ── Config ──
PROJECT_DIR = os.path.dirname(os.path.abspath('__file__'))
MD_DIR = os.path.join(PROJECT_DIR, 'trajectories')
RECEPTOR_PDB = os.path.join(PROJECT_DIR, 'receptor', 'SMVT_prepared.pdb')
DOCKING_DIR = os.path.join(PROJECT_DIR, 'ligands')

NAME = 'HYDROMORPHONE'
SMILES = 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5'
PROD_NS = 100

os.makedirs(MD_DIR, exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.FileHandler(f'{MD_DIR}/hydromorphone_100ns.log', mode='w'),
              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

log.info(f'HYDROMORPHONE 100ns MD — Start')

# ── 1. Prepare protein ──
log.info('Step 1/9: Preparing protein...')
protein_clean = f'{MD_DIR}/protein_prepared.pdb'
if not os.path.exists(protein_clean):
    fixer = PDBFixer(filename=RECEPTOR_PDB)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(protein_clean, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
log.info(f'  Protein saved: {protein_clean}')

# ── 2. Extract Vina pose ──
log.info('Step 2/9: Loading Vina docking pose...')
pdbqt = f'{DOCKING_DIR}/{NAME}.pdbqt'
vina_pos = []
with open(pdbqt) as f:
    for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            try:
                vina_pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except: pass
        if len(vina_pos) > 0 and not (line.startswith('HETATM') or line.startswith('ATOM')):
            break
log.info(f'  Loaded {len(vina_pos)} atoms from Vina pose')

# ── 3. Parameterize ligand ──
log.info('Step 3/9: Parameterizing ligand (OpenFF Sage + Gasteiger)...')
off_mol = OFFMolecule.from_smiles(SMILES, allow_undefined_stereo=True)
off_mol.generate_conformers(n_conformers=1)
conf = off_mol.conformers[0]
for i in range(min(off_mol.n_atoms, len(vina_pos))):
    conf[i] = np.array(vina_pos[i]) * off_unit.angstrom
off_mol.assign_partial_charges(partial_charge_method='gasteiger')

# Save ligand PDB for Modeller
lig_pdb_path = f'{MD_DIR}/ligand.pdb'
rd_mol = Chem.MolFromSmiles(SMILES)
rd_mol = Chem.AddHs(rd_mol)
AllChem.EmbedMolecule(rd_mol, randomSeed=42)
rd_conf = rd_mol.GetConformer()
for i, pos in enumerate(vina_pos):
    if i < rd_mol.GetNumAtoms():
        rd_conf.SetAtomPosition(i, (pos[0]/10.0, pos[1]/10.0, pos[2]/10.0))
Chem.MolToPDBFile(rd_mol, lig_pdb_path)

# ── 4. Build system ──
log.info('Step 4/9: Building solvated system (ff14SB + TIP3P + 0.15M NaCl, 8A padding)...')
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='openff-2.1.0',
    molecules=[off_mol],
    cache=f'{MD_DIR}/openff_cache.json'
)
prot_pdb = app.PDBFile(protein_clean)
lig_pdb = app.PDBFile(lig_pdb_path)
modeller = app.Modeller(prot_pdb.topology, prot_pdb.positions)
modeller.add(lig_pdb.topology, lig_pdb.positions)
modeller.addSolvent(system_generator.forcefield, model='tip3p',
                    padding=0.8*unit.nanometers, ionicStrength=0.15*unit.molar, neutralize=True)
system = system_generator.create_system(modeller.topology, molecules=[off_mol])
n_particles = system.getNumParticles()
n_water = (n_particles - len(vina_pos)) / 3
log.info(f'  System: {n_particles} particles (~{n_water:.0f} waters)')

# ── 5. Setup simulation ──
log.info('Step 5/9: Setting up integrator + OpenCL platform...')
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
platform = mm.Platform.getPlatformByName('OpenCL')
simulation = app.Simulation(modeller.topology, system, integrator, platform, {'Precision': 'mixed'})
simulation.context.setPositions(modeller.positions)
log.info(f'  Platform: {platform.getName()}')

# ── 6. Minimization ──
log.info('Step 6/9: Energy minimization (5000 steps)...')
simulation.minimizeEnergy(maxIterations=5000)
state = simulation.context.getState(getEnergy=True)
log.info(f'  Final potential: {state.getPotentialEnergy()}')

# ── 7. NVT heating ──
log.info('Step 7/9: NVT equilibration (50K→300K gradual heating)...')
simulation.integrator.setStepSize(1.0*unit.femtoseconds)
for t in [50, 150, 300]:
    simulation.integrator.setTemperature(t*unit.kelvin)
    simulation.step(10000)  # 10ps each
simulation.integrator.setTemperature(300*unit.kelvin)
simulation.integrator.setStepSize(2.0*unit.femtoseconds)
simulation.step(35000)  # remaining 70ps
log.info('  NVT done')

# ── 8. NPT ──
log.info('Step 8/9: NPT equilibration (100ps)...')
system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300*unit.kelvin))
simulation.context.reinitialize(preserveState=True)
simulation.step(50000)
log.info('  NPT done')

# ── 9. Production 100ns ──
log.info(f'Step 9/9: Production {PROD_NS}ns...')
out_dir = f'{MD_DIR}/{NAME}_100ns'
os.makedirs(out_dir, exist_ok=True)
simulation.reporters.append(app.DCDReporter(f'{out_dir}/{NAME}_prod.dcd', 25000))
simulation.reporters.append(app.StateDataReporter(
    f'{out_dir}/{NAME}_prod.csv', 25000,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    temperature=True, volume=True, density=True
))
simulation.reporters.append(app.CheckpointReporter(f'{out_dir}/{NAME}_prod.chk', 250000))

prod_steps = int(PROD_NS * 500000)  # 2fs timestep, 1ns = 500000 steps
log.info(f'  Total steps: {prod_steps}')

t0 = time.time()
n_reports = 20
for i in range(n_reports):
    chunk = prod_steps // n_reports
    simulation.step(chunk)
    elapsed = (time.time() - t0) / 3600
    done_pct = (i+1)/n_reports*100
    ns_done = (i+1)/n_reports * PROD_NS
    rate_ns_day = ns_done / elapsed * 24 if elapsed > 0.01 else 999
    eta_h = elapsed / done_pct * 100 - elapsed
    log.info(f'  [{done_pct:.0f}%] {ns_done:.1f}ns | {elapsed:.1f}h elapsed | ~{rate_ns_day:.1f} ns/day | ETA {eta_h:.1f}h')

# Save final frame
final = simulation.context.getState(getPositions=True)
with open(f'{out_dir}/{NAME}_final.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, final.getPositions(), f)
simulation.saveCheckpoint(f'{out_dir}/{NAME}_prod.chk')

total_h = (time.time() - t0) / 3600
log.info(f'')
log.info(f'========================================')
log.info(f' DONE! HYDROMORPHONE 100ns completed')
log.info(f' Total time: {total_h:.1f} hours')
log.info(f' Rate: {PROD_NS/total_h:.1f} ns/h = {PROD_NS/total_h*24:.1f} ns/day')
log.info(f' Output: {out_dir}/')
log.info(f'========================================')
"
