"""
SMVT MD v6 — Diagnostic: protein-only NVT test.
If protein alone passes 200K, the issue is the ligand/Vina pose.
If protein alone also crashes, the issue is the system setup or OpenCL.
"""
import os, time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

plat = None
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        break
    except: continue
print(f"OpenMM {mm.__version__} | Platform: {plat.getName()}")

DATA_DIR = "/content"
PADDING = 1.0 * unit.nanometers
TEMP = 300.0 * unit.kelvin
DT = 2.0 * unit.femtoseconds

print(f"\n{'='*55}")
print("v6 DIAGNOSTIC: Protein-only NVT test")
print(f"{'='*55}")

# Prepare protein
print("PDBFixer...")
fixer = PDBFixer(filename=f"{DATA_DIR}/AF-Q9Y289-F1.pdb")
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)

# ForceField (no ligand template needed!)
print("ForceField: amber14 + tip3p (no ligand)")
forcefield = app.ForceField("amber14-all.xml", "tip3p.xml")

# Build system
modeller = app.Modeller(fixer.topology, fixer.positions)
modeller.addSolvent(forcefield, model="tip3p", padding=PADDING,
                    ionicStrength=0.0 * unit.molar, neutralize=False)
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds, ignoreExternalBonds=True)
print(f"System: {system.getNumParticles()} particles")

# Restraints
k_rest = 5.0 * unit.kilocalories_per_mole / (unit.angstroms ** 2)
force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addGlobalParameter("k", k_rest)
force.addPerParticleParameter("x0"); force.addPerParticleParameter("y0"); force.addPerParticleParameter("z0")
for atom in modeller.topology.atoms():
    if atom.residue.name not in ("HOH", "WAT", "SOL", "NA", "CL"):
        if atom.name in ("CA", "C", "N"):
            p = modeller.positions[atom.index]
            force.addParticle(atom.index, [p.x, p.y, p.z])
system.addForce(force)

# Simulation
integrator = mm.LangevinMiddleIntegrator(TEMP, 1.0 / unit.picosecond, DT)
simulation = app.Simulation(modeller.topology, system, integrator, plat,
                            {"Precision": "mixed"})
simulation.context.setPositions(modeller.positions)

def check(label):
    pe = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    pe_val = pe.value_in_unit(unit.kilocalories_per_mole)
    ok = not (np.isnan(pe_val) or np.isinf(pe_val))
    print(f"  {label}: PE={pe_val:.0f} kcal/mol {'✓' if ok else '❌ NAN!'}")
    return ok

# Minimize
print("Minimize 5000...")
simulation.minimizeEnergy(maxIterations=5000)
if not check("minimize"): exit(1)

# NVT
print("NVT heating 50→300K...")
simulation.integrator.setStepSize(1.0 * unit.femtoseconds)
for t in [50, 100, 150, 200, 250, 300]:
    simulation.integrator.setTemperature(t * unit.kelvin)
    simulation.step(10_000)
    if not check(f"NVT {t}K"):
        print(f"❌ PROTEIN ITSELF CRASHES AT {t}K!")
        print("→ Problem is system setup or OpenCL, not the ligand.")
        exit(1)

print("\n✓ PROTEIN PASSES NVT 300K! Problem is the ligand/Vina pose.")
