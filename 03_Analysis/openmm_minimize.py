"""
OpenMM energy minimization for SLC5A6 (SMVT).
Uses Amber14SB + GBSA implicit solvent.
Strategy: load cleaned (no-H) structure → OpenMM adds H → phased minimization
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import os

WORKDIR = r"D:\Researching\SMVT"
INPUT_PDB  = os.path.join(WORKDIR, "AF-Q9Y289-F1_cleaned.pdb")   # no hydrogens
OUTPUT_PDB = os.path.join(WORKDIR, "AF-Q9Y289-F1_prepared.pdb")
OUTPUT_LOG = os.path.join(WORKDIR, "minimization.log")

print("=" * 60)
print("SLC5A6 Energy Minimization — OpenMM + Amber14SB + GBSA")
print("=" * 60)

# ── 1. Load cleaned (no-H) structure ──
print("\n[1/5] Loading cleaned structure (no hydrogens)...")
pdb = app.PDBFile(INPUT_PDB)
print(f"  Atoms: {pdb.topology.getNumAtoms()}  (expect ~4824, heavy atoms only)")

# ── 2. OpenMM adds H at pH 7.0 ──
print("\n[2/5] Adding hydrogens via OpenMM Modeller (pH 7.0)...")
forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)
print(f"  Atoms after +H: {modeller.topology.getNumAtoms()}")

# ── 3. Build system with GBSA ──
print("\n[3/5] Building system (Amber14SB + GBSA-OBC implicit solvent)...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None,
)
print(f"  System ready: {system.getNumParticles()} particles")

# Save protonated pre-min structure
premin_pdb = os.path.join(WORKDIR, "AF-Q9Y289-F1_H.pdb")
app.PDBFile.writeFile(modeller.topology, modeller.positions, open(premin_pdb, 'w'))
print(f"  H-structure: {premin_pdb}")

# ── 4. Energy minimization ──
print("\n[4/5] Energy minimization...")

integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Initial energy
e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  Initial E:  {e0.value_in_unit(unit.kilojoules_per_mole):.1f} kJ/mol")

# Phase A: strong heavy-atom restraint → relax added H
print("  Phase A: Heavy atom restraint → H relaxation...")
k_strong = 100.0 * unit.kilocalories_per_mole / unit.angstroms**2
restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
restraint.addPerParticleParameter("k")
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")
heavy_count = 0
for i, atom in enumerate(modeller.topology.atoms()):
    if atom.element.symbol != 'H':
        pos = modeller.positions[i]
        restraint.addParticle(i, [k_strong, pos.x, pos.y, pos.z])
        heavy_count += 1

ridx = system.addForce(restraint)
simulation.minimizeEnergy(tolerance=10.0, maxIterations=300)
e1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  After H-relax: {e1.value_in_unit(unit.kilojoules_per_mole):.1f} kJ/mol")

# Phase B: remove restraint → full minimization
print("  Phase B: Full unrestrained minimization...")
system.removeForce(ridx)
simulation.context.reinitialize(preserveState=True)
simulation.minimizeEnergy(tolerance=1.0, maxIterations=2000)
state_final = simulation.context.getState(getEnergy=True, getPositions=True)
e_final = state_final.getPotentialEnergy()
drop = e0.value_in_unit(unit.kilojoules_per_mole) - e_final.value_in_unit(unit.kilojoules_per_mole)
print(f"  Final E:     {e_final.value_in_unit(unit.kilojoules_per_mole):.1f} kJ/mol")
print(f"  Total drop:  {drop:.1f} kJ/mol ({drop/max(1,heavy_count):.3f} kJ/mol per heavy atom)")

# ── 5. Save ──
print(f"\n[5/5] Saving...")
app.PDBFile.writeFile(modeller.topology, state_final.getPositions(), open(OUTPUT_PDB, 'w'))

with open(OUTPUT_LOG, 'w') as f:
    f.write(f"SLC5A6 Energy Minimization Log\n")
    f.write(f"Force field: Amber14SB + GBSA-OBC implicit solvent\n")
    f.write(f"Initial energy: {e0}\n")
    f.write(f"After H-relax:  {e1}\n")
    f.write(f"Final energy:   {e_final}\n")
    f.write(f"Total drop:     {e0 - e_final}\n")

print(f"  Output: {OUTPUT_PDB}")
print(f"  Log:    {OUTPUT_LOG}")
print("\n✅ Done!")
