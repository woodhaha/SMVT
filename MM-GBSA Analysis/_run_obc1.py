"""Run OBC1 MM-GBSA for all 8 SMVT compounds. Saves obc1_results.json."""
import os, sys, json, time, warnings
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
warnings.filterwarnings("ignore")
KB = unit.kilojoule_per_mole
CUTOFF = 2.0 * unit.nanometer
TEMP = 300.0 * unit.kelvin

ANALYSIS = r"D:\Researching\SMVT\MM-GBSA Analysis"
LIGANDS = r"D:\Researching\SMVT\SMVT_MD\ligands"
GAFF_XML = "C:/anaconda3/envs/smvt-md/Lib/site-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]

def run_obc1(compound):
    complex_pdb = os.path.join(ANALYSIS, f"{compound}_complex.pdb")
    lig_tmpl = os.path.join(LIGANDS, f"{compound}_template.xml")
    for f in [complex_pdb, lig_tmpl, GAFF_XML]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing: {f}")
    pdb = app.PDBFile(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "implicit/obc1.xml", GAFF_XML, lig_tmpl)
    system = ff.createSystem(pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF, constraints=app.HBonds)

    caf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    caf.addGlobalParameter("k", 500.0 * KB / (unit.nanometer**2))
    caf.addPerParticleParameter("x0"); caf.addPerParticleParameter("y0"); caf.addPerParticleParameter("z0")
    for a in pdb.topology.atoms():
        if a.name.strip() == "CA":
            x,y,z = pdb.positions[a.index].value_in_unit(unit.nanometer)
            caf.addParticle(a.index, [x,y,z])
    system.addForce(caf)
    cpu = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx.setPositions(pdb.positions)
    mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=40, tolerance=100.0)
    system.removeForce(system.getNumForces()-1)
    ctx.reinitialize(preserveState=True)
    pos_min = ctx.getState(getPositions=True).getPositions()
    ec = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)

    mod_r = app.Modeller(pdb.topology, pos_min)
    to_del = [a for a in mod_r.topology.atoms() if a.residue.name.strip() == "LIG"]
    mod_r.delete(to_del)
    ff_rec = app.ForceField("amber14-all.xml", "implicit/obc1.xml")
    sys_rec = ff_rec.createSystem(mod_r.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
    ctx_r = mm.Context(sys_rec, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx_r.setPositions(mod_r.positions)
    er = ctx_r.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)

    mod_l = app.Modeller(pdb.topology, pos_min)
    to_del_l = [a for a in mod_l.topology.atoms() if a.residue.name.strip() != "LIG"]
    mod_l.delete(to_del_l)
    ff_lig = app.ForceField("amber14-all.xml", "implicit/obc1.xml", GAFF_XML, lig_tmpl)
    sys_lig = ff_lig.createSystem(mod_l.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
    ctx_l = mm.Context(sys_lig, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx_l.setPositions(mod_l.positions)
    el = ctx_l.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)

    dg_kcal = (ec - er - el) / 4.184
    return round(dg_kcal, 2)

results = {}
for c in COMPS:
    t0 = time.time()
    try:
        dg = run_obc1(c)
        results[c] = dg
        print(f"{c:28s}  OBC1={dg:+.2f} kcal/mol  ({time.time()-t0:.0f}s)", flush=True)
    except Exception as e:
        print(f"{c:28s}  ERROR: {e}", flush=True)

out = os.path.join(ANALYSIS, "obc1_results.json")
with open(out, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSaved: {out}")
print("\nOBC1 Summary:")
for c, v in sorted(results.items(), key=lambda x: x[1]):
    print(f"  {c:28s}  {v:+.2f}")
