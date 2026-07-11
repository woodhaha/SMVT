"""
MM-GBSA — SMVT-8 compounds.
OpenMM amber14 + GB(OBC2) + GAFF-2.11.

Protocol:
1. Load complex PDB (protein+LIG from _final.pdb, waters stripped)
2. Build FF system, restrained minimization
3. Build receptor/ligand subsystems from minimized complex
4. ΔG = E_complex - E_receptor - E_ligand

Usage:
    python _run_mmgbsa.py [COMPOUND ...]
    No args = run all 8.
"""
import os, sys, json, time, warnings
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
warnings.filterwarnings("ignore")
KB = unit.kilojoule_per_mole

SMVT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ANALYSIS_DIR = os.path.join(SMVT_ROOT, "MM-GBSA Analysis")
LIGANDS_DIR = os.path.join(SMVT_ROOT, "SMVT_MD", "ligands")
TRAJ_DIR = os.path.join(SMVT_ROOT, "SMVT_MD", "trajectories")
SYNC_PATH = os.path.join(SMVT_ROOT, "SMVT-MD-Analysis", "_SYNC_PENDING.md")
GAFF_XML = "C:/anaconda3/envs/smvt-md/Lib/site-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
TEMP = 300.0 * unit.kelvin
CUTOFF = 2.0 * unit.nanometer


def build(compound):
    """Load protein+LIG complex PDB, assign FF."""
    complex_pdb = os.path.join(ANALYSIS_DIR, f"{compound}_complex.pdb")
    lig_tmpl = os.path.join(LIGANDS_DIR, f"{compound}_template.xml")
    for f in [complex_pdb, lig_tmpl, GAFF_XML]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing: {f}")

    pdb = app.PDBFile(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl)
    system = ff.createSystem(pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)

    lig_atoms = set(a.index for a in pdb.topology.atoms()
                    if a.residue.name.strip() == "LIG")
    prot_atoms = set(range(pdb.topology.getNumAtoms())) - lig_atoms
    return system, pdb, lig_atoms, prot_atoms


def minimize(ctx, system, topology, positions):
    """Cα + ligand restrained minimization."""
    caf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    caf.addGlobalParameter("k", 500.0 * KB / (unit.nanometer**2))
    caf.addPerParticleParameter("x0")
    caf.addPerParticleParameter("y0")
    caf.addPerParticleParameter("z0")
    for a in topology.atoms():
        if a.name.strip() == "CA":
            x, y, z = positions[a.index].value_in_unit(unit.nanometer)
            caf.addParticle(a.index, [x, y, z])
    system.addForce(caf)
    ctx.reinitialize(preserveState=True)
    ctx.setPositions(positions)
    try:
        mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=40, tolerance=100.0)
    except Exception:
        pass
    system.removeForce(system.getNumForces() - 1)
    ctx.reinitialize(preserveState=True)
    return ctx.getState(getPositions=True).getPositions()


def subsystem(system, topology, positions, keep_resnames, ff_files):
    """Extract subsystem by residue name(s), rebuild FF."""
    mod = app.Modeller(topology, positions)
    to_del = [a for a in mod.topology.atoms()
              if a.residue.name.strip() not in keep_resnames]
    mod.delete(to_del)
    ff = app.ForceField(*ff_files)
    sys = ff.createSystem(mod.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)
    return sys, mod.positions


def run(compound):
    print(f"\n{'='*60}", flush=True)
    print(f"  {compound}", flush=True)
    t0 = time.time()

    sys, pdb, lig_a, prot_a = build(compound)
    natoms = pdb.topology.getNumAtoms()
    print(f"  Complex: {natoms} atoms ({len(lig_a)} lig) | "
          f"built in {time.time()-t0:.1f}s", flush=True)

    cpu = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(sys, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    pos_min = minimize(ctx, sys, pdb.topology, pdb.positions)
    ctx.setPositions(pos_min)
    e_min = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
    print(f"  E_min: {e_min:.0f} kJ/mol ({time.time()-t0:.1f}s)", flush=True)

    lig_tmpl = os.path.join(LIGANDS_DIR, f"{compound}_template.xml")
    sys_rec, pos_rec = subsystem(sys, pdb.topology, pos_min,
        {r.name.strip() for r in pdb.topology.residues() if r.name.strip() != "LIG"},
        ["amber14-all.xml", "implicit/obc2.xml"])
    sys_lig, pos_lig = subsystem(sys, pdb.topology, pos_min, {"LIG"},
        ["amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl])
    print(f"  Receptor: {sys_rec.getNumParticles()} | Ligand: {sys_lig.getNumParticles()} atoms", flush=True)

    ctx_r = mm.Context(sys_rec, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx_l = mm.Context(sys_lig, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx.setPositions(pos_min)
    ctx_r.setPositions(pos_rec)
    ctx_l.setPositions(pos_lig)

    ec = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
    er = ctx_r.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
    el = ctx_l.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
    dg = ec - er - el
    dg_kcal = dg / 4.184

    print(f"  E_complex: {ec:+.1f} | E_receptor: {er:+.1f} | E_ligand: {el:+.1f}")
    print(f"  ΔG = {dg:+.1f} kJ/mol = {dg_kcal:+.2f} kcal/mol", flush=True)

    out_dir = os.path.join(ANALYSIS_DIR, compound)
    os.makedirs(out_dir, exist_ok=True)
    result = {
        "compound": compound,
        "dG_kJ": round(dg, 2), "dG_kcal": round(dg_kcal, 2),
        "E_complex_kJ": round(ec, 2), "E_receptor_kJ": round(er, 2),
        "E_ligand_kJ": round(el, 2), "E_min_kJ": round(e_min, 2),
        "n_atoms": natoms, "n_lig_atoms": len(lig_a),
    }
    with open(os.path.join(out_dir, "mmgbsa_results.json"), "w") as f:
        json.dump(result, f, indent=2)
    print(f"  -> {out_dir}/ in {time.time()-t0:.0f}s", flush=True)
    return result


def main():
    targets = sys.argv[1:] if len(sys.argv) > 1 else [
        "NAFTAZONE", "BIOTIN", "ESKETAMINE", "FUROSEMIDE",
        "GABAPENTIN_ENACARBIL", "HYDROMORPHONE", "PHENOBARBITAL", "RIBOFLAVIN"]
    results = {}
    for c in targets:
        try:
            r = run(c.upper())
            if r: results[c.upper()] = r
        except Exception as e:
            print(f"  ERROR {c}: {e}", flush=True)
            import traceback; traceback.print_exc()

    if len(results) > 1:
        print(f"\n{'='*60}")
        print(f"  MM-GBSA SUMMARY")
        print(f"{'='*60}")
        print(f"  {'Compound':28s}  {'ΔG(kJ/mol)':>10s}  {'ΔG(kcal/mol)':>10s}")
        print(f"  {'-'*54}")
        for c, r in sorted(results.items()):
            print(f"  {c:28s}  {r['dG_kJ']:8.1f} kJ  {r['dG_kcal']:8.2f} kcal")
        print(f"\n  Results: {os.path.join(ANALYSIS_DIR, '<COMPOUND>/mmgbsa_results.json')}")

if __name__ == "__main__":
    main()
