"""
Energy Decomposition for MM-GBSA (SMVT-8).
Force groups: NonbondedForce (Coulomb+LJ) → group 0
             CustomGBForce (GB+SA)     → group 1
             All other (bonds, angles) → group 2
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
GAFF_XML = "C:/anaconda3/envs/smvt-md/Lib/site-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
TEMP = 300.0 * unit.kelvin
CUTOFF = 2.0 * unit.nanometer


def assign_force_groups(system):
    """Assign force groups for energy decomposition."""
    nb_group = 0   # NonbondedForce: Coulomb + LJ
    gb_group = 1   # CustomGBForce: GB + SA
    other_group = 2

    for i in range(system.getNumForces()):
        f = system.getForce(i)
        name = f.__class__.__name__
        if name == 'NonbondedForce':
            f.setForceGroup(nb_group)
        elif name == 'CustomGBForce':
            f.setForceGroup(gb_group)
        else:
            f.setForceGroup(other_group)


def build(compound):
    """Load complex PDB, assign FF with force groups."""
    complex_pdb = os.path.join(ANALYSIS_DIR, f"{compound}_complex.pdb")
    lig_tmpl = os.path.join(LIGANDS_DIR, f"{compound}_template.xml")
    pdb = app.PDBFile(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl)
    system = ff.createSystem(pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)
    assign_force_groups(system)

    lig_atoms = set(a.index for a in pdb.topology.atoms()
                    if a.residue.name.strip() == "LIG")
    return system, pdb, lig_atoms


def minimize(system, pdb):
    """Calpha + ligand restrained minimization (matching _run_mmgbsa.py)."""
    caf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    caf.addGlobalParameter("k", 500.0 * KB / (unit.nanometer**2))
    caf.addPerParticleParameter("x0")
    caf.addPerParticleParameter("y0")
    caf.addPerParticleParameter("z0")
    for a in pdb.topology.atoms():
        if a.name.strip() == "CA":
            x, y, z = pdb.positions[a.index].value_in_unit(unit.nanometer)
            caf.addParticle(a.index, [x, y, z])
    system.addForce(caf)
    caf.setForceGroup(3)  # restraint group

    cpu = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx.setPositions(pdb.positions)

    try:
        mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=40, tolerance=100.0)
    except Exception:
        pass

    system.removeForce(system.getNumForces() - 1)
    ctx.reinitialize(preserveState=True)

    positions = ctx.getState(getPositions=True).getPositions()
    return ctx, positions


def decompose(system, topology, positions, lig_atoms):
    """Group energy by force type for a given system state."""
    ctx = mm.Context(system, mm.LangevinIntegrator(TEMP, 1/unit.picosecond, 1*unit.femtosecond),
                     mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(positions)

    # Total
    s = ctx.getState(getEnergy=True)
    total = s.getPotentialEnergy().value_in_unit(KB)

    # By force group
    e_nb = ctx.getState(getEnergy=True, groups={0}).getPotentialEnergy().value_in_unit(KB)   # Coulomb+LJ
    e_gb = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(KB)   # GB+SA
    e_other = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(KB)  # bonded
    del ctx

    return total, e_nb, e_gb, e_other


def run(compound):
    print(f"\n{'='*60}")
    print(f"  ENERGY DECOMPOSITION: {compound}")
    t0 = time.time()

    system, pdb, lig_atoms = build(compound)
    natoms = pdb.topology.getNumAtoms()
    prot_atoms = set(range(natoms)) - lig_atoms

    # Minimize
    ctx, pos_min = minimize(system, pdb)
    e_min = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
    del ctx

    # Decompose complex
    ec_total, ec_nb, ec_gb, ec_other = decompose(system, pdb.topology, pos_min, lig_atoms)

    # Receptor-only subsystem (remove LIG)
    mod_rec = app.Modeller(pdb.topology, pos_min)
    to_del_rec = [a for a in mod_rec.topology.atoms() if a.residue.name.strip() == "LIG"]
    mod_rec.delete(to_del_rec)

    ff_rec = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
    sys_rec = ff_rec.createSystem(mod_rec.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)
    assign_force_groups(sys_rec)
    er_total, er_nb, er_gb, er_other = decompose(sys_rec, mod_rec.topology, mod_rec.positions, set())

    # Ligand-only subsystem (keep only LIG)
    mod_lig = app.Modeller(pdb.topology, pos_min)
    to_del_lig = [a for a in mod_lig.topology.atoms() if a.residue.name.strip() != "LIG"]
    mod_lig.delete(to_del_lig)

    lig_tmpl = os.path.join(LIGANDS_DIR, f"{compound}_template.xml")
    ff_lig = app.ForceField("amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl)
    sys_lig = ff_lig.createSystem(mod_lig.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)
    assign_force_groups(sys_lig)
    el_total, el_nb, el_gb, el_other = decompose(sys_lig, mod_lig.topology, mod_lig.positions, set())

    # Decomposition
    dg_total = ec_total - er_total - el_total
    dg_nb = ec_nb - er_nb - el_nb       # ΔE_coulomb + ΔE_LJ
    dg_gb = ec_gb - er_gb - el_gb       # ΔG_GB
    dg_other = ec_other - er_other - el_other  # bonded (should ~0)
    dg_sa = dg_total - dg_nb - dg_gb - dg_other  # surface area within GB

    print(f"  {'Component':>15s}  {'Complex':>10s}  {'Receptor':>10s}  {'Ligand':>10s}  {'Delta':>10s}")
    print(f"  {'-'*55}")
    print(f"  {'Total (kJ/mol)':15s}  {ec_total:>10.1f}  {er_total:>10.1f}  {el_total:>10.1f}  {dg_total:>+10.1f}")
    print(f"  {'Coulomb+LJ':15s}  {ec_nb:>10.1f}  {er_nb:>10.1f}  {el_nb:>10.1f}  {dg_nb:>+10.1f}")
    print(f"  {'GB polar':15s}  {ec_gb:>10.1f}  {er_gb:>10.1f}  {el_gb:>10.1f}  {dg_gb:>+10.1f}")
    print(f"  {'SA (implied)':15s}  {'':>10s}  {'':>10s}  {'':>10s}  {dg_sa:>+10.1f}")
    print(f"  {'Bonded/other':15s}  {ec_other:>10.1f}  {er_other:>10.1f}  {el_other:>10.1f}  {dg_other:>+10.1f}")
    print(f"  {'Check sum':15s}  {'':>10s}  {'':>10s}  {'':>10s}  {dg_nb+dg_gb+dg_sa+dg_other:>+10.1f} kJ  =  {(dg_nb+dg_gb+dg_sa+dg_other)/4.184:>+8.2f} kcal")
    print(f"  Time: {time.time()-t0:.0f}s", flush=True)

    return {
        "compound": compound,
        "dG_total": round(dg_total, 2), "dG_kcal": round(dg_total/4.184, 2),
        "dE_coulomb_LJ": round(dg_nb, 2),
        "dG_GB_polar": round(dg_gb, 2),
        "dG_SA_implied": round(dg_sa, 2),
        "dE_bonded": round(dg_other, 2),
        "E_complex_total": round(ec_total, 1),
        "E_receptor_total": round(er_total, 1),
        "E_ligand_total": round(el_total, 1),
        "n_atoms": natoms, "n_lig_atoms": len(lig_atoms),
    }


def main():
    targets = sys.argv[1:] if len(sys.argv) > 1 else [
        "NAFTAZONE", "BIOTIN", "ESKETAMINE", "FUROSEMIDE",
        "GABAPENTIN_ENACARBIL", "HYDROMORPHONE", "PHENOBARBITAL", "RIBOFLAVIN"]

    all_results = {}
    for c in targets:
        try:
            r = run(c.upper())
            if r: all_results[c.upper()] = r
        except Exception as e:
            print(f"  ERROR {c}: {e}")
            import traceback; traceback.print_exc()

    # Summary table
    if all_results:
        print(f"\n{'='*72}")
        print(f"  ENERGY DECOMPOSITION SUMMARY (kJ/mol)")
        print(f"{'='*72}")
        print(f"  {'Compound':28s}  {'dG_total':>9s}  {'dE_CLJ':>9s}  {'dG_GB':>9s}  {'dG_SA':>9s}  {'dBond':>8s}")
        print(f"  {'-'*68}")
        for c, r in sorted(all_results.items()):
            print(f"  {c:28s}  {r['dG_total']:>+8.1f}  {r['dE_coulomb_LJ']:>+8.1f}  "
                  f"{r['dG_GB_polar']:>+8.1f}  {r['dG_SA_implied']:>+8.1f}  {r['dE_bonded']:>+8.1f}")

        # Save
        out = os.path.join(ANALYSIS_DIR, "energy_decomposition.json")
        with open(out, "w") as f:
            json.dump(all_results, f, indent=2)
        print(f"\n  Saved: {out}")
        print(f"  Note: dG_SA_implied = dG_total - dE_CLJ - dG_GB - dBond (since OBC2 combines GB+SA)")


if __name__ == "__main__":
    main()
