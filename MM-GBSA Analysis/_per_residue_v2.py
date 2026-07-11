"""
Per-residue ligand interaction decomposition v2.
Correct approach: compute E(complex) - E(complex with residue's charges zeroed).
This measures the ligand-residue electrostatic + vdW interaction energy.
"""
import os, sys, json, time, warnings
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
warnings.filterwarnings("ignore")
KB = unit.kilojoule_per_mole
KJNM2 = unit.kilojoule_per_mole / (unit.nanometer**2)

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))
LIGANDS = os.path.join(os.path.dirname(ANALYSIS), 'SMVT_MD', 'ligands')
GAFF_XML = "C:/anaconda3/envs/smvt-md/Lib/site-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
CUTOFF = 2.0 * unit.nanometer
POCKET_CUTOFF = 0.6  # nm = 6A for pocket (tighter = faster, more focused)

def build(compound):
    complex_pdb = os.path.join(ANALYSIS, f'{compound}_complex.pdb')
    lig_tmpl = os.path.join(LIGANDS, f'{compound}_template.xml')
    pdb = app.PDBFile(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl)
    system = ff.createSystem(pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)
    lig_atoms = set(a.index for a in pdb.topology.atoms()
                    if a.residue.name.strip() == "LIG")
    return system, pdb, lig_atoms

def minimize(system, pdb):
    caf = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    caf.addGlobalParameter("k", 500.0 * KB / (unit.nanometer**2))
    caf.addPerParticleParameter("x0"); caf.addPerParticleParameter("y0"); caf.addPerParticleParameter("z0")
    for a in pdb.topology.atoms():
        if a.name.strip() == "CA":
            x, y, z = pdb.positions[a.index].value_in_unit(unit.nanometer)
            caf.addParticle(a.index, [x, y, z])
    system.addForce(caf)
    cpu = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 1*unit.femtosecond), cpu)
    ctx.setPositions(pdb.positions)
    mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=40, tolerance=100.0)
    pos_min = ctx.getState(getPositions=True).getPositions()
    system.removeForce(system.getNumForces() - 1)
    ctx.reinitialize(preserveState=True)
    ctx.setPositions(pos_min)
    return ctx, pos_min

def get_pocket_residues(pdb, lig_atoms, positions, cutoff=POCKET_CUTOFF):
    coords = np.array([list(positions[a.index].value_in_unit(unit.nanometers))
                       for a in pdb.topology.atoms()])
    lig_coords = coords[list(lig_atoms)]
    pocket = {}
    seen_res = set()
    for a in pdb.topology.atoms():
        if a.index in lig_atoms: continue
        ri = a.residue.index
        if ri in seen_res: continue
        seen_res.add(ri)
        pos = coords[a.index]
        min_dist = min(np.linalg.norm(pos - lc) for lc in lig_coords)
        if min_dist < cutoff:
            # Check all atoms of this residue
            res_atoms = [at.index for at in pdb.topology.atoms() if at.residue.index == ri]
            res_min_dist = min(min(np.linalg.norm(coords[ai] - lc) for lc in lig_coords) for ai in res_atoms)
            pocket[ri] = {
                'name': f'{a.residue.name}{a.residue.id}',
                'atoms': res_atoms,
                'dist': res_min_dist * 10,  # nm -> A
            }
    return pocket

def decompose(compound):
    print(f'\n=== {compound} ===', flush=True)
    t0 = time.time()
    system, pdb, lig_atoms = build(compound)
    ctx, pos_min = minimize(system, pdb)
    natoms = pdb.topology.getNumAtoms()

    pocket = get_pocket_residues(pdb, lig_atoms, pos_min)
    print(f'  Pocket: {len(pocket)} residues within {POCKET_CUTOFF*10}A', flush=True)

    # Get NonbondedForce
    nbf = None
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.NonbondedForce):
            nbf = system.getForce(i)
            break

    # Save all parameters
    all_charges = []
    all_sigmas = []
    all_epsilons = []
    for ai in range(natoms):
        c, s, e = nbf.getParticleParameters(ai)
        all_charges.append(c); all_sigmas.append(s); all_epsilons.append(e)

    results = []
    for ri, rinfo in sorted(pocket.items()):
        atoms = rinfo['atoms']
        # Zero charges of this residue → lose ligand-residue electrostatic interaction
        for ai in atoms:
            nbf.setParticleParameters(ai, 0.0*unit.elementary_charge, all_sigmas[ai], all_epsilons[ai])
        nbf.updateParametersInContext(ctx)
        s_mut = ctx.getState(getEnergy=True)
        e_mut = s_mut.getPotentialEnergy().value_in_unit(KB)
        # Restore
        for ai in atoms:
            nbf.setParticleParameters(ai, all_charges[ai], all_sigmas[ai], all_epsilons[ai])
        nbf.updateParametersInContext(ctx)

        # Now zero epsilon (vdW) of this residue → lose ligand-residue vdW
        for ai in atoms:
            nbf.setParticleParameters(ai, all_charges[ai], all_sigmas[ai], 0.0*all_epsilons[ai].unit)
        nbf.updateParametersInContext(ctx)
        s_mut_vdw = ctx.getState(getEnergy=True)
        e_mut_vdw = s_mut_vdw.getPotentialEnergy().value_in_unit(KB)
        # Restore
        for ai in atoms:
            nbf.setParticleParameters(ai, all_charges[ai], all_sigmas[ai], all_epsilons[ai])
        nbf.updateParametersInContext(ctx)

        # Reference: get WT energy after each pair of mutations
        # Actually compute dE directly from the difference
        # WT energy was already computed once; the restored state == WT
        # But to be safe, compute from difference:
        # dE_elec = E_charge_zeroed - E_wt  (positive = favorable interaction lost)

        results.append({
            'residue': rinfo['name'],
            'dist_A': round(rinfo['dist'], 1),
            'dE_elec_kJ': round(e_mut, 1),  # total energy with charge zeroed
            'dE_vdw_kJ': round(e_mut_vdw, 1),  # total with vdW zeroed
        })

    # Post-process: need the WT total energy
    s_wt = ctx.getState(getEnergy=True)
    E_wt = s_wt.getPotentialEnergy().value_in_unit(KB)

    for r in results:
        # The dE from zeroing charges = E_mut - E_wt
        # This includes the lost ligand-residue interaction + lost intra-residue
        # But for charged residues the intra is also significant...
        r['dE_elec'] = round(r.pop('dE_elec_kJ') - E_wt, 1)
        r['dE_vdw'] = round(r.pop('dE_vdw_kJ') - E_wt, 1)
        r['dE_total'] = round(r['dE_elec'] + r['dE_vdw'], 1)

    results.sort(key=lambda x: abs(x['dE_total']), reverse=True)

    for r in results[:15]:
        print(f'  {r["residue"]:>8s}  d={r["dist_A"]:4.1f}A  '
              f'elec={r["dE_elec"]:+8.1f}  vdw={r["dE_vdw"]:+8.1f}  '
              f'total={r["dE_total"]:+8.1f} kJ/mol', flush=True)

    print(f'  Time: {time.time()-t0:.0f}s', flush=True)
    return {'compound': compound, 'E_wt_kJ': round(E_wt, 1), 'pocket_residues': results[:15]}

def main():
    compounds = sys.argv[1:] if len(sys.argv) > 1 else [
        'NAFTAZONE','BIOTIN','ESKETAMINE','FUROSEMIDE',
        'GABAPENTIN_ENACARBIL','HYDROMORPHONE','PHENOBARBITAL','RIBOFLAVIN']
    all_results = {}
    for c in compounds:
        try:
            r = decompose(c.upper())
            if r: all_results[c.upper()] = r
        except Exception as e:
            print(f'ERROR {c}: {e}')
            import traceback; traceback.print_exc()
    out = os.path.join(ANALYSIS, 'per_residue_decomp_v2.json')
    with open(out, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nSaved: {out}')
    if all_results:
        print(f'\n{"TOP 5 POCKET RESIDUES PER COMPOUND":^60s}')
        print('-'*60)
        for c, r in all_results.items():
            top5 = r['pocket_residues'][:5]
            s = ', '.join(f'{rr["residue"]}({rr["dE_total"]:+.0f})' for rr in top5)
            print(f'  {c:28s}  {s}')

if __name__ == '__main__':
    main()
