"""
Per-residue energy decomposition: ligand interaction with each pocket residue.
Uses OpenMM energy-difference approach (zero-out charges for each residue).
Runs ~1 min per compound.
"""
import os, sys, json, time, warnings
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
warnings.filterwarnings("ignore")
KB = unit.kilojoule_per_mole

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))
LIGANDS = os.path.join(os.path.dirname(ANALYSIS), 'SMVT_MD', 'ligands')
GAFF_XML = "C:/anaconda3/envs/smvt-md/Lib/site-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"
CUTOFF = 2.0 * unit.nanometer
CUTOFF_DIST = 0.8  # nm = 8A for pocket definition

def get_pocket_residues(pdb, lig_atoms, cutoff=CUTOFF_DIST):
    """Find residues with any atom within cutoff of any ligand atom."""
    coords = np.array([list(pdb.positions[a.index].value_in_unit(unit.nanometers))
                       for a in pdb.topology.atoms()])
    lig_coords = coords[list(lig_atoms)]

    pocket = {}
    for a in pdb.topology.atoms():
        if a.index in lig_atoms:
            continue
        res = a.residue
        ri = res.index
        if ri not in pocket:
            pocket[ri] = {
                'name': f'{res.name}{res.id}',
                'chain': res.chain.index,
                'atoms': [],
                'dist': None,
            }
        pos = coords[a.index]
        min_dist = min(np.linalg.norm(pos - lc) for lc in lig_coords)
        pocket[ri]['atoms'].append(a.index)
        if pocket[ri]['dist'] is None or min_dist < pocket[ri]['dist']:
            pocket[ri]['dist'] = min_dist

    # Filter by cutoff
    pocket = {k: v for k, v in pocket.items() if v['dist'] < cutoff}
    return pocket


def decompose_one(compound):
    """Per-residue decomposition for one compound."""
    print(f'\n=== {compound} ===', flush=True)
    t0 = time.time()

    complex_pdb = os.path.join(ANALYSIS, f'{compound}_complex.pdb')
    lig_tmpl = os.path.join(LIGANDS, f'{compound}_template.xml')

    pdb = app.PDBFile(complex_pdb)
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml", GAFF_XML, lig_tmpl)
    system = ff.createSystem(pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=CUTOFF,
        constraints=app.HBonds)

    lig_atoms = set(a.index for a in pdb.topology.atoms()
                    if a.residue.name.strip() == "LIG")

    # Restrained minimization (same as _run_mmgbsa.py)
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

    # Get minimized positions & system (remove restraint force)
    pos_min = ctx.getState(getPositions=True).getPositions()
    system.removeForce(system.getNumForces() - 1)
    ctx.reinitialize(preserveState=True)
    ctx.setPositions(pos_min)

    # Full complex energy
    state_wt = ctx.getState(getEnergy=True)
    e_wt = state_wt.getPotentialEnergy().value_in_unit(KB)
    print(f'  E_wt: {e_wt:.0f} kJ/mol', flush=True)

    # Pocket residues
    pocket = get_pocket_residues(pdb, lig_atoms)
    print(f'  Pocket: {len(pocket)} residues within 8A of ligand', flush=True)

    # Get the NonbondedForce
    nbf = None
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.NonbondedForce):
            nbf = system.getForce(i)
            break

    if nbf is None:
        print('  ERROR: No NonbondedForce found')
        return None

    # Per-residue scan: zero out charges
    results = []
    for ri, rinfo in sorted(pocket.items()):
        res_name = rinfo['name']
        atom_indices = rinfo['atoms']
        dist = rinfo['dist']

        # Save wildtype params
        saved = []
        for ai in atom_indices:
            charge, sigma, epsilon = nbf.getParticleParameters(ai)
            saved.append((charge, sigma, epsilon))
            # Zero the charge
            nbf.setParticleParameters(ai, 0.0*charge.unit, sigma, epsilon)

        nbf.updateParametersInContext(ctx)
        e_mut = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
        dE_elec = e_wt - e_mut

        # Restore charges
        for ai, (charge, sigma, epsilon) in zip(atom_indices, saved):
            nbf.setParticleParameters(ai, charge, sigma, epsilon)
        nbf.updateParametersInContext(ctx)

        # Now zero vdW (epsilon) for the same residue
        saved_vdw = []
        for ai in atom_indices:
            charge, sigma, epsilon = nbf.getParticleParameters(ai)
            saved_vdw.append((charge, sigma, epsilon))
            nbf.setParticleParameters(ai, charge, sigma, 0.0*epsilon.unit)
        nbf.updateParametersInContext(ctx)
        e_mut_vdw = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KB)
        dE_vdw = e_wt - e_mut_vdw

        # Restore
        for ai, (charge, sigma, epsilon) in zip(atom_indices, saved_vdw):
            nbf.setParticleParameters(ai, charge, sigma, epsilon)
        nbf.updateParametersInContext(ctx)

        results.append({
            'residue': res_name,
            'chain': rinfo['chain'],
            'dist_A': dist * 10,
            'dE_elec': round(dE_elec, 1),
            'dE_vdw': round(dE_vdw, 1),
            'dE_total': round(dE_elec + dE_vdw, 1),
        })

        print(f'  {res_name:>8s}  d={dist*10:.1f}A  '
              f'elec={dE_elec:+7.1f}  vdw={dE_vdw:+7.1f}  '
              f'total={dE_elec+dE_vdw:+7.1f} kJ/mol',
              flush=True)

    # Sort by absolute total contribution
    results.sort(key=lambda x: abs(x['dE_total']), reverse=True)
    print(f'  Time: {time.time()-t0:.0f}s', flush=True)

    return {'compound': compound, 'pocket_residues': results[:15]}


def main():
    compounds = sys.argv[1:] if len(sys.argv) > 1 else [
        'BIOTIN', 'ESKETAMINE', 'FUROSEMIDE', 'GABAPENTIN_ENACARBIL',
        'HYDROMORPHONE', 'NAFTAZONE', 'PHENOBARBITAL', 'RIBOFLAVIN']

    all_results = {}
    for c in compounds:
        try:
            r = decompose_one(c.upper())
            if r:
                all_results[c.upper()] = r
        except Exception as e:
            print(f'ERROR {c}: {e}')
            import traceback; traceback.print_exc()

    out = os.path.join(ANALYSIS, 'per_residue_decomp.json')
    with open(out, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nSaved: {out}')

    # Summary: top residues per compound
    print(f'\n{"="*60}')
    print(f'  TOP POCKET RESIDUES')
    print(f'{"="*60}')
    for c, r in all_results.items():
        top5 = r['pocket_residues'][:5]
        top_str = ', '.join(f'{rr["residue"]}({rr["dE_total"]:+.0f})' for rr in top5)
        print(f'  {c:28s}  {top_str}')


if __name__ == '__main__':
    main()
