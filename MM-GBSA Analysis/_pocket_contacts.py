"""
Per-residue ligand contact analysis: identify pocket residues + measure distances.
Simple geometry-based: no charge perturbation issues.
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

def analyze_pocket(compound):
    """Identify pocket residues and compute interaction features."""
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

    # Restrained minimization (same as original)
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
    try: mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=40, tolerance=100.0)
    except: pass
    pos_min = ctx.getState(getPositions=True).getPositions()
    system.removeForce(system.getNumForces() - 1)
    ctx.reinitialize(preserveState=True)
    ctx.setPositions(pos_min)

    # Get coords in nm
    coords = np.array([list(pos_min[a.index].value_in_unit(unit.nanometers)) for a in pdb.topology.atoms()])
    lig_coords = coords[list(lig_atoms)]
    lig_center = lig_coords.mean(axis=0)

    # Analyze each residue within 6A of ligand
    residues = {}
    for a in pdb.topology.atoms():
        if a.index in lig_atoms: continue
        ri = a.residue.index
        pos = coords[a.index]
        min_dist = min(np.linalg.norm(pos - lc) for lc in lig_coords)
        if min_dist < 0.6:  # 6A cutoff
            if ri not in residues:
                residues[ri] = {
                    'name': f'{a.residue.name}{a.residue.id}',
                    'chain': a.residue.chain.index,
                    'min_dist': min_dist,
                    'has_hbond': False,
                }
            residues[ri]['min_dist'] = min(residues[ri]['min_dist'], min_dist)

    # Identify HBond donors/acceptors (N/O within 3.5A of ligand N/O)
    for a in pdb.topology.atoms():
        if a.index in lig_atoms:
            if a.element.symbol in ('N', 'O'):
                for b in pdb.topology.atoms():
                    if b.index in lig_atoms: continue
                    if b.element.symbol in ('N', 'O'):
                        d = np.linalg.norm(coords[a.index] - coords[b.index])
                        if d < 0.35:  # 3.5A
                            ri = b.residue.index
                            if ri in residues:
                                residues[ri]['has_hbond'] = True

    # Format output
    pocket_list = []
    for ri, r in sorted(residues.items(), key=lambda x: x[1]['min_dist']):
        pocket_list.append({
            'residue': r['name'],
            'dist_A': round(r['min_dist']*10, 2),
            'hbond': r['has_hbond'],
        })

    natoms = pdb.topology.getNumAtoms()
    print(f'  {compound}: {len(pocket_list)} contact residues (6A), '
          f'{sum(1 for r in pocket_list if r["hbond"])} H-bond capable, '
          f'time={time.time()-t0:.0f}s', flush=True)
    return {
        'compound': compound,
        'n_atoms': natoms,
        'n_lig_atoms': len(lig_atoms),
        'pocket_residues': pocket_list,
    }

def main():
    compounds = sys.argv[1:] if len(sys.argv) > 1 else [
        'NAFTAZONE','BIOTIN','ESKETAMINE','FUROSEMIDE',
        'GABAPENTIN_ENACARBIL','HYDROMORPHONE','PHENOBARBITAL','RIBOFLAVIN']
    all_results = {}
    for c in compounds:
        try: all_results[c.upper()] = analyze_pocket(c.upper())
        except Exception as e:
            print(f'ERROR {c}: {e}')
            import traceback; traceback.print_exc()

    out = os.path.join(ANALYSIS, 'pocket_contacts.json')
    with open(out, 'w') as f:
        json.dump(all_results, f, indent=2)

    # Summary: consensus pocket residues
    print(f'\n{"POCKET CONTACT CONSENSUS":^60s}')
    print(f'{"Residue":>8s}  ', end='')
    for c in compounds:
        print(f'{c[:5]:>6s}', end='')
    print(f'  {"Count":>5s}  {"H-bonds":>7s}')

    # All residue names across compounds
    all_res = sorted(set(r['residue'] for res in all_results.values() for r in res['pocket_residues']))
    consensus = {}
    for res_name in all_res:
        count = 0
        hbonds = 0
        for c in compounds:
            if c not in all_results: continue
            for r in all_results[c]['pocket_residues']:
                if r['residue'] == res_name:
                    count += 1
                    if r['hbond']: hbonds += 1
                    break
        consensus[res_name] = {'count': count, 'hbonds': hbonds}

    for res_name in sorted(consensus.keys(), key=lambda x: -consensus[x]['count']):
        c = consensus[res_name]
        if c['count'] >= 4:
            print(f'{res_name:>8s}  ', end='')
            for comp in compounds:
                if comp not in all_results: print(f'{"":>6s}', end='')
                else:
                    found = [r for r in all_results[comp]['pocket_residues'] if r['residue'] == res_name]
                    if found: print(f'{found[0]["dist_A"]:>5.1f}A', end=' ')
                    else: print(f'{"":>6s}', end='')
            print(f'  {c["count"]:>4d}/8  {c["hbonds"]:>3d}/8')

    print(f'\nSaved: {out}')
    print(f'Note: pocket = residues with any atom within 6A of ligand')

if __name__ == '__main__':
    main()
