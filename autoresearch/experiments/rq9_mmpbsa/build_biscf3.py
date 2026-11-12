#!/usr/bin/env python
"""Build COOH-bisCF3 (carboxylate, -1) SMVT membrane system — reuse biotin pipeline.
Mirrors /root/md/smvt_build_ligands.py but for the single COOH-bisCF3 ligand."""
import subprocess, os, collections
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

MD_ENV = "/usr/local/miniconda3/envs/md_env"
os.environ["AMBERHOME"] = MD_ENV
os.environ["PATH"] = f"{MD_ENV}/bin:" + os.environ["PATH"]

BASE = "/root/md/build"
MEM_PACK = "/root/md/manual/mem_pack.pdb"
ZHEN = "/root/md/zhen2026.pdb"
NA = [(27.555, 28.768, 1.41), (25.991, 32.037, 7.32)]

NAME = "BISCF3"
SMILES = "[O-]C(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2"
CHARGE = -1


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def cryo_centroid():
    pts = []
    for l in open(ZHEN):
        if l.startswith("HETATM") and l[17:20].strip() == "BTN":
            pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(pts).mean(0) + np.array([0, 0, -70.95])


def rdkit_pdb(smiles, out):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    ok = AllChem.EmbedMolecule(mol, randomSeed=42)
    if ok != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    block = Chem.MolToPDBBlock(mol)
    lines = []
    for l in block.splitlines():
        if l.startswith(("HETATM", "ATOM")):
            lines.append(l[:17] + "UNL" + l[20:])
        elif l.startswith("END"):
            lines.append(l)
    open(out, "w").write("\n".join(lines) + "\n")
    return mol


def read_pdb_coords(pdb):
    pts, names, elems = [], [], []
    for l in open(pdb):
        if l.startswith(("HETATM", "ATOM")):
            pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            names.append(l[12:16].strip())
            elems.append((l[76:78].strip() or l[12:14].strip()))
    return np.array(pts), names, elems


def write_combined(lig_pdb, shift, out):
    lines = []
    for l in open(MEM_PACK):
        if l.startswith("END"):
            continue
        lines.append(l.rstrip("\n"))
    pts, names, elems = read_pdb_coords(lig_pdb)
    for i, (p, n, e) in enumerate(zip(pts, names, elems)):
        x, y, z = p + shift
        lines.append(f"HETATM{i+1:5d} {n:>4s} UNL   501    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {e:>2s}")
    for j, (x, y, z) in enumerate(NA):
        lines.append(f"HETATM{len(pts)+1+j:5d}  Na+ Na+  {502+j:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          Na")
    lines.append("END")
    open(out, "w").write("\n".join(lines) + "\n")
    return len(pts)


def verify(prmtop, inpcrd, name):
    import openmm as mm, openmm.app as app
    import openmm.unit as unit
    prm = app.AmberPrmtopFile(prmtop)
    crd = app.AmberInpcrdFile(inpcrd)
    cnt = collections.Counter(r.name for r in prm.topology.residues())
    print(f"  [{name}] atoms={prm.topology.getNumAtoms()} UNL={cnt.get('UNL')} Na+={cnt.get('Na+')} Cl-={cnt.get('Cl-')} PA={cnt.get('PA')}", flush=True)
    sys_ = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
    plat = mm.Platform.getPlatformByName("OpenCL")
    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)
    sim = app.Simulation(prm.topology, sys_, integ, plat, {"OpenCLPrecision": "mixed"})
    sim.context.setPositions(crd.positions)
    sim.context.setPeriodicBoxVectors(*crd.boxVectors)
    pe0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    sim.minimizeEnergy(maxIterations=20000, tolerance=1)
    pe1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"  [{name}] PE before={pe0:.1f} after min={pe1:.1f} kcal/mol", flush=True)
    import parmed as pmd
    s = pmd.load_file(prmtop)
    for res in s.residues:
        if res.name == "UNL":
            q = sum(a.charge for a in res.atoms)
            print(f"  [{name}] UNL natoms={len(res.atoms)} charge={q:.3f}", flush=True)
    return pe1


def main():
    cryo = cryo_centroid()
    print(f"cryo biotin centroid: {np.round(cryo, 2)}", flush=True)
    d = f"{BASE}/{NAME}"
    os.makedirs(d, exist_ok=True)
    print(f"\n=== {NAME} (charge {CHARGE}) ===", flush=True)
    lig_pdb = f"{d}/{NAME}_h.pdb"
    mol = rdkit_pdb(SMILES, lig_pdb)
    natoms = mol.GetNumAtoms()
    print(f"  RDKit atoms={natoms} formula={Chem.rdMolDescriptors.CalcMolFormula(mol)}", flush=True)
    rc, o, e = sh(f"cd {d} && antechamber -i {NAME}_h.pdb -fi pdb -o {NAME}.prepin -fo prepi -c bcc -nc {CHARGE} -at gaff2 -s 2 2>&1 | tail -3")
    print(f"  antechamber rc={rc} {o.strip()} {e.strip()}", flush=True)
    rc, o, e = sh(f"cd {d} && parmchk2 -i {NAME}.prepin -f prepi -o {NAME}.frcmod -s gaff2 2>&1 | tail -2")
    print(f"  parmchk2 rc={rc} {o.strip()}", flush=True)
    pts, _, _ = read_pdb_coords(lig_pdb)
    shift = cryo - pts.mean(0)
    combined = f"{d}/combined.pdb"
    n = write_combined(lig_pdb, shift, combined)
    print(f"  combined.pdb ligand atoms={n} shift={np.round(shift, 2)}", flush=True)
    tleap_in = f"{d}/tleap.in"
    with open(tleap_in, "w") as f:
        f.write(f"""source leaprc.protein.ff14SB
source leaprc.lipid21
source leaprc.water.tip3p
source leaprc.gaff2
loadamberprep {NAME}.prepin
loadamberparams {NAME}.frcmod
x = loadpdb combined.pdb
setBox x vdw
addIons x Cl- 0
addIons x Na+ 90 Cl- 90
saveamberparm x system_full.prmtop system_full.inpcrd
savepdb x system_full.pdb
quit
""")
    rc, o, e = sh(f"cd {d} && tleap -f tleap.in > tleap.log 2>&1")
    rc2, o2, _ = sh(f"cd {d} && grep -iE 'error|unknown residue|FATAL' tleap.log | head -20")
    print(f"  tleap rc={rc}; problems: {o2.strip() or '(none)'}", flush=True)
    if not os.path.exists(f"{d}/system_full.prmtop"):
        _, tail, _ = sh(f"cd {d} && tail -15 tleap.log")
        print("  !! prmtop MISSING — tail tleap.log:\n" + tail, flush=True)
        return
    try:
        verify(f"{d}/system_full.prmtop", f"{d}/system_full.inpcrd", NAME)
    except Exception as ex:
        print(f"  !! verify FAILED: {type(ex).__name__}: {ex}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
