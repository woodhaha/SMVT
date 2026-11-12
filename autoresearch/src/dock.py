"""Reusable SMVT docking utility — RDKit 3D + meeko prep + AutoDock Vina binary.

Used by all autoresearch experiments. Single function `dock_smiles` + batch helper.
Receptor/config fixed to the SMVT binding pocket (cryo-EM 26va-aligned AF2 pocket).
"""
import os, subprocess, tempfile, json, re
from rdkit import Chem
from rdkit.Chem import AllChem

VINA = r"D:\Researching\SMVT\03_Analysis\tools\vina.exe"
RECEPTOR = r"D:\Researching\SMVT\03_Analysis\docking\SMVT_receptor.pdbqt"  # meeko-prepared AF2 receptor (vina-compatible)
MK_PREP_LIG = r"C:\anaconda3\Scripts\mk_prepare_ligand.exe"

# SMVT pocket (from 03_Analysis/data/docking_config.json)
CENTER = (5.03, 1.5, -2.44)
SIZE = (25.0, 25.0, 25.0)
EXHAUSTIVENESS = 16
NUM_MODES = 9
ENERGY_RANGE = 3


def _to_pdbqt(smiles: str, workdir: str) -> str | None:
    """SMILES -> meeko PDBQT via the proven mk_prepare_ligand.exe CLI.

    RDKit (ETKDG + MMFF) -> SDF -> mk_prepare_ligand (pH 7.4, GAFF) -> PDBQT.
    """
    os.makedirs(workdir, exist_ok=True)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    ok = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    if ok != 0:
        ok = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=200, useRandomCoords=True)
    if ok != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    sdf = os.path.join(workdir, "ligand.sdf")
    Chem.MolToMolFile(mol, sdf)
    pdbqt = os.path.join(workdir, "ligand.pdbqt")
    r = subprocess.run([MK_PREP_LIG, "-i", sdf, "-o", pdbqt],
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0 or not os.path.exists(pdbqt):
        return None
    return pdbqt


def _parse_vina(stdout: str) -> list[float]:
    """Parse Vina stdout affinity table -> list of mode energies (kcal/mol).
    Matches only the 'mode | affinity' rows (e.g. '   1    -6.881    0    0')."""
    out = []
    for line in stdout.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s", line)
        if m:
            out.append(float(m.group(2)))
    return out


def dock_smiles(smiles: str, name: str, workdir: str = None, seed: int = None) -> dict:
    """Dock one SMILES into SMVT pocket. Returns dict(name, best_affinity, all_affinities, ok, error).
    `seed` fixes the Vina random seed (default: random, for reproducibility set 0-1000000)."""
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix="smvt_dock_")
    os.makedirs(workdir, exist_ok=True)
    ligand = _to_pdbqt(smiles, workdir)
    if ligand is None:
        return {"name": name, "best_affinity": 99.0, "all_affinities": [],
                "ok": False, "error": "ligand prep failed"}
    out_pdbqt = os.path.join(workdir, "docked.pdbqt")
    cmd = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
           "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
           "--out", out_pdbqt, "--exhaustiveness", str(EXHAUSTIVENESS),
           "--num_modes", str(NUM_MODES), "--energy_range", str(ENERGY_RANGE)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return {"name": name, "best_affinity": 99.0, "all_affinities": [],
                "ok": False, "error": "timeout"}
    affs = _parse_vina(proc.stdout)
    if not affs:
        return {"name": name, "best_affinity": 99.0, "all_affinities": [],
                "ok": False, "error": "no affinity parsed: " + proc.stderr[:120]}
    return {"name": name, "best_affinity": min(affs), "all_affinities": affs,
            "ok": True, "error": ""}


def dock_batch(ligands: list[tuple[str, str]], workdir: str) -> list[dict]:
    """Dock a list of (name, SMILES). Returns list of result dicts sorted by best_affinity."""
    results = []
    for name, smi in ligands:
        r = dock_smiles(smi, name, os.path.join(workdir, name))
        results.append(r)
    results.sort(key=lambda r: r["best_affinity"])
    return results


if __name__ == "__main__":
    # self-check: biotin should dock ~ -6.8 (reproduce baseline)
    import sys
    smi = sys.argv[1] if len(sys.argv) > 1 else "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"
    name = sys.argv[2] if len(sys.argv) > 2 else "biotin"
    r = dock_smiles(smi, name)
    print(json.dumps(r, indent=2))
