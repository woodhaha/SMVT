"""
GAFF 2.11 template generator — PROPER atom types via openff-toolkit SMIRNOFF vdW labeling.

Uses the SMIRNOFF force field's vdW handler to assign per-atom chemical types
from SMIRKS patterns, cross-walked to GAFF 2.11 atom types.

Charges: RDKit MMFF94 (better than Gasteiger for organic ligands).
"""

import os, sys, re
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField


# ── SMIRKS pattern → GAFF 2.11 type ──
# Derived from SMIRNOFF vdW handler SMIRKS. Each pattern encodes the
# element, hybridization, and chemical environment.

def _smirks_to_gaff(smirks: str, aromatic: bool = False,
                    is_carbonyl: bool = False, is_amide_n: bool = False) -> str:
    """
    Parse a SMIRNOFF vdW SMIRKS pattern and return the GAFF 2.11 atom type.

    Uses molecule-level aromaticity/is_carbonyl as fallback when SMIRKS
    patterns are ambiguous (e.g. generic [#6:1] can't distinguish
    aromatic 'ca' from carbonyl 'c').

    Examples:
      [#6X4:1]           → c3   (sp3 carbon)
      [#6X3:1]           → c2   (sp2 carbon, alkene)
      [#6:1]  + aromatic → ca   (aromatic)
      [#6:1]  + carbonyl → c    (carbonyl)
      [#6a:1]            → ca   (aromatic carbon, explicit)
      [#8:1]             → o    (carbonyl oxygen)
      [#8X2:1]           → os   (ether/ester oxygen)
      [#7:1]  + amide    → n    (amide nitrogen)
      [#7X3:1]           → n3   (amine nitrogen)
      [#1:1]-[#6X4]      → hc   (H on sp3 carbon)
      [#1:1]-[#6X3]      → ha   (H on sp2/aromatic carbon)
      [#17:1]            → cl   (chlorine)
    """
    # Extract element number from [#N...]
    elem_match = re.match(r'\[#(\d+)', smirks)
    if not elem_match:
        return "X"
    atomic_num = int(elem_match.group(1))

    # Check for explicit aromatic marker [#6a:1]
    is_aromatic_explicit = bool(re.search(r'\[#\d+a', smirks))

    # Extract coordination number Xn
    coord_match = re.search(r'X(\d)', smirks)
    coord = int(coord_match.group(1)) if coord_match else None

    # ── Element-specific mapping ──
    if atomic_num == 1:  # Hydrogen
        if '[#6X4]' in smirks: return "hc"
        if '[#6X3]' in smirks or '[#6a]' in smirks: return "ha"
        if '[#7]' in smirks: return "hn"
        if '[#8]' in smirks: return "ho"
        if '[#16]' in smirks: return "hs"
        return "hc"

    elif atomic_num == 6:  # Carbon
        if is_aromatic_explicit or aromatic:
            return "ca"
        if coord == 4:
            return "c3"
        if coord == 3:
            return "c2"
        if coord == 1:
            return "c1"
        # coord is None: distinguish carbonyl (c) vs aromatic (ca)
        if is_carbonyl:
            return "c"
        if aromatic:
            return "ca"
        # Truly ambiguous — default to sp3 as safest fallback
        return "c3"

    elif atomic_num == 7:  # Nitrogen
        if is_aromatic_explicit or aromatic:
            return "nb"
        if coord == 4:
            return "n4"  # ammonium
        if coord == 3:
            return "n3"
        if coord == 2:
            return "n2"
        if coord == 1:
            return "n1"
        if is_amide_n:
            return "n"
        return "n3"  # default: amine (most common for drug ligands)

    elif atomic_num == 8:  # Oxygen
        if coord == 2:
            return "os"  # ether/ester
        if is_carbonyl:
            return "o"
        # Default: hydroxyl (oh) — safer than generic o
        return "oh"

    elif atomic_num == 9:   return "f"
    elif atomic_num == 15:  return "p5"
    elif atomic_num == 16:  # Sulfur
        if coord == 4: return "s4"
        if coord == 6: return "s6"
        return "s"
    elif atomic_num == 17:  return "cl"
    elif atomic_num == 35:  return "br"
    elif atomic_num == 53:  return "i"

    return "X"


def _detect_carbonyl(rd_mol) -> set:
    """Return set of atom indices that are carbonyl C or O (C=O double bond)."""
    carbonyl_atoms = set()
    for bond in rd_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8:
                carbonyl_atoms.add(bond.GetBeginAtomIdx())  # C
                carbonyl_atoms.add(bond.GetEndAtomIdx())    # O
            elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8:
                carbonyl_atoms.add(bond.GetEndAtomIdx())    # C
                carbonyl_atoms.add(bond.GetBeginAtomIdx())  # O
    return carbonyl_atoms


def _detect_amide_n(rd_mol) -> set:
    """Return set of N atom indices directly bonded to a carbonyl C."""
    amide_n = set()
    for bond in rd_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            for (n_idx, c_candidate) in [(bond.GetBeginAtomIdx(), a2),
                                          (bond.GetEndAtomIdx(), a1)]:
                atom = rd_mol.GetAtomWithIdx(n_idx)
                if atom.GetAtomicNum() != 7:
                    continue
                if c_candidate.GetAtomicNum() != 6:
                    continue
                # Check if this C is a carbonyl (double-bonded to O)
                for nb in c_candidate.GetNeighbors():
                    if nb.GetAtomicNum() == 8:
                        nb_bond = rd_mol.GetBondBetweenAtoms(
                            c_candidate.GetIdx(), nb.GetIdx())
                        if nb_bond and nb_bond.GetBondType() == Chem.BondType.DOUBLE:
                            amide_n.add(n_idx)
                            break
    return amide_n


# ── Main ──

def generate_template(smiles: str, residue_name: str = "LIG",
                      charge_method: str = "mmff94") -> str:
    """
    Generate an OpenMM ForceField XML residue template with correct GAFF types.

    Uses SMIRNOFF vdW SMIRKS patterns cross-walked to GAFF 2.11,
    with RDKit aromaticity/carbonyl info as fallback for ambiguous patterns.
    """
    # RDKit: primary structure analysis (aromaticity, hybridization, carbonyl)
    rd_mol = Chem.MolFromSmiles(smiles)
    if rd_mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    rd_mol = Chem.AddHs(rd_mol)
    AllChem.EmbedMolecule(rd_mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(rd_mol)
    mmff = AllChem.MMFFGetMoleculeProperties(rd_mol)

    # Detect chemical features
    carbonyl_set = _detect_carbonyl(rd_mol)
    amide_n_set = _detect_amide_n(rd_mol)

    # Per-atom aromaticity from RDKit (heavy atoms only)
    rd_aromatic = [atom.GetIsAromatic() for atom in rd_mol.GetAtoms()]

    # Load with openff-toolkit for SMIRNOFF vdW labeling
    off_mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    off_mol.generate_conformers(n_conformers=1)
    off_mol.name = residue_name

    ff = ForceField("openff-2.1.0.offxml")
    labels = ff.label_molecules(off_mol.to_topology())[0]

    # Extract per-atom SMIRKS pattern from vdW handler
    n_atoms = off_mol.n_atoms
    smirks_list = _extract_smirks(labels, n_atoms)

    # Verify atom count match between RDKit and openff
    n_rd = rd_mol.GetNumAtoms()
    if n_rd != n_atoms:
        print(f"  ⚠ Atom count mismatch: RDKit={n_rd}, openff={n_atoms}")

    # Atom names (RDKit order)
    counters = {}
    names = []
    for atom in rd_mol.GetAtoms():
        elem = atom.GetSymbol()
        counters[elem] = counters.get(elem, 0) + 1
        names.append(f"{elem.upper()}{counters[elem]}")

    # Build XML
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<ForceField>',
        f'  <!-- GAFF 2.11 template: {residue_name} -->',
        f'  <!-- SMILES: {smiles} -->',
        f'  <!-- Charges: {charge_method.upper()} -->',
        f'  <Residues>',
        f'    <Residue name="{residue_name}">',
    ]

    for i in range(n_atoms):
        aromatic_i = rd_aromatic[i] if i < n_rd else False
        is_carbonyl = i in carbonyl_set
        is_amide = i in amide_n_set

        gaff = _smirks_to_gaff(smirks_list[i],
                               aromatic=aromatic_i,
                               is_carbonyl=is_carbonyl,
                               is_amide_n=is_amide)
        chg = mmff.GetMMFFPartialCharge(i) if mmff else 0.0
        lines.append(
            f'      <Atom name="{names[i]}" type="{gaff}" charge="{chg:.6f}"/>'
        )

    # Bonds
    seen = set()
    for bond in rd_mol.GetBonds():
        u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (u, v) in seen or (v, u) in seen:
            continue
        seen.add((u, v))
        lines.append(f'      <Bond from="{u}" to="{v}"/>')

    lines.append('    </Residue>')
    lines.append('  </Residues>')
    lines.append('</ForceField>')

    return "\n".join(lines), smirks_list


def _extract_smirks(labels, n_atoms: int) -> list:
    """
    Extract per-atom SMIRKS patterns from the vdW handler labels.

    The vdW handler has entries like:
      (0,): <vdWType with smirks: [#6X4:1]  epsilon: ...  id: n16  ...>

    We parse the SMIRKS string for each atom index.
    """
    smirks_by_atom = [None] * n_atoms

    if "vdW" in labels:
        for (atom_idx,), vdw_type in labels["vdW"].items():
            smirks_str = getattr(vdw_type, "smirks", "")
            if smirks_str:
                smirks_by_atom[atom_idx] = smirks_str

    # Fill any gaps with element-based fallbacks
    for i in range(n_atoms):
        if smirks_by_atom[i] is None:
            smirks_by_atom[i] = f"[#X:1]"  # generic fallback

    return smirks_by_atom


# ═══ CLI ═══
if __name__ == "__main__":
    TOP3 = [
        ("NAFTAZONE", "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N"),
        ("PHENOBARBITAL", "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2"),
        ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"),
    ]

    output_dir = sys.argv[1] if len(sys.argv) > 1 else "ligand_params_v3"
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("GAFF Template Generator v3 — RDKit aromaticity + SMIRNOFF vdW")
    print("=" * 60)
    print("Fix: aromatic C → ca (not c), amide N → n (not n3)")
    print("Uses RDKit aromaticity/carbonyl as fallback for ambiguous SMIRKS")
    print()

    for name, smiles in TOP3:
        print(f"[{name}] {smiles}")
        xml, smirks_list = generate_template(smiles, residue_name="LIG")

        xml_path = os.path.join(output_dir, f"{name}_template.xml")
        with open(xml_path, "w") as f:
            f.write(xml)

        # Parse XML back for type summary
        type_counts = {}
        for line in xml.split("\n"):
            if 'type="' in line and '<Atom ' in line:
                m = re.search(r'type="(\S+)"', line)
                if m:
                    t = m.group(1)
                    type_counts[t] = type_counts.get(t, 0) + 1

        print(f"  Atoms: {len(smirks_list)} | Bonds: {xml.count('<Bond ')}")
        print(f"  Types: {dict(sorted(type_counts.items()))}")
        print()

    print(f"Done → {output_dir}/")
