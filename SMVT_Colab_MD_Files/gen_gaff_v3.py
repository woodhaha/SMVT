"""
GAFF 2.11 template generator v3 — RDKit atom typing + MMFF94 charges.

Combines the best of both approaches:
- v1 (gen_gaff_templates.py): correct RDKit aromaticity/hybridization → GAFF types
- v2 (gen_gaff_proper.py): MMFF94 charges + clean XML with SMILES header

No openff-toolkit dependency — works immediately.
MMFF94 charges are significantly better than Gasteiger for MD simulations.

Usage:
    python gen_gaff_v3.py [output_dir]
"""

import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem


# ═══════════════════════════════════════════════════════════════════════
# GAFF 2.11 Atom Type Assignment (RDKit-based)
# ═══════════════════════════════════════════════════════════════════════

def assign_gaff_type(atom, mol):
    """
    Assign GAFF 2.11 atom type based on RDKit atom properties.

    Priority: aromaticity > carbonyl detection > hybridization > defaults.
    """
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHybridization()
    degree = atom.GetDegree()
    is_aromatic = atom.GetIsAromatic()
    formal_charge = atom.GetFormalCharge()

    if atomic_num == 1:  # Hydrogen
        if degree >= 1:
            neighbor = atom.GetNeighbors()[0]
            n_elem = neighbor.GetAtomicNum()
            if n_elem == 6:
                if neighbor.GetIsAromatic() or neighbor.GetHybridization() in (
                    Chem.HybridizationType.SP2, Chem.HybridizationType.SP
                ):
                    return "ha"
                return "hc"
            elif n_elem == 7:
                return "hn"
            elif n_elem == 8:
                return "ho"
            elif n_elem == 16:
                return "hs"
        return "hc"

    elif atomic_num == 6:  # Carbon
        if is_aromatic:
            return "ca"
        # Carbonyl C: sp2 carbon double-bonded to O
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
                if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return "c"  # carbonyl
        if hyb == Chem.HybridizationType.SP:
            return "c1"
        if hyb == Chem.HybridizationType.SP2:
            return "c2"
        return "c3"

    elif atomic_num == 7:  # Nitrogen
        if is_aromatic:
            return "nb"
        # Amide N: N directly bonded to carbonyl C
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 6:
                for nn in n.GetNeighbors():
                    if nn.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            return "n"  # amide
        if hyb == Chem.HybridizationType.SP:
            return "n1"
        if hyb == Chem.HybridizationType.SP2:
            return "n2"
        if formal_charge > 0:
            return "n4"
        return "n3"

    elif atomic_num == 8:  # Oxygen
        # Carbonyl O: double-bonded to C
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
                if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return "o"  # carbonyl oxygen
        # Hydroxyl O: bonded to H
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 1:
                return "oh"
        return "os"  # ether/ester

    elif atomic_num == 9:   return "f"
    elif atomic_num == 15:  return "p5"
    elif atomic_num == 16:  # Sulfur
        if is_aromatic: return "ss"
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 1: return "sh"
        return "s"
    elif atomic_num == 17:  return "cl"
    elif atomic_num == 35:  return "br"
    elif atomic_num == 53:  return "i"

    return "X"


# ═══════════════════════════════════════════════════════════════════════
# Template Generator
# ═══════════════════════════════════════════════════════════════════════

def generate_template(smiles: str, residue_name: str = "LIG") -> str:
    """
    Generate OpenMM ForceField XML residue template.

    Atom types: GAFF 2.11 via RDKit
    Charges: MMFF94 (superior to Gasteiger for organic drug molecules)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)

    # 3D conformer for charge assignment
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    # MMFF94 charges
    mmff = AllChem.MMFFGetMoleculeProperties(mol)
    if mmff is None:
        print(f"  ⚠ MMFF94 failed, falling back to Gasteiger")
        AllChem.ComputeGasteigerCharges(mol)

    atoms = list(mol.GetAtoms())
    n_atoms = len(atoms)

    # Build XML
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<ForceField>',
        f'  <!-- GAFF 2.11 template: {residue_name} -->',
        f'  <!-- SMILES: {smiles} -->',
        f'  <!-- Generator: gen_gaff_v3.py (RDKit + MMFF94) -->',
        f'  <!-- Charge model: MMFF94 -->',
        f'  <Residues>',
        f'    <Residue name="{residue_name}">',
    ]

    type_counts = {}
    for i, atom in enumerate(atoms):
        gaff_type = assign_gaff_type(atom, mol)
        elem = atom.GetSymbol()
        type_counts[gaff_type] = type_counts.get(gaff_type, 0) + 1

        if mmff:
            charge = mmff.GetMMFFPartialCharge(i)
        elif atom.HasProp('_GasteigerCharge'):
            charge = float(atom.GetProp('_GasteigerCharge'))
        else:
            charge = 0.0

        atom_name = f"{elem.upper()}{i + 1}"
        lines.append(
            f'      <Atom name="{atom_name}" type="{gaff_type}" charge="{charge:.6f}"/>'
        )

    # Bonds (each bond once)
    seen = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (i, j) in seen or (j, i) in seen:
            continue
        seen.add((i, j))
        lines.append(f'      <Bond from="{i}" to="{j}"/>')

    lines.append('    </Residue>')
    lines.append('  </Residues>')
    lines.append('</ForceField>')

    return "\n".join(lines), n_atoms, type_counts


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    TOP3 = [
        ("NAFTAZONE", "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N"),
        ("PHENOBARBITAL", "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2"),
        ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"),
    ]

    output_dir = sys.argv[1] if len(sys.argv) > 1 else "ligand_params_v3"
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("GAFF Template Generator v3 — RDKit + MMFF94")
    print("=" * 60)
    print("Key: ca=aromatic C, c=carbonyl C, c3=sp3 C, c2=sp2 C")
    print("     n=amide N, n3=amine N, nb=aromatic N")
    print("     o=carbonyl O, oh=hydroxyl O, os=ether O")
    print("     cl=Cl, ha=H-aromatic, hc=H-sp3C, hn=H-N, ho=H-O")
    print()

    for name, smiles in TOP3:
        print(f"[{name}] {smiles}")
        xml, n_atoms, type_counts = generate_template(smiles, residue_name="LIG")

        xml_path = os.path.join(output_dir, f"{name}_template.xml")
        with open(xml_path, "w") as f:
            f.write(xml)

        n_bonds = xml.count('<Bond ')
        print(f"  Atoms: {n_atoms} | Bonds: {n_bonds}")
        print(f"  Types: {dict(sorted(type_counts.items()))}")
        print()

    print(f"Done → {output_dir}/")
    print("Next: upload to Colab → colab exec notebook → verify NVT no NaN")
