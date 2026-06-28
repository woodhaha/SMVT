# GAFF atom type assigner using ONLY RDKit — no openff-toolkit needed
# Maps RDKit atom properties → GAFF atom types for OpenMM residue templates

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# GAFF atom type mapping based on element, atomic number, hybridization, neighbors
def assign_gaff_type(atom, mol):
    """Assign GAFF atom type to an RDKit atom. Returns (type, element)."""
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHybridization()
    degree = atom.GetDegree()
    total_valence = atom.GetTotalValence()
    is_in_ring = atom.IsInRing()
    is_aromatic = atom.GetIsAromatic()
    formal_charge = atom.GetFormalCharge()
    elem = atom.GetSymbol()

    if atomic_num == 1:  # Hydrogen
        # Check what it's bonded to
        if degree >= 1:
            neighbor = atom.GetNeighbors()[0]
            n_elem = neighbor.GetAtomicNum()
            if n_elem == 6:  # H-C
                if neighbor.GetIsAromatic():
                    return ("ha", "H")
                elif neighbor.GetHybridization() == Chem.HybridizationType.SP2:
                    return ("ha", "H")
                elif neighbor.GetHybridization() == Chem.HybridizationType.SP:
                    return ("ha", "H")
                else:
                    return ("hc", "H")
            elif n_elem == 7:  # H-N
                return ("hn", "H")
            elif n_elem == 8:  # H-O
                return ("ho", "H")
            elif n_elem == 16:  # H-S
                return ("hs", "H")
        return ("hc", "H")

    elif atomic_num == 6:  # Carbon
        if is_aromatic:
            return ("ca", "C")
        if hyb == Chem.HybridizationType.SP:
            return ("c1", "C")
        if hyb == Chem.HybridizationType.SP2:
            # Check if carbonyl (C=O)
            for n in atom.GetNeighbors():
                if n.GetAtomicNum() == 8 and n.GetHybridization() == Chem.HybridizationType.SP2:
                    return ("c", "C")  # carbonyl
            return ("c2", "C")
        return ("c3", "C")

    elif atomic_num == 7:  # Nitrogen
        if is_aromatic:
            return ("nb", "N")
        if hyb == Chem.HybridizationType.SP2:
            return ("n", "N")
        if hyb == Chem.HybridizationType.SP:
            return ("n1", "N")
        if formal_charge > 0:
            return ("n4", "N")
        return ("n3", "N")

    elif atomic_num == 8:  # Oxygen
        if is_aromatic:
            return ("os", "O")
        # Check if carbonyl (double bonded to C)
        is_carbonyl = False
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            return ("o", "O")
        # Check if hydroxyl (bonded to H)
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 1:
                return ("oh", "O")
        return ("os", "O")

    elif atomic_num == 16:  # Sulfur
        if is_aromatic:
            return ("ss", "S")
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 1:
                return ("sh", "S")
        return ("s", "S")

    elif atomic_num == 9:  # Fluorine
        return ("f", "F")

    elif atomic_num == 17:  # Chlorine
        return ("cl", "Cl")

    elif atomic_num == 35:  # Bromine
        return ("br", "Br")

    elif atomic_num == 15:  # Phosphorus
        return ("p5", "P")

    return ("X", elem)


def generate_residue_template_xml(smiles, residue_name="LIG"):
    """
    Generate an OpenMM ForceField XML residue template for a ligand.
    Uses RDKit for structure + heuristic GAFF atom typing + Gasteiger charges.
    Returns XML string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    # Compute Gasteiger partial charges
    AllChem.ComputeGasteigerCharges(mol)
    atoms = mol.GetAtoms()

    # Build XML
    lines = []
    lines.append('<ForceField>')
    lines.append(f' <Residues>')
    lines.append(f'  <Residue name="{residue_name}">')

    for atom in atoms:
        gaff_type, elem = assign_gaff_type(atom, mol)
        atom_name = f"{elem.upper()}{atom.GetIdx() + 1}"
        charge = float(atom.GetProp('_GasteigerCharge')) if atom.HasProp('_GasteigerCharge') else 0.0
        lines.append(f'   <Atom name="{atom_name}" type="{gaff_type}" charge="{charge:.6f}"/>')

    # Bonds — one per bond
    seen = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (i, j) in seen or (j, i) in seen:
            continue
        seen.add((i, j))
        lines.append(f'   <Bond from="{i}" to="{j}"/>')

    lines.append(f'  </Residue>')
    lines.append(f' </Residues>')
    lines.append(f'</ForceField>')

    return "\n".join(lines)


def generate_ligand_pdb(smiles, vina_pos, out_path, residue_name="LIG"):
    """Generate ligand PDB from SMILES with Vina docking pose coordinates."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    # Set positions from Vina pose (Angstrom → Angstrom, same units in RDKit PDB)
    n_atoms = mol.GetNumAtoms()
    for i in range(min(n_atoms, len(vina_pos))):
        conf.SetAtomPosition(i, (vina_pos[i][0], vina_pos[i][1], vina_pos[i][2]))

    # Write PDB and rename residue
    pdb_str = Chem.MolToPDBBlock(mol)
    pdb_str = pdb_str.replace("UNL", f"{residue_name:>3}")
    with open(out_path, "w") as f:
        f.write(pdb_str)
    return out_path


# ═══ Main: generate for all 3 ligands ═══
if __name__ == "__main__":
    import os, sys

    TOP3 = [
        ("NAFTAZONE", "C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N"),
        ("PHENOBARBITAL", "CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2"),
        ("ESKETAMINE", "CNC1(C2=CC=CC=C2Cl)CCCCC1=O"),
    ]

    output_dir = sys.argv[1] if len(sys.argv) > 1 else "ligand_params"
    os.makedirs(output_dir, exist_ok=True)

    for name, smiles in TOP3:
        xml = generate_residue_template_xml(smiles, residue_name="LIG")
        xml_path = os.path.join(output_dir, f"{name}_template.xml")
        with open(xml_path, "w") as f:
            f.write(xml)
        print(f"Generated {xml_path}")
        print(f"  SMILES: {smiles}")

    print(f"\nDone! {len(TOP3)} templates in {output_dir}/")
