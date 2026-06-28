#!/usr/bin/env python3
"""
01_prepare_ligands.py
Prepare known SMVT interactors and FDA-approved drug dataset for docking.

Creates:
1. Known SMVT interactors as 3D SDF files (validation set)
2. FDA-approved drug library as 3D SDF files
3. All ligands converted to PDBQT format
"""

import os, sys, json, csv
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors

BASE = Path(r"D:\Researching\SMVT")
DATA_DIR = BASE / "03_Analysis" / "data"
LIGANDS_DIR = DATA_DIR / "ligands"
PDBQT_DIR = DATA_DIR / "pdbqt_ligands"
os.makedirs(LIGANDS_DIR, exist_ok=True)
os.makedirs(PDBQT_DIR, exist_ok=True)

# ============================================================
# PART 1: Known SMVT interactors
# ============================================================
KNOWN_INTERACTORS = {
    "Biotin": {
        "SMILES": "O=C1N[C@@H]2[C@@H](CS1)[C@@H](CCCCC(=O)O)S2",  # PubChem CID 171548
        "PubChem_CID": 171548, "Type": "Natural substrate", "Known_Kd_uM": 2.8,
    },
    "Lipoic_acid": {
        "SMILES": "OC(=O)CCCCC1CCSS1",
        "PubChem_CID": 864, "Type": "Natural substrate", "Known_Kd_uM": 12.0,
    },
    "Pantothenic_acid": {
        "SMILES": "OC(=O)CCNC(=O)C(O)C(C)(C)CO",
        "PubChem_CID": 6613, "Type": "Natural substrate", "Known_Kd_uM": 5.0,
    },
    "Gabapentin_enacarbil": {
        "SMILES": "CCCCCOC(=O)CC1(CC(C)C1)CNC(=O)OCC(C)OC",
        "PubChem_CID": 9852637, "Type": "Prodrug (targeted)", "Known_Kd_uM": 8.5,
    },
    "Indomethacin": {
        "SMILES": "CC1=C(C2=C(N1CC(=O)O)C=CC(=C2)Cl)OC3=CC=C(C=C3)C(=O)O",
        "PubChem_CID": 3715, "Type": "NSAID (known interactor)", "Known_Kd_uM": 45.0,
    },
    "Ibuprofen": {
        "SMILES": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "PubChem_CID": 3672, "Type": "NSAID (known interactor)", "Known_Kd_uM": 120.0,
    },
    "Diclofenac": {
        "SMILES": "OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl",
        "PubChem_CID": 3033, "Type": "NSAID (known interactor)", "Known_Kd_uM": 85.0,
    },
    "Ketoprofen": {
        "SMILES": "CC(C(=O)O)C1=CC=CC=C1C(=O)C2=CC=CC=C2",
        "PubChem_CID": 3825, "Type": "NSAID (known interactor)", "Known_Kd_uM": 95.0,
    },
}

# ============================================================
# PART 2: FDA-approved drug library
# ============================================================
FDA_DRUGS = {
    # Cardiovascular
    "Atorvastatin": "CC(C)C1=CC(=C(C=C1)C2=CC=C(C=C2)C3(CC(CC3=O)O)O)NC4=CC=C(C=C4)F",
    "Lisinopril": "C1CCNC(C1)C(=O)NCCCC2=CC=CC=C2C(=O)O",
    "Metoprolol": "CC(C)NCC(COC1=CC=C(C=C1)CCO)O",
    "Losartan": "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)Cl)C(=O)O",
    "Amlodipine": "CCOC(=O)C1=C(CC(C)C(C)=C(C1c2ccccc2Cl)C(=O)OC)CN",
    "Valsartan": "CCCCC(=O)N(CC1=CC=C(C=C1)C2=CC=CC=C2C3=NNN=N3)C(C(=O)O)C(C)C",
    "Simvastatin": "CC(C)C(=O)OC1CC(C=C2C1C(C)C(=O)C2)(C)CCC3CC(CC3=O)O",
    "Warfarin": "CC(=O)CC(C1=CC=CC=C1)C2=C(O)OC3=CC=CC=C3C2=O",
    "Furosemide": "NS(=O)(=O)c1ccc(Cl)c(NCc2ccco2)c1C(=O)O",
    "Hydrochlorothiazide": "C1=C2C(=CC(=C1Cl)S(=O)(=O)N)S(=O)(=O)NCN2",
    # CNS
    "Fluoxetine": "CNC(C)COC1=CC=C(C=C1)C(C2=CC=CC=C2)C(F)(F)F",
    "Sertraline": "CNC1CCC(C2=C1C3=CC=CC=C3C4=CC=CC=C24)Cl",
    "Diazepam": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
    "Olanzapine": "CN1CCN(CC1)C2=C(C3=CC=CC=C3S2)N4C=CC=CC4=N",
    "Risperidone": "CC1=C(C(=O)N2CCCC2N1)CCC3=CC=C4N3C(=O)C5=C4C=CC=C5F",
    "Haloperidol": "C1CN(CCC1C(=O)C2=CC=C(F)C=C2)CCCC3=CC(=O)C4=C(C3)OC=C4",
    "Clozapine": "CN1CCN(CC1)C2=C3C=CC=CC3=NC4=C2C5=CC=CC5N4",
    "Mirtazapine": "CN1CCN2C(C1)C3=CC=CC=C3C4=C2C=CC=N4",
    "Venlafaxine": "CN(C)CC(C1=CC=C(C=C1)OC)C2(CCCC2)O",
    "Levodopa": "NC(CC1=CC(O)=C(O)C=C1)C(=O)O",
    # Diabetes
    "Metformin": "CN(C)C(=N)N=C(N)N",
    "Pioglitazone": "CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3",
    "Sitagliptin": "CC1(NC(=O)C(C1)N2CC3=C(C2=NN=N3)C4=CC(=C(C=C4)F)F)C(F)(F)F",
    "Empagliflozin": "C1CC(OC1(C2=CC=C(C=C2)Cl)C3=CC=C(C=C3)CC4C(C(C(C(O4)O)O)O)O)O",
    "Dapagliflozin": "CCOC1=CC=C(C=C1)CC2=CC=C(C=C2)C3C(C(C(C(O3)O)O)O)O",
    # Oncology
    "Imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    "Nilotinib": "CN1CCN(CC1)C2=CC=C(C=C2)C3=NC=C(C=N3)NC4=CC(=C(C=C4)C(F)(F)F)NC(=O)C5=CC=CC=C5",
    "Dasatinib": "CC1=C(C(=NN1C2=CC(=C(C=C2)Cl)C(=O)N)NC3=NC(=CS3)C4=CC=CC=C4)CN5CCN(CC5)CCO",
    "Tamoxifen": "CCCC1=CC=C(C=C1)C(=C(C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3)",
    "Letrozole": "C1=CC(=CC=C1C#N)C(CN2C=NC=N2)(C3=CC=C(C=C3)C#N)C#N",
    # Anti-infective
    "Amoxicillin": "CC1(C(=O)NC2(C1SC2(C)C)C(=O)O)NC(=O)C3=CC(=C(C=C3)O)O",
    "Ciprofloxacin": "C1CN(CC1)C2=C(C=C3C(=C2)C(=O)C(=CN3C4CC4)C(=O)O)F",
    "Acyclovir": "C1=NC2=C(N1COCCO)N=C(NC2=O)N",
    # Anti-inflammatory
    "Prednisolone": "CC12CC(C3C(C1CC2C(=O)CO)C4=CC(=O)C=CC4(C3)O)O",
    "Dexamethasone": "CC1CC2C3CCC4=CC(=O)C=C(C4(C3(F)C(CC2(C1(C(=O)CO)O)C)O)C)F",
    "Methotrexate": "CN(CC1=CN=C2N=C(N)N=C(C2=N1)N)C3=CC=C(C=C3)C(=O)O",
    # Respiratory
    "Montelukast": "CC(C)(C(=O)NCC1=CC=C(C=C1)C=C2C3=CC=CC=C3C(=C2)SCCC4=CC=CC=C4)O",
    "Salbutamol": "CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O",
    "Loratadine": "CCOC(=O)N1CCC(=C2C3=CC=CC=C3CCC4=C2N=CC=C4)CC1",
    # Anticoagulant
    "Rivaroxaban": "CC1=CC=C(C=C1)C(=O)N2CCN(C2=O)C3=CC=C(C=C3)NC(=O)NC4=CC=CS4",
    "Apixaban": "CN1C(=O)C2=C(N=C1C3=CC=CC=C3)N(C(=O)N2C4=CC=CC=C4)C5=CC=CC=C5",
    # Endocrine
    "Levothyroxine": "NC(CC1=CC(I)=C(OC2=CC(I)=C(O)C(I)=C2)C(I)=C1)C(=O)O",
    "Spironolactone": "CC12CCC3C4CCC5=CC(=O)CCC5(C4C(=O)C3(C1)CC2)C(=O)SC",
    # GI
    "Omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=CC=C3N2",
    "Ondansetron": "CC1=C(C(=O)N2C1CC3=C2C=CC4=C3CC4)N5CCN(CC5)C",
    # Additional important FDA drugs
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Acetaminophen": "CC(=O)Nc1ccc(O)cc1",
    "Ivermectin_B1a": "CC1C(CCC2(O1)CC3CC(O2)CCC4C3(C(C(C4(C)C)O)OC5CC(C(C(O5)C)OC6CC(C(C(O6)C)O)OC7CC(C(C(O7)C)O)OC(=O)C8COC9C8C(C=C(C)C9)C)C)C)OC",  # noqa: simplified
    "Famotidine": "NC1=NC=C(S1)CS(=O)(=O)N=C(N)N",
    "Ranitidine": "CN(C)Cc1cccs1CNOCCS(=O)(=O)CC#C",
    "Fentanyl": "CCC(=O)N(C1CCN(CC1)Cc2ccccc2)c3ccccc3",
    "Morphine": "Cn1ccC24c5c3c(O)ccc5OC1C4C(C2)NC3",
}


def prepare_3d_molecule(name, smiles, properties=None):
    """Generate 3D conformer from SMILES, write SDF, return record."""
    if properties is None:
        properties = {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol.SetProp("_Name", name)
    for k, v in properties.items():
        mol.SetProp(k, str(v))

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        params = AllChem.ETKDG()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
    if result == 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            pass

    safe_name = name.replace("/", "_").replace(" ", "_").replace(",", "")
    sdf_path = LIGANDS_DIR / f"{safe_name}.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    return {
        "name": safe_name,
        "sdf": str(sdf_path),
        "smiles": smiles,
        "mwt": round(MolWt(mol), 2),
        "logp": round(MolLogP(mol), 2),
        "hbd": NumHDonors(mol),
        "hba": NumHAcceptors(mol),
        "rotatable": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def convert_to_pdbqt_meeko(record):
    """Convert SDF to PDBQT using Meeko."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    sdf_path = record["sdf"]
    name = record["name"]
    pdbqt_path = PDBQT_DIR / f"{name}.pdbqt"

    if pdbqt_path.exists() and os.path.getsize(pdbqt_path) > 100:
        record["pdbqt"] = str(pdbqt_path)
        return True

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(iter(suppl), None)
    if mol is None:
        print(f"  SKIP {name}: cannot read SDF")
        return False

    # Ensure explicit hydrogens
    mol = Chem.AddHs(mol)

    try:
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            print(f"  SKIP {name}: meeko prepare returned empty")
            return False
        pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            print(f"  SKIP {name}: meeko write: {err_msg}")
            return False
        with open(pdbqt_path, "w") as f:
            f.write(pdbqt_string)
        record["pdbqt"] = str(pdbqt_path)
        return True
    except Exception as e:
        print(f"  SKIP {name}: {e}")
        return False


def convert_to_pdbqt_obabel(record):
    """Convert SDF to PDBQT using OpenBabel."""
    sdf_path = record["sdf"]
    name = record["name"]
    pdbqt_path = PDBQT_DIR / f"{name}.pdbqt"

    if pdbqt_path.exists() and os.path.getsize(pdbqt_path) > 100:
        record["pdbqt"] = str(pdbqt_path)
        return True

    obabel = r"C:\anaconda3\Scripts\obabel.exe"
    cmd = f'"{obabel}" "{sdf_path}" -O "{pdbqt_path}" 2>nul'
    ret = os.system(cmd)
    if ret == 0 and pdbqt_path.exists() and os.path.getsize(pdbqt_path) > 50:
        record["pdbqt"] = str(pdbqt_path)
        return True
    return False


def convert_to_pdbqt_manual(record):
    """Fallback: write PDBQT manually from RDKit molecule."""
    from rdkit import Chem

    sdf_path = record["sdf"]
    name = record["name"]
    pdbqt_path = PDBQT_DIR / f"{name}.pdbqt"

    if pdbqt_path.exists() and os.path.getsize(pdbqt_path) > 100:
        record["pdbqt"] = str(pdbqt_path)
        return True

    suppl = Chem.SDMolSupplier(sdf_path)
    mol = next(iter(suppl), None)
    if mol is None:
        return False

    conf = mol.GetConformer()
    lines = [f"REMARK  GENERATED BY RDKit (manual PDBQT writer)\n"]
    for i, atom in enumerate(mol.GetAtoms(), 1):
        pos = conf.GetAtomPosition(atom.GetIdx())
        elem = atom.GetSymbol()
        ad4 = elem if elem != 'H' else 'H'  # simplified typing
        if elem == 'C': ad4 = 'C'
        elif elem == 'N': ad4 = 'N'
        elif elem == 'O': ad4 = 'O'
        elif elem == 'S': ad4 = 'S'
        elif elem == 'F': ad4 = 'F'
        elif elem == 'Cl': ad4 = 'Cl'
        elif elem == 'Br': ad4 = 'Br'
        elif elem == 'I': ad4 = 'I'
        elif elem == 'P': ad4 = 'P'

        # RDKit atoms store their names via Prop: "computed::Name" or _Name
        atom_name = mol.GetAtomWithIdx(atom.GetIdx()).GetProp("_Name") if mol.GetAtomWithIdx(atom.GetIdx()).HasProp("_Name") else atom.GetSymbol()
        pdb_name = (atom.GetSymbol() + str(i)).ljust(4)[:4]
        line = f"ATOM  {i:>5d} {pdb_name:<4s} UNK     1    {pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  0.00          {ad4:<2s}\n"
        lines.append(line)
    lines.append("END\n")

    with open(pdbqt_path, "w") as f:
        f.writelines(lines)
    record["pdbqt"] = str(pdbqt_path)
    return True


def main():
    print("=" * 60)
    print("STEP 1: Prepare known SMVT interactors")
    print("=" * 60)
    known_records = []
    for name, info in KNOWN_INTERACTORS.items():
        rec = prepare_3d_molecule(name, info["SMILES"], {
            "PUBCHEM_CID": str(info.get("PubChem_CID", "")),
            "Type": info.get("Type", ""),
            "Known_Kd_uM": str(info.get("Known_Kd_uM", "")),
        })
        if rec:
            rec["type"] = info.get("Type", "")
            rec["known_kd_uM"] = info.get("Known_Kd_uM", "")
            known_records.append(rec)
            print(f"  {name}: MW={rec['mwt']}, LogP={rec['logp']}")
        else:
            print(f"  FAILED: {name}")
    print(f"  Total known interactors: {len(known_records)}")

    print("\n" + "=" * 60)
    print("STEP 2: Build FDA-approved drug library")
    print("=" * 60)
    fda_records = []
    for name, smiles in FDA_DRUGS.items():
        rec = prepare_3d_molecule(name, smiles, {"Type": "FDA-approved"})
        if rec:
            rec["type"] = "FDA-approved"
            rec["known_kd_uM"] = ""
            fda_records.append(rec)
        else:
            print(f"  FAILED: {name}")
    print(f"  Total FDA drugs prepared: {len(fda_records)}")

    all_records = known_records + fda_records
    print(f"\nTotal ligands: {len(all_records)}")

    # Test meeko first
    meeko_ok = False
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        meeko_ok = True
    except ImportError:
        pass

    # Convert to PDBQT
    print("\n" + "=" * 60)
    print("STEP 3: Convert all ligands to PDBQT")
    print("=" * 60)

    meeko_success = 0
    obabel_success = 0
    manual_success = 0

    for rec in all_records:
        converted = False
        if meeko_ok:
            if convert_to_pdbqt_meeko(rec):
                meeko_success += 1
                converted = True
        if not converted and convert_to_pdbqt_obabel(rec):
            obabel_success += 1
            converted = True
        if not converted and convert_to_pdbqt_manual(rec):
            manual_success += 1
            converted = True
        if not converted:
            print(f"  FAILED PDBQT: {rec['name']}")

    total_pdbqt = meeko_success + obabel_success + manual_success
    print(f"  PDBQT conversions: {total_pdbqt}/{len(all_records)}")
    print(f"    Meeko:  {meeko_success}")
    print(f"    Obabel: {obabel_success}")
    print(f"    Manual: {manual_success}")

    # Save ligand database
    with open(DATA_DIR / "ligands_database.json", "w") as f:
        json.dump(all_records, f, indent=2)
    print(f"\nLigand database saved: {DATA_DIR / 'ligands_database.json'}")
    print(f"  Total: {len(all_records)} ({len(known_records)} known + {len(fda_records)} FDA)")


if __name__ == "__main__":
    main()
