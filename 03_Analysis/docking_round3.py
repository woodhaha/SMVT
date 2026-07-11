#!/usr/bin/env python3
"""SMVT Docking Round 3 — SAR-guided expansion (fenamates + small acids + biotin analogs)"""
import os, sys, subprocess, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
vina_bin = "os.environ.get("VINA_BIN", "vina")"
mk_prep_rec = "C:/anaconda3/Scripts/mk_prepare_receptor.exe"
mk_prep_lig = "C:/anaconda3/Scripts/mk_prepare_ligand.exe"
obabel_bin = "C:/anaconda3/Scripts/obabel"
center = [-2.5, 1.0, -1.0]; box_size = [22, 22, 22]

# Load existing results to skip
existing = set()
for f in ["03_Analysis/docking/docking_expanded_results.csv"]:
    if os.path.exists(f):
        df = pd.read_csv(f)
        existing.update(df['name'].tolist())

# ═══ Round 3: SAR-guided expansion ═══
compounds = [
    # ── More fenamates (the #1 chemical family) ──
    ("Tolfenamic_Acid", "CC1=CC=CC(=C1NC2=CC=CC=C2C(=O)O)Cl", "nsaid"),
    ("Niflumic_Acid", "FC(F)(F)C1=CC=NC(=C1)NC2=CC=CC=C2C(=O)O", "nsaid"),
    ("Etofenamate", "CCOCCOC(=O)C1=CC=CC=C1NC2=CC=CC(=C2)C(F)(F)F", "nsaid"),
    ("Meclofenamic_Acid", "CC1=CC=C(C(=C1)NC2=CC=CC=C2C(=O)O)Cl", "nsaid"),
    ("Flufenamic_Acid", "FC(F)(F)C1=CC=CC(=C1)NC2=CC=CC=C2C(=O)O", "nsaid"),

    # ── Other NSAID subclasses ──
    ("Nabumetone", "COC1=CC2=C(C=C1)C=C(C=C2)CCC(=O)C", "nsaid"),
    ("Oxaprozin", "C1=CC=C(C=C1)C2=CC(=NO2)CCC(=O)O", "nsaid"),
    ("Etodolac", "CCC1=C2C(=CC=C1)OC3=C(C2=CC(=C3)CC(=O)O)CC", "nsaid"),
    ("Fenbufen", "C1=CC=C(C=C1)C(=O)CCC2=CC=C(C=C2)CC(=O)O", "nsaid"),
    ("Tiaprofenic_Acid", "CC(C1=CC=CC=C1C(=O)C2=CC=CS2)C(=O)O", "nsaid"),

    # ── Small carboxylic acids (< 200 Da) ──
    ("Benzoic_Acid", "C1=CC=CC=C1C(=O)O", "fda"),
    ("Phenylacetic_Acid", "C1=CC=CC=C1CC(=O)O", "fda"),
    ("Hippuric_Acid", "C1=CC=CC=C1C(=O)NCC(=O)O", "fda"),
    ("Succinic_Acid", "C(CC(=O)O)C(=O)O", "fda"),
    ("Glutaric_Acid", "C(CC(=O)O)CC(=O)O", "fda"),
    ("Adipic_Acid", "C(CCC(=O)O)CC(=O)O", "fda"),
    ("Pimelic_Acid", "C(CCCC(=O)O)CC(=O)O", "fda"),
    ("Suberic_Acid", "C(CCCCC(=O)O)CC(=O)O", "fda"),
    ("Azelaic_Acid", "C(CCCCCCC(=O)O)CC(=O)O", "fda"),
    ("Sebacic_Acid", "C(CCCCCCCC(=O)O)CC(=O)O", "fda"),

    # ── Biotin analogs / derivatives ──
    ("Biotin_Methyl_Ester", "C1C2C(C(S1)CCCCC(=O)OC)NC(=O)N2", "substrate"),
    ("Biotin_Sulfoxide", "C1C2C(C(S1(=O))CCCCC(=O)O)NC(=O)N2", "substrate"),
    ("Biotin_Sulfone", "C1C2C(C(S1(=O)=O)CCCCC(=O)O)NC(=O)N2", "substrate"),
    ("Norbiotin", "C1C2C(C(S1)CCCC(=O)O)NC(=O)N2", "substrate"),
    ("Homobiotin", "C1C2C(C(S1)CCCCCC(=O)O)NC(=O)N2", "substrate"),

    # ── Aromatic amino acids (all have -COOH) ──
    ("L-Phenylalanine", "C1=CC=C(C=C1)CC(C(=O)O)N", "fda"),
    ("L-Tyrosine", "C1=CC(=CC=C1CC(C(=O)O)N)O", "fda"),
    ("L-Tryptophan", "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N", "fda"),
    ("5-Hydroxytryptophan", "C1=CC2=C(C=C1O)C(=CN2)CC(C(=O)O)N", "fda"),

    # ── Bile acids (cholanic acid family) ──
    ("Cholic_Acid", "CC(CCC(=O)O)C1CCC2C1(C(CC3C2C(CC4C3(CCC(C4)O)C)O)O)C", "fda"),
    ("Deoxycholic_Acid", "CC(CCC(=O)O)C1CCC2C1(C(CC3C2CC(C4C3(CCC(C4)O)C)O)O)C", "fda"),
    ("Ursodeoxycholic_Acid", "CC(CCC(=O)O)C1CCC2C1(CCC3C2C(CC4C3(CCC(C4)O)C)O)C", "fda"),

    # ── Other FDA drugs with -COOH (diverse scaffolds) ──
    ("Gemfibrozil", "CC1=CC(=C(C=C1)C)OCCCC(C)(C)C(=O)O", "fda"),
    ("Bezafibrate", "CC(C)(C(=O)O)OC1=CC=C(C=C1)CCNC(=O)C2=CC=C(C=C2)Cl", "fda"),
    ("Ciprofibrate", "CC(C)(C(=O)O)OC1=CC=C(C=C1)C2CC2", "fda"),
    ("Montelukast", "CC(C)(C1=CC=CC=C1CC[C@H](C2=CC=CC(=C2)CC=C(C)C)SCC3(CC3)CC(=O)O)O", "fda"),
    ("Zafirlukast", "CC1=CC(=C(C=C1)C2=CC=C(C=C2)C(=O)NS(=O)(=O)C3=CC=CC=C3C(=O)O)N", "fda"),
]

# Filter out already-docked
new_compounds = [(n,s,t) for n,s,t in compounds if n not in existing]
print(f"Round 3: {len(new_compounds)} new compounds (skipping {len(compounds)-len(new_compounds)} already docked)")

receptor_pdbqt = "03_Analysis/docking/SMVT_receptor.pdbqt"
results = []

for i, (name, smiles, ctype) in enumerate(new_compounds):
    ligand_sdf = f"03_Analysis/docking/{name}.sdf"
    ligand_pdbqt = f"03_Analysis/docking/{name}.pdbqt"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [{i+1}/{len(new_compounds)}] {name}: SKIP (SMILES)")
        continue
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print(f"  [{i+1}/{len(new_compounds)}] {name}: SKIP (3D)")
        continue
    Chem.MolToMolFile(mol, ligand_sdf)

    r = subprocess.run([mk_prep_lig, "-i", ligand_sdf, "-o", ligand_pdbqt], capture_output=True, text=True)
    if not os.path.exists(ligand_pdbqt):
        subprocess.run([obabel_bin, ligand_sdf, "-O", ligand_pdbqt, "--gen3d"], capture_output=True)

    output = f"03_Analysis/docking/{name}_docked.pdbqt"
    cmd = [vina_bin, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(box_size[0]), "--size_y", str(box_size[1]), "--size_z", str(box_size[2]),
           "--out", output, "--exhaustiveness", "16", "--num_modes", "5"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        affs = []
        for line in proc.stdout.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in '12345':
                try: affs.append(float(parts[1]))
                except: pass
        best = min(affs) if affs else 99.0
        results.append({'name': name, 'type': ctype, 'best_affinity': best})
        hit = " ***" if best < -7 else " **" if best < -6.5 else ""
        print(f"  [{i+1}/{len(new_compounds)}] {name:25s}  {best:7.1f}{hit}")
    except:
        print(f"  [{i+1}/{len(new_compounds)}] {name}: FAILED")

# Merge with previous results
prev = pd.read_csv("03_Analysis/docking/docking_expanded_results.csv")
merged = pd.concat([prev, pd.DataFrame(results)], ignore_index=True).sort_values('best_affinity')
merged.to_csv("03_Analysis/docking/docking_expanded_results.csv", index=False)
merged.to_csv("03_Analysis/outputs/docking_expanded_results.csv", index=False)

print(f"\n{'='*60}")
print(f"ROUND 3 COMPLETE — Total: {len(merged)} compounds")
print(f"{'='*60}")
print(f"Hits (< -7): {len(merged[merged.best_affinity < -7])}")
print(f"Strong (< -6.5): {len(merged[merged.best_affinity < -6.5])}")

# Top new finds
new_df = pd.DataFrame(results).sort_values('best_affinity')
print(f"\nNew this round (top 10):")
for _, r in new_df.head(10).iterrows():
    print(f"  {r['name']:25s} {r['best_affinity']:6.1f} ({r['type']})")

print("\nDONE")
