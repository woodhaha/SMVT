#!/usr/bin/env python3
"""SMVT Deep Docking — Expanded library (50+ compounds)"""
import os, sys, subprocess
import pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

os.chdir("D:/Researching/SMVT")
vina_bin = "C:/Users/woodh/bin/vina"
mk_prep_rec = "C:/anaconda3/Scripts/mk_prepare_receptor.exe"
mk_prep_lig = "C:/anaconda3/Scripts/mk_prepare_ligand.exe"
obabel_bin = "C:/anaconda3/Scripts/obabel"

center = [-2.5, 1.0, -1.0]
box_size = [22, 22, 22]

os.makedirs("03_Analysis/docking", exist_ok=True)

# ═══ Extended compound library ═══
compounds = [
    # ── Round 1 (existing, re-dock with higher exhaustiveness) ──
    ("Biotin", "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2", "substrate"),
    ("Lipoic_Acid", "C1CSSC1CCCCC(=O)O", "substrate"),
    ("Pantothenic_Acid", "CC(C)(CO)C(O)C(=O)NCCC(=O)O", "substrate"),
    ("Gabapentin_enacarbil", "CC(C)(C(=O)O)OC(=O)NCC1(CCCCC1)CC(=O)O", "substrate"),
    ("Desthiobiotin", "CC1C(NC(=O)N1)CCCCCC(=O)O", "substrate"),
    ("Indomethacin", "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=C(C=C2)OC)CC(=O)O", "nsaid"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "nsaid"),
    ("Diclofenac", "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl", "nsaid"),
    ("Ketoprofen", "CC(C1=CC=CC(=C1)C(=O)C2=CC=CC=C2)C(=O)O", "nsaid"),
    ("Flurbiprofen", "CC1=CC(=CC=C1)C2=CC=C(C=C2)C(C)C(=O)O", "nsaid"),
    ("Phenylbutazone", "CCCCC1C(=O)N(N(C1=O)C2=CC=CC=C2)C3=CC=CC=C3", "nsaid"),
    ("Naproxen", "CC(C1=CC2=C(C=C1)C=CC(=C2)OC)C(=O)O", "nsaid"),
    ("Celecoxib", "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F", "nsaid"),
    ("Piroxicam", "CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3", "nsaid"),
    ("Meloxicam", "CC1=CC=CS(=O)(=O)N1C(=O)NC2=C(C3=NC=C(S3)C)O2", "nsaid"),
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "nsaid"),
    ("Sulindac", "CC1=C(C2=C(C=CC(=C2)F)C(=C1)CC(=O)O)C3=CC=C(C=C3)S(=O)C", "nsaid"),
    ("Mefenamic_Acid", "CC1=CC=CC(=C1NC2=CC=CC=C2C(=O)O)C", "nsaid"),

    # ── Statins (all have carboxylic acid moiety) ──
    ("Simvastatin", "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(O)(CC(=O)O3)C)C", "fda"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "fda"),
    ("Rosuvastatin", "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "fda"),
    ("Pravastatin", "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(O)(CC(=O)O3)C)C", "fda"),
    ("Fluvastatin", "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "fda"),

    # ── ACE inhibitors (carboxylic acid) ──
    ("Enalaprilat", "CC(C(=O)O)NC(C)C(=O)N1CC2=CC=CC=C2CC1C(=O)O", "fda"),
    ("Lisinopril", "CC(C(=O)O)NC(CCC1=CC=CC=C1)C(=O)N2CCCC2C(=O)O", "fda"),
    ("Captopril", "CC(CS)C(=O)N1CCCC1C(=O)O", "fda"),

    # ── Other carboxylic acid drugs ──
    ("Valproic_Acid", "CCCC(CCC)C(=O)O", "fda"),
    ("Methotrexate", "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O", "fda"),
    ("Levodopa", "C1=CC(=C(C=C1CC(C(=O)O)N)O)O", "fda"),
    ("Tranexamic_Acid", "C1CC(CCC1CN)C(=O)O", "fda"),

    # ── Other FDA drugs ──
    ("Metformin", "CN(C)C(=N)NC(=N)N", "fda"),
    ("Omeprazole", "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC", "fda"),
    ("Metronidazole", "CC1=NC=C(N1CCO)[N+](=O)[O-]", "fda"),
    ("Fluoxetine", "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F", "fda"),
    ("Sertraline", "CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl", "fda"),
    ("Warfarin", "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O", "fda"),
    ("Furosemide", "C1=CC(=C(C=C1C(=O)O)S(=O)(=O)N)NC2=NC=CC(=C2)Cl", "fda"),
    ("Hydrochlorothiazide", "C1=NC2=C(CC(NS(=O)(=O)C3=CC=CC=C3Cl)S2)C(=O)N1", "fda"),
    ("Acetazolamide", "CC(=O)NC1=NC(=C(S1)S(=O)(=O)N)S(=O)(=O)N", "fda"),
    ("Probenecid", "CCCN(CCC)S(=O)(=O)C1=CC=C(C=C1)C(=O)O", "fda"),

    # ── Vitamins / cofactors ──
    ("Thiamine", "CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO", "vitamin"),
    ("Riboflavin", "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)CC(C(C(CO)O)O)O", "vitamin"),
    ("Ascorbic_Acid", "C(C(C1C(=C(C(=O)O1)O)O)O)O", "vitamin"),
    ("Folic_Acid", "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)N=C(N3)N", "vitamin"),
    ("Pyridoxine", "CC1=NC=C(C(=C1O)CO)CO", "vitamin"),
    ("Nicotinic_Acid", "C1=CC(=CN=C1)C(=O)O", "vitamin"),
    ("Retinoic_Acid", "CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC(=O)O)C)C", "vitamin"),

    # ── Additional diverse compounds ──
    ("Salicylic_Acid", "C1=CC=C(C(=C1)C(=O)O)O", "fda"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "fda"),
    ("Acetaminophen", "CC(=O)NC1=CC=C(C=C1)O", "fda"),
    ("Isoniazid", "C1=CN=CC=C1C(=O)NN", "fda"),
    ("Allopurinol", "C1=NC(=O)C2=C(N1)N=CN2", "fda"),
]

print(f"Expanded library: {len(compounds)} compounds")
print(f"  Substrates: {sum(1 for _,_,t in compounds if t=='substrate')}")
print(f"  NSAIDs: {sum(1 for _,_,t in compounds if t=='nsaid')}")
print(f"  FDA drugs: {sum(1 for _,_,t in compounds if t=='fda')}")
print(f"  Vitamins: {sum(1 for _,_,t in compounds if t=='vitamin')}")

# ═══ Prepare receptor (once) ═══
receptor_pdbqt = "03_Analysis/docking/SMVT_receptor.pdbqt"
if not os.path.exists(receptor_pdbqt):
    print("\nPreparing receptor...")
    subprocess.run([mk_prep_rec, "-i", "02_Data/cleaned/AF-Q9Y289-F1_prepared.pdb",
                    "-o", receptor_pdbqt], capture_output=True)

# ═══ Dock all ═══
results = []
existing = set()  # skip already-docked from round 1
for f in os.listdir("03_Analysis/docking"):
    if f.endswith("_docked.pdbqt"):
        existing.add(f.replace("_docked.pdbqt", ""))

print(f"\nDocking {len(compounds)} compounds (skipping {len(existing)} already docked)...")
for i, (name, smiles, ctype) in enumerate(compounds):
    if name in existing:
        # Load previous result
        prev_df = pd.read_csv("03_Analysis/docking/docking_results.csv") if os.path.exists("03_Analysis/docking/docking_results.csv") else None
        if prev_df is not None and name in prev_df['name'].values:
            row = prev_df[prev_df['name'] == name].iloc[0]
            results.append({'name': name, 'type': ctype, 'best_affinity': row['best_affinity']})
            continue

    ligand_sdf = f"03_Analysis/docking/{name}.sdf"
    ligand_pdbqt = f"03_Analysis/docking/{name}.pdbqt"

    # Generate 3D
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [{i+1}/{len(compounds)}] {name}: SKIP (SMILES)")
        continue
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print(f"  [{i+1}/{len(compounds)}] {name}: SKIP (3D gen)")
        continue
    Chem.MolToMolFile(mol, ligand_sdf)

    # Prepare
    r = subprocess.run([mk_prep_lig, "-i", ligand_sdf, "-o", ligand_pdbqt], capture_output=True, text=True)
    if not os.path.exists(ligand_pdbqt):
        subprocess.run([obabel_bin, ligand_sdf, "-O", ligand_pdbqt, "--gen3d"], capture_output=True)

    # Dock
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
        hit = " ***" if best < -7 else " **" if best < -6 else ""
        print(f"  [{i+1}/{len(compounds)}] {name:20s}  {best:7.1f}{hit}")
    except subprocess.TimeoutExpired:
        print(f"  [{i+1}/{len(compounds)}] {name}: TIMEOUT")
    except Exception as e:
        print(f"  [{i+1}/{len(compounds)}] {name}: ERROR")

# ═══ Save results ═══
df = pd.DataFrame(results).sort_values('best_affinity')
df.to_csv("03_Analysis/docking/docking_expanded_results.csv", index=False)
df.to_csv("03_Analysis/outputs/docking_expanded_results.csv", index=False)

# Summary
print(f"\n{'='*60}")
print(f"EXPANDED DOCKING RESULTS ({len(df)} compounds)")
print(f"{'='*60}")
for _, r in df.head(20).iterrows():
    m = " ***" if r['best_affinity'] < -7 else " **" if r['best_affinity'] < -6 else ""
    print(f"  {r['name']:22s} {r['type']:10s} {r['best_affinity']:7.1f}{m}")

hits = df[df['best_affinity'] < -7]
novel = df[(df['best_affinity'] < -6.5) & (df['type'].isin(['nsaid','fda']))]
print(f"\nHits (< -7): {len(hits)}")
print(f"Novel FDA/NSAID hits (< -6.5): {len(novel)}")
if len(novel) > 0:
    print("Novel hits:", ', '.join(novel['name'].tolist()))

print("\nDONE")
