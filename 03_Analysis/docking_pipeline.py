#!/usr/bin/env python3
"""SMVT Molecular Docking — vina binary + meeko + rdkit"""
import os, sys, subprocess, shutil
import numpy as np
import pandas as pd
from pathlib import Path

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
vina_bin = "os.environ.get("VINA_BIN", "vina")"
mk_prep_rec = "C:/anaconda3/Scripts/mk_prepare_receptor.exe"
mk_prep_lig = "C:/anaconda3/Scripts/mk_prepare_ligand.exe"
obabel_bin = "C:/anaconda3/Scripts/obabel"
receptor_pdb = "02_Data/cleaned/AF-Q9Y289-F1_prepared.pdb"

print("="*60)
print("SMVT (SLC5A6) Molecular Docking — FDA Drug Virtual Screening")
print("="*60)

# ═══ Step 1: Prepare receptor ══════════════════════════════════
print("\n[1/5] Preparing SMVT receptor...")

os.makedirs("03_Analysis/docking", exist_ok=True)
receptor_pdbqt = "03_Analysis/docking/SMVT_receptor.pdbqt"

# Use meeko to prepare protein
result = subprocess.run([
    mk_prep_rec, "-i", receptor_pdb, "-o", receptor_pdbqt
], capture_output=True, text=True)
if result.returncode != 0:
    print(f"  mk_prepare_receptor failed: {result.stderr[:200]}")
    # Fallback: use obabel
    print("  Falling back to obabel...")
    result = subprocess.run([
        obabel_bin, receptor_pdb, "-O", receptor_pdbqt, "-xr"
    ], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  obabel also failed: {result.stderr[:200]}")
        sys.exit(1)
print(f"  Receptor prepared: {receptor_pdbqt}")

# ═══ Step 2: Define docking box ═════════════════════════════════
# SMVT central cavity — approximate center from AlphaFold structure
# Protein center roughly at (0,0,0) in AF coordinates
# Box covers the central substrate-binding cavity
center = [-2.5, 1.0, -1.0]  # x, y, z approximate
box_size = [22, 22, 22]  # 22A cube covers the entire TM cavity
print(f"\n[2/5] Docking box: center={center}, size={box_size}A")

# ═══ Step 3: Build FDA drug library ═══════════════════════════
print("\n[3/5] Building ligand library...")

# Known SMVT interactors + FDA drugs with structural diversity
ligands = [
    # ── SMVT substrates (positive controls, should dock well) ──
    ("Biotin", "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2", "DB00121"),
    ("Lipoic_Acid", "C1CSSC1CCCCC(=O)O", "DB00166"),
    ("Pantothenic_Acid", "CC(C)(CO)C(O)C(=O)NCCC(=O)O", "DB01783"),
    ("Gabapentin_enacarbil", "CC(C)(C(=O)O)OC(=O)NCC1(CCCCC1)CC(=O)O", "DB08848"),

    # ── NSAID inhibitors (known SMVT inhibitors from literature) ──
    ("Indomethacin", "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=C(C=C2)OC)CC(=O)O", "DB00328"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "DB01050"),
    ("Diclofenac", "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl", "DB00586"),
    ("Ketoprofen", "CC(C1=CC=CC(=C1)C(=O)C2=CC=CC=C2)C(=O)O", "DB01009"),
    ("Flurbiprofen", "CC1=CC(=CC=C1)C2=CC=C(C=C2)C(C)C(=O)O", "DB00712"),

    # ── Other FDA drugs from DrugBank ──
    ("Metformin", "CN(C)C(=N)NC(=N)N", "DB00331"),
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "DB00945"),
    ("Simvastatin", "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(O)(CC(=O)O3)C)C", "DB00641"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "DB01076"),
    ("Omeprazole", "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC", "DB00338"),
    ("Metronidazole", "CC1=NC=C(N1CCO)[N+](=O)[O-]", "DB00916"),
    ("Fluoxetine", "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F", "DB00472"),
    ("Sertraline", "CN[C@H]1CC[C@H](C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl", "DB01104"),
    ("Warfarin", "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O", "DB00682"),
    ("Furosemide", "C1=CC(=C(C=C1C(=O)O)S(=O)(=O)N)NC2=NC=CC(=C2)Cl", "DB00695"),

    # ── Vitamin/cofactor analogs ──
    ("Thiamine", "CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO", "DB00152"),
    ("Riboflavin", "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)CC(C(C(CO)O)O)O", "DB00140"),
    ("Ascorbic_Acid", "C(C(C1C(=C(C(=O)O1)O)O)O)O", "DB00126"),
    ("Folic_Acid", "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)N=C(N3)N", "DB00158"),
]

print(f"  Library: {len(ligands)} compounds")
print(f"  Positive controls (SMVT substrates): Biotin, Lipoic_Acid, Pantothenic_Acid, Gabapentin_enacarbil")
print(f"  Literature inhibitors (NSAIDs): Indomethacin, Ibuprofen, Diclofenac, Ketoprofen, Flurbiprofen")

# ═══ Step 4: Prepare & dock all ligands ═══════════════════════
print("\n[4/5] Docking...")
results = []

for name, smiles, db_id in ligands:
    ligand_sdf = f"03_Analysis/docking/{name}.sdf"
    ligand_pdbqt = f"03_Analysis/docking/{name}.pdbqt"

    # Generate 3D structure from SMILES
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  {name}: SMILES parse FAILED")
        continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, ligand_sdf)

    # Convert to PDBQT with meeko
    result = subprocess.run([
        mk_prep_lig, "-i", ligand_sdf, "-o", ligand_pdbqt
    ], capture_output=True, text=True)
    if result.returncode != 0:
        # Fallback: obabel
        result = subprocess.run([
            obabel_bin, ligand_sdf, "-O", ligand_pdbqt, "--gen3d"
        ], capture_output=True, text=True)

    if not os.path.exists(ligand_pdbqt):
        print(f"  {name}: PDBQT preparation FAILED")
        continue

    # Run vina
    output_pdbqt = f"03_Analysis/docking/{name}_docked.pdbqt"
    cmd = [
        vina_bin,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--out", output_pdbqt,
        "--exhaustiveness", "8",
        "--num_modes", "5",
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        # Parse affinity from output
        affinities = []
        for line in proc.stdout.split('\n'):
            if line.strip().startswith(('1','2','3','4','5')):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        aff = float(parts[1])
                        affinities.append(aff)
                    except ValueError:
                        pass

        best_aff = min(affinities) if affinities else 99.0
        results.append({
            'name': name,
            'drugbank_id': db_id,
            'best_affinity': best_aff,
            'all_affinities': str(affinities),
            'type': 'substrate' if name in ['Biotin','Lipoic_Acid','Pantothenic_Acid','Gabapentin_enacarbil']
                    else 'NSAID_inhibitor' if name in ['Indomethacin','Ibuprofen','Diclofenac','Ketoprofen','Flurbiprofen']
                    else 'FDA_drug',
        })
        print(f"  {name:20s}: best={best_aff:.1f} kcal/mol")

    except subprocess.TimeoutExpired:
        print(f"  {name}: TIMEOUT")
    except Exception as e:
        print(f"  {name}: ERROR — {e}")

# ═══ Step 5: Analyze results ═══════════════════════════════════
print(f"\n[5/5] Results")
results_df = pd.DataFrame(results).sort_values('best_affinity')

# Save
results_df.to_csv("03_Analysis/docking/docking_results.csv", index=False)

# Print ranked results
print(f"\n{'='*60}")
print("DOCKING RESULTS — Ranked by binding affinity")
print(f"{'='*60}")
print(f"{'Rank':<5} {'Compound':<22} {'Type':<18} {'Affinity':>8}")
print("-"*55)
hit_count = 0
for i, (_, r) in enumerate(results_df.iterrows()):
    marker = " ***" if r['best_affinity'] < -7 else ""
    if r['best_affinity'] < -7:
        hit_count += 1
    print(f"{i+1:<5} {r['name']:<22} {r['type']:<18} {r['best_affinity']:>7.1f}{marker}")

print(f"\nHits (affinity < -7 kcal/mol): {hit_count}/{len(results_df)}")

# Check if positive controls rank well (validation)
substrates = results_df[results_df['type'] == 'substrate']
if len(substrates) > 0:
    best_substrate = substrates.iloc[0]
    print(f"Best substrate (positive control): {best_substrate['name']} — {best_substrate['best_affinity']:.1f} kcal/mol")
    print(f"Substrate mean rank: {results_df.index[results_df['type']=='substrate'].mean()+1:.1f} / {len(results_df)}")

results_df.to_csv("03_Analysis/outputs/docking_results.csv", index=False)
print(f"\nSaved: 03_Analysis/docking/docking_results.csv")
print(f"Saved: 03_Analysis/outputs/docking_results.csv")
print("DONE")
