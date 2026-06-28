# SMVT MD Simulation Package — Top 4 Hits + Controls

> **Target**: SMVT (SLC5A6) Na⁺-dependent multivitamin transporter
> **Receptor**: AlphaFold AF-Q9Y289-F1 (energy-minimized)
> **Protocol**: Minim → NVT 100ps → NPT 200ps → Production 100ns
> **Date**: 2026-06-27

## Compounds (7 total)

### Test Compounds (Top 4 from virtual screening)

| # | Compound | ΔG (kcal/mol) | Class | Role |
|:---:|----------|:---:|-------|------|
| 1 | **Hydromorphone** | −8.58 | Opioid analgesic | Best hit overall |
| 2 | **Furosemide** | −8.36 | Loop diuretic | Sulfonamide class |
| 3 | **Naftazone** | −8.34 | Naphthoquinone | Phase A best |
| 4 | **Phenobarbital** | −8.30 | Barbiturate | Scaffold representative |

### Controls

| # | Compound | ΔG (kcal/mol) | Class | Role |
|:---:|----------|:---:|-------|------|
| 5 | **Biotin** | −6.76 | Vitamin B7 | **Natural substrate (reference)** |
| 6 | **Gabapentin enacarbil** | −6.63 | Prodrug | **Positive control (FDA-approved SMVT drug)** |
| 7 | **Riboflavin** | −0.01 | Vitamin B2 | **Negative control (vitamin, no SMVT binding)** |

## SMILES

```
HYDROMORPHONE:       CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5
FUROSEMIDE:          NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl
NAFTAZONE:           C1=CC=C2C(=C1)C=CC(=O)C2=NNC(=N)N
PHENOBARBITAL:       CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2
BIOTIN:              C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2
GABAPENTIN_ENACARBIL: CC(OC(=O)NCC1(CC(=O)O)CCCCC1)OC(=O)C(C)C
RIBOFLAVIN:          CC1=CC2=C(C=C1C)N(C[C@@H]([C@@H]([C@@H](CO)O)O)O)C3=NC(=O)NC(=O)C3=N2
```

## MD Protocol

| Phase | Duration | Ensemble | dt | Purpose |
|-------|:---:|------|:---:|------|
| Minimization | convergence | NVE | — | Remove steric clashes |
| NVT Equilibration | 100 ps | NVT | 2 fs | Thermalize to 300K |
| NPT Equilibration | 200 ps | NPT | 2 fs | Equilibrate density (1 atm) |
| **Production** | **100 ns** | NPT | 2 fs | **Binding stability assessment** |

### Parameters
- Force field: AMBER ff14SB (protein) + GAFF2 (ligand) + TIP3P (water)
- Temperature: 300 K (Langevin, γ=1 ps⁻¹)
- Pressure: 1 atm (Monte Carlo barostat)
- Cutoff: 10 Å (PME for long-range electrostatics)
- Padding: 12 Å (0.15 M NaCl)

## Analysis Outputs

| Metric | Description |
|--------|-------------|
| RMSD (protein Cα) | Receptor stability |
| RMSD (ligand) | Ligand pose stability |
| RMSF (per-residue) | Binding site flexibility |
| H-bond occupancy | Key interactions persistence |
| Contact distance | Ligand-receptor distance over time |
| Binding free energy | MM/GBSA (last 20ns) |
| SASA | Solvent accessible surface area changes |

## Expected Outcomes

- **Hits**: RMSD < 3 Å → stable binding, confirms docking predictions
- **Biotin**: Moderate RMSD (~2-4 Å) → natural substrate behavior
- **Gabapentin enacarbil**: Higher RMSD → transient substrate interaction
- **Riboflavin**: Unstable binding, ligand dissociation → validates screening specificity

## Usage

```bash
# 1. Prepare ligands
python scripts/prepare_ligands.py

# 2. Run MD (submit one per compound, each 100ns ~2-4 days on GPU)
python scripts/run_md.py --compound HYDROMORPHONE --ns 100
python scripts/run_md.py --compound FUROSEMIDE --ns 100
# ... etc

# 3. Analyze trajectories
python scripts/analyze_md.py --all

# 4. Generate figures
python scripts/visualize_md.py --all
```

## Receptor

```
02_Data/raw/AF-Q9Y289-F1.pdb  — AlphaFold v6, pLDDT=79.4
02_Data/cleaned/SMVT_prepared.pdb  — energy-minimized (UFF → OpenMM)
```

## Dependencies

```
conda create -n smvt-md python=3.11
conda install -c conda-forge openmm openff-toolkit pdbfixer
pip install mdtraj nglview matplotlib seaborn
```
