# SMVT (SLC5A6) — Pan-cancer analysis & structure-based virtual screening

Companion code for *Structure-Based Discovery of Carboxyl-Bearing Ligands of the Sodium-Dependent Multivitamin Transporter (SMVT/SLC5A6)*.

## Pipeline

| Module | Scripts |
|--------|---------|
| 1. Pan-cancer expression & survival | `03_Analysis/run_pancancer.py`, `multivariate_survival.py`, `survival_gdc.py`, `survival_analysis_real.py` |
| 2. Mutation landscape + ESM-2 pathogenicity | `03_Analysis/extract_slc5a6_maf.py`, `fig_mutation_final.py`, `esm_pathogenicity_mc3.py`, `esm650_crossval.py` |
| 3. PPI & pathway enrichment | `03_Analysis/visualize_string_network.py`, `pathlinkR_analysis.R` |
| 4. Virtual KO | `03_Analysis/run_scTenifoldKnk_light.R`, `pc_regression_GRN.py` |
| 5. Virtual screening (ML + docking + SAR) | `03_Analysis/pharmacophore_ml_screen.py`, `docking_batch_screen.py`, `docking_expanded.py`, `merge_anticancer.py`, `sar_deep_analysis.py`, `pathc_v2.py`, `glutarimide_hunt.py` |
| 6. MD + MM-GBSA | `SMVT-MD-Analysis/`, `MM-GBSA Analysis/` |

## Requirements

- Python 3.x: RDKit, AutoDock Vina, OpenMM ≥ 8.5, mdtraj, matplotlib, numpy/pandas, meeko (ligand preparation)
- R: pathlinkR, scTenifoldKnk

## Data

Input datasets and access URLs are listed in [`data_sources.md`](data_sources.md). Docking results are under `03_Analysis/outputs/docking_final/`.

## Key results

- Virtual screening: 11,823 compounds → 20 hits with ΔG < −8.0 kcal/mol (three chemotypes: adamantane-carboxylic acids, barbiturates/glutarimides, NMDA-axis ligands).
- MD: 8 compounds × 100 ns, all stable (RMSD < 3 Å).
- Best head group: adamantane-carboxylic acid (ΔG −7.73, cryo-EM validated).
