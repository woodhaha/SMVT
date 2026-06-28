# SMVT Top-3 Hit MD Simulation — Colab Package

## Files

| File | Description |
|------|-------------|
| `AF-Q9Y289-F1.pdb` | AlphaFold SMVT receptor structure |
| `NAFTAZONE_docked.pdbqt` | Top hit (-8.34 kcal/mol) Vina pose |
| `PHENOBARBITAL_docked.pdbqt` | #2 hit (-8.30 kcal/mol) Vina pose |
| `ESKETAMINE_docked.pdbqt` | Novel hit (-7.58 kcal/mol) Vina pose |
| `colab_md_top3.ipynb` | Colab GPU MD notebook |

## Instructions

### 1. Upload to Google Drive
- Go to [drive.google.com](https://drive.google.com)
- Create folder `SMVT_MD` in root
- Upload ALL 5 files above into `SMVT_MD/`

### 2. Run Colab Notebook
- Go to [colab.research.google.com](https://colab.research.google.com)
- Upload `colab_md_top3.ipynb`
- Runtime → Change runtime type → T4 GPU
- Run cells top to bottom

### 3. Expected Output
- Each compound: 50ns production MD with TIP3P water
- Trajectories (DCD format), energy logs (CSV), final frames (PDB)
- Saved to `SMVT_MD/trajectories/` on your Drive
- Approximate runtime: 2-4 hours per compound on T4

### 4. After MD Completes
- Download trajectories back to `03_Analysis/md/trajectories/`
- Run `analyze_md_results.py` locally for RMSD/RMSF/contact analysis

## Top 3 Compounds

| Compound | Docking ΔG | Type | Clinical Use |
|----------|:---:|------|-------------|
| Naftazone | -8.34 | Naphthoquinone | Hemostatic |
| Phenobarbital | -8.30 | Barbiturate | Anticonvulsant (WHO) |
| Esketamine | -7.58 | Arylcyclohexylamine | Antidepressant (Spravato) |

## Methods Summary
- Force field: AMBER ff14SB (protein) + OpenFF Sage 2.1.0 (ligand) + TIP3P (water)
- Protocol: Minim 5000 → NVT 50→300K → NPT 100ps → Production 50ns
- Timestep: 2fs with SHAKE constraints
- GPU: NVIDIA CUDA (T4/V100/A100)
