# Data Sources

All data are publicly available. No new data were generated for this study.

| Resource | Use | Access |
|----------|-----|--------|
| TCGA (NCI Genomic Data Commons) | RNA-seq (TOIL FPKM) + clinical + MC3 somatic mutations (33 cancer types, 10,967 samples) | https://portal.gdc.cancer.gov |
| GEO GSE178341 | Colorectal cancer scRNA-seq (371,000 cells) for virtual-KO | https://www.ncbi.nlm.nih.gov/geo/ |
| gnomAD v4 | pLI / LOEUF / %HI constraint metrics | https://gnomad.broadinstitute.org |
| STRING v12 | Protein–protein interaction network | https://string-db.org |
| ChEMBL v34 | 3,311 approved drugs (initial screening set) | https://www.ebi.ac.uk/chembl |
| Drug Repurposing Hub | 7,423 compounds (expanded screening) | https://repo-hub.broadinstitute.org |
| Selleck L3000 | 4,247 compounds (expanded screening) | https://www.selleckchem.com |
| WHO ATC L01 | Antineoplastic agent list | https://www.whocc.no |
| NCI oncology library | ~3,500 compounds (expanded screening) | https://dtp.cancer.gov |
| AlphaFold AF-Q9Y289-F1 | SMVT structure for docking | https://alphafold.ebi.ac.uk/entry/Q9Y289 |
| RCSB PDB 26va–26ve | Cryo-EM SMVT structures (occluded/outward-open) | https://www.rcsb.org |
| Hugging Face esm2_t12_35M_UR50D / esm2_t33_650M_UR50D | ESM-2 pathogenicity scoring | https://huggingface.co/facebook |

Compound libraries (SDF/annotation files) are large and re-downloadable from the sources above; curated candidate CSVs are under `03_Analysis/data/`.
