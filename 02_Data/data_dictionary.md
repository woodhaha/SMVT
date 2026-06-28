# Data Dictionary — SMVT Project

> 最后更新: 2026-06-23

---

## 1. TCGA Pan-Cancer Expression

| Variable | Type | Source | Description |
|----------|------|--------|-------------|
| `project` | character | TCGAbiolinks | TCGA project ID (e.g., TCGA-BRCA) |
| `cancer_type` | character | TCGAbiolinks | Cancer type abbreviation |
| `sample` | character | TCGAbiolinks | Sample barcode |
| `sample_type` | character | TCGAbiolinks | Tumor / Normal / Metastatic |
| `SLC5A6_expr` | numeric | TCGAbiolinks | SLC5A6 expression (log2(FPKM+1)) |
| `expression_category` | character | derived | High / Low (median split) |

---

## 2. TCGA Mutation Data

| Variable | Type | Source | Description |
|----------|------|--------|-------------|
| `sample` | character | TCGAbiolinks | Sample barcode |
| `mutation_type` | character | TCGA MAF | Missense_Mutation / Silent / etc. |
| `protein_change` | character | TCGA MAF | Amino acid change (e.g., R142H) |
| `position` | integer | TCGA MAF | Amino acid position in SMVT |
| `domain` | character | UniProt | Transmembrane / Loop / N-term / C-term |
| `is_mutant` | binary | derived | 0 = WT, 1 = mutant |

---

## 3. STRING Network

| Variable | Type | Source | Description |
|----------|------|--------|-------------|
| `node_gene` | character | STRING v12 | Gene symbol |
| `interaction_score` | numeric | STRING v12 | Combined score (0–1) |
| `evidence_type` | character | STRING v12 | Experiment / Database / Coexpression / etc. |

---

## 4. Enrichment (pathlinkR)

| File | Source | Description |
|------|--------|-------------|
| `SMVT_GO_enrichment.csv` | pathlinkR | Gene Ontology enrichment |
| `SMVT_KEGG_enrichment.csv` | pathlinkR | KEGG pathway enrichment |
| `SMVT_Reactome_enrichment.csv` | pathlinkR | Reactome pathway enrichment |

---

## 5. Structural Data

| File | Type | Source | Description |
|------|------|--------|-------------|
| `AF-Q9Y289-F1.pdb` | PDB | AlphaFold DB | Raw structure (read-only) |
| `AF-Q9Y289-F1_cleaned.pdb` | PDB | derived | Hydrogens removed, non-standard residues stripped |
| `AF-Q9Y289-F1_prepared.pdb` | PDB | derived | Prepared for docking (H added, charges assigned) |
| `AF-Q9Y289-F1_H.pdb` | PDB | derived | Hydrogens added |

---

## Column Value Encodings

| Variable | Encoding |
|----------|----------|
| `sample_type` | `Primary Tumor`, `Solid Tissue Normal`, `Metastatic` |
| `mutation_type` | `Missense_Mutation`, `Silent`, `Nonsense_Mutation`, `Splice_Site` |
| `expression_category` | `High` (>median), `Low` (≤median) |
| `domain` | `TM1-12` (transmembrane helix), `N-term`, `C-term`, `EL1-6` (extracellular loop), `IL1-5` (intracellular loop) |
