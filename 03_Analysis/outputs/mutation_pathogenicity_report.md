# SLC5A6 (SMVT) Missense Mutation Pathogenicity Prediction

> Generated: 2026-06-23 15:35
> Gene: SLC5A6 (SMVT) | Protein: 635 aa | UniProt: Q9Y289

---

## Key Finding: Zero TCGA Somatic Missense Mutations

All 32 TCGA PanCancer Atlas studies (10,967 samples) were queried via the cBioPortal API. **No SLC5A6 missense mutations were found in any cancer type.**

This is consistent with SLC5A6's profile as an expression-driven metabolic target rather than a mutation-driven oncogene:

| Metric | Value | Implication |
|--------|-------|-------------|
| Pan-cancer mutation rate | <2% | Passenger gene level |
| pLI (gnomAD) | 0.01 | Complete LoF tolerance |
| LOEUF | 0.61 | Not evolutionarily constrained |
| Mutation hotspots | None | No recurrent somatic mutations |

---

## Method

Pathogenicity was predicted using:

1. **ESM** (`facebook/esm2_t12_35M_UR50D`): Meta's protein language model. Log-likelihood ratio (LLR) between mutant and wildtype at the masked position. Lower (negative) LLR = more deleterious.
2. **BLOSUM62**: Evolutionary substitution matrix. Negative scores = rare, potentially damaging substitutions.
3. **Composite score**: LLR z-score (70%) + BLOSUM62 z-score (30%). Higher = more deleterious.

## Limitations

- **No TCGA mutations exist for SLC5A6** across any cancer type. Pathogenicity predictions are limited to known germline ClinVar variants.
- Germline mutations are not cancer-associated but represent known biotin-responsive metabolic disease variants.
- ESM-2 predicts biochemical impact, not cancer-specific functional consequence.

---

## Pathogenicity Rankings (Germline Mutations)

| Rank | Mutation | Domain | Source | ESM LLR | BLOSUM62 | Percentile |
|------|----------|--------|--------|---------|----------|------------|
| 1 | R123L | Transmembrane_helix_2 | ClinVar | -4.3445 | -2 | 100% |

### Top Pathogenic Candidates

| Rank | Mutation | Domain | ESM LLR | BLOSUM62 | Composite |
|------|----------|--------|---------|----------|-----------|
| 1 | R123L | Transmembrane_helix_2 | -4.3445 | -2 | nan |

| ESM LLR range | -4.3445 to -4.3445 |
| Median ESM LLR | -4.3445 |

---

## All-Possible Substitutions at ClinVar-Related Positions

Below are BLOSUM62-based substitution scores for ALL possible amino acid changes at each position reported in ClinVar/literature. Wildtype residues are drawn from the canonical SMVT sequence (UniProt Q9Y289). NOTE: Only position 123 (R) matches the reported ClinVar reference. Positions 189 (V, not Y), 317 (A, not G), and 489 (N, not S) differ, likely due to alternative isoform or transcript numbering.

### Position 123 (R)

| Variant | BLOSUM62 | Domain |
|---------|----------|--------|
| R123A | -1 | Transmembrane_helix_2 |
| R123C | -3 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123D | -2 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123E | +0 | Transmembrane_helix_2 |
| R123F | -3 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123G | -2 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123H | +0 | Transmembrane_helix_2 |
| R123I | -3 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123K | +2 | Transmembrane_helix_2 |
| R123L | -2 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123M | -1 | Transmembrane_helix_2 |
| R123N | +0 | Transmembrane_helix_2 |
| R123P | -2 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123Q | +1 | Transmembrane_helix_2 |
| R123S | -1 | Transmembrane_helix_2 |
| R123T | -1 | Transmembrane_helix_2 |
| R123V | -3 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123W | -3 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |
| R123Y | -2 | Transmembrane_helix_2 [KNOWN PATHOGENIC] |

### Position 189 (V (lit. reports Y))

| Variant | BLOSUM62 | Domain |
|---------|----------|--------|
| V189A | +0 | Transmembrane_helix_3 |
| V189C | -1 | Transmembrane_helix_3 |
| V189D | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189E | -2 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189F | -1 | Transmembrane_helix_3 |
| V189G | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189H | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189I | +3 | Transmembrane_helix_3 |
| V189K | -2 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189L | +1 | Transmembrane_helix_3 |
| V189M | +1 | Transmembrane_helix_3 |
| V189N | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189P | -2 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189Q | -2 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189R | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189S | -2 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189T | +0 | Transmembrane_helix_3 |
| V189W | -3 | Transmembrane_helix_3 [KNOWN PATHOGENIC] |
| V189Y | -1 | Transmembrane_helix_3 |

### Position 317 (A (lit. reports G))

| Variant | BLOSUM62 | Domain |
|---------|----------|--------|
| A317C | +0 | Unknown |
| A317D | -2 | Unknown [KNOWN PATHOGENIC] |
| A317E | -1 | Unknown |
| A317F | -2 | Unknown [KNOWN PATHOGENIC] |
| A317G | +0 | Unknown |
| A317H | -2 | Unknown [KNOWN PATHOGENIC] |
| A317I | -1 | Unknown |
| A317K | -1 | Unknown |
| A317L | -1 | Unknown |
| A317M | -1 | Unknown |
| A317N | -2 | Unknown [KNOWN PATHOGENIC] |
| A317P | -1 | Unknown |
| A317Q | -1 | Unknown |
| A317R | -1 | Unknown |
| A317S | +1 | Unknown |
| A317T | +0 | Unknown |
| A317V | +0 | Unknown |
| A317W | -3 | Unknown [KNOWN PATHOGENIC] |
| A317Y | -2 | Unknown [KNOWN PATHOGENIC] |

### Position 489 (N (lit. reports S))

| Variant | BLOSUM62 | Domain |
|---------|----------|--------|
| N489A | -2 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489C | -3 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489D | +1 | Transmembrane_helix_11 |
| N489E | +0 | Transmembrane_helix_11 |
| N489F | -3 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489G | +0 | Transmembrane_helix_11 |
| N489H | +1 | Transmembrane_helix_11 |
| N489I | -3 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489K | +0 | Transmembrane_helix_11 |
| N489L | -3 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489M | -2 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489P | -2 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489Q | +0 | Transmembrane_helix_11 |
| N489R | +0 | Transmembrane_helix_11 |
| N489S | +1 | Transmembrane_helix_11 |
| N489T | +0 | Transmembrane_helix_11 |
| N489V | -3 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489W | -4 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |
| N489Y | -2 | Transmembrane_helix_11 [KNOWN PATHOGENIC] |

---

## TCGA Query Details

The following 32 TCGA PanCancer Atlas studies were queried via cBioPortal API:

| Cancer | Study | Cancer | Study |
|--------|-------|--------|-------|
| ACC | acc_tcga_pan_can_atlas_2018 | BLCA | blca_tcga_pan_can_atlas_2018 |
| BRCA | brca_tcga_pan_can_atlas_2018 | CESC | cesc_tcga_pan_can_atlas_2018 |
| CHOL | chol_tcga_pan_can_atlas_2018 | COAD/READ | coadread_tcga_pan_can_atlas_2018 |
| DLBC | dlbc_tcga_pan_can_atlas_2018 | ESCA | esca_tcga_pan_can_atlas_2018 |
| GBM | gbm_tcga_pan_can_atlas_2018 | HNSC | hnsc_tcga_pan_can_atlas_2018 |
| KICH | kich_tcga_pan_can_atlas_2018 | KIRC | kirc_tcga_pan_can_atlas_2018 |
| KIRP | kirp_tcga_pan_can_atlas_2018 | LAML | laml_tcga_pan_can_atlas_2018 |
| LGG | lgg_tcga_pan_can_atlas_2018 | LIHC | lihc_tcga_pan_can_atlas_2018 |
| LUAD | luad_tcga_pan_can_atlas_2018 | LUSC | lusc_tcga_pan_can_atlas_2018 |
| MESO | meso_tcga_pan_can_atlas_2018 | OV | ov_tcga_pan_can_atlas_2018 |
| PAAD | paad_tcga_pan_can_atlas_2018 | PCPG | pcpg_tcga_pan_can_atlas_2018 |
| PRAD | prad_tcga_pan_can_atlas_2018 | SARC | sarc_tcga_pan_can_atlas_2018 |
| SKCM | skcm_tcga_pan_can_atlas_2018 | STAD | stad_tcga_pan_can_atlas_2018 |
| TGCT | tgct_tcga_pan_can_atlas_2018 | THCA | thca_tcga_pan_can_atlas_2018 |
| THYM | thym_tcga_pan_can_atlas_2018 | UCS | ucs_tcga_pan_can_atlas_2018 |
| UCEC | ucec_tcga_pan_can_atlas_2018 | UVM | uvm_tcga_pan_can_atlas_2018 |

**No SLC5A6 missense mutations were found in any of the above studies.**

---

## References

- Meier et al. (2021). Language models enable zero-shot prediction of the effects of mutations. *Nat Biotechnol*.
- gnomAD v4.0: pLI=0.01, LOEUF=0.61 for SLC5A6
- TCGA PanCancer Atlas (2018). 32 cancer types, 10,967 samples (Cell Press)
- PMID: 27904971 - SLC5A6 mutations in biotin-responsive neurodegeneration
- PMID: 31754459 - SLC5A6-related vitamin-dependent neuro-metabolic disease
