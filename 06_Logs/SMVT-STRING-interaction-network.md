# SLC5A6 (SMVT) — STRING Protein-Protein Interaction Network

> Query: 2026-06-23 | Database: STRING v12 | Species: Homo sapiens (9606)

---

## 1. Network Overview

SLC5A6 interacts with **10 proteins** at medium confidence (combined score ≥ 0.4).
Only **1** high-confidence partner (≥ 0.7).

### Interaction Partners (ranked by combined score)

| # | Gene | Protein Name | Score | Evidence Channels | Functional Link |
|---|------|-------------|-------|-------------------|-----------------|
| 1 | **PDZD11** | PDZ domain-containing 11 | **0.969** | DB(0.900), TXT(0.609), EXP(0.292) | Apical membrane scaffold |
| 2 | **HLCS** | Holocarboxylase synthetase | 0.635 | TXT(0.625), COE(0.067), NGH(0.042) | Biotin metabolism |
| 3 | **SLC5A7** | Na+/choline cotransporter | 0.655 | TXT(0.581), PHY(0.212) | Solute carrier family member |
| 4 | **SLC22A12** | URAT1 — urate transporter | 0.631 | TXT(0.629), EXP(0.045), COE(0.042) | Renal solute transport |
| 5 | **SLC26A4** | Pendrin | 0.578 | TXT(0.574), COE(0.050) | Anion transporter |
| 6 | **BTD** | Biotinidase | 0.555 | TXT(0.551), COE(0.049) | Biotin recycling |
| 7 | **SLC5A3** | Na+/myo-inositol cotransporter | 0.511 | TXT(0.316), PHY(0.228), EXP(0.137), COE(0.054) | Solute carrier family |
| 8 | **SLC23A1** | Na+/vitamin C cotransporter 1 | 0.493 | TXT(0.454), COE(0.110), FUS(0.003) | Vitamin transport |
| 9 | **DPH2** | Diphthamide biosynthesis 2 | 0.471 | COE(0.471) | Co-expression only |
| 10 | **SLC19A2** | Thiamine transporter 1 | 0.469 | TXT(0.459), COE(0.059) | Vitamin B1 transport |

> Evidence codes: EXP=Experimental, DB=Database(annotated), TXT=Text mining, COE=Co-expression, PHY=Phylogenetic co-occurrence, NGH=Neighborhood(gene), FUS=Gene fusion

---

## 2. Network Topology

```
                    ┌── SLC5A7 (choline transporter)
                    │
        ┌── SLC5A3 (inositol transporter)
        │           │
        │           └── SLC22A12 (URAT1, urate)
        │
SLC5A6 ──★── PDZD11 (scaffold, score=0.969)  ← ANCHOR
  (SMVT) │
        ├── HLCS ── BTD (biotin metabolism hub)
        │
        ├── SLC23A1 (vitamin C) ── SLC19A2 (vitamin B1)
        │
        ├── SLC26A4 (pendrin, anion)
        │
        └── DPH2 (co-expression only, unclear functional link)
```

**Network characteristics**:
- **Hub-bottleneck**: PDZD11 is the sole high-confidence anchor — an apical membrane scaffold that may organize SMVT localization
- **Biotin module**: HLCS (biotin ligase) + BTD (biotin recycling) form a metabolic sub-network
- **Solute carrier cluster**: 7 of 10 partners are SLC family transporters
- **Co-expression dominant**: Most edges rely on text mining + co-expression, few on direct experimental evidence

---

## 3. Functional Module Analysis

### Module A: Biotin / Vitamin Metabolism
| Gene | Function | Link to SMVT |
|------|----------|-------------|
| **HLCS** | Biotinylates histones + carboxylases; feedback-regulates SLC5A6 promoter | Direct regulatory loop |
| **BTD** | Liberates biotin from dietary proteins; recycles endogenous biotin | Shared substrate pool |

**Key regulatory circuit**:  
HLCS senses biotin → biotinylates H4K12 at SLC5A6 promoter → silences SLC5A6. When biotin is low, this repression is relieved → SLC5A6 upregulated. This is the only known transcriptional feedback loop for SLC5A6.

### Module B: Solute Carrier Co-Transport Systems
| Gene | Substrate | Co-transport ion | Tissue |
|------|----------|-----------------|--------|
| SLC5A6 | Biotin, pantothenate, lipoate | Na+ (2:1) | Intestine, kidney, placenta |
| SLC5A7 | Choline | Na+, Cl- | Neurons |
| SLC5A3 | Myo-inositol | Na+ (2:1) | Kidney, brain |
| SLC22A12 | Urate | — (exchanger) | Kidney |
| SLC23A1 | Vitamin C (ascorbate) | Na+ (2:1) | Intestine, kidney |
| SLC19A2 | Thiamine (B1) | — (facilitated) | Ubiquitous |
| SLC26A4 | I-, Cl-, HCO3- | — (exchanger) | Thyroid, inner ear |
| SLC5A7 | Choline | Na+, Cl- | Cholinergic neurons |

**Pattern**: SMVT clusters with other Na+-coupled vitamin/cofactor transporters, suggesting co-regulation at the transcriptional level and potential functional redundancy in vitamin uptake across epithelia.

### Module C: PDZD11 — The Anchor
PDZD11 (PDZ domain-containing 11) is an apical membrane scaffold protein that:
- Localizes to the brush border of polarized epithelia
- Contains a PDZ domain that binds C-terminal motifs of transmembrane proteins
- SLC5A6 C-terminus: `...S-E-R-T-L` — **class I PDZ-binding motif** (x-S/T-x-Φ)
- **Implication**: PDZD11 likely anchors SMVT at the apical membrane of intestinal/renal epithelial cells

---

## 4. Cancer Relevance of the Network

| Gene | TCGA Expression | Cancer Link |
|------|----------------|-------------|
| SLC5A6 | Up in LUAD/LUSC/COAD/BLCA/STAD/ESCA | SMVT→biotin→ACC→FASN axis |
| PDZD11 | Variable | Scaffold — could be rate-limiting for SMVT membrane localization |
| HLCS | Down in some cancers | Biotin-sensing broken → SLC5A6 de-repressed |
| SLC22A12 | Up in renal cancer | URAT1 — shared renal expression with SLC5A6 |
| SLC23A1 | Variable | Vitamin C uptake — potential antioxidant role in tumor |
| SLC19A2 | Up in breast cancer | Thiamine — pentose phosphate pathway → nucleotide synthesis |

---

## 5. Experimental Validation Status

| Interaction | Evidence Level | Suggested Experiment |
|-------------|---------------|---------------------|
| SLC5A6—PDZD11 | **Strongest** (DB + EXP + TXT) | Co-IP / PLA in Caco-2 or HEK293 |
| SLC5A6—HLCS | Medium (TXT dominant) | ChIP for H4K12bio at SLC5A6 promoter |
| SLC5A6—SLC22A12 | Weak (TXT dominant) | Co-expression validation in kidney IHC |
| Others | Text mining only | Need co-expression validation |

**Bottom line**: Only the SLC5A6—PDZD11 interaction has database-level support. Most other links derive from text mining of co-mentioned genes in literature abstracts. The HLCS regulatory loop is mechanistically the most interesting but lacks direct PPI evidence (it's a transcriptional interaction, not a physical binding).

---

## 6. Data Source

- STRING v12: `https://string-db.org/`  
- Query: `SLC5A6` + `species=9606` + `required_score=400`
- High-confidence subset: `required_score=700` → PDZD11 only
