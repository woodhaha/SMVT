# Cover Letter — PLOS Computational Biology

**Date**: 2026-11-08
**To**: The Editors, PLOS Computational Biology
**Re**: Submission of "Structure-Based Discovery of Carboxyl-Bearing Ligands for the Multivitamin Transporter SLC5A6, a Pan-Cancer Metabolic Vulnerability"

---

Dear Editor,

We are pleased to submit our manuscript for consideration as a research article in *PLOS Computational Biology*.

## Summary

The sodium-dependent multivitamin transporter (SMVT/SLC5A6) sits at the interface of vitamin uptake and lipogenic metabolism, transporting biotin, pantothenate, and lipoate — cofactors required for fatty-acid synthesis and TCA-cycle function. Although the cryo-EM structure of SMVT was recently solved, the clinical landscape of SLC5A6 and its druggability have remained unexplored. We address this gap with an integrated computational workflow spanning three scales:

1. **Genome scale** — a pan-cancer analysis of 10,967 TCGA samples (33 cancer types) integrating expression, survival, mutation constraint, and protein–protein interaction networks, with a machine-learning pharmacophore screen (AUC = 0.888) over 3,311 approved drugs and an expanded 11,823-compound library.
2. **Structure scale** — structure-based virtual screening against the SMVT binding pocket, followed by 100-ns molecular dynamics simulations and cross-structure re-docking against the experimental cryo-EM conformations.
3. **Mechanism scale** — a network-perturbation and pathway-enrichment analysis linking SMVT transport to the biotin–ACC–SREBP–FASN lipogenic axis.

## Key Findings

1. **SLC5A6 is an expression-driven — not mutation-driven — metabolic vulnerability.** SLC5A6 is significantly upregulated in 14 of 21 cancer types with matched normal tissue, yet harbors only passenger-level somatic mutations (1.5% of samples; no positive selection by constraint metrics or ESM-2 pathogenicity profiling). High expression associates with worse overall survival in six cancer types (hazard ratios 1.8–4.0), remaining significant after multivariable adjustment for age and sex (and stage, where available).

2. **Structure-based screening identifies barbiturates and glutarimides as predicted SMVT ligands.** All eight barbiturate compounds in the screened library docked favourably (phenobarbital −8.30 kcal/mol vs. biotin −6.76), rationalized by the malonylurea core mimicking biotin's ureido ring. An expanded glutarimide screen reaches −9.36 kcal/mol. We frame these as *predicted* ligands: the barbiturate scaffold lacks the substrate carboxyl group, and recent experimental work (biotin conjugates lacking the carboxyl are not transported) cautions against over-interpreting this result.

3. **Carboxyl-bearing chemotypes — 5-hydroxytryptophan, NSAIDs, and adamantane carboxylic acids — are the more plausible ligands.** 5-HTP (−7.66 kcal/mol) docks in the biotin site (pose overlap 0.37 Å; stable over an 11.2-ns MD window), and multiple NSAIDs reproduce the experimentally reported SMVT inhibition at the chemotype level. An expanded 11,823-compound screen identifies adamantane/diamondoid carboxylic acids (best −8.84 kcal/mol) as a novel carboxylate pharmacophore. These findings motivate a testable SMVT–serotonin metabolic-axis hypothesis.

4. **Cross-structure validation against cryo-EM.** Re-docking against the recently solved SMVT structures (occluded/outward-open/inward-open) confirms naftazone as the most conformation-robust ligand and reveals that hydromorphone's score is highly receptor-conformation dependent — a caution that generalizes to AF2-based screening.

## Fit for PLOS Computational Biology

This study documents a clear biological conclusion — that the multivitamin transporter SLC5A6 is a pan-cancer metabolic vulnerability whose druggable site discriminates carboxyl-bearing from non-carboxyl ligands — reached strictly through computational methods that are themselves a contribution to the field. The work integrates three complementary computational disciplines: **machine learning** (a pharmacophore screen that prioritizes 500 compounds from a 3,311-drug library), **bioinformatics** (pan-cancer expression, survival, and mutation-constraint analysis of 10,967 TCGA samples), and **molecular simulation** (AutoDock Vina screening plus 100-ns OpenMM MD with MM-GBSA). A key methodological advance is our benchmarking of the AlphaFold-based screen against the experimental cryo-EM conformations, which quantifies the receptor-conformation dependence of docking scores — a result of general relevance to any structure-based virtual-screening campaign. The pipeline is fully reproducible from public resources (TCGA/gnomAD/ChEMBL/AlphaFold/RCSB) and yields concrete, falsifiable hypotheses (barbiturate/glutarimide scaffolds for medicinal chemistry and 5-HTP competition-uptake assays).

## Transparency and Availability

All analyses use public data (TCGA, gnomAD, STRING, ChEMBL, AlphaFold, RCSB). Analysis code will be deposited in a public repository (GitHub/Zenodo) and is available to reviewers upon request. The manuscript acknowledges its limitations explicitly: all ligand predictions are computational and pending experimental validation, and the MD was performed in a simplified water-box setup without the Na⁺ co-substrate or a lipid bilayer. As the immediate next step, we propose biotin-competition uptake assays in SLC5A6-overexpressing cells followed by SPR/ITC to test the top carboxyl-bearing candidates (5-HTP, NSAIDs, and adamantane carboxylic acids), and evaluation of adamantane-carboxylate-functionalised nanoparticles for SMVT-targeted uptake in SLC5A6-overexpressing versus knockdown cells.

We confirm this manuscript is not under consideration elsewhere and that all authors have approved the submission. No competing interests.

Thank you for your consideration.

Sincerely,

Zhuha Zhou (Corresponding author)
The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, China
zhouzhuha@wmu.edu.cn
