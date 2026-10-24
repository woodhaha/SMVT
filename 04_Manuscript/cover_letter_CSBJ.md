# Cover Letter — Computational and Structural Biotechnology Journal (CSBJ)

**Date**: 2026-10-24
**To**: The Editors, Computational and Structural Biotechnology Journal
**Re**: Submission of "Structure-Based Virtual Screening of FDA-Approved Drugs and Molecular Dynamics Simulation Identifies Predicted Ligands of the Pan-Cancer Multivitamin Transporter SLC5A6 (SMVT)"

---

Dear Editor,

We are pleased to submit our manuscript for consideration as a research article in *Computational and Structural Biotechnology Journal*.

## Summary

The sodium-dependent multivitamin transporter (SMVT/SLC5A6) sits at the interface of vitamin uptake and lipogenic metabolism, transporting biotin, pantothenate, and lipoate — cofactors required for fatty-acid synthesis and TCA-cycle function. Although the cryo-EM structure of SMVT was recently solved, the clinical landscape of SLC5A6 and its druggability have remained unexplored. We address this gap with an integrated computational workflow spanning three scales:

1. **Genome scale** — a pan-cancer analysis of 10,967 TCGA samples (33 cancer types) integrating expression, survival, mutation constraint, and protein–protein interaction networks.
2. **Structure scale** — structure-based virtual screening of 3,311 approved drugs against the SMVT binding pocket, followed by 100-ns molecular dynamics simulations and cross-structure re-docking against the experimental cryo-EM conformations.
3. **Mechanism scale** — a network-perturbation and pathway-enrichment analysis linking SMVT transport to the biotin–ACC–SREBP–FASN lipogenic axis.

## Key Findings

1. **SLC5A6 is an expression-driven — not mutation-driven — metabolic vulnerability.** SLC5A6 is significantly upregulated in 14 of 21 cancer types with matched normal tissue, yet harbors only passenger-level somatic mutations (1.5% of samples; no positive selection by constraint metrics or ESM-2 pathogenicity profiling). High expression associates with worse overall survival in six cancer types (hazard ratios 1.8–4.0), and these associations **remain significant after multivariable adjustment for age, sex, and stage** — an analysis we emphasize for its robustness.

2. **Structure-based screening identifies barbiturates and glutarimides as predicted SMVT ligands.** Barbiturates show a 100% hit rate in docking (phenobarbital −8.30 kcal/mol vs. biotin −6.76), rationalized by the malonylurea core mimicking biotin's ureido ring. An expanded glutarimide screen reaches −9.36 kcal/mol, with a tentative CF₃-substitution structure–activity trend. We frame these as *predicted* ligands: the barbiturate scaffold lacks the substrate carboxyl group, and we explicitly note that recent experimental work (biotin conjugates lacking the carboxyl are not transported) cautions against over-interpreting this result.

3. **Carboxyl-bearing chemotypes — 5-hydroxytryptophan and NSAIDs — are the more plausible ligands.** 5-HTP (−7.66 kcal/mol) docks in the biotin site (pose overlap 0.37 Å; stable over an 11.2-ns MD window), and multiple NSAIDs reproduce the experimentally reported SMVT inhibition at the chemotype level. These findings motivate a testable SMVT–serotonin metabolic-axis hypothesis.

4. **Cross-structure validation against cryo-EM.** Re-docking against the recently solved SMVT structures (occluded/outward-open/inward-open) confirms naftazone as the most conformation-robust ligand and reveals that hydromorphone's score is highly receptor-conformation dependent — a caution that generalizes to AF2-based screening.

## Fit for CSBJ

This study aligns closely with the journal's emphasis on *computational chemistry and drug discovery* and on *functional/mechanistic understanding through computational methods*. The work is fully computational, delivers a reproducible docking–MD pipeline, and produces concrete, testable hypotheses (barbiturate/glutarimide scaffolds for medicinal chemistry; 5-HTP competition-uptake experiments). No structural data are required to contextualize our findings, but we integrate the available cryo-EM structures to benchmark our AF2-based screen — a methodological contribution of general interest.

## Transparency and Availability

All analyses use public data (TCGA, gnomAD, STRING, ChEMBL, AlphaFold, RCSB). Analysis code will be deposited in a public repository (Zenodo/Figshare) upon acceptance and is available to reviewers upon request. The manuscript acknowledges its limitations explicitly: all ligand predictions are computational and pending experimental validation, and the MD was performed in a simplified water-box setup without the Na⁺ co-substrate or a lipid bilayer.

We confirm this manuscript is not under consideration elsewhere and that all authors have approved the submission. No competing interests.

Thank you for your consideration.

Sincerely,

Zhouzhuha (Woodhaha)
Wenzhou Medical University
zhouzhuha@wmu.edu.cn
