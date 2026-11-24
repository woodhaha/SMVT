# Cover Letter — npj Systems Biology and Applications

**Date**: 2026-11-24
**To**: The Editors, npj Systems Biology and Applications
**Re**: Submission of "Structure-Based Discovery of Carboxyl-Bearing Ligands for the Multivitamin Transporter SLC5A6, a Tumor-Intrinsic Pan-Cancer Metabolic Vulnerability"

---

Dear Editor,

We are pleased to submit our manuscript for consideration as an original research article in *npj Systems Biology and Applications*.

## Summary

The sodium-dependent multivitamin transporter (SMVT/SLC5A6) sits at the interface of vitamin uptake and lipogenic metabolism, transporting biotin, pantothenate, and lipoate — cofactors required for fatty-acid synthesis and TCA-cycle function. Although the cryo-EM structure of SMVT was recently solved, the systems-level landscape of SLC5A6 and its druggability have remained unexplored. We address this gap with an integrated computational workflow spanning three scales:

1. **Genome scale** — a pan-cancer analysis of 10,967 TCGA samples (33 cancer types) integrating expression, survival, mutation constraint, and protein–protein interaction networks, together with a machine-learning pharmacophore screen (AUC = 0.888) over 3,311 approved drugs and an expanded 11,823-compound library.
2. **Structure scale** — structure-based virtual screening against the SMVT binding pocket, followed by 100-ns molecular dynamics simulations (including membrane-embedded Na⁺-bound MD of the carboxyl-bearing ligands) and cross-structure re-docking against the experimental cryo-EM conformations.
3. **Mechanism scale** — network-perturbation (virtual knockout), pathway-enrichment, and single-cell cell-type-localization analyses linking SMVT transport to the biotin–ACC–FASN lipogenic axis.

## Key Findings

1. **SLC5A6 is an expression-driven — not mutation-driven, and not cell-autonomously essential — nutrient-transport vulnerability.** DepMap CRISPR screens show SLC5A6 is not a common-essential gene (mean gene effect +0.015 across 1,150 lines), so we frame the vulnerability as nutrient-conditional and propose testing it in nutrient-defined medium rather than standard culture. SLC5A6 is significantly upregulated in 14 of 21 cancer types with matched normal tissue, yet harbors only passenger-level somatic mutations (1.5% of samples; no positive selection by constraint metrics or ESM-2 pathogenicity profiling). Single-cell analysis of a 370,115-cell colorectal cancer atlas (GSE178341) shows the transcript is epithelial-enriched (21.3% of epithelial cells positive vs. 0.9–6.1% across immune/stromal populations) and markedly higher in tumor than normal epithelium (28.3% vs. 6.95% positive), confirming the bulk upregulation is tumor-cell-intrinsic. High expression associates with worse overall survival in six cancer types (hazard ratios 1.8–4.0), remaining significant after multivariable adjustment.

2. **Structure-based screening identifies barbiturates and glutarimides as predicted SMVT ligands.** All eight barbiturate compounds in the screened library docked favourably (phenobarbital −8.30 kcal/mol vs. biotin −6.76), rationalized by the malonylurea core mimicking biotin's ureido ring. We frame these as *predicted* ligands: the barbiturate scaffold lacks the substrate carboxyl group required for SMVT-mediated transport, and recent experimental work cautions against over-interpreting this result.

3. **Carboxyl-bearing chemotypes — 5-hydroxytryptophan, NSAIDs, and adamantane carboxylic acids — are the more plausible ligands.** 5-HTP (−7.66 kcal/mol) docks in the biotin site, and multiple NSAIDs reproduce the experimentally reported SMVT inhibition at the chemotype level. We further designed a lead compound, 1-carboxy-3,5-bis(trifluoromethyl)adamantane (BISCF3), and validated a carboxyl-discrimination principle by membrane-embedded MM-GBSA: the two carboxylate ligands bind most favourably (biotin −56.7, BISCF3 −54.6 kcal/mol), followed by zwitterionic 5-HTP (−36.8) and neutral phenobarbital (−22.8) — a charge-driven gradient consistent with carboxyl recognition.

4. **Independent protein-level corroboration, and its limits.** CPTAC mass spectrometry of matched tumour/adjacent normal tissue reproduces the SLC5A6 tumour bias in colon adenocarcinoma (log2FC +1.05, p_adj = 5.6e-17) and lung squamous carcinoma (+0.50, p_adj = 8.4e-7) — two of the cancer types our TCGA analysis flags. We also report what the same resources do *not* show: SLC5A6 protein is quantified in 0 of 20 normal tissues by MS (versus 12–20 for the downstream axis enzymes) and never exceeds a "Medium" IHC call, so the normal-tissue baseline is constrained by a detection floor rather than confirmed. UMAP localisation of all 370,115 atlas cells (no subsampling) places the transcript in the epithelial compartment at 28.2% tumour versus 6.9% normal epithelial positivity.

5. **Cross-structure validation against cryo-EM.** Re-docking against the recently solved SMVT occluded-state structure confirms naftazone as the most conformation-robust ligand, and re-docking of biotin reproduced the experimental co-crystal pose with a heavy-atom RMSD of 1.33 Å (top-ranked pose), validating the docking protocol against experimental binding geometry.

## Fit for npj Systems Biology and Applications

This study documents a clear systems-level conclusion — that the multivitamin transporter SLC5A6 is a tumor-intrinsic pan-cancer metabolic vulnerability whose druggable site discriminates carboxyl-bearing from non-carboxyl ligands — reached strictly through computational methods that are themselves a contribution to the field. The work integrates three complementary computational disciplines aligned with the journal's scope on **systems biology and drug discovery**: machine learning (a pharmacophore screen prioritizing 500 compounds from a 3,311-drug library), bioinformatics (pan-cancer expression, survival, mutation-constraint, and single-cell analysis of 10,967 TCGA samples plus a 370,115-cell atlas), and molecular simulation (AutoDock Vina screening plus 100-ns OpenMM MD with MM-GBSA). A key methodological advance is our benchmarking of the AlphaFold-based screen against the experimental cryo-EM conformations, which quantifies the receptor-conformation dependence of docking scores — a result of general relevance to any structure-based virtual-screening campaign. The pipeline is fully reproducible from public resources (TCGA/gnomAD/ChEMBL/AlphaFold/RCSB) and yields concrete, falsifiable hypotheses (barbiturate/glutarimide scaffolds for medicinal chemistry and 5-HTP competition-uptake assays).

## Transparency and Availability

All analyses use public data (TCGA, gnomAD, STRING, ChEMBL, AlphaFold, RCSB, GEO). Analysis code is publicly available at https://github.com/woodhaha/SMVT. The manuscript acknowledges its limitations explicitly: all ligand predictions are computational and pending experimental validation. We include membrane-embedded, Na⁺-bound MD of the carboxyl-bearing ligands (biotin, phenobarbital, 5-HTP), which confirms stable substrate binding. As the immediate next step, we propose biotin-competition uptake assays in SLC5A6-overexpressing cells followed by SPR/ITC to test the top carboxyl-bearing candidates.

We confirm this manuscript is not under consideration elsewhere and that all authors have approved the submission. No competing interests. During the preparation of this work, we used Claude (Anthropic), an AI-assisted tool, for drafting portions of the text and generating analysis code; all content was reviewed and verified by the authors, who take full responsibility for the work.

Thank you for your consideration.

Sincerely,

Zhuha Zhou (Corresponding author)
The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, China
zhouzhuha@wmu.edu.cn
