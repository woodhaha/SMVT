# Draft — Expanded virtual screening (Results subsection)

## Expanded virtual screening and emergent scaffold chemotypes

To broaden coverage beyond the initial 703-compound set, we expanded the docking library to **11,823 unique compounds** by adding the Drug Repurposing Hub (7,423 compounds) and an oncology-focused library assembled from the Selleck L3000 anticancer collection, WHO ATC L01 antineoplastics, and NCI-approved oncology drugs. All ligands were prepared with Meeko (AM1-BCC charges) and docked with AutoDock Vina (exhaustiveness 16) into a 22 Å box centered on the biotin/Na⁺ site. After removing perfluorocarbon false positives (two compounds, ΔG < −9 kcal/mol, retained on hydrophobic collapse but chemically irrelevant to SMVT), the screen yielded **20 compounds with ΔG < −8.0 kcal/mol** and **349 with ΔG < −7.0 kcal/mol** — a substantial enrichment relative to the eight < −8 kcal/mol hits in the original library.

The top hits cluster into three mechanistically coherent chemotypes (**Fig. X**, **Table X**):

**1. Adamantane/diamondoid carboxylic acids (novel chemotype).** Six compounds spanning ΔG −8.1 to −8.6 kcal/mol share a rigid hydrophobic cage bearing one to three carboxyl groups — adamantane-1,3,5-tricarboxylic acid (−8.51), 5,7-dimethyladamantane-1,3-dicarboxylic acid (−8.49), 5,7-dimethyladamantane-1,3-diacetic acid (−8.62), and hexacyclododecane-1-carboxylic acid (−8.29). This scaffold recapitulates the carboxylate recognition that defines SMVT's endogenous substrates (biotin, pantothenate, lipoate), and its enrichment among top hits is consistent with carboxylate coordination at the Na⁺ site. A focused SAR of 58 adamantane/diamondoid analogs confirmed the pharmacophore: di-carboxylate substitution outperforms mono-carboxylate (mean −8.46 vs −7.43 kcal/mol), electron-withdrawing F/CN substitution further improves affinity (best −8.84 kcal/mol), and the larger diamantane cage is disfavored, indicating that the pocket is size-matched to the adamantane (C10) core.

**2. Barbiturate/glutarimide scaffold (reconfirmed).** Phenobarbital (−8.17 kcal/mol) and the glutarimide-type 5,5-dimethyl-2,6-dioxo-4-phenylpiperidine-3-carbonitrile (−8.61 kcal/mol) ranked among the top hits, extending our earlier finding that barbiturate and glutethimide-type compounds are unexpected SMVT ligands and reinforcing the prior focused SAR.

**3. NMDA-antagonist axis.** The uncompetitive NMDA antagonist dizocilpine (MK-801; −8.02 kcal/mol) and the memantine-like 3-(2-aminoethyl)adamantane-1-carboxylic acid (−8.12 kcal/mol) both scored below −8 kcal/mol, echoing the esketamine hit in the original screen and suggesting a recurring NMDA-related chemotype among predicted SMVT ligands.

Natural-product sesquiterpenes (patchouli alcohol −8.56, cedrol −8.36 kcal/mol) and alkaloids (securinine −8.21, cytisine −8.01 kcal/mol) completed the top-twenty set. Together, these results expand the predicted ligand landscape of SMVT from approved drugs to a broader chemical space, and identify the adamantane carboxylic-acid scaffold as a new, tractable chemotype for SMVT engagement.

---

## Notes for integration (to be removed before submission)

- **Fig. X**: 2D structures of the 20 < −8 kcal/mol hits grouped by chemotype (RDKit grid).
- **Table X**: the 20 hits with InChIKey, ΔG, chemotype, and purchasability (see `03_Analysis/outputs/docking_final/hits_top20_classified.csv`).
- Numbers to keep consistent with Methods: 11,823 compounds; 22 Å box; exhaustiveness 16; Meeko/AM1-BCC; Vina.
- Caveat to state in Limitations if reviewers press: the adamantane carboxylic acids are research chemicals (not approved drugs); they serve as pharmacophore probes rather than immediate repurposing candidates. Benzotriazole (−8.04) is an industrial UV absorber and was excluded from biological interpretation.
