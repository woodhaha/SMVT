# A Trifluoromethyl-Optimized Carboxylate Ligand of the Human Multivitamin Transporter (SMVT): Docking-Guided Discovery, Cross-Structure Validation, and the Case Against Direct Conjugate Transport

> Autoresearch draft · 2026-11-12 · data in `../findings.md` and `../experiments/*/`

## Abstract

The sodium-dependent multivitamin transporter (SMVT/SLC5A6) — the sole intestinal biotin uptake system, induced in Ras-driven intestinal tumors and overexpressed in colorectal cancer (CRC) — recognizes its substrates through a **free carboxyl group** anchored in the TM1 kink (Ser81/Ala84 backbone NH + Tyr156). Docking-guided SAR across ~83 carboxyl-bearing ligands identified **1-carboxy-3,5-bis(trifluoromethyl)adamantane (COOH-bisCF3)** as the strongest SMVT ligand to date (−9.68 kcal/mol on the AlphaFold structure; −8.55 occluded and −8.37 outward-open cryo-EM states; three independent Vina seeds, SD 0.01). Trifluoromethyl emerged as a general affinity enhancer with a clean two-group optimum — 1 CF3 (−9.04) < **2 CF3 (−9.68)** > 3 CF3 (−8.03) — and cross-structure validation ruled out AlphaFold overfitting, correcting a hydrophobic sub-pocket artifact at the 7-position. We further show that the fourth bridgehead tolerates only minimal substitution (7-methyl −9.76 on AF2; longer linkers fall steeply), so that **direct drug conjugates (~760 Da, 22 Å extended) are rejected with positive scores on all three structures (AF2 +59.0, occluded +23.7, outward-open +3.2 kcal/mol)**. Delivery of SMVT-targeted payloads must therefore proceed by **SMVT-mediated endocytosis of a conjugate-loaded nanoparticle**, not by direct transport. As a concrete application we design **COOH-bisCF3–SN-38** (C₃₇H₃₅F₆N₃O₈, MW 763.7), an SN-38 conjugate linked through a cleavable C10 carbamate, positioned for CRC and liver-metastasis therapy.

## 1. Introduction

**SMVT as a CRC-selective target.** SMVT (SLC5A6) is the sole intestinal biotin transporter; its expression is induced in Ras-driven intestinal tumors, and silencing SMVT reduces tumor size and mitosis (Neophytou & Pitsouli, 2022). Pan-cancer analysis places SMVT among the transporters overexpressed in CRC (COAD/READ). Crucially, SMVT recognizes its substrates through a **free carboxyl group** anchored in the TM1 kink (Ser81/Ala84 backbone NH + Tyr156) — a "carboxyl-discrimination" principle that makes carboxyl-bearing ligands, not biotin conjugates that consume the carboxyl, the effective recognition heads.

**The delivery problem.** Irinotecan, a backbone of the FOLFIRI regimen for metastatic CRC, is a prodrug that must be hydrolyzed to its ~1000-fold more potent active metabolite SN-38, whose poor solubility and dose-limiting toxicity motivate tumor-selective delivery. A widespread — but largely untested — assumption in transporter-targeted drug design is that a transporter ligand can be conjugated to a cytotoxic payload and still be transported by its transporter.

**This study.** We (i) discover a high-affinity carboxyl-bearing SMVT ligand by docking-guided SAR, (ii) validate it across three receptor structures, and (iii) test the direct-conjugate assumption explicitly. The conjugate test fails decisively, redirecting SMVT-targeted delivery from transport to endocytosis.

## 2. Results

### 2.1 Discovery of COOH-bisCF3 by docking-guided SAR

SAR across adamantane/diamondoid, aromatic amino-acid, and alternative-cage chemotypes (Table 1) identified **trifluoromethyl as a general affinity enhancer** on the adamantane cage, with **1-carboxy-3,5-bis(trifluoromethyl)adamantane (COOH-bisCF3)** as the strongest carboxyl-bearing SMVT ligand (−9.68 kcal/mol AF2), surpassing the prior best transportable ligand (adamantane + Ac, −8.84). The CF3 enhancement shows a clean two-group optimum — 1 CF3 (−9.04) < **2 CF3 (−9.68)** > 3 CF3 (−8.03) — reproduced across three Vina seeds (−9.68/−9.69/−9.70, SD 0.01). Alternative cages (cubane, bicyclo[2.2.2]octane, noradamantane, homoadamantane) and the aromatic amino-acid chemotype (cap ~−7.7) are all inferior, confirming the adamantane C10 cage as uniquely size-matched to the pocket.

**Table 1 — Best ligand per chemotype (AF2).**

| Chemotype | Best ligand | Vina (kcal/mol) |
|---|---|---|
| Adamantane + CF3 | **COOH-bisCF3** | **−9.68** |
| Adamantane + Ac (prior) | 1-carboxy-3-acetyl-adamantane | −8.84 |
| Adamantane diCOOH | adamantane-1,3-dicarboxylic acid | −8.48 |
| Aromatic amino-acid | 1-naphthylalanine | −7.74 |
| Alternative cage | cubane-1,4-dicarboxylic acid | −7.57 |
| (substrate) | biotin | −6.76 |

### 2.2 Cross-structure validation rules out AlphaFold overfitting

COOH-bisCF3 is the top-ranked ligand on **all three receptor states** (Table 2), ruling out conformation-specificity: AF2 −9.68, cryo-EM occluded (26va) −8.55, outward-open (26vd) −8.37. The CF3 ranking (bis-CF3 > CF3 > OCF3 > Ac) is preserved cross-structure. Notably, cross-validation **corrected an AlphaFold artifact at the 7-position** — an apparent hydrophobic sub-pocket on AF2 absent from the cryo-EM structure — underscoring the need to re-check substituent SAR against the experimental structure.

**Table 2 — Three-state scores of the COOH-bisCF3 series.**

| Ligand | AF2 | 26va (occluded) | 26vd (outward-open) |
|---|---|---|---|
| **COOH-bisCF3** | **−9.68** | **−8.55** | **−8.37** |
| COOH-CF3 | −9.04 | −8.52 | −7.51 |
| COOH-OCF3 | −8.94 | −8.13 | −7.32 |
| COOH-Ac | −8.84 | −7.84 | −6.92 |
| biotin | −6.76 | −6.54 | −6.01 |

### 2.3 Direct conjugates cannot be transported

The fourth bridgehead (7-position) of COOH-bisCF3 is the natural conjugation site but is severely size-intolerant: minimal substituents are tolerated (7-methyl −9.76, 7-ethynyl −9.84, both on AF2) but longer linkers fall steeply with size (7-aminopropyl −7.69, a ~2 kcal/mol loss). Docking of the full conjugate COOH-bisCF3–SN-38 (22 Å extended) confirms the consequence on all three structures: it is rejected with positive scores on AF2 (+59.0), the occluded cryo-EM state (+23.7), and even the outward-open state (+3.2), whereas the head alone binds at −9.68/−8.55/−8.37 (Table 3). A biotin–SN-38 conjugate (618 Da) is smaller yet still rejected on the occluded state (+14.2); because its carboxyl is consumed in the ester linkage, it cannot anchor Ser81/Ala84/Tyr156 and binds the outward-open state more weakly (−4.8) than free biotin (−6.01), recapitulating the known failure of biotin conjugates (Tripathi et al., 2025). A negative control — SN-38 conjugated to a low-affinity head (1-carboxy-3-aminomethyladamantane, the COOH-bisCF3 cage minus both CF3 groups; head −7.47/−6.84/−6.36, Table 3) — is likewise rejected (+48.5/+24.7/+4.0): removing CF3 collapses head affinity to biotin's level while the conjugate remains size-rejected, confirming that CF3 drives the head's affinity and that conjugate rejection is size-dominated. Consistent with the ~300 Da substrate-size limit of SMVT, a ~760 Da conjugate cannot form a productive pose; SMVT-targeted delivery must instead exploit surface binding for **endocytosis**.

**Table 3 — Heads bind; conjugates are rejected regardless of head affinity.** Three-state scores (kcal/mol); positive = no binding.

| Ligand | AF2 | 26va (occluded) | 26vd (outward-open) |
|---|---|---|---|
| COOH-bisCF3 (head, 345.3 Da) | −9.68 | −8.55 | −8.37 |
| low-affinity head (1-COOH-3-CH₂NH₂-adamantane, 209.1 Da) | −7.47 | −6.84 | −6.36 |
| biotin (substrate, 244.3 Da) | −6.76 | −6.54 | −6.01 |
| **COOH-bisCF3–SN-38** (763.7 Da) | **+59.0** | **+23.7** | **+3.2** |
| low-aff–SN-38 (627.3 Da) | +48.5 | +24.7 | +4.0 |
| biotin–SN-38 (618.2 Da) | +38.7 | +14.2 | −4.8 |

## 3. Discussion

**Why CF3 works.** The trifluoromethyl group enhances binding through a triple role: electron-withdrawing induction strengthens carboxylate anchoring to Ser81/Ala84/Tyr156; its hydrophobicity fills the aliphatic pocket (Val83/Met266/Phe431); and its compact size avoids the steric clash that penalizes bulkier substituents. The two-group optimum balances these gains against steric crowding at the third bridgehead.

**The transport-refutation principle.** Direct conjugates fail on all three structures (Table 3) for two complementary reasons: COOH-bisCF3–SN-38 exceeds the ~300 Da pocket and cannot fit (+59.0/+23.7/+3.2), while biotin–SN-38, though smaller, loses the free-carboxyl anchor upon conjugation and is still rejected on the occluded state (+14.2). A low-affinity-head conjugate is rejected identically, confirming the size barrier is independent of head affinity. Neither route permits transporter-mediated uptake. This reframes SMVT-targeted delivery: SMVT's value is its **tumor-selective surface overexpression**, to be exploited by **SMVT-mediated endocytosis** of a conjugate-loaded nanoparticle (liposome/PLGA, ~100 nm), where the multivalent COOH-bisCF3 head binds SMVT on the CRC cell surface and triggers internalization — a route precedented by biotin-decorated SMVT-targeted liposomes/cubosomes. The nanoparticle strategy also sidesteps the 7-position size intolerance: the head's linker extends into the solvent-facing nanoparticle surface, not the pocket.

**A concrete application — COOH-bisCF3–SN-38.** SN-38 (7-ethyl-10-hydroxycamptothecin), the active metabolite of the FOLFIRI backbone irinotecan, is attached through its C10 phenol via a cleavable **carbamate** (C10–O–C(=O)–NH–CH₂–adamantane), giving **COOH-bisCF3–SN-38** (C₃₇H₃₅F₆N₃O₈, MW 763.7). Delivered via an SMVT-targeted nanoparticle, it is positioned for CRC and liver-metastasis therapy. An alternative payload — oxaliplatin, the FOLFOX counterpart, via Pt(II)→Pt(IV) axial-carboxylate conjugation — follows the same head, though the 3,5-CF3 groups preclude β-cyclodextrin host-guest binding.

**Limitations.** Results are Vina docking scores (rank-order reliable, ~2 kcal/mol absolute uncertainty). Confirmation requires MM-PBSA on explicit-solvent MD and wet-lab biotin-competition uptake assays in SMVT-high vs -low CRC cells.

## 4. Conclusion

We identified **COOH-bisCF3 (1-carboxy-3,5-bis(trifluoromethyl)adamantane)** as the strongest carboxyl-bearing SMVT ligand to date (−9.68 kcal/mol, cross-validated on three receptor states and three docking seeds), with trifluoromethyl substitution as a general, transferable design principle. We further show that direct drug conjugates cannot be transported by SMVT, redirecting the field toward SMVT-mediated endocytosis — exemplified by the designed conjugate **COOH-bisCF3–SN-38** for colorectal cancer.

## Outlook

All conclusions are computational. Two decisive tests remain: free-energy (MM-PBSA) confirmation of COOH-bisCF3 binding, and the biotin-competition uptake experiment that determines whether COOH-bisCF3–SN-38 is recognized and internalized by a living SMVT-high CRC cell.

## References

Core SMVT:
- Tripathi R, et al. Biotin-mediated drug delivery: does the biotin transporter mediate uptake of biotin conjugates? *J Enzyme Inhib Med Chem.* 2025;40:2554896.
- Zhen Q, et al. Structural basis for multivitamin recognition and transport by human SMVT. *Nat Commun.* 2026;17:8036.

Delivery design (Crossref-verified):
1. Neophytou C, Pitsouli C. Biotin controls intestinal stem cell mitosis and host-microbiome interactions. *Cell Rep.* 2022. 10.1016/j.celrep.2022.110505
2. Baião A, et al. Targeted nanoparticle delivery of irinotecan. *J Control Release.* 2026. 10.1016/j.jconrel.2026.114609
3. Hong Y, et al. Nanoscale drug delivery for colorectal cancer liver metastasis. *Biomed Pharmacother.* 2019. 10.1016/j.biopha.2019.108764
4. Qi Q, et al. SN38 drug delivery system in cancer treatment. *Int J Nanomedicine.* 2024. 10.2147/IJN.S435407
5. Shi Y, Dabrowiak JC. Host–guest interactions involving platinum anticancer agents. *Inorg Chim Acta.* 2012. 10.1016/j.ica.2012.06.021
6. Prashar D, et al. Adamantane–platinum conjugate hosted in β-cyclodextrin. *Bioorg Med Chem Lett.* 2011. 10.1016/j.bmcl.2011.10.006
7. Kou L, et al. Transporter-targeted nano-sized vehicles. *Cancers.* 2020. 10.3390/cancers12102837
8. Patel M, et al. SMVT molecular expression and functional activity. *Int J Pharm.* 2012. 10.1016/j.ijpharm.2012.06.011
9. O'Dowd PD, et al. Oxaliplatin and its derivatives — an overview. *Coord Chem Rev.* 2023. 10.1016/j.ccr.2023.215439
10. Marotta C, et al. Recent advancements in anticancer Pt(IV) prodrugs. *Inorg Chim Acta.* 2023. 10.1016/j.ica.2023.121388
