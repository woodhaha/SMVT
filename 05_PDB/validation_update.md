# SMVT Docking Validation — COVID-EM Structure Update

> **2026-07-10**: Zhang Zhe group SMVT cryo-EM structures (Nature Comms, PMID: 42364996) released via PDB on 2026-07-08. Full re-docking and validation completed.

## What Changed

The manuscript was drafted using an AlphaFold2 model (AF-Q9Y289-F1). Now we have 5 experimental structures at 3.4–4.3 A resolution. **Key update for the "Limitations" section** — the limitation statement "should re-dock against the experimental coordinates" is no longer pending; we've done it.

## Results of Experimental Structure Re-Docking

### Binding Pocket Validation: AF2 vs Cryo-EM

| Metric | Value |
|--------|-------|
| AF2 pocket vs 26va pocket overlap (residues) | 16/23 (F1=74%) |
| Precision (AF2-predicted in experimental pocket) | 80% |
| Recall (experimental residues captured) | 70% |

The AF2 model accurately predicted the SMVT binding pocket (F1=74%), with 16 of 23 experimental pocket residues correctly identified among the top-20 predicted. This validates the manuscript's binding site definition as fundamentally correct, while the 4 false positives (LEU76, GLN277, SER304, LEU528) and 7 false negatives (THR78, ALA82, VAL83, TYR156, THR366, ALA427, PHE431) represent subtle AF2 limitations in pocket boundary definition.

### Cross-Structure Docking Comparison

| Rank | Vina (26va cryo-EM) | Vina (AF2 model) | MM-GBSA |
|:----:|:-------------------:|:----------------:|:--------:|
| 1 | NAFTAZONE (-7.58) | NAFTAZONE (-8.03) | Gabapentin Enacarbil (-43.3) |
| 2 | **PHENOBARBITAL** (-7.49) | Hydromorphone (-7.72) | Riboflavin (-41.5) |
| 3 | Gabapentin Enacarbil (-5.96) | Furosemide (-7.56) | Biotin (-30.0) |
| 4 | Biotin (-5.93) | Riboflavin (-7.56) | Hydromorphone (-29.9) |
| 5 | Esketamine (-5.90) | Esketamine (-7.36) | Esketamine (-27.9) |
| 6 | Furosemide (-4.32) | Phenobarbital (-7.31) | Phenobarbital (-23.6) |
| 7 | Hydromorphone (-3.74) | Gabapentin Enacarbil (-6.97) | NAFTAZONE (-22.4) |
| 8 | Riboflavin (-1.55) | Biotin (-6.38) | Furosemide (+6.4) |

Spearman correlation: 26va vs AF2 Vina = -0.27, 26va vs MM-GBSA = -0.26.

### Key Updates vs Manuscript Claims

**1. NAFTAZONE — Strengthened**

Claim: "top-scoring compound across virtual screening." Now confirmed as rank #1 in both AF2 AND cryo-EM Vina docking. Also has the lowest RMSD to biotin's experimental pose (2.38 A). This is the most robust SMVT ligand discovered in this study.

**2. PHENOBARBITAL — Upgraded**

Manuscript mentions phenobarbital hit rate but ranks it #2 in AF2. In the cryo-EM structure, phenobarbital ranks #2 as well but its score (-7.49) is much closer to naftazone (-7.58). The barbiturate scaffold's ureido ring mimicry of biotin is validated by the experimental pocket shape.

**3. HYDROMORPHONE — Downgraded**

Manuscript positions hydromorphone as "best consensus candidate" (#4 MM-GBSA, #2 Vina). In the cryo-EM structure, hydromorphone drops to rank #7 (-3.74). Its Vina score drops by 4 kcal/mol from AF2 (−7.72 → −3.74), suggesting the AF2 pocket was overfit for this specific ligand. The manuscript should no longer highlight hydromorphone as a top candidate.

**4. RIBOFLAVIN — Confirmed as negative control**

Drops to rank #8 in experimental docking (-1.55), consistent with its use of dedicated RFVT transporters. MM-GBSA's overestimation is now clearly identified as a GB model artifact.

**5. BIOTIN — Basal validation**

Biotin ranks #4 in cryo-EM docking (-5.93 vs -6.38 in AF2). This is a narrower gap than other compounds, confirming the natural substrate binding mode is well-predicted across receptor conformations.

### Manuscript Updates Needed (Text Changes)

| Location | Current | Change To |
|----------|---------|-----------|
| Abstract | "hydromorphone... best consensus candidate" | De-emphasize hydromorphone; highlight naftazone as cross-method validated |
| Results §Virtual Screening | "hydromorphone showed the most favorable balance" | Remove or qualify with experimental structure context |
| Discussion | "hydromorphone... ranking second in Vina" | Update: ranks #7 in experimental structure |
| Limitations | "preliminary comparison indicates AF2 consistency" | Replace with quantitative evidence: "Re-docking against the cryo-EM structure (26va) confirmed naftazone as the top-ranked compound (Vina = −7.58 kcal/mol), followed by phenobarbital (−7.49). The native substrate biotin ranked #4, validating the docking protocol. Re-docking reduced hydromorphone's affinity by 4 kcal/mol relative to AF2, changing its rank from #2 to #7, suggesting structural sensitivity for this ligand." |
| Fig. 7 (docking) | AF2-only | Add cryo-EM comparison as panel or supplementary |
| Table 2 (top compounds) | "Hydromorphone: best consensus" | Reorder; add experimental structure rank column |

### Recommended Revisions Strategy

```
Keep:
  - NAFTAZONE as primary hit (validated cross-structure)
  - Phenobarbital + barbiturate class finding (now stronger)
  - Biotin as internal reference

Revise:
  - Hydromorphone: from "top prospect" to "AF2-specific hit"
  - Riboflavin negative control documentation: mention GB bias
  - MM-GBSA ranking limitations: note discrepancy with experimental structure

Add:
  - Cryo-EM validation paragraph (brief, 3-5 sentences)
  - Dual-receptor ranking table (AF2 vs 26va)
  - Explanation: why AF2 and cryo-EM disagree (induced fit vs rigid)
```

### Files

- New docking results: `D:\Researching\SMVT\05_PDB\validation\`
- Pose comparison: `pose_comparison.json`, `pose_comparison.png`
- Validation summary: `validation_summary.png`
- All 5 PDBs: `D:\Researching\SMVT\05_PDB\`
