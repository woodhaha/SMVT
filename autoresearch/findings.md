# SMVT Drug Design — Findings (autoresearch)

## Current Understanding (bootstrap, 2026-11-12)

**The carboxyl-discrimination principle** is the established foundation:
- Free carboxyl = the SMVT recognition "lock" (Ser81/Ala84 backbone NH + Tyr156 in TM1 kink, cryo-EM 26va).
- tripathi2025 experiments: conjugates that consume the carboxyl are NOT transported.
- So any drug-design target must preserve a free carboxyl.

**Scaffold landscape so far** (from the manuscript + prior SAR):
| Scaffold | Best score | Carboxyl? | Transportable? |
|---|---|---|---|
| Glutarimide (3-benzyl-3-chloro-piperidine-2,6-dione) | −9.36 | ❌ | No (scaffold only) |
| Adamantane/diamondoid carboxylic acid (F/CN) | **−8.84** | ✅ | **Yes (best transportable)** |
| Barbiturate | −8.79 | ❌ | No |
| Adamantane-1,3,5-tricarboxylic acid | −8.51 | ✅ | Yes |
| 5-HTP (aromatic amino acid) | −7.66 | ✅ | Yes |
| Biotin (natural substrate) | −6.76/−6.82 | ✅ | Yes |

**Conjugate delivery is refuted**: SMVT substrate recognition requires "carboxyl in TM1 kink + precise pocket geometry" (~300 Da limit). Large-payload conjugates (SN-38/doxorubicin) can't anchor the carboxyl → non-specific. Direction shifted to nanoparticle-mediated endocytosis (adamantane-carboxylate = stronger targeting head than biotin).

## Open question (this project)

Can we find a **carboxyl-bearing** ligand that beats −8.84 (the current transportable best) — i.e., a better SMVT substrate/inhibitor head? The glutarimide −9.36 is carboxyl-free and therefore not transportable; a carboxyl-bearing equivalent would be a materially stronger result.

## Patterns and Insights (prior SAR rules to exploit)

1. **Di-/tri-carboxylate > mono-carboxylate** (adamantane: mean −8.46 di- vs −7.43 mono).
2. **Electron-withdrawing F/CN substitution enhances** (adamantane SAR → −8.84).
3. **Pocket is size-matched to adamantane (C10)**; diamantane (C14) too big.
4. **Hydrophobic cage + polar carboxyl** is the winning combination (matches Val83/Met266/Phe431 hydrophobic pocket).

## Lessons and Constraints

- Docking infra: use `03_Analysis/docking/SMVT_receptor.pdbqt` (meeko AF2, vina-compatible). The `03_Analysis/data/receptor.pdbqt` is BROKEN (old AD4 format, parse error). Vina 1.2.7.
- Ligand prep: RDKit SMILES → ETKDG+MMFF → `mk_prepare_ligand.exe` (meeko CLI) → PDBQT. (meeko Python `PDBQTWriterLegacy` string output also failed to parse — use the CLI.)
- Vina scores have ~2 kcal/mol stochastic variance; rank-order within a series, don't over-read single-mode numbers.
- Biotin baseline reproduces at −6.5 to −6.9 on this receptor.

## Open Questions

- Does H1 (EWG) generalize from adamantane cages to aromatic amino-acid scaffolds?
- Does H3 (ring extension) push 5-HTP (−7.66) below −8.0?
- Is there a carboxyl-bearing scaffold that reaches < −9.0 (the glutarimide bar, but transportable)?

## E1 results (2026-11-12) — aromatic amino-acid chemotype

28 compounds docked. Baseline reproduced: 5-HTP −7.61, biotin −6.86.

Top scorers:
| Compound | Score | Note |
|---|---|---|
| 1-naphthylalanine | **−7.74** | best new; ring extension works |
| 5-F-tryptophan | −7.63 | small F slightly helps |
| 5-HTP (baseline) | −7.61 | reproduced |
| kynurenic acid | −7.56 | carboxyl + quinolone |
| 7-F-tryptophan | −7.50 | |
| indole-3-acetic acid | −7.39 | amino group NOT required |

### Outer-loop synthesis

1. **H3 (ring extension) CONFIRMED**: 1-naphthylalanine −7.74 > 5-HTP −7.61 > tryptophan −7.34.
   Position is decisive: 1-naphthyl fits, 2-naphthyl (−6.62) does NOT.
2. **H1 (EWG) NUANCED**: only *small* F at 5-/7-position helps (5-F −7.63, 7-F −7.50).
   Bulky EWG hurt via steric clash: 5-Cl −6.58, 5-Br −5.95, 5-CF3 −6.36, 5-CN −6.99.
   → "EWG helps" is wrong as a general rule; it's "small substituents at tolerant positions".
3. **New insight — amino group optional**: indole-3-acetic acid (−7.39) ≈ tryptophan (−7.34).
   Carboxyl + aromatic is the recognition core; the amine is dispensable.
4. **H2 (multi-carboxyl) does NOT generalize**: xanthurenic acid (−6.40, 2 COOH + OH) worse
   than single-carboxyl indoles. Multi-carboxyl anchoring is adamantane-cage-specific.

### Direction decision → PIVOT to diamondoid cages (H4)

The aromatic amino-acid chemotype **caps at ~−7.7** (naphthylalanine); it cannot reach the
−8.84 adamantane bar. The hydrophobic cage + carboxyl combination (adamantane) is where
the affinity is. E2 = enumerate diamondoid/bicyclic carboxyl cages + small EWG to beat −8.84.

## Lessons and Constraints (updated)

- Aromatic amino acids: amino group optional, small-F tolerant, bulky substituents steric-fail.
- 1-naphthyl > 2-naphthyl (pocket geometry prefers linear extension at the 1-position).
- **Adamantane cage SMILES gotcha**: the 1-COOH-3-X template is `OC(=O)C12CC3CC(CC(C3)(C1)X)C2`
  (substituent X inside the parens at C3 bridgehead). Appending X after the final `C2`
  puts it at a non-bridgehead CH2 and produces an un-embeddable strained molecule
  (EmbedMolecule rc=-1 → "ligand prep failed"). Also use `useRandomCoords=True` in
  `EmbedMolecule` for rigid cages.

## E2 results (2026-11-12) — cage carboxylates

Calibration reproduced: adamantane-1,3,5-triCOOH −8.35 (manuscript −8.51).

Top scorers:
| Compound | Score | Note |
|---|---|---|
| adamantane-1,3-diCOOH | **−8.48** | di-carboxyl confirms (prior −8.46) |
| adamantane-1,3,5-triCOOH | −8.35 | calibration ✓ |
| 3,5-dimethyl-adam-1-COOH | −8.20 | **NEW: methyl on cage helps** |
| 3-methyl-adam-1-COOH | −7.95 | |
| 1-naphthylalanine | −7.77 | E1 winner reproduced |
| cubane-1,4-diCOOH | −7.57 | |
| bicyclo[2.2.2]octane-1,4-diCOOH | −7.05 | |
| noradamantane / homoadamantane / cubane-mono | −6.4 to −6.9 | all inferior |

### Outer-loop synthesis (E2)

1. **H4 (new cages) REFUTED**: cubane, bicyclo[2.2.2]octane, noradamantane, homoadamantane
   ALL inferior to adamantane. The C10 adamantane cage is uniquely size-matched to the pocket.
2. **Methyl enhancement (new rule)**: 3,5-dimethyl-adam-1-COOH −8.20 > 3-methyl −7.95 >
   1-adamantanecarboxylic −7.42. Small hydrophobic methyls on the cage help (consistent with
   E1 "small substituents help, bulky hurt").
3. **di-COOH ≈ tri-COOH** (−8.48 vs −8.35): two carboxyls suffice; third adds nothing.
4. **Still < −8.84**. The prior focused SAR already found adamantane-COOH-Ac −8.84 (carboxyl +
   acetyl, extra ketone H-bond acceptor). This is the ceiling.

### Direction → E3 (final push)

Combine the two best elements not yet co-tested: carboxyl + acetyl + small group, and
carboxyl + iPr/cPr (not in the prior Me/Et/F/CN/Ac/OH partner set). Target: beat −8.84.
If E3 fails, the answer is "adamantane-COOH-Ac −8.84 is the pocket ceiling".

## E3 + E4 results (2026-11-12) — CF3 BREAKTHROUGH

**E3**: COOH-CF3 (adamantane-1-carboxyl-3-trifluoromethyl) = **−9.02**, first to beat −8.84.
COOH-Ac reproduced −8.82; iPr −8.60, cPr −8.54, CN −8.25.

**E4 (confirmation + expansion)**:
- COOH-CF3 confirmed across 3 seeds: −9.04 / −9.05 / −9.02 → reproducible, not noise.
- **COOH-bisCF3 (two CF3) = −9.69** — best carboxyl-bearing (transportable) ligand yet,
  beating even the carboxyl-free glutarimide (−9.36).
- diCOOH-CF3 −9.16; COOH-CH2CF3 −9.11; COOH-OCF3 −8.94; COOH-CHF2 −8.50.

### CF3 SAR (headline finding)

**Trifluoromethyl on the adamantane cage is a strong affinity enhancer** — the single
biggest discovery of this project. CF3 combines EWG (strengthens carboxyl anchoring) +
hydrophobicity (fills the pocket) + small size (no steric clash). Rank: bis-CF3 (−9.69)
> diCOOH-CF3 (−9.16) > CF3 (−9.04) > OCF3 (−8.94) > CHF2 (−8.50). More CF3 = better.

**This directly beats the manuscript's −8.84 ceiling and the −9.36 glutarimide (carboxyl-free)
bar with a carboxyl-bearing (transportable) ligand.**

## E5 results (2026-11-12) — SAR closed

- **COOH-bisCF3 confirmed across 3 seeds: −9.70 / −9.68 / −9.69 → −9.68 ± 0.01 (robust).**
- triCF3 (3 CF3) −8.03 — steric crowding overcomes EWG benefit.
- diCOOH-bisCF3 −9.16 — adding a 2nd COOH to bisCF3 hurts (geometry disruption).

### FINAL SAR — headline

**1-carboxy-3,5-bis(trifluoromethyl)adamantane = −9.68 kcal/mol** is the strongest
carboxyl-bearing (transportable) SMVT ligand discovered, confirmed across 3 Vina seeds.

CF3 count optimum: 1 (−9.04) < **2 (−9.68)** > 3 (−8.03).

## Conclusion — research question answered

Q: Can we find a carboxyl-bearing SMVT ligand beating −8.84?
A: **Yes.** COOH-bisCF3 (−9.68) beats the −8.84 ceiling by ~0.84 kcal/mol AND the
carboxyl-free glutarimide (−9.36), while retaining the free carboxyl required for
transport. Mechanism: CF3 = electron-withdrawing (strengthens carboxyl anchoring to
Ser81/Ala84/Tyr156) + hydrophobic (fills the Val83/Met266/Phe431 pocket) + small
(no steric clash); two CF3 at the 3,5-bridgeheads is the geometric optimum.

### Full 5-experiment trajectory (~75 compounds)

| Experiment | Chemotype | Best | Verdict |
|---|---|---|---|
| E1 | aromatic amino-acid | 1-naphthylalanine −7.74 | secondary lead; caps ~−7.7 |
| E2 | new cages | adamantane diCOOH −8.48 | alt cages all inferior |
| E3 | carboxyl + small partner | COOH-CF3 −9.02 | first to beat −8.84 |
| E4 | CF3 variants | bisCF3 −9.69 | breakthrough |
| E5 | confirm + push | bisCF3 −9.68 (3 seeds) | confirmed optimal |

### Next steps (validation, not more docking)
1. Re-dock COOH-bisCF3 against cryo-EM 26va (occluded) to confirm cross-structure.
2. MM-PBSA/MD on the SMVT–bisCF3 complex (does it stay bound? carboxyl anchored?).
3. Wet-lab: biotin-competition uptake assay (COOH-bisCF3 is a plausible SMVT substrate).

## RQ2 results (2026-11-12) — cross-structure CONFIRMED

Re-docked CF3 hits against cryo-EM 26va (occluded, biotin-stripped, obabel -xr receptor,
box = biotin centroid 109.09/109.19/107.32). Biotin calibration −6.54 (manuscript −6.49 ✓).

| Compound | AF2 | 26va |
|---|---|---|
| COOH-bisCF3 | −9.68 | **−8.55** |
| COOH-CF3 | −9.04 | −8.52 |
| COOH-CH2CF3 | −9.11 | −8.22 |
| COOH-OCF3 | −8.94 | −8.13 |
| COOH-Ac | −8.84 | −7.84 |

**The CF3 enhancement is NOT AF2 overfitting** — bisCF3 is #1 on both structures, and the
CF3 ranking (bis > CF3 > OCF3 > Ac) is preserved cross-structure. 26va scores are uniformly
less favorable (occluded = tighter pocket), but the relative order holds. COOH-bisCF3 (−8.55
on 26va) also beats the manuscript's 26va naftazone (−7.58) and phenobarbital (−7.49).

## 26va receptor prep (for reuse)

- Source: 05_PDB/26va.pdb → strip BTN (26va_nolig.pdb) → `obabel -xr` → strip non-ATOM/
  HETATM/REMARK lines (obabel keeps HEADER/TITLE/COMPND which vina rejects).
- File: autoresearch/data/26va_receptor_clean.pdbqt. Box = biotin centroid, 25^3.
- Gotchas: meeko mk_prepare_receptor needs `prody` (not installed); obabel leaves HEADER
  lines that must be stripped; box center differs from AF2 (109 vs 5).

## RQ3 results (2026-11-12) — THREE-STATE validation complete

Re-docked CF3 hits against 26vd (outward-open, box = pocket-residue centroid
104.99/103.44/105.33). COOH-bisCF3 is **#1 on all three receptor states**:

| Compound | AF2 | 26va (occluded) | 26vd (outward-open) |
|---|---|---|---|
| COOH-bisCF3 | −9.68 | −8.55 | **−8.37** |
| COOH-CF3 | −9.04 | −8.52 | −7.51 |
| COOH-OCF3 | −8.94 | −8.13 | −7.32 |
| COOH-Ac | −8.84 | −7.84 | −6.92 |
| biotin | — | −6.54 | −6.01 |

**The CF3 enhancement is robust across receptor conformations** — not AF2 overfitting,
not occluded-state-specific. COOH-bisCF3 is the strongest carboxyl-bearing SMVT ligand
across all three states, beating the −8.84 ceiling and the carboxyl-free glutarimide.

## FINAL CONCLUSION

**1-carboxy-3,5-bis(trifluoromethyl)adamantane (COOH-bisCF3) is the strongest
carboxyl-bearing (transportable) SMVT ligand discovered** — #1 on AF2 (−9.68), 26va
occluded (−8.55), and 26vd outward-open (−8.37), with the CF3 ranking preserved across
states. This beats the manuscript's −8.84 adamantane ceiling and the −9.36 carboxyl-free
glutarimide, while retaining the free carboxyl required for SMVT transport.

Remaining validation (for the paper, not more docking):
1. MM-PBSA / membrane MD free energy on the SMVT–bisCF3 complex.
2. Wet-lab biotin-competition uptake assay (bisCF3 is a plausible transportable substrate).

## RQ5 results (2026-11-12) — 7-position linker design (nanoparticle targeting head)

Docked COOH-bisCF3 with a linker at the 4th bridgehead (7-position), AF2 + 26va.

**AF2 only (misleading):** propargyl −9.84, methyl −9.76 both > base −9.68, suggesting a
hydrophobic sub-pocket at position 7.

**Cross-validation (26va) CORRECTS this — the "enhancement" is AF2 overfitting:**

| 7-substituent | AF2 | 26va |
|---|---|---|
| base bisCF3 | −9.68 | **−8.55** (optimal) |
| methyl | −9.76 | −8.04 |
| CH2OH | −9.32 | −6.90 |
| CH2NH2 | −9.12 | −6.59 |
| propargyl | −9.84 | −6.50 (worst!) |

**Corrected conclusion:** on the experimental cryo-EM structure, ANY 7-position substituent
degrades binding (larger = worse). The AF2 pocket has a spurious hydrophobic sub-pocket at
position 7 that the cryo-EM structure lacks. The **robust targeting head is the parent
COOH-bisCF3 (−8.55 on 26va)**; a conjugation handle must be as small as possible (methyl
−8.04 is tolerable). Design implication: (a) minimal covalent linker (accept ~0.5 kcal/mol
loss), or (b) **cyclodextrin–adamantane host-guest non-covalent attachment** (classic
β-CD⊃adamantane interaction) that needs NO 7-position covalent linker. **Lesson: always
cross-validate substituent SAR on 26va — AF2 alone over-predicts pocket tolerance.**
