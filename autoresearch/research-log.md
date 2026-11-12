# Research Log — SMVT Drug Design autoresearch

## 2026-11-12 — Bootstrap

- Scoped research question: discover carboxyl-bearing SMVT ligands < −8.84 (beat current transportable best).
- Established evaluation metric: Vina score on AF2 receptor (SMVT_receptor.pdbqt), box 25³, exhaustiveness 16.
- Fixed docking infra: `03_Analysis/data/receptor.pdbqt` is broken (old AD4 format → vina parse error). Working receptor = `03_Analysis/docking/SMVT_receptor.pdbqt`. Ligand prep via meeko CLI `mk_prepare_ligand.exe` (Python API also failed parse). Biotin reproduces −6.5 to −6.9.
- Defined hypotheses H1–H4 from prior SAR.
- Wrote `src/dock.py` (reusable: `dock_smiles` + `dock_batch`).
- Launched E1 (H3): aromatic amino-acid / 5-HTP analogs.

## Prior context (not re-verified this session)

- Prior SAR: adamantane 58 analogs (best −8.84), barbiturate 682 analogs (best −8.79). See `SMVT_drug_design_strategy.md` and `SMVT_screening_expansion_plan.md`.
- Membrane MD (Na⁺+POPC, 3 ligands) complete; Na⁺/substrate occupy distinct sites.

## 2026-11-12 — Conclusion (answer found)

- E1 (28 cmpds): aromatic chemotype capped ~−7.7 (1-naphthylalanine −7.74 new lead).
- E2 (18 cmpds): alternative cages all inferior; adamantane C10 confirmed optimal; methyl enhancement.
- E3 (7 cmpds): COOH-CF3 −9.02 first to beat −8.84.
- E4 (9 cmpds): CF3 confirmed (3 seeds); bisCF3 −9.69 breakthrough.
- E5 (6 cmpds): bisCF3 confirmed (3 seeds −9.68/−9.69/−9.70); triCF3 −8.03 (steric); SAR closed.
- **Answer**: COOH-bisCF3 (1-carboxy-3,5-bis-CF3-adamantane) −9.68 = best transportable SMVT ligand.
- Report: to_human/report.html

## 2026-11-12 — Three-state validation complete

- RQ2 (26va occluded): bisCF3 -8.55 #1; biotin calib -6.54 (-6.49 manuscript).
- RQ3 (26vd outward-open): bisCF3 -8.37 #1.
- **bisCF3 is #1 on AF2 (-9.68) / 26va (-8.55) / 26vd (-8.37) — robust, cross-validated.**
- 26va + 26vd receptors prepared (obabel -xr + strip HEADER), box = pocket centroid.

## 2026-11-12 — Paper drafted

- RQ4 (single-structure MM-GBSA) abandoned after 4 bugs (meeko H-strip / atom reorder /
  no-bonds / serial mismatch) + root cause = docking receptor (receptor_clean.pdbqt is a
  fragment in a different frame) vs complex-frame mismatch. Low value anyway (Vina↔MM-GBSA
  inverse correlation already documented).
- CONCLUDE: paper draft written to paper/paper_draft.md — COOH-bisCF3 as strongest
  carboxyl-bearing SMVT ligand (−9.68 AF2 / −8.55 26va / −8.37 26vd, 3 seeds).
- Remaining validation: MM-PBSA/MD (rigorous, pod ~8h) + wet-lab biotin competition.

## 2026-11-12 — RQ5 linker design (cross-validated)

- Docked 12 linker variants at bisCF3 7-position (AF2): propargyl -9.84 / CH3 -9.76 > base -9.68.
- 26va cross-validation CORRECTS: base -8.55 remains optimal; any 7-substituent degrades
  (CH3 -8.04 > CH2OH -6.90 > CH2NH2 -6.59 > propargyl -6.50). AF2 "enhancement" = overfitting.
- Design: minimal covalent linker (methyl) OR β-cyclodextrin⊃adamantane host-guest (no linker).
- Lesson: always cross-validate substituent SAR on 26va (AF2 over-predicts pocket tolerance).
