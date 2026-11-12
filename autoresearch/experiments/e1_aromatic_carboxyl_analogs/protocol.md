# E1 — Aromatic carboxyl-bearing analogs (H3 + partial H1)

## What
Dock ~28 aromatic amino-acid / indole-carboxyl analogs (all free-carboxyl, MW ~160–270)
into the SMVT AF2 pocket.

## Why
5-HTP (−7.66) is a known carboxyl-bearing hit but not optimized. The adamantane cage
(−8.84) is a different chemotype. H3 predicts hydrophobic ring extension + carboxyl
improves affinity; H1 predicts EWG substitution helps. This tests whether the aromatic
amino-acid chemotype can be pushed toward −8.0.

## Prediction
- F/CF3/CN-substituted tryptophans beat 5-HTP (H1).
- Naphthalene/benzofuran-extended analogs beat tryptophan (H3).
- Some approach −8.0; likely none beat −8.84 (adamantane remains best).

## Metric
Vina best-affinity (kcal/mol) on AF2 receptor, exhaustiveness 16, 9 modes.
Baseline: 5-HTP −7.66, biotin −6.76.

## Status
Running (background). Results → e1_results.csv.
