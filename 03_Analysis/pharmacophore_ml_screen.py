#!/usr/bin/env python3
"""
Phase A — Pharmacophore-Guided ML Pre-Screening for SMVT Virtual Screening
=======================================================================
1. Load 84 docked compounds + SMILES + affinity
2. Generate ECFP4 fingerprints (2048-bit)
3. Train Random Forest regressor → predict binding affinity
4. Train RF classifier → predict hit/non-hit (hit = ΔG < −6.5 kcal/mol)
5. Feature importance → map key fingerprint bits back to pharmacophore features
6. Save model for DrugBank screening

Output:
  - models/smvt_rf_regressor.pkl       (affinity predictor)
  - models/smvt_rf_classifier.pkl      (hit/non-hit predictor)
  - models/pharmacophore_bits.csv      (top predictive fingerprint bits)
  - pharmacophore_ml_screen.log        (diagnostics)
"""

import os, sys, pickle, logging, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw, rdFingerprintGenerator
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Avalon import pyAvalonTools

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import cross_val_score, LeaveOneOut, StratifiedKFold
from sklearn.metrics import r2_score, mean_absolute_error, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

os.chdir("D:/Researching/SMVT")
os.makedirs("03_Analysis/models", exist_ok=True)

# ═══ Logging ═══
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("03_Analysis/models/pharmacophore_ml.log", mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════
# 1. LOAD DATA — SMILES for all 84 docked compounds
# ═══════════════════════════════════════════════════════════════

# Full compound library with SMILES (from docking_expanded.py + docking_round3.py)
COMPOUND_SMILES = {
    # ── Natural substrates ──
    "Biotin":                  "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2",
    "Lipoic_Acid":             "C1CSSC1CCCCC(=O)O",
    "Pantothenic_Acid":        "CC(C)(CO)C(O)C(=O)NCCC(=O)O",
    "Gabapentin_enacarbil":    "CC(C)(C(=O)O)OC(=O)NCC1(CCCCC1)CC(=O)O",
    "Desthiobiotin":           "CC1C(NC(=O)N1)CCCCCC(=O)O",
    "Biotin_Methyl_Ester":     "C1C2C(C(S1)CCCCC(=O)OC)NC(=O)N2",
    "Biotin_Sulfoxide":        "C1C2C(C(S1(=O))CCCCC(=O)O)NC(=O)N2",
    "Biotin_Sulfone":          "C1C2C(C(S1(=O)=O)CCCCC(=O)O)NC(=O)N2",
    "Norbiotin":               "C1C2C(C(S1)CCCC(=O)O)NC(=O)N2",
    "Homobiotin":              "C1C2C(C(S1)CCCCCC(=O)O)NC(=O)N2",

    # ── NSAIDs ──
    "Indomethacin":   "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=C(C=C2)OC)CC(=O)O",
    "Ibuprofen":      "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Diclofenac":     "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl",
    "Ketoprofen":     "CC(C1=CC=CC(=C1)C(=O)C2=CC=CC=C2)C(=O)O",
    "Flurbiprofen":   "CC1=CC(=CC=C1)C2=CC=C(C=C2)C(C)C(=O)O",
    "Phenylbutazone": "CCCCC1C(=O)N(N(C1=O)C2=CC=CC=C2)C3=CC=CC=C3",
    "Naproxen":       "CC(C1=CC2=C(C=C1)C=CC(=C2)OC)C(=O)O",
    "Celecoxib":      "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",
    "Piroxicam":      "CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3",
    "Meloxicam":      "CC1=CC=CS(=O)(=O)N1C(=O)NC2=C(C3=NC=C(S3)C)O2",
    "Aspirin":        "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Sulindac":       "CC1=C(C2=C(C=CC(=C2)F)C(=C1)CC(=O)O)C3=CC=C(C=C3)S(=O)C",
    "Mefenamic_Acid": "CC1=CC=CC(=C1NC2=CC=CC=C2C(=O)O)C",
    "Tolfenamic_Acid": "CC1=CC=CC(=C1NC2=CC=CC=C2C(=O)O)Cl",
    "Niflumic_Acid":   "FC(F)(F)C1=CC=NC(=C1)NC2=CC=CC=C2C(=O)O",
    "Etofenamate":     "CCOCCOC(=O)C1=CC=CC=C1NC2=CC=CC(=C2)C(F)(F)F",
    "Meclofenamic_Acid": "CC1=CC=C(C(=C1)NC2=CC=CC=C2C(=O)O)Cl",
    "Flufenamic_Acid":   "FC(F)(F)C1=CC=CC(=C1)NC2=CC=CC=C2C(=O)O",
    "Nabumetone":      "COC1=CC2=C(C=C1)C=C(C=C2)CCC(=O)C",
    "Oxaprozin":       "C1=CC=C(C=C1)C2=CC(=NO2)CCC(=O)O",
    "Etodolac":        "CCC1=C2C(=CC=C1)OC3=C(C2=CC(=C3)CC(=O)O)CC",
    "Fenbufen":        "C1=CC=C(C=C1)C(=O)CCC2=CC=C(C=C2)CC(=O)O",
    "Tiaprofenic_Acid": "CC(C1=CC=CC=C1C(=O)C2=CC=CS2)C(=O)O",

    # ── Statins ──
    "Simvastatin":     "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(O)(CC(=O)O3)C)C",
    "Atorvastatin":    "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
    "Rosuvastatin":    "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
    "Pravastatin":     "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(O)(CC(=O)O3)C)C",
    "Fluvastatin":     "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",

    # ── ACE inhibitors ──
    "Enalaprilat":     "CC(C(=O)O)NC(C)C(=O)N1CC2=CC=CC=C2CC1C(=O)O",
    "Lisinopril":      "CC(C(=O)O)NC(CCC1=CC=CC=C1)C(=O)N2CCCC2C(=O)O",
    "Captopril":       "CC(CS)C(=O)N1CCCC1C(=O)O",

    # ── Other FDA drugs ──
    "Valproic_Acid":   "CCCC(CCC)C(=O)O",
    "Methotrexate":    "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O",
    "Levodopa":        "C1=CC(=C(C=C1CC(C(=O)O)N)O)O",
    "Tranexamic_Acid": "C1CC(CCC1CN)C(=O)O",
    "Metformin":       "CN(C)C(=N)NC(=N)N",
    "Omeprazole":      "CC1=CN=C(C(=C1OC)C)CSC2=NC=C(C(=C2C)OC)C",
    "Furosemide":      "C1=CC(=C(C=C1C(=O)O)S(=O)(=O)N)Cl",
    "Sertraline":      "CNC(C1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl)C",
    "Hydrochlorothiazide": "C1NC2=CC(=C(C=C2S(=O)(=O)N1)S(=O)(=O)N)Cl",
    "Metronidazole":   "CC1=NC=C(N1CCO)[N+](=O)[O-]",
    "Fluoxetine":      "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
    "Probenecid":      "CCCN(CCC)S(=O)(=O)C1=CC=C(C=C1)C(=O)O",
    "Acetaminophen":   "CC(=O)NC1=CC=C(C=C1)O",
    "Allopurinol":     "C1=NC2=C(N1)C(=O)NC(=O)N2",
    "Isoniazid":       "C1=CC=NC=C1C(=O)NN",
    "Caffeine":        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Acetazolamide":   "CC(=O)NC1=NN=C(S1)S(=O)(=O)N",
    "Warfarin":        "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",

    # ── Small acids ──
    "Benzoic_Acid":       "C1=CC=CC=C1C(=O)O",
    "Salicylic_Acid":     "C1=CC=CC(=C1C(=O)O)O",
    "Phenylacetic_Acid":  "C1=CC=CC=C1CC(=O)O",
    "Hippuric_Acid":      "C1=CC=CC=C1C(=O)NCC(=O)O",
    "Succinic_Acid":      "C(CC(=O)O)C(=O)O",
    "Glutaric_Acid":      "C(CC(=O)O)CC(=O)O",
    "Adipic_Acid":        "C(CCC(=O)O)CC(=O)O",
    "Pimelic_Acid":       "C(CCCC(=O)O)CC(=O)O",
    "Suberic_Acid":       "C(CCCCC(=O)O)CC(=O)O",
    "Azelaic_Acid":       "C(CCCCCCC(=O)O)CC(=O)O",
    "Sebacic_Acid":       "C(CCCCCCCC(=O)O)CC(=O)O",

    # ── Aromatic amino acids ──
    "L-Phenylalanine":     "C1=CC=C(C=C1)CC(C(=O)O)N",
    "L-Tyrosine":          "C1=CC(=CC=C1CC(C(=O)O)N)O",
    "L-Tryptophan":        "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N",
    "5-Hydroxytryptophan": "C1=CC2=C(C=C1O)C(=CN2)CC(C(=O)O)N",

    # ── Bile acids ──
    "Cholic_Acid":          "CC(CCC(=O)O)C1CCC2C1(C(CC3C2C(CC4C3(CCC(C4)O)C)O)O)C",
    "Deoxycholic_Acid":     "CC(CCC(=O)O)C1CCC2C1(C(CC3C2CC(C4C3(CCC(C4)O)C)O)O)C",
    "Ursodeoxycholic_Acid": "CC(CCC(=O)O)C1CCC2C1(CCC3C2C(CC4C3(CCC(C4)O)C)O)C",

    # ── Fibrates ──
    "Gemfibrozil":   "CC1=CC(=C(C=C1)C)OCCCC(C)(C)C(=O)O",
    "Bezafibrate":   "CC(C)(C(=O)O)OC1=CC=C(C=C1)CCNC(=O)C2=CC=C(C=C2)Cl",
    "Ciprofibrate":  "CC(C)(C(=O)O)OC1=CC=C(C=C1)C2CC2",

    # ── Leukotriene antagonists ──
    "Montelukast":  "CC(C)(C1=CC=CC=C1CC[C@H](C2=CC=CC(=C2)CC=C(C)C)SCC3(CC3)CC(=O)O)O",
    "Zafirlukast":  "CC1=CC(=C(C=C1)C2=CC=C(C=C2)C(=O)NS(=O)(=O)C3=CC=CC=C3C(=O)O)N",

    # ── Vitamins ──
    "Ascorbic_Acid":  "C(C(C1C(=C(C(=O)O1)O)O)O)O",
    "Nicotinic_Acid": "C1=CC=CN=C1C(=O)O",
    "Pyridoxine":     "CC1=NC=C(C(=C1O)CO)CO",
    "Thiamine":       "CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO",
    "Riboflavin":     "CC1=C(C=C(C(=C1)C)N=C2C(=O)NC(=O)N=C2N3CC(C(C3O)O)O)O",
    "Folic_Acid":     "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)N=C(N3)N",
    "Retinoic_Acid":  "CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC(=O)O)C)C",
}

# Load docking results
df_results = pd.read_csv("03_Analysis/outputs/docking_expanded_results.csv")
log.info(f"Loaded {len(df_results)} docking results")

# Merge SMILES
df_results["smiles"] = df_results["name"].map(COMPOUND_SMILES)
missing_smiles = df_results[df_results["smiles"].isna()]["name"].tolist()
if missing_smiles:
    log.warning(f"Missing SMILES for {len(missing_smiles)} compounds: {missing_smiles}")
    df_results = df_results[df_results["smiles"].notna()].copy()

log.info(f"Compounds with SMILES: {len(df_results)}")

# ═══════════════════════════════════════════════════════════════
# 2. GENERATE FINGERPRINTS
# ═══════════════════════════════════════════════════════════════

# Morgan/ECFP4 fingerprints (radius=2, 2048 bits) — gold standard for drug-like molecules
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fp(smi):
    """Convert SMILES to ECFP4 fingerprint, return None on failure."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return mfpgen.GetFingerprint(mol)

fps = []
valid_idx = []
for i, row in df_results.iterrows():
    fp = smiles_to_fp(row["smiles"])
    if fp is not None:
        fps.append(fp)
        valid_idx.append(i)

df_valid = df_results.loc[valid_idx].copy()
fp_array = np.array([list(fp) for fp in fps])  # 2048-bit binary

log.info(f"Valid fingerprints: {len(fp_array)} ({len(df_results) - len(fp_array)} failed)")

# ═══════════════════════════════════════════════════════════════
# 3. FEATURE ENGINEERING — add simple molecular descriptors
# ═══════════════════════════════════════════════════════════════

mols = [Chem.MolFromSmiles(s) for s in df_valid["smiles"]]

descriptor_features = pd.DataFrame({
    "MW":           [Descriptors.MolWt(m) for m in mols],
    "LogP":         [Descriptors.MolLogP(m) for m in mols],
    "HBA":          [Descriptors.NumHAcceptors(m) for m in mols],
    "HBD":          [Descriptors.NumHDonors(m) for m in mols],
    "RotBonds":     [Descriptors.NumRotatableBonds(m) for m in mols],
    "TPSA":         [Descriptors.TPSA(m) for m in mols],
    "RingCount":    [Descriptors.RingCount(m) for m in mols],
    "AromaticRing": [Descriptors.NumAromaticRings(m) for m in mols],
    "Carboxyl":     [s.count("C(=O)O") for s in df_valid["smiles"]],
    "FractionCsp3": [Descriptors.FractionCSP3(m) for m in mols],
    "HeavyAtom":    [Descriptors.HeavyAtomCount(m) for m in mols],
})

# Normalize descriptors
scaler = StandardScaler()
desc_scaled = scaler.fit_transform(descriptor_features)

# Combine fingerprints + descriptors
X = np.hstack([fp_array, desc_scaled])
y = df_valid["best_affinity"].values

log.info(f"Feature matrix: {X.shape} (ECFP4=2048 + descriptors={desc_scaled.shape[1]})")
log.info(f"Affinity range: {y.min():.2f} to {y.max():.2f} kcal/mol")

# ═══════════════════════════════════════════════════════════════
# 4. TRAIN RANDOM FOREST REGRESSOR
# ═══════════════════════════════════════════════════════════════

rf_reg = RandomForestRegressor(
    n_estimators=500,
    max_depth=8,
    min_samples_leaf=3,
    max_features="sqrt",
    random_state=42,
    n_jobs=-1
)

# Leave-One-Out cross-validation (appropriate for N=84)
loo = LeaveOneOut()
cv_scores = cross_val_score(rf_reg, X, y, cv=loo, scoring="neg_mean_absolute_error")
mae_cv = -cv_scores.mean()
r2_cv = cross_val_score(rf_reg, X, y, cv=loo, scoring="r2").mean()

log.info(f"RF Regressor LOO-CV: MAE = {mae_cv:.3f} kcal/mol, R² = {r2_cv:.3f}")

# Train on full dataset
rf_reg.fit(X, y)
y_pred = rf_reg.predict(X)
r2_full = r2_score(y, y_pred)
mae_full = mean_absolute_error(y, y_pred)
log.info(f"RF Regressor Full: R² = {r2_full:.3f}, MAE = {mae_full:.3f} kcal/mol")

# ═══════════════════════════════════════════════════════════════
# 5. TRAIN RANDOM FOREST CLASSIFIER (hit vs non-hit)
# ═══════════════════════════════════════════════════════════════

# Hit definition: ΔG < −6.5 kcal/mol (~top 25% of current library, better than avg substrate)
HIT_THRESHOLD = -6.5
y_class = (y < HIT_THRESHOLD).astype(int)
n_hits = y_class.sum()
log.info(f"Hit threshold: ΔG < {HIT_THRESHOLD} kcal/mol → {n_hits}/{len(y)} compounds ({100*n_hits/len(y):.1f}%)")

rf_clf = RandomForestClassifier(
    n_estimators=500,
    max_depth=8,
    min_samples_leaf=3,
    max_features="sqrt",
    class_weight="balanced",
    random_state=42,
    n_jobs=-1
)

# Stratified CV
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_auc = cross_val_score(rf_clf, X, y_class, cv=skf, scoring="roc_auc").mean()
cv_mcc = cross_val_score(rf_clf, X, y_class, cv=skf, scoring="matthews_corrcoef").mean()

log.info(f"RF Classifier 5-fold CV: AUC = {cv_auc:.3f}, MCC = {cv_mcc:.3f}")

rf_clf.fit(X, y_class)
y_prob = rf_clf.predict_proba(X)[:, 1]
auc_full = roc_auc_score(y_class, y_prob)
log.info(f"RF Classifier Full: AUC = {auc_full:.3f}")

# ═══════════════════════════════════════════════════════════════
# 6. FEATURE IMPORTANCE — Pharmacophore Mapping
# ═══════════════════════════════════════════════════════════════

# Combined feature importance
fp_importance = rf_reg.feature_importances_[:2048]
desc_importance = rf_reg.feature_importances_[2048:]

# Top fingerprint bits
top_n = 50
top_bits = np.argsort(fp_importance)[::-1][:top_n]

# Map bits back to chemical meaning using representative compounds
bit_info = []
for bit_idx in top_bits:
    imp = fp_importance[bit_idx]
    # Find compounds where this bit is set
    active_compounds = []
    for i, fp in enumerate(fps):
        if fp[int(bit_idx)]:
            active_compounds.append((df_valid.iloc[i]["name"], y[i]))

    # Sort by affinity (best binders first)
    active_compounds.sort(key=lambda x: x[1])

    bit_info.append({
        "bit": bit_idx,
        "importance": imp,
        "n_active": len(active_compounds),
        "top_compounds": ", ".join([c for c, a in active_compounds[:8]]),
        "mean_affinity": np.mean([a for _, a in active_compounds]) if active_compounds else 0,
        "hit_enrichment": sum(1 for _, a in active_compounds if a < HIT_THRESHOLD) / max(1, len(active_compounds))
    })

df_bits = pd.DataFrame(bit_info)
df_bits.to_csv("03_Analysis/models/pharmacophore_bits.csv", index=False)
log.info(f"Top {top_n} pharmacophore bits saved")

# Key pharmacophore features summary
hit_enriched_bits = df_bits[df_bits["hit_enrichment"] > 0.4].sort_values("importance", ascending=False)
log.info(f"\n{'='*60}")
log.info(f"PHARMACOPHORE FEATURES (bits enriched in hits > 40%)")
log.info(f"{'='*60}")
for _, row in hit_enriched_bits.head(15).iterrows():
    log.info(f"  Bit {row['bit']:4d} | imp={row['importance']:.4f} | "
             f"hit%={row['hit_enrichment']:.1%} | n={row['n_active']:2d} | "
             f"top: {row['top_compounds'][:100]}")

# Descriptor importance
desc_names = list(descriptor_features.columns)
desc_ranked = sorted(zip(desc_names, desc_importance), key=lambda x: -x[1])
log.info(f"\n{'='*60}")
log.info(f"DESCRIPTOR IMPORTANCE")
log.info(f"{'='*60}")
for name, imp in desc_ranked:
    log.info(f"  {name:20s}: {imp:.4f}")

# ═══════════════════════════════════════════════════════════════
# 7. SAVE MODELS + METADATA
# ═══════════════════════════════════════════════════════════════

model_artifacts = {
    "rf_regressor": rf_reg,
    "rf_classifier": rf_clf,
    "scaler": scaler,
    "descriptor_names": list(descriptor_features.columns),
    "fp_generator": "Morgan_ECFP4_r2_n2048",
    "hit_threshold": HIT_THRESHOLD,
    "feature_names": [f"ECFP4_{i}" for i in range(2048)] + [f"DESC_{n}" for n in descriptor_features.columns],
    "training_size": len(X),
    "cv_mae": mae_cv,
    "cv_r2": r2_cv,
    "cv_auc": cv_auc,
    "cv_mcc": cv_mcc,
    "top_pharmacophore_bits": df_bits.to_dict("records"),
}

with open("03_Analysis/models/smvt_ml_screen.pkl", "wb") as f:
    pickle.dump(model_artifacts, f)

log.info(f"\n{'='*60}")
log.info(f"MODEL SAVED → 03_Analysis/models/smvt_ml_screen.pkl")
log.info(f"{'='*60}")
log.info(f"  CV MAE:          {mae_cv:.3f} kcal/mol")
log.info(f"  CV R²:           {r2_cv:.3f}")
log.info(f"  CV AUC:          {cv_auc:.3f}")
log.info(f"  CV MCC:          {cv_mcc:.3f}")
log.info(f"  Hit threshold:   < {HIT_THRESHOLD} kcal/mol")
log.info(f"  Training set:    {len(X)} compounds")

# ═══════════════════════════════════════════════════════════════
# 8. PREDICT ON TRAINING SET — diagnostic table
# ═══════════════════════════════════════════════════════════════

df_valid_copy = df_valid.copy()
df_valid_copy["predicted_affinity"] = y_pred
df_valid_copy["predicted_hit_prob"] = y_prob
df_valid_copy["residual"] = y - y_pred
df_valid_copy["abs_error"] = np.abs(y - y_pred)

# Top predicted hits
log.info(f"\n{'='*60}")
log.info(f"TOP 20 PREDICTED HITS (from training set, for validation)")
log.info(f"{'='*60}")
df_top = df_valid_copy.nsmallest(20, "predicted_affinity")
for i, (_, row) in enumerate(df_top.iterrows()):
    flag = "← MISMATCH" if (row["best_affinity"] < HIT_THRESHOLD) != (row["predicted_affinity"] < HIT_THRESHOLD) else ""
    log.info(f"  {i+1:2d}. {row['name']:25s} | true={row['best_affinity']:.2f} | pred={row['predicted_affinity']:.2f} | "
             f"p(hit)={row['predicted_hit_prob']:.2f} | err={row['abs_error']:.2f} {flag}")

# Worst predictions
log.info(f"\nLARGEST PREDICTION ERRORS:")
df_worst = df_valid_copy.nlargest(10, "abs_error")
for _, row in df_worst.iterrows():
    log.info(f"  {row['name']:25s} | true={row['best_affinity']:.2f} | pred={row['predicted_affinity']:.2f} | "
             f"err={row['abs_error']:.2f} kcal/mol")

log.info("\nDone. Model ready for DrugBank screening.")
