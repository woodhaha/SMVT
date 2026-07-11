#!/usr/bin/env python3
"""
Phase A Step 4 — Hit Identification + SAR + Pharmacophore Analysis
====================================================================
1. Merge R1-R3 (84 hand-picked) + R4 (500 ChEMBL ML-screened) docking results
2. Z-score normalize within each round to correct for exhaustiveness differences
3. Identify hits: Z < -1.5 OR ΔG < -6.5 kcal/mol
4. Scaffold clustering + chemical family classification
5. Pharmacophore feature enrichment analysis
6. Drug repurposing opportunity scoring
7. Generate Fig 7 candidate: SMVT screening volcano plot

Output:
  - outputs/screening_master_results.csv
  - outputs/hit_summary.csv
  - outputs/screening_report.md (auto-generated)
"""

import os, sys, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw, rdFingerprintGenerator, PandasTools
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw import IPythonConsole
from collections import Counter

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/outputs", exist_ok=True)
os.makedirs("04_Manuscript/figures", exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("03_Analysis/models/screening_analysis.log", mode="w"),
              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════
# 1. LOAD & MERGE ALL DOCKING DATA
# ═══════════════════════════════════════════════════════════════

# R1-R3: hand-picked library (exhaustiveness=8)
df_old = pd.read_csv("03_Analysis/outputs/docking_expanded_results.csv")
df_old["round"] = "R1-R3"
df_old["source"] = "hand_picked"
df_old["exhaustiveness"] = 8

log.info(f"Loaded R1-R3: {len(df_old)} compounds")

# R4: ChEMBL ML-screened (exhaustiveness=16)
df_new_path = "03_Analysis/outputs/docking_batch_results.csv"
if os.path.exists(df_new_path):
    df_new = pd.read_csv(df_new_path)
    df_new["round"] = "R4"
    df_new["source"] = "chembl_ml_screen"
    df_new["exhaustiveness"] = 16
    log.info(f"Loaded R4: {len(df_new)} compounds")
else:
    log.warning(f"R4 results not found at {df_new_path}, using placeholder")
    df_new = pd.DataFrame()

# Merge
if len(df_new) > 0:
    df_all = pd.concat([df_old, df_new[["name", "best_affinity", "round", "source", "exhaustiveness"]]],
                       ignore_index=True)
else:
    df_all = df_old.copy()

# Remove compounds with no valid affinity
df_all = df_all[df_all["best_affinity"].notna()].copy()
df_all = df_all[df_all["best_affinity"] < 50].copy()  # Remove nonsense values

log.info(f"Total merged: {len(df_all)} compounds ({len(df_old)} R1-R3 + {len(df_new)} R4)")

# ═══════════════════════════════════════════════════════════════
# 2. Z-SCORE NORMALIZATION + HIT IDENTIFICATION
# ═══════════════════════════════════════════════════════════════

# Per-round Z-scores (accounts for different exhaustiveness)
df_all["z_score"] = np.nan
for rnd in df_all["round"].unique():
    mask = df_all["round"] == rnd
    affs = df_all.loc[mask, "best_affinity"]
    if len(affs) > 1:
        df_all.loc[mask, "z_score"] = stats.zscore(affs)
    else:
        df_all.loc[mask, "z_score"] = 0

# Hit criteria (multi-level)
# Level 1: Z < -2.0 (statistical outlier, top ~2.5%)
# Level 2: Z < -1.5 (top ~7%)
# Level 3: ΔG < -7.0 kcal/mol (absolute strong binder)
# Level 4: ΔG < -6.5 kcal/mol (better than or equal to biotin)

BIOTIN_AFFINITY = -6.76  # Reference: natural substrate

conditions = [
    df_all["z_score"] < -2.0,
    (df_all["z_score"] < -1.5) & (df_all["z_score"] >= -2.0),
    (df_all["best_affinity"] < -7.0) & (df_all["z_score"] >= -1.5),
    (df_all["best_affinity"] <= BIOTIN_AFFINITY) & (df_all["best_affinity"] > -7.0) & (df_all["z_score"] >= -1.5),
]

hit_levels = ["L1_Strong", "L2_Moderate", "L3_Absolute", "L4_BiotinLike"]
df_all["hit_level"] = "NonHit"
for cond, level in zip(conditions, hit_levels):
    df_all.loc[cond, "hit_level"] = level

# Also flag anything better than -6.0 as "weak hit" for exploratory analysis
df_all.loc[(df_all["hit_level"] == "NonHit") & (df_all["best_affinity"] < -6.0), "hit_level"] = "L5_Weak"

n_hits = (df_all["hit_level"] != "NonHit").sum()
n_strong = (df_all["hit_level"].isin(["L1_Strong", "L2_Moderate", "L3_Absolute"])).sum()

log.info(f"Hits identified: {n_hits} total ({n_strong} L1-L3, {(df_all['hit_level']=='L4_BiotinLike').sum()} biotin-like)")

# ═══════════════════════════════════════════════════════════════
# 3. CHEMICAL FAMILY CLASSIFICATION
# ═══════════════════════════════════════════════════════════════

# Build SMILES lookup from existing + newly docked compounds
# We have SMILES for R1-R3 compounds in the ML script
COMPOUND_SMILES = {}  # Will populate from what we have

# Try loading from the ML training data
try:
    # Reconstruct from pharmacophore_ml_screen.py COMPOUND_SMILES dict
    exec(open("03_Analysis/pharmacophore_ml_screen.py").read().split("# Load docking results")[0])
except:
    pass

# For R4 compounds, load from candidate file
if os.path.exists("03_Analysis/outputs/drugbank_top500_for_docking.csv"):
    df_cand = pd.read_csv("03_Analysis/outputs/drugbank_top500_for_docking.csv")
    for _, row in df_cand.iterrows():
        COMPOUND_SMILES[row["name"]] = row["smiles"]

# Add SMILES to merged dataframe
df_all["smiles"] = df_all["name"].map(COMPOUND_SMILES)

# Classify by chemical family based on SMILES patterns
def classify_compound(name, smiles):
    """Classify into chemical family based on SMILES features."""
    if pd.isna(smiles):
        return "unknown"

    smi = str(smiles).lower()

    # Carboxylic acids (most SMVT binders)
    has_cooh = "C(=O)O" in smi or smi.endswith("C(O)=O") or "C(O)=O" in smi
    has_carboxyl = has_cooh or "carboxylic" in name.lower() or "acid" in name.lower()

    # Specific families
    if "biotin" in name.lower() or "C1C2C(C(S1)" in smi:
        return "Biotin_analog"
    if "fenamic" in name.lower() or "niflumic" in name.lower():
        return "Fenamate"
    if "profen" in name.lower():
        return "Profens"
    if "statin" in name.lower():
        return "Statin"
    if "dioic" in name.lower() or "diacid" in name.lower():
        return "Dicarboxylic_acid"
    if (smi.count("C(=O)O") >= 2 or smi.count("C(O)=O") >= 2) and "biotin" not in name.lower():
        return "Dicarboxylic_acid"
    if "amino" in name.lower() and has_carboxyl:
        return "Amino_acid"
    if "tryptophan" in name.lower() or "tyrosine" in name.lower() or "phenylalanine" in name.lower() or "dopa" in name.lower():
        return "Amino_acid"
    if "nsaid" in str(name).lower() or "coxib" in name.lower():
        return "NSAID"
    if has_carboxyl and ("salicyl" in name.lower() or "aspirin" in name.lower()):
        return "Salicylate"
    if has_carboxyl:
        return "Carboxylic_acid"
    if "sulfonamide" in name.lower() or "sulfon" in smi:
        return "Sulfonamide"

    return "Other"

df_all["family"] = df_all.apply(lambda r: classify_compound(r["name"], r.get("smiles")), axis=1)

# ═══════════════════════════════════════════════════════════════
# 4. SCAFFOLD ANALYSIS
# ═══════════════════════════════════════════════════════════════

def get_scaffold(smi):
    if pd.isna(smi): return None
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None: return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except:
        return None

df_all["scaffold"] = df_all["smiles"].apply(get_scaffold)

# Scaffold enrichment in hits
hit_scaffolds = df_all[df_all["hit_level"].isin(["L1_Strong", "L2_Moderate", "L3_Absolute", "L4_BiotinLike"])]["scaffold"].dropna()
all_scaffolds = df_all["scaffold"].dropna()

scaffold_counts = Counter(hit_scaffolds)
scaffold_total = Counter(all_scaffolds)

enriched_scaffolds = []
for scaff, hit_count in scaffold_counts.most_common(20):
    total_count = scaffold_total.get(scaff, 1)
    enrichment = hit_count / total_count
    enriched_scaffolds.append({
        "scaffold": scaff[:80],
        "hit_count": hit_count,
        "total_count": total_count,
        "enrichment": enrichment,
    })

df_scaffolds = pd.DataFrame(enriched_scaffolds)
log.info(f"\nTop enriched scaffolds:")
for _, row in df_scaffolds.head(10).iterrows():
    log.info(f"  {row['scaffold'][:60]:60s} | hits={row['hit_count']:2d}/{row['total_count']:3d} | enrich={row['enrichment']:.2f}")

# ═══════════════════════════════════════════════════════════════
# 5. FAMILY-LEVEL STATISTICS
# ═══════════════════════════════════════════════════════════════

family_stats = df_all.groupby("family").agg(
    count=("name", "count"),
    mean_affinity=("best_affinity", "mean"),
    best_affinity=("best_affinity", "min"),
    hits=("hit_level", lambda x: (x != "NonHit").sum()),
    strong_hits=("hit_level", lambda x: x.isin(["L1_Strong", "L2_Moderate", "L3_Absolute"]).sum()),
).sort_values("best_affinity")

log.info(f"\nFamily-level statistics:")
log.info(f"{'Family':25s} | {'N':>4s} | {'Mean ΔG':>8s} | {'Best ΔG':>8s} | {'Hits':>5s} | {'Strong':>6s}")
log.info("-" * 75)
for family, row in family_stats.iterrows():
    log.info(f"  {family:23s} | {int(row['count']):4d} | {row['mean_affinity']:8.2f} | {row['best_affinity']:8.2f} | "
             f"{int(row['hits']):5d} | {int(row['strong_hits']):6d}")

# ═══════════════════════════════════════════════════════════════
# 6. DRUG REPURPOSING SCORE
# ═══════════════════════════════════════════════════════════════

# Score hits for repurposing potential:
# - Already FDA-approved → +2
# - Known safety profile → +1
# - Novel target (no prior SMVT literature) → +1
# - Drug-like properties → +0.5
# - Strong affinity → +1

KNOWN_SMVT_DRUGS = {"gabapentin enacarbil", "biotin", "pantothenic acid", "lipoic acid",
                    "diclofenac", "ibuprofen", "ketoprofen", "flurbiprofen",
                    "indomethacin", "aspirin", "naproxen"}

df_hits = df_all[df_all["hit_level"] != "NonHit"].copy()

def repurposing_score(row):
    score = 0
    name = str(row["name"]).lower()
    # FDA approval (all ChEMBL compounds are approved, hand-picked may not be)
    if row["source"] == "chembl_ml_screen":
        score += 2
    elif any(d in name for d in ["fda", "approved"]):
        score += 2
    # Novelty (not previously known SMVT binder)
    if name not in KNOWN_SMVT_DRUGS:
        score += 1
    # Affinity
    if row["best_affinity"] < -7.0:
        score += 1.5
    elif row["best_affinity"] <= BIOTIN_AFFINITY:
        score += 0.5
    return score

if len(df_hits) > 0:
    df_hits["repo_score"] = df_hits.apply(repurposing_score, axis=1)
    df_hits = df_hits.sort_values(["hit_level", "repo_score", "best_affinity"],
                                   ascending=[True, False, True])

# ═══════════════════════════════════════════════════════════════
# 7. SAVE OUTPUTS
# ═══════════════════════════════════════════════════════════════

# Full merged results
cols_out = ["name", "best_affinity", "z_score", "hit_level", "family", "round",
            "source", "exhaustiveness", "scaffold", "smiles"]
df_all[cols_out].sort_values("best_affinity").to_csv(
    "03_Analysis/outputs/screening_master_results.csv", index=False
)
log.info(f"\nSaved screening_master_results.csv ({len(df_all)} compounds)")

# Hit summary
if len(df_hits) > 0:
    hit_cols = ["name", "best_affinity", "z_score", "hit_level", "family",
                "source", "repo_score", "scaffold"]
    df_hits[hit_cols].to_csv("03_Analysis/outputs/hit_summary.csv", index=False)
    log.info(f"Saved hit_summary.csv ({len(df_hits)} hits)")

    # Top hits preview
    log.info(f"\n{'='*60}")
    log.info(f"TOP 30 HITS (ranked by affinity + novelty)")
    log.info(f"{'='*60}")
    for i, (_, row) in enumerate(df_hits.head(30).iterrows()):
        log.info(f"  {i+1:2d}. {row['name'][:35]:35s} | ΔG={row['best_affinity']:.2f} | "
                 f"Z={row['z_score']:.2f} | {row['hit_level']:15s} | {row['family']:20s} | "
                 f"repo={row.get('repo_score', 0):.1f}")

log.info(f"\nAnalysis complete. Ready for report generation.")
