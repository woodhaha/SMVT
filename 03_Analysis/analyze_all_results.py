#!/usr/bin/env python3
"""
Post-docking analysis — merge all results + generate updated report
Runs after dock_parallel.py completes.
"""
import os, sys, json, logging, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
from collections import Counter

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/outputs", exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# ── 1. LOAD ALL DATA ──────────────────────────────────────────

# R1-R3: hand-picked 84 compounds (ex=8)
df_old = pd.read_csv("03_Analysis/outputs/docking_expanded_results.csv")
df_old["round"] = "R1-R3"
df_old["exhaustiveness"] = 8

# R4: ChEMBL ML-screened 356 compounds (ex=16)
df_r4 = pd.read_csv("03_Analysis/outputs/docking_batch_results.csv")
df_r4["round"] = "R4"
df_r4["exhaustiveness"] = 16

# FDA leftover: 318+ compounds (ex=16)
df_fda = pd.read_csv("03_Analysis/outputs/docking_fda_leftover_results.csv")
df_fda["round"] = "R5-FDA"
df_fda["exhaustiveness"] = 16

log.info(f"Loaded R1-R3: {len(df_old)} | R4: {len(df_r4)} | FDA: {len(df_fda)}")

# ── 2. MERGE ────────────────────────────────────────────────

df_all = pd.concat([
    df_old[["name", "best_affinity", "round", "exhaustiveness"]],
    df_r4[["name", "best_affinity", "round", "exhaustiveness"]],
    df_fda[["name", "best_affinity", "round", "exhaustiveness"]],
], ignore_index=True)

df_all = df_all[df_all["best_affinity"].notna()]
df_all = df_all[df_all["best_affinity"] < 0]  # Keep valid negative values
log.info(f"Total merged: {len(df_all)} compounds")

# ── 3. Z-SCORE NORMALIZATION ──────────────────────────────────

BIOTIN = -6.76
df_all["z_score"] = np.nan
for rnd in df_all["round"].unique():
    mask = df_all["round"] == rnd
    affs = df_all.loc[mask, "best_affinity"]
    if len(affs) > 1:
        df_all.loc[mask, "z_score"] = stats.zscore(affs)

# Hit levels
df_all["hit_level"] = "NonHit"
df_all.loc[df_all["z_score"] < -2.0, "hit_level"] = "L1_Strong"
df_all.loc[(df_all["z_score"] < -1.5) & (df_all["z_score"] >= -2.0), "hit_level"] = "L2_Moderate"
df_all.loc[(df_all["best_affinity"] < -7.0) & (df_all["z_score"] >= -1.5), "hit_level"] = "L3_Absolute"
df_all.loc[(df_all["best_affinity"] <= BIOTIN) & (df_all["best_affinity"] > -7.0) & (df_all["z_score"] >= -1.5), "hit_level"] = "L4_BiotinLike"
df_all.loc[(df_all["hit_level"] == "NonHit") & (df_all["best_affinity"] < -6.0), "hit_level"] = "L5_Weak"

# ── 4. RANKING ──────────────────────────────────────────────

df_all = df_all.sort_values("best_affinity").reset_index(drop=True)

n_elite = (df_all["best_affinity"] < -8.0).sum()
n_hits = (df_all["best_affinity"] < -7.0).sum()
n_biotin = (df_all["best_affinity"] <= BIOTIN).sum()
n_total_hits = (df_all["hit_level"] != "NonHit").sum()

log.info(f"Elite (<-8.0): {n_elite} | Hits (<-7.0): {n_hits} | Biotin-level (<={BIOTIN}): {n_biotin} | Total hits: {n_total_hits}")

# ── 5. TOP 40 ──────────────────────────────────────────────

print("\n" + "="*80)
print("TOP 40 HITS — SMVT Virtual Screening (All Rounds Merged)")
print("="*80)
print(f"{'Rank':<5} {'Compound':<30} {'DG':>7} {'Z':>7} {'Level':<14} {'Round':<8}")
print("-"*80)

for i, (_, row) in enumerate(df_all.head(40).iterrows()):
    print(f"{i+1:<5} {row['name'][:29]:<30} {row['best_affinity']:>7.2f} {row['z_score']:>7.2f} {row['hit_level']:<14} {row['round']:<8}")

# ── 6. BY ROUND SUMMARY ──────────────────────────────────────

print("\n" + "="*80)
print("BY ROUND SUMMARY")
print("="*80)
for rnd in ["R1-R3", "R4", "R5-FDA"]:
    subset = df_all[df_all["round"] == rnd]
    if len(subset) == 0: continue
    print(f"\n{rnd}: {len(subset)} compounds")
    print(f"  Best: {subset['best_affinity'].min():.2f} ({subset.iloc[subset['best_affinity'].argmin()]['name']})")
    print(f"  Elite (<-8): {(subset['best_affinity']<-8.0).sum()}")
    print(f"  Hits (<-7): {(subset['best_affinity']<-7.0).sum()}")
    print(f"  Mean DG: {subset['best_affinity'].mean():.2f} ± {subset['best_affinity'].std():.2f}")

# ── 7. SAVE ──────────────────────────────────────────────

df_all.to_csv("03_Analysis/outputs/screening_master_results_ALL.csv", index=False)
log.info("Saved: screening_master_results_ALL.csv")

# Top hits only
df_hits = df_all[df_all["hit_level"] != "NonHit"].copy()
df_hits.to_csv("03_Analysis/outputs/hit_summary_ALL.csv", index=False)
log.info(f"Saved: hit_summary_ALL.csv ({len(df_hits)} hits)")

print("\nDone.")
