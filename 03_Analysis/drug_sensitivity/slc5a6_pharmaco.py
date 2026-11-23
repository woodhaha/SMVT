"""SLC5A6 pharmacogenomics layer: essentiality + unbiased expression-response screen.

Context that shaped this script (verified 2026-11-23 against the local files):

* The planned headline test - "do SLC5A6-high lines respond differently to the
  predicted SMVT ligands?" - is NOT answerable with local data. GDSC2 here holds
  286 oncology/kinase-focused compounds and screens none of the manuscript's
  ligands (diclofenac, carprofen, phenobarbital, pentobarbital, butalbital,
  furosemide, naftazone, hydromorphone, biotin, 5-HTP, gabapentin enacarbil all
  return zero matches in DRUG_NAME). We report that as the coverage statement
  rather than quietly substituting a proxy drug.
* Two things ARE answerable, and both bear on the manuscript's central claim
  that SLC5A6 is an "expression-driven metabolic vulnerability":
    (1) CRISPR essentiality of SLC5A6 itself, per lineage. If SLC5A6 is not a
        dependency anywhere, "vulnerability" has to be softened.
    (2) An unbiased Spearman screen of SLC5A6 expression against every screened
        compound's AUC. Hypothesis-free by design; we then look for mechanism in
        whatever comes out. One mechanistically interesting compound is present:
        CPI-637 (devimistat) acts on the lipoamide/lipoate axis, and lipoate is
        an SMVT substrate - so it earns a pre-registered look, not a fishing run.

Inputs are local, read by column index so the 461 MB / 19,194-column DepMap
matrices never enter memory in full:
  DOCTOR/Paper/depmap_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv  ("SLC5A6 (8884)")
  DOCTOR/Paper/depmap_data/CRISPRGeneEffect.csv                          ("SLC5A6 (8884)")
  DOCTOR/Paper/depmap_data/Model.csv                                     (ModelID, OncotreeLineage, COSMICID)
  ImmunoTargets/data/GDSC_temp/GDSC2_fitted_dose_response.csv             (COSMIC_ID, DRUG_NAME, AUC)

No downsampling: every cell line with data enters every statistic.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DEPMAP = r"D:\Researching\DOCTOR\Paper\depmap_data"
GDSC = r"D:\Researching\ImmunoTargets\data\GDSC_temp\GDSC2_fitted_dose_response.csv"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
GENE_COL = "SLC5A6 (8884)"
# Chronos/gene-effect convention: <= -0.5 common essential, <= -0.7 highly, ~0 neutral
ESSENTIAL_CUT = -0.5

os.makedirs(OUT, exist_ok=True)


def read_depmap_column(path, column_name):
    """Pull one gene column out of a huge wide DepMap CSV by header index."""
    header = pd.read_csv(path, nrows=0)
    if column_name not in header.columns:
        raise SystemExit(f"FATAL: {column_name!r} absent from {os.path.basename(path)}"
                         f" ({len(header.columns)} columns). Not guessing an alternative.")
    idx = list(header.columns).index(column_name)
    df = pd.read_csv(path, usecols=[0, idx], index_col=0)
    df.columns = [column_name]
    print(f"  {os.path.basename(path)}: {len(df)} models, column {idx} -> {column_name}")
    return df[column_name]


def main():
    model = pd.read_csv(os.path.join(DEPMAP, "Model.csv"),
                        usecols=["ModelID", "StrippedCellLineName", "OncotreeLineage",
                                 "SangerModelID", "COSMICID"])
    expr = read_depmap_column(os.path.join(DEPMAP, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), GENE_COL)
    crispr = read_depmap_column(os.path.join(DEPMAP, "CRISPRGeneEffect.csv"), GENE_COL)

    # ---------- 1. essentiality ----------
    ess = pd.concat([model.set_index("ModelID"),
                    crispr.rename("gene_effect").to_frame()], axis=1, join="inner")
    ess = ess[ess["gene_effect"].notna()]
    by_lin = (ess.groupby("OncotreeLineage")["gene_effect"]
                .agg(n="size", mean="mean", median="median")
                .assign(pct_essential=lambda d: (100 * (d["mean"] <= ESSENTIAL_CUT)))
                )
    # fraction of lines within each lineage that are "essential"
    frac = (ess.assign(flag=ess["gene_effect"] <= ESSENTIAL_CUT)
               .groupby("OncotreeLineage")["flag"].mean() * 100)
    by_lin = by_lin.join(frac.rename("pct_lines_essential")).sort_values("mean")
    print(f"\n[1] SLC5A6 CRISPR essentiality: n={len(ess)} models, "
          f"overall mean gene effect {ess['gene_effect'].mean():.3f}, "
          f"median {ess['gene_effect'].median():.3f}, "
          f"{(ess['gene_effect'] <= ESSENTIAL_CUT).sum()} of {len(ess)} lines at/below {ESSENTIAL_CUT} "
          f"({100*(ess['gene_effect'] <= ESSENTIAL_CUT).mean():.1f}%)")
    print(by_lin.head(12).round(3).to_string())
    by_lin.to_csv(os.path.join(OUT, "slc5a6_crispr_by_lineage.csv"))

    # ---------- 2. expression-response screen ----------
    # Crosswalk: join on Sanger IDs, not COSMIC. In Model.csv, COSMICID parses as
    # float ("905933.0") and is populated for only 977/1959 models, while
    # SangerModelID ("SIDM00105") matches GDSC2's SANGER_MODEL_ID for 1216 models.
    gd = pd.read_csv(GDSC, usecols=["SANGER_MODEL_ID", "DRUG_NAME", "AUC",
                                    "PUTATIVE_TARGET", "PATHWAY_NAME"])
    xwalk = model[["ModelID", "SangerModelID", "OncotreeLineage"]].dropna(subset=["SangerModelID"])
    gd = gd.merge(xwalk, left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")
    n_unmapped = gd["ModelID"].isna().sum()
    mapped_lines = gd["SANGER_MODEL_ID"].nunique() - gd.loc[gd["ModelID"].isna(), "SANGER_MODEL_ID"].nunique()
    gd = gd.dropna(subset=["ModelID", "AUC"])
    if gd.empty:
        raise SystemExit("FATAL: no GDSC cell lines mapped to DepMap IDs. "
                         "Check the SANGER_MODEL_ID <-> SangerModelID key before trusting any screen.")
    e = expr.reset_index().rename(columns={"index": "ModelID"})
    if "ModelID" not in e.columns:            # name it if the CSV had no header for col 0
        e.columns = ["ModelID"] + list(e.columns[1:])
    gd = gd.merge(e, on="ModelID", how="inner")
    print(f"\n[2] GDSC2 rows: {len(pd.read_csv(GDSC, usecols=['DRUG_NAME']))} read, "
          f"{n_unmapped} rows with no COSMIC->DepMap mapping, "
          f"{len(gd)} rows usable (drug x line with expression), "
          f"{gd['DRUG_NAME'].nunique()} distinct drugs")

    rows = []
    for drug, sub in gd.groupby("DRUG_NAME"):
        sub = sub.dropna(subset=["SLC5A6 (8884)", "AUC"])
        n = len(sub)
        if n < 30:                    # stated floor; compounds below it are listed, not tested
            rows.append({"drug": drug, "n": n, "rho": np.nan, "p": np.nan,
                         "note": "n<30 floor, not tested"})
            continue
        rho, p = stats.spearmanr(sub["SLC5A6 (8884)"], sub["AUC"])
        rows.append({"drug": drug, "n": n, "rho": rho, "p": p,
                     "targets": ";".join(sorted(set(sub["PUTATIVE_TARGET"].dropna().astype(str)))[:4]),
                     "pathway": ";".join(sorted(set(sub["PATHWAY_NAME"].dropna().astype(str)))[:3]),
                     "note": ""})
    res = pd.DataFrame(rows)
    tested = res["rho"].notna()
    res.loc[tested, "fdr"] = stats.false_discovery_control(res.loc[tested, "p"].values)
    res = res.sort_values("rho")

    # Higher AUC = less sensitive. So rho<0 => SLC5A6-high lines are MORE sensitive.
    res["direction"] = np.where(res["rho"] < 0, "SLC5A6-high = more sensitive",
                                np.where(res["rho"] > 0, "SLC5A6-high = more resistant", ""))
    res.to_csv(os.path.join(OUT, "slc5a6_drug_response_screen.csv"), index=False)
    sig = res[tested & (res["fdr"] < 0.05)]
    print(f"\n    tested {tested.sum()} drugs, BH-FDR<0.05: {len(sig)} "
          f"(more-sensitive end {len(sig[sig['rho']<0])}, resistant end {len(sig[sig['rho']>0])})")
    print(res.head(12)[["drug", "n", "rho", "p", "fdr", "pathway"]].round(4).to_string(index=False))
    print("    ...")
    print(res.tail(12)[["drug", "n", "rho", "p", "fdr", "pathway"]].round(4).to_string(index=False))

    # ---------- 3. coverage statement for the manuscript's ligands ----------
    ligands = {  # name -> what the manuscript claims about it
        "diclofenac": "NSAID hit, dG -8.07", "carprofen": "NSAID hit, dG -8.04",
        "phenobarbital": "barbiturate hit, dG -8.30", "pentobarbital": "barbiturate hit, dG -8.18",
        "butalbital": "barbiturate hit, dG -8.02", "furosemide": "sulfonamide hit, dG -8.36",
        "naftazone": "naphthoquinone hit, dG -8.34", "hydromorphone": "opioid top hit, dG -8.58",
        "biotin": "natural substrate", "5-hydroxytryptophan": "5-HTP, dG -7.66",
        "gabapentin enacarbil": "FDA SMVT prodrug positive control",
    }
    present = sorted(set(gd["DRUG_NAME"].astype(str).str.lower()))
    cov = []
    for lig, claim in ligands.items():
        hits = [d for d in gd["DRUG_NAME"].astype(str).unique() if lig.split()[0].lower() in d.lower()]
        cov.append({"manuscript_ligand": lig, "claim": claim,
                    "gdsc2_match": "; ".join(hits) if hits else "NOT SCREENED"})
    coverage = pd.DataFrame(cov)
    coverage.to_csv(os.path.join(OUT, "ligand_coverage_GDSC2.csv"), index=False)
    print("\n[3] coverage of the manuscript's predicted ligands in GDSC2:")
    print(coverage.to_string(index=False))

    # lipoate-axis compound present? pre-registered look, not post-hoc
    cpi = res[res["drug"].str.contains("CPI-637", case=False, na=False)]
    if len(cpi):
        print("\n[4] pre-registered lipoate-axis look (CPI-637/devimistat):")
        print(cpi[["drug", "n", "rho", "p", "fdr", "direction"]].round(4).to_string(index=False))
    return res, by_lin, coverage


if __name__ == "__main__":
    main()
