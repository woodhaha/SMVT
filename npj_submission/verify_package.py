"""Machine verifier for the npj submission package: inventory, numbers, hygiene.

Deterministic checks only - each one recomputes from the source files rather than
trusting the manuscript text, so a stale number in the tex fails here instead of
being re-read by another human/LLM pass.

Run from npj_submission/:  python verify_package.py
"""
import hashlib
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, "npj_manuscript.tex")
FIGDIR = os.path.join(HERE, "figures")
MS_FIGDIR = os.path.join(HERE, os.pardir, "04_Manuscript", "figures")
ANALYSIS = os.path.join(HERE, os.pardir, "03_Analysis")
DATA = os.path.join(HERE, os.pardir, "02_Data", "HPA", "extract")

tex = open(TEX, encoding="utf-8").read()
fails, notes = [], []


def check(label, ok, detail=""):
    (notes if ok else fails).append(f"{'PASS' if ok else 'FAIL'}  {label}" + (f"  [{detail}]" if detail else ""))


# ---------- 1. figure inventory ----------
inc = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex)
missing = [p for p in inc if not os.path.exists(os.path.join(HERE, p))]
check(f"includegraphics resolve ({len(inc)} refs)", not missing, "; ".join(missing))

env_captions = re.findall(r"\\caption\{", tex)
legend_main = sorted({int(m) for m in re.findall(r"\\textbf\{Fig\.~(\d+)\}", tex)})
legend_si = sorted({int(m) for m in re.findall(r"\\textbf\{Fig\.~S(\d+)\}", tex)})
cited_si = sorted({int(m) for m in re.findall(r"Fig\.~S(\d+)", tex)})
cited_main = sorted({int(m) for m in re.findall(r"Fig\.~(\d+)[^0-9]", tex)})
check(f"legend list covers Fig 1-{max(legend_main)} contiguously",
      legend_main == list(range(1, max(legend_main) + 1)), str(legend_main))
check(f"legend list covers Fig S1-S{max(legend_si)} contiguously",
      legend_si == list(range(1, max(legend_si) + 1)),
      f"missing {sorted(set(range(1, max(legend_si)+1)) - set(legend_si))}")
orphan = sorted(set(cited_si) - set(legend_si))
check("every SI figure cited in body has a legend", not orphan, f"orphans: {orphan}")
decl = re.search(r"Supplementary Figures S1--S(\d+)", tex)
check("declared SI range matches the legends",
      decl and int(decl.group(1)) == max(legend_si),
      f"declared S1--{decl.group(1) if decl else '?'}, legends to S{max(legend_si)}")

# ---------- 2. numbers recomputed from source ----------
lin = pd.read_csv(os.path.join(ANALYSIS, "drug_sensitivity", "outputs",
                               "slc5a6_crispr_by_lineage.csv"))
def stat(lineage):
    r = lin[lin.OncotreeLineage == lineage].iloc[0]
    return round(float(r["mean"]), 3), int(r["n"])

for lineage, want_mean, want_n in [("Testis", -0.134, 4), ("Bowel", -0.017, 61),
                                    ("Lung", -0.006, 125)]:
    got = stat(lineage)
    check(f"DepMap {lineage.lower()} mean {want_mean:+.3f} n={want_n}",
          got == (want_mean, want_n), f"got {got}")

scr = pd.read_csv(os.path.join(ANALYSIS, "drug_sensitivity", "outputs",
                               "slc5a6_drug_response_screen.csv"))
tested = scr[scr.rho.notna()]
sig = tested[tested.fdr < 0.05]
check("GDSC2 286 compounds tested", len(tested) == 286, f"{len(tested)}")
check("GDSC2 27 at FDR<0.05", len(sig) == 27, f"{len(sig)}")
check("all significant |rho| <= 0.32", sig.rho.abs().max() <= 0.32,
      f"max {sig.rho.abs().max():.3f}")
cpi = tested[tested.drug == "CPI-637"].iloc[0]
check("CPI-637 rho +0.018 p 0.66", round(cpi.rho, 3) == 0.018 and round(cpi.p, 2) == 0.66,
      f"rho {cpi.rho:.4f}, p {cpi.p:.3f}")
cov = pd.read_csv(os.path.join(ANALYSIS, "drug_sensitivity", "outputs",
                               "ligand_coverage_GDSC2.csv"))
unscreened = int((cov.gdsc2_match == "NOT SCREENED").sum())
check("coverage table = 11 ligands, all NOT SCREENED",
      len(cov) == 11 and unscreened == 11, f"{len(cov)} rows, {unscreened} NOT SCREENED")

cp = pd.read_csv(os.path.join(DATA, "cancer_cptac_axis.tsv"), sep="\t")
cp["pval_adj"] = cp.pval_adj.astype(float); cp["logFC"] = cp.logFC.astype(float)
slc = cp[cp.Gene == "ENSG00000138074"]
col = slc[slc.Cancer == "Colon AC"].iloc[0]
lsc = slc[slc.Cancer == "Lung SQCC"].iloc[0]
check("CPTAC colon logFC +1.05 p 5.6e-17",
      round(col.logFC, 2) == 1.05 and abs(col.pval_adj - 5.64e-17) < 1e-18,
      f"{col.logFC:.4f}, {col.pval_adj:.2e}")
check("CPTAC lung squamous +0.50 p 8.4e-7",
      round(lsc.logFC, 2) == 0.50 and abs(lsc.pval_adj - 8.42e-7) < 1e-8,
      f"{lsc.logFC:.4f}, {lsc.pval_adj:.2e}")

ihc = pd.read_csv(os.path.join(DATA, "normal_ihc_data_axis.tsv"), sep="\t")
s = ihc[ihc.Gene == "ENSG00000138074"].copy()
rank = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
s["r"] = s.Level.map(rank)
best = s.loc[s.groupby("Tissue")["r"].idxmax(), "Level"].value_counts()
counts = {k: int(best.get(k, 0)) for k in ("High", "Medium", "Low", "Not detected")}
check("IHC 45 tissues: Medium 8 / Low 8 / ND 29 / High 0",
      s.Tissue.nunique() == 45 and counts["Medium"] == 8 and counts["Low"] == 8
      and counts["Not detected"] == 29 and counts["High"] == 0, str(counts))

ms = pd.read_csv(os.path.join(DATA, "ms_tissue_axis.tsv"), sep="\t")
ms["det"] = ms.Intensity.notna()
tab = ms.groupby("Gene name")["det"].agg(["sum", "size"])
check("MS SLC5A6 0/20 while axis enzymes 12-20/20",
      tab.loc["SLC5A6", "sum"] == 0 and tab.drop(index="SLC5A6")["sum"].min() >= 12,
      str(tab["sum"].to_dict()))

um = pd.read_csv(os.path.join(ANALYSIS, "outputs", "umap_gse178341",
                              "umap_coordinates.csv.gz"))
te = 100 * um.loc[um.group == "Tumor epithelial", "slc5a6_pos"].mean()
ne = 100 * um.loc[um.group == "Normal epithelial", "slc5a6_pos"].mean()
epi = 100 * um.loc[um.clTopLevel == "Epi", "slc5a6_pos"].mean()
check("UMAP 370,115 cells", len(um) == 370115, f"{len(um)}")
check("UMAP 28.2 / 6.9 / 21.3 percent",
      round(te, 1) == 28.2 and round(ne, 1) == 6.9 and round(epi, 1) == 21.3,
      f"{te:.2f} / {ne:.2f} / {epi:.2f}")

# ---------- 3. packaging hygiene ----------
# Only live lines count: the tex carries an explanatory comment that names
# \bibliographystyle precisely so nobody re-adds it, and a naive substring test
# would report that as a defect.
live = "\n".join(ln for ln in tex.splitlines() if not ln.lstrip().startswith("%"))
for bad, why in [("\\begin{abstract}", "PLOS abstract env left over"),
                 ("authblk", "PLOS author package left over"),
                 ("Author Summary", "PLOS-only section left over"),
                 ("\\bibliographystyle", "doubles the class's \\bibstyle and breaks BibTeX")]:
    check(f"no {bad} ({why})", bad not in live)

m = re.search(r"\\abstract\{(.*?)\n\}", tex, re.S)
n_words = len(m.group(1).split()) if m else -1
check("abstract <= 150 words", 0 < n_words <= 150, f"{n_words} words")

for sec in ["Acknowledgements", "Author contributions", "Competing interests",
            "Data availability", "Code availability", "Ethics",
            "Use of Generative AI"]:
    check(f"Nature-required statement present: {sec}", sec in tex)

stale = []
for f in os.listdir(FIGDIR):
    if not f.endswith(".pdf"):
        continue
    p = os.path.join(FIGDIR, f)
    q = os.path.join(MS_FIGDIR, f)
    if os.path.exists(q):
        h1 = hashlib.md5(open(p, "rb").read()).hexdigest()
        h2 = hashlib.md5(open(q, "rb").read()).hexdigest()
        if h1 != h2:
            stale.append(f)
check("submission figures match the canonical renders (MD5)", not stale, ", ".join(stale))

orcid = re.findall(r"000[09]-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9Xx]", open(
    os.path.join(HERE, "authors.md"), encoding="utf-8").read())
check("corresponding ORCID is the Springer-bound 0009 id",
      "0009-0005-1818-4789" in orcid, str(sorted(set(orcid))))

# ---------- report ----------
print("\n".join(notes))
print("-" * 70)
if fails:
    print("\n".join(fails))
    print(f"\nVERDICT: {len(fails)} FAILURE(S) of {len(fails) + len(notes)} checks")
    sys.exit(1)
print(f"VERDICT: all {len(notes)} checks pass")
