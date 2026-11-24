"""Stage the npj SBA submission files where SNAPP can actually take them.

Two constraints this exists to satisfy (both learned the hard way on prior SNAPP
submissions): the figure upload field accepts only .jpeg/.tiff/.eps - not the PDFs
the manuscript compiles against - and the browser automation can only read files
under the workspace root, so nothing under D:\\Researching can be uploaded directly.

Produces C:\\Users\\woodh\\Documents\\_snapp_npj\\ :
  manuscript_source.zip  tex + references.bib + .bbl + sn-jnl.cls + figures/*.pdf
  main_figures/Fig1..9.tiff   600 dpi LZW masters for the figure field
  cover_letter.txt       plain text (pandoc chokes on LaTeX math in .md -> .pdf)
  MANIFEST.txt           file list + sizes + MD5, so the upload can be verified
                         afterwards rather than eyeballed
"""
import hashlib
import os
import shutil
import zipfile

SRC = r"D:\Researching\SMVT\npj_submission"
FIGS = os.path.join(SRC, "figures")
MS_FIGS = r"D:\Researching\SMVT\04_Manuscript\figures"
CLS = r"C:\texlive\2026\texmf-dist\tex\latex\sn-jnl\sn-jnl.cls"
OUT = r"C:\Users\woodh\Documents\_snapp_npj"

os.makedirs(OUT, exist_ok=True)
for sub in ("main_figures", "_tex"):
    d = os.path.join(OUT, sub)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)

manifest = []


def note(path, label=""):
    h = hashlib.md5(open(path, "rb").read()).hexdigest()[:12]
    sz = os.path.getsize(path)
    manifest.append(f"{label or os.path.basename(path):44s} {sz:10,d} B  md5:{h}")


# 1. LaTeX source bundle - SNAPP compiles this, so it needs the class file and the
#    pre-built .bbl (the server does not run BibTeX).
tex_files = ["npj_manuscript.tex", "references.bib", "npj_manuscript.bbl"]
stage = os.path.join(OUT, "_tex")
for f in tex_files:
    shutil.copy2(os.path.join(SRC, f), stage)
shutil.copy2(CLS, os.path.join(stage, "sn-jnl.cls"))
pdfs = sorted(x for x in os.listdir(FIGS) if x.endswith(".pdf"))
fig_stage = os.path.join(stage, "figures")
os.makedirs(fig_stage, exist_ok=True)
for f in pdfs:
    shutil.copy2(os.path.join(FIGS, f), fig_stage)

zip_path = os.path.join(OUT, "manuscript_source.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(stage):
        for f in files:
            full = os.path.join(root, f)
            z.write(full, os.path.relpath(full, stage))
note(zip_path, "manuscript_source.zip")
manifest.append(f"    -> contains {len(tex_files)+1} source files + {len(pdfs)} figure PDFs")

# 2. Main figures as 600 dpi LZW TIFF (the only raster formats the field accepts).
main_stems = ["Fig1_SMVT_TCGA_pan_cancer", "Fig2_SMVT_mechanism",
              "Fig3_SMVT_mutation_landscape", "Fig4_SMVT_mutation_vs_expression",
              "Fig5_SMVT_STRING_network", "Fig6_SMVT_STRING_evidence_heatmap",
              "Fig7_SMVT_virtual_screening", "Fig8_SMVT_barbiturate_pharmacophore",
              "Fig9_expanded_screening_hits"]
missing = []
for stem in main_stems:
    tiff = os.path.join(MS_FIGS, stem + ".tiff")
    if not os.path.exists(tiff):
        missing.append(stem)
        continue
    dst = os.path.join(OUT, "main_figures", stem + ".tiff")
    shutil.copy2(tiff, dst)
    note(dst)
if missing:
    manifest.append("!! NO TIFF FOR: " + ", ".join(missing))

# 3. Cover letter as plain text.
md = open(os.path.join(SRC, "cover_letter_npj.md"), encoding="utf-8").read()
lines = []
for l in md.splitlines():
    l = l.replace("**", "").replace("`", "")
    s = l.strip()
    if s and set(s) <= {"-", "="}:            # markdown rule -> plain divider
        l = "-" * 66
    elif s.startswith("#"):                   # heading marker must not reach the editor
        l = s.lstrip("#").strip()
    lines.append(l)
txt = "\n".join(lines)
cl = os.path.join(OUT, "cover_letter.txt")
open(cl, "w", encoding="utf-8").write(txt)
note(cl, "cover_letter.txt")

open(os.path.join(OUT, "MANIFEST.txt"), "w", encoding="utf-8").write(
    "SNAPP staging for npj Systems Biology and Applications - SMVT/SLC5A6\n"
    + "-" * 66 + "\n" + "\n".join(manifest) + "\n")
print("\n".join(manifest))
print(f"\nstaged at {OUT}")
print("TIFF coverage:", f"{len(main_stems)-len(missing)}/{len(main_stems)}",
      "| missing: " + (", ".join(missing) if missing else "none"))
