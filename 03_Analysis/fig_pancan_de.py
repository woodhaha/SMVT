"""Fig1: bar + volcano(swapped) + heatmap — FIXED mapping."""
import pandas as pd, numpy as np, json, csv, os, shutil
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
from fig_style import apply_style, PALETTE
apply_style()

OUT = r"D:/Researching/SMVT/03_Analysis/figures"
FINAL = r"D:/Researching/SMVT/04_Manuscript/figures"
de = pd.read_csv(r"D:/Researching/SMVT/03_Analysis/outputs/pancan_de.csv")
de = de.sort_values("log2FC", ascending=True); n = len(de)
colors = ["#c62828" if r["sig"]=="*" and r["log2FC"]>0 else "#1565c0" if r["sig"]=="*" and r["log2FC"]<0 else "#bdbdbd" for _, r in de.iterrows()]
ylabels = [r["cancer"]+" ("+str(r["n_t"])+"/"+str(r["n_n"])+")"+(" ***" if r["pval"]<1e-30 else " **" if r["pval"]<1e-10 else " *" if r["sig"]=="*" else "") for _, r in de.iterrows()]
cancers = list(de["cancer"])

t2a = {"Acute Myeloid Leukemia":"LAML","Adrenocortical Cancer":"ACC","Bladder Urothelial Carcinoma":"BLCA",
    "Brain Lower Grade Glioma":"LGG","Breast Invasive Carcinoma":"BRCA","Cervical & Endocervical Cancer":"CESC",
    "Cholangiocarcinoma":"CHOL","Colon Adenocarcinoma":"COAD","Diffuse Large B-Cell Lymphoma":"DLBC",
    "Esophageal Carcinoma":"ESCA","Glioblastoma Multiforme":"GBM","Head & Neck Squamous Cell Carcinoma":"HNSC",
    "Kidney Chromophobe":"KICH","Kidney Clear Cell Carcinoma":"KIRC","Kidney Papillary Cell Carcinoma":"KIRP",
    "Liver Hepatocellular Carcinoma":"LIHC","Lung Adenocarcinoma":"LUAD","Lung Squamous Cell Carcinoma":"LUSC",
    "Mesothelioma":"MESO","Ovarian Serous Cystadenocarcinoma":"OV","Pancreatic Adenocarcinoma":"PAAD",
    "Pheochromocytoma & Paraganglioma":"PCPG","Prostate Adenocarcinoma":"PRAD","Rectum Adenocarcinoma":"READ",
    "Sarcoma":"SARC","Skin Cutaneous Melanoma":"SKCM","Stomach Adenocarcinoma":"STAD",
    "Testicular Germ Cell Tumor":"TGCT","Thymoma":"THYM","Thyroid Carcinoma":"THCA",
    "Uterine Carcinosarcoma":"UCS","Uterine Corpus Endometrioid Carcinoma":"UCEC","Uveal Melanoma":"UVM"}

cache = r"D:/Researching/SMVT/03_Analysis/data/slc5a6_tn.json"
if os.path.exists(cache):
    T = json.load(open(cache)); tumors, normals = T["tumors"], T["normals"]
else:
    CRC = r"D:/Researching/CRCproject/02_Data/raw"
    pheno = pd.read_csv(CRC+"/TcgaTargetGTEX_phenotype.txt", sep="	", encoding="ISO-8859-1")
    s2c = dict(zip(pheno["sample"], pheno["primary disease or tissue"]))
    with open(CRC+"/TcgaTargetGtex_rsem_gene_fpkm", "r", encoding="ISO-8859-1") as f:
        reader = csv.reader(f, delimiter="	"); hdr = next(reader)
        for row in reader:
            if row[0].startswith("ENSG00000138074"):
                vals = {hdr[i]: float(v) for i,v in enumerate(row[1:],1) if v and float(v)>0}; break
    tumors, normals = {}, {}
    for sid, v in vals.items():
        if not sid.startswith("TCGA"): continue
        full = s2c.get(sid,""); ab = t2a.get(full,"")
        if not ab: continue
        lv = float(np.log2(v+1))
        if sid.endswith(("-01","-03")): tumors.setdefault(ab,[]).append(lv)
        elif sid.endswith(("-11","-10")): normals.setdefault(ab,[]).append(lv)
    json.dump({"tumors":tumors,"normals":normals}, open(cache,"w"))

fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.7], height_ratios=[0.5, 0.5], hspace=0.28, wspace=0.18)

ax1 = fig.add_subplot(gs[0, 0])
ax1.barh(range(n), de["log2FC"].values, color=colors, ec="white", lw=0.4, height=0.6)
ax1.axvline(0, color="#333", lw=0.8); ax1.axvline(0.2, color="#888", ls="--", lw=0.5, alpha=0.5)
ax1.set_yticks(range(n)); ax1.set_yticklabels(ylabels, fontsize=7.5)
ax1.set_xlabel("Log Fold Change", fontsize=9)
ax1.set_title("a  SLC5A6 Differential Expression", fontsize=10, fontweight="bold", loc="left")
ax1.set_xlim(de["log2FC"].min()-0.03, de["log2FC"].max()+0.18)
ax1.grid(axis="x", alpha=0.1); ax1.set_axisbelow(True); ax1.tick_params(axis="y", left=False, labelsize=7.5)
ax1.legend(handles=[Patch(color="#c62828",label="Up"),Patch(color="#1565c0",label="Down"),Patch(color="#bdbdbd",label="NS")], fontsize=6.5, loc="lower right")

ax2 = fig.add_subplot(gs[0, 1])
for _, r in de.iterrows():
    c = "#c62828" if r["sig"]=="*" and r["log2FC"]>0 else "#1565c0" if r["sig"]=="*" else "#bdbdbd"
    ax2.scatter(-np.log10(r["pval"]), r["log2FC"], s=15 if r["sig"]=="*" else 8, color=c, alpha=0.8)
for _, r in de.iterrows():
    if r["pval"] < 1e-15:
        ax2.annotate(r["cancer"], (-np.log10(r["pval"]), r["log2FC"]), fontsize=6, ha="right", va="center", xytext=(-3,0), textcoords="offset points")
ax2.axvline(-np.log10(0.05), color="#888", ls="--", lw=0.5, alpha=0.5)
ax2.set_xlabel("-log10(P-value)", fontsize=9); ax2.set_ylabel("Log Fold Change", fontsize=9)
ax2.set_title("b  Volcano Plot", fontsize=10, fontweight="bold", loc="left")
ax2.grid(alpha=0.1); ax2.axhline(0, color="#ccc", lw=0.3)

ax3 = fig.add_subplot(gs[1, :])
M, G = 60, 5; W = M+G+M
mat = np.full((len(cancers), W), np.nan)
tc, nc = [], []
for i, c in enumerate(cancers):
    tv = sorted(tumors.get(c,[]), reverse=True)[:M]
    nv = sorted(normals.get(c,[]), reverse=True)[:M]
    mat[i, :len(tv)] = tv; mat[i, M+G:M+G+len(nv)] = nv
    tc.append(len(tv)); nc.append(len(nv))
for s in [0, M+G]:
    b = mat[:, s:s+M]; o = np.argsort(np.nanmedian(b, axis=0))
    mat[:, s:s+M] = b[:, o]

cmp = plt.cm.YlOrRd.copy(); cmp.set_bad("#f0f0f0")
ax3.imshow(mat, aspect="auto", cmap=cmp, vmin=1.5, vmax=3.0)
ax3.axvline(M-0.5, color="#333", lw=1.5)
ax3.text(M/2, -0.8, "Tumor", ha="center", fontsize=9, fontweight="bold")
ax3.text(M+G+M/2, -0.8, "Normal", ha="center", fontsize=9, fontweight="bold")
ax3.set_yticks(range(len(cancers))); ax3.set_yticklabels(cancers, fontsize=7.5)
ax3.set_xticks([]); ax3.set_xlim(-0.5, W+2.5)
ax3.set_title("c  SLC5A6 Expression", fontsize=10, fontweight="bold", loc="left")
ax3.tick_params(axis="y", left=False, labelsize=7.5)

for i, c in enumerate(cancers):
    ax3.plot(W+0.3, i, "o", color="#c62828", ms=4)
    ax3.plot(W+1.0, i, "o", color="#1565c0", ms=4)
    ax3.text(W+1.8, i, "n="+str(tc[i])+"/"+str(nc[i]), va="center", fontsize=5.5, color="#555")

plt.colorbar(ax3.images[0], ax=ax3, shrink=0.4, pad=0.06).set_label("log2(FPKM+1)", fontsize=8)
ax3.legend(handles=[Line2D([0],[0],marker="o",color="w",mfc="#c62828",ms=5,label="Tumor mean"),Line2D([0],[0],marker="o",color="w",mfc="#1565c0",ms=5,label="Normal mean")], fontsize=6, loc="upper right")

plt.tight_layout()
fig.savefig(OUT+"/Fig1_SMVT_TCGA_pan_cancer.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT+"/Fig1_SMVT_TCGA_pan_cancer.png", dpi=300, bbox_inches="tight")
shutil.copy(OUT+"/Fig1_SMVT_TCGA_pan_cancer.pdf", FINAL+"/Fig1_SMVT_TCGA_pan_cancer.pdf")
shutil.copy(OUT+"/Fig1_SMVT_TCGA_pan_cancer.png", FINAL+"/Fig1_SMVT_TCGA_pan_cancer.png")
plt.close(); print("Fig1 done")
