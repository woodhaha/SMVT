"""Compute UMAP for GSE178341 (370,115 CRC cells) and cache coordinates.

Why: the GEO metadata for this atlas carries clinical fields (metatables) and cell-type
labels (cluster.csv) but NO embedding coordinates, so the UMAP has to be derived from
the 10x matrix. Coordinates are cached to a small CSV so figure work never re-runs the
~1-hour embedding.

All 370,115 cells enter the embedding - no subsampling. Memory is kept in check by
discarding the full expression matrix right after PCA; SLC5A6 values come from the
already-verified per-cell table (outputs/slc5a6_percell_gse178341.csv.gz), so the
embedding and the FigS23 statistics cannot drift apart.
"""
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.n_jobs = 8

RAW = r"D:\Researching\SMVT\02_Data\raw"
H5 = os.path.join(RAW, "GSE178341_crc10x_full_c295v4_submit.h5")
CLUSTER = os.path.join(RAW, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz")
PERCELL = r"D:\Researching\SMVT\03_Analysis\outputs\slc5a6_percell_gse178341.csv.gz"
OUT = r"D:\Researching\SMVT\03_Analysis\outputs\umap_gse178341"
N_PCS, N_NEIGH, N_TOP = 50, 15, 2000
os.makedirs(OUT, exist_ok=True)
t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


log(f"reading 10x h5: {os.path.basename(H5)}")
adata = sc.read_10x_h5(H5)
adata.var_names_make_unique()
log(f"loaded {adata.shape[0]} cells x {adata.shape[1]} genes")

cl = pd.read_csv(CLUSTER, index_col=0)
cl.index.name = "barcode"
meta = pd.read_csv(PERCELL, usecols=["barcode", "cpm", "log1p_cpm", "slc5a6_pos",
                                     "clTopLevel", "SPECIMEN_TYPE", "group"],
                   index_col="barcode")
log(f"cluster rows={len(cl)}  per-cell rows={len(meta)}")
common = adata.obs_names.intersection(meta.index)
if len(common) != len(meta):
    log(f"WARN barcode overlap {len(common)} of {len(meta)} per-cell rows; "
        f"restricting to the intersection")
adata = adata[adata.obs_names.isin(common)].copy()
meta = meta.loc[meta.index.intersection(adata.obs_names)]
cl = cl.reindex(adata.obs_names)
log(f"aligned to {adata.n_obs} cells")

adata.obs["group"] = meta["group"].values
adata.obs["clTopLevel"] = meta["clTopLevel"].values
adata.obs["SPECIMEN_TYPE"] = meta["SPECIMEN_TYPE"].values
adata.obs["slc5a6_cpm"] = meta["cpm"].values
adata.obs["slc5a6_log1p"] = meta["log1p_cpm"].values
adata.obs["slc5a6_pos"] = meta["slc5a6_pos"].values
if "cl295v11SubFull" in cl.columns:
    adata.obs["subtype"] = cl["cl295v11SubFull"].values

log("normalize_total -> log1p -> HVG")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP, subset=False, flavor="seurat")
hvg = adata.var["highly_variable"].values
log(f"HVG selected: {int(hvg.sum())}")

log("PCA (this is the heavy sparse step)")
# sc.pp.pca on an AnnData, not sc.tl.pca on a raw array: tl.pca's keyword is
# `return_info`, and guessing against the 1.12 signature cost one full 63 s load.
ahvg = adata[:, hvg].copy()
sc.pp.pca(ahvg, n_comps=N_PCS, svd_solver="randomized", random_state=0, zero_center=True)
adata.obsm["X_pca"] = ahvg.obsm["X_pca"]
evr = ahvg.uns.get("pca", {}).get("variance_ratio")
if evr is not None:
    log(f"PCA explained variance top5 = {np.round(np.asarray(evr)[:5], 4).tolist()}")
del ahvg
log(f"PCA stored: {adata.obsm['X_pca'].shape}")
log("neighbours (pynndescent)")
sc.pp.neighbors(adata, n_neighbors=N_NEIGH, n_pcs=N_PCS, random_state=0, use_rep="X_pca")
log("UMAP layout")
sc.tl.umap(adata, random_state=0)
emb = adata.obsm["X_umap"]
log(f"UMAP done; range x [{emb[:,0].min():.2f}, {emb[:,0].max():.2f}] "
    f"y [{emb[:,1].min():.2f}, {emb[:,1].max():.2f}]")

coords = pd.DataFrame({"barcode": adata.obs_names,
                       "UMAP_1": emb[:, 0], "UMAP_2": emb[:, 1]})
for c in ("group", "clTopLevel", "SPECIMEN_TYPE", "slc5a6_cpm", "slc5a6_log1p",
          "slc5a6_pos", "subtype"):
    if c in adata.obs:
        coords[c] = adata.obs[c].values
out_csv = os.path.join(OUT, "umap_coordinates.csv.gz")
coords.to_csv(out_csv, index=False, compression="gzip")
np.save(os.path.join(OUT, "umap_pca.npy"), adata.obsm["X_pca"].astype("float32"))
log(f"wrote {out_csv} ({os.path.getsize(out_csv)/1e6:.1f} MB), rows={len(coords)}")

# consistency guard: the cached table must reproduce the FigS23 numbers
chk = coords.groupby("group")["slc5a6_pos"].agg(["size", "mean"])
print("\npositivity by group (must match FigS23 / manuscript):")
print((chk.assign(pct=(100 * chk["mean"]).round(2))[["size", "pct"]]).to_string())
te = coords.loc[coords.group == "Tumor epithelial", "slc5a6_pos"].mean() * 100
ne = coords.loc[coords.group == "Normal epithelial", "slc5a6_pos"].mean() * 100
epi = coords.loc[coords.clTopLevel == "Epi", "slc5a6_pos"].mean() * 100
print(f"\ntumor epi {te:.2f}% | normal epi {ne:.2f}% | all epi {epi:.2f}% "
      f"(manuscript: 28.3 / 6.95 / 21.3)")
assert abs(te - 28.25) < 0.05 and abs(ne - 6.95) < 0.05 and abs(epi - 21.28) < 0.05, \
    "UMAP table does not reproduce FigS23 numbers - investigate before plotting"
log("consistency check PASSED")
