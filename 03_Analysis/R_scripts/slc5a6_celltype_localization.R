# SLC5A6 cell-type localization from GSE178341 (Pelka 2021 CRC scRNA-seq)
# Extracts SLC5A6 expression per cell, maps to cell type + specimen, summarizes.
# Outputs: per-cell CSV + summary CSVs for figure generation.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR  <- if (length(args) >= 1) args[1] else "02_Data"
OUT_DIR   <- if (length(args) >= 2) args[2] else "03_Analysis/outputs"

H5_FILE      <- file.path(DATA_DIR, "raw", "GSE178341_crc10x_full_c295v4_submit.h5")
CLUSTER_FILE <- file.path(DATA_DIR, "GSE178341", "cluster.csv.gz")
META_FILE    <- file.path(DATA_DIR, "GSE178341", "metatables.csv.gz")

cat("=== Reading H5 (Read10X_h5) ===\n")
mat <- Read10X_h5(H5_FILE)
cat(sprintf("Matrix: %d genes x %d cells\n", nrow(mat), ncol(mat)))
if (!"SLC5A6" %in% rownames(mat)) stop("SLC5A6 not found in rownames")

barcodes <- colnames(mat)
slc5a6 <- as.numeric(mat["SLC5A6", , drop = FALSE])
total  <- as.numeric(Matrix::colSums(mat))
rm(mat); gc()
names(slc5a6) <- barcodes
names(total)  <- barcodes

# per-cell table
d <- data.table(barcode = barcodes, slc5a6 = slc5a6, total = total)
d[, cpm := slc5a6 / total * 1e4]
d[, log1p_cpm := log1p(cpm)]
d[, slc5a6_pos := slc5a6 > 0]

cat("=== Reading annotations ===\n")
clu <- fread(CLUSTER_FILE)
met <- fread(META_FILE)

cat(sprintf("cluster rows=%d cols=%s\n", nrow(clu), paste(names(clu), collapse=",")))
cat(sprintf("meta rows=%d cols=%s\n", nrow(met), paste(names(met)[1:6], collapse=",")))

# map barcode -> cell type (clTopLevel) and specimen type
clu_map <- clu[, .(barcode = sampleID, clTopLevel, clMidwayPr)]
met_map <- met[, .(barcode = cellID, SPECIMEN_TYPE)]
d <- merge(d, clu_map, by = "barcode", all.x = TRUE)
d <- merge(d, met_map, by = "barcode", all.x = TRUE)

# tumor vs normal epithelial
d[, group := ifelse(clTopLevel == "Epi" & SPECIMEN_TYPE == "T", "Tumor epithelial",
             ifelse(clTopLevel == "Epi" & SPECIMEN_TYPE == "N", "Normal epithelial",
             ifelse(clTopLevel == "Epi", "Other epithelial", clTopLevel)))]

# --- summary 1: per clTopLevel ---
sum1 <- d[, .(
  n_cells = .N,
  pct_positive = round(100 * mean(slc5a6_pos), 2),
  mean_cpm = round(mean(cpm), 4),
  mean_log1p_cpm = round(mean(log1p_cpm), 4),
  median_cpm = round(median(cpm), 4)
), by = clTopLevel][order(-mean_log1p_cpm)]

# --- summary 2: tumor vs normal epithelial ---
sum2 <- d[clTopLevel == "Epi", .(
  n_cells = .N,
  pct_positive = round(100 * mean(slc5a6_pos), 2),
  mean_cpm = round(mean(cpm), 4),
  mean_log1p_cpm = round(mean(log1p_cpm), 4)
), by = SPECIMEN_TYPE][order(SPECIMEN_TYPE)]

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(d,    file.path(OUT_DIR, "slc5a6_percell_gse178341.csv.gz"))
fwrite(sum1, file.path(OUT_DIR, "slc5a6_celltype_summary.csv"))
fwrite(sum2, file.path(OUT_DIR, "slc5a6_epi_tumor_normal_summary.csv"))

cat("\n=== per-clTopLevel summary ===\n"); print(sum1, n = 50)
cat("\n=== tumor vs normal epithelial ===\n"); print(sum2)
cat("\nDone. Per-cell CSV saved (", nrow(d), "cells).\n", sep = "")
