# Lightweight scTenifoldKnk — smaller subset for faster completion
# Load cached matrix, subsample aggressively, run KO

library(future)
library(data.table)
library(Matrix)

plan(multisession, workers = 8)  # fewer workers, less memory
setDTthreads(8)

message("=== LOADING RDS ===")
mat <- readRDS("03_Analysis/outputs/scTenifoldKnk_input_matrix.rds")
message(sprintf("Loaded: %d genes x %d cells", nrow(mat), ncol(mat)))

# Aggressive subsample: 2000 cells + 1500 genes (much faster)
set.seed(42)
if (ncol(mat) > 2000) {
  keep_cells <- sample(colnames(mat), 2000)
  mat <- mat[, keep_cells]
}
# Keep top variable genes
vars <- apply(mat, 1, var)
keep_genes <- names(sort(vars, decreasing = TRUE)[1:1500])
keep_genes <- union(keep_genes, "SLC5A6")
mat <- mat[rownames(mat) %in% keep_genes, ]

message(sprintf("Subset: %d genes x %d cells", nrow(mat), ncol(mat)))

# Save light matrix
saveRDS(mat, "03_Analysis/outputs/scTenifoldKnk_light_matrix.rds", compress = FALSE)

message("=== RUNNING scTenifoldKnk ===")
library(scTenifoldKnk)

result <- scTenifoldKnk(
  countMatrix = mat,
  gKO = "SLC5A6",
  qc = TRUE,
  nCores = 8
)

message("=== SAVING RESULTS ===")
saveRDS(result, "03_Analysis/outputs/scTenifoldKnk_result.rds")

# Extract DRGs
drg <- result$diffRegulation
drg <- drg[order(-abs(drg$distance)), ]
fwrite(drg, "03_Analysis/outputs/scTenifoldKnk_DRGs.csv")

# Top 100 DRGs
top100 <- head(drg, 100)

# Check partner genes
partners <- c("PDZD11","HLCS","BTD","SLC5A7","SLC19A2","SLC22A12","SLC26A4","SLC5A3","SLC23A1","DPH2")
partner_hits <- drg[drg$gene %in% partners, ]
fwrite(partner_hits, "03_Analysis/outputs/scTenifoldKnk_partner_hits.csv")

# Enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

drg_genes <- head(drg$gene, 300)
ego <- enrichGO(drg_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
fwrite(as.data.frame(ego), "03_Analysis/outputs/scTenifoldKnk_GO_enrichment.csv")

# Quick report
sink("03_Analysis/outputs/scTenifoldKnk_report.md")
cat("# scTenifoldKnk SMVT Virtual KO — Results\n\n")
cat(sprintf("**Cells**: %d | **Genes**: %d | **Cores**: 8\n\n", ncol(mat), nrow(mat)))
cat("## Top 20 DRGs\n\n")
knitr::kable(head(drg, 20))
cat("\n## Partner Gene Effects\n\n")
knitr::kable(partner_hits)
cat(sprintf("\n## SLC5A6 rank in DRGs: %d / %d\n", which(drg$gene == "SLC5A6"), nrow(drg)))
cat("\n## Top GO Terms\n\n")
if (nrow(ego) > 0) knitr::kable(head(as.data.frame(ego), 10))
sink()

message("=== DONE ===")
