###############################################################################
# scTenifoldKnk Virtual Knockout of SMVT (SLC5A6) on GSE178341
#
# Dataset: CRC 10X Full ~295K cells (GSE178341)
# Method: PC regression GRN -> virtual KO -> DRG ranking -> pathway enrichment
# Acceleration: data.table + future multisession (15 workers)
#
# Loads pre-saved 6.4GB input matrix RDS, subsets to
# 10K cells x 2000 HVG, then runs scTenifoldKnk.
###############################################################################

# -- 0. Setup -----------------------------------------------------------------
library(Seurat)
library(data.table)
library(future)
library(furrr)
library(scTenifoldKnk)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

plan(multisession, workers = 15)
setDTthreads(15)
set.seed(42)

BASE <- getwd()
OUTPUT_DIR <- file.path(BASE, "03_Analysis", "outputs")
FIGURE_DIR <- file.path(BASE, "03_Analysis", "figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

INTERIM_MATRIX <- file.path(OUTPUT_DIR, "scTenifoldKnk_input_matrix.rds")
RESULT_FILE   <- file.path(OUTPUT_DIR, "scTenifoldKnk_result.rds")
DRG_CSV       <- file.path(OUTPUT_DIR, "scTenifoldKnk_DRGs.csv")
PARTNER_CSV   <- file.path(OUTPUT_DIR, "scTenifoldKnk_partner_hits.csv")
GO_CSV        <- file.path(OUTPUT_DIR, "scTenifoldKnk_GO_enrichment.csv")
KEGG_CSV      <- file.path(OUTPUT_DIR, "scTenifoldKnk_KEGG_enrichment.csv")

NCELLS <- 10000
NHVG   <- 2000

Sys.setenv(OMP_NUM_THREADS = "2")
message("=== scTenifoldKnk SMVT Virtual KO -- GSE178341 ===")

# -- 1. Load matrix and subset -------------------------------------------------
message("(1/4) Loading pre-saved count matrix (RDS) and subsetting...")
counts_matrix <- readRDS(INTERIM_MATRIX)
message("  Full matrix: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " cells")

all_cells <- seq_len(ncol(counts_matrix))
keep_cells <- sample(all_cells, min(NCELLS, length(all_cells)))
counts_matrix <- counts_matrix[, keep_cells]
message("  Subsampled cells: ", ncol(counts_matrix))

keep_genes <- rowSums(counts_matrix > 0) >= 3
counts_matrix <- counts_matrix[keep_genes, ]
message("  Genes expressed >=3 cells: ", nrow(counts_matrix))

gene_vars <- apply(counts_matrix, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(NHVG, length(gene_vars))])
if (!"SLC5A6" %in% top_genes) {
  top_genes <- c("SLC5A6", top_genes[1:(length(top_genes)-1)])
}
counts_matrix <- counts_matrix[top_genes, ]
message("  Final matrix: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " cells")
remove(keep_genes, keep_cells, gene_vars, top_genes, all_cells)
gc()

# -- 2. Virtual KO ------------------------------------------------------------
message("(2/4) Running scTenifoldKnk virtual KO of SLC5A6...")

ko_result <- tryCatch({
  scTenifoldKnk::scTenifoldKnk(
    countMatrix = counts_matrix,
    gKO = "SLC5A6",
    nCores = 15,
    qc = TRUE
  )
}, error = function(e) {
  message("  scTenifoldKnk() with qc=TRUE failed: ", e$message)
  message("  Retrying with qc=FALSE...")
  scTenifoldKnk::scTenifoldKnk(
    countMatrix = counts_matrix,
    gKO = "SLC5A6",
    nCores = 15,
    qc = FALSE
  )
})

saveRDS(ko_result, RESULT_FILE)
message("  Saved result to: ", RESULT_FILE)

# -- 3. Extract DRGs and Partners ---------------------------------------------
message("(3/4) Extracting differentially regulated genes...")

drg_table <- as.data.table(ko_result$diffRegulation)
drg_table <- drg_table[order(-abs(distance))]

message("  Total DRGs: ", nrow(drg_table))
message("  DRGs with p < 0.05: ", sum(drg_table$p.value < 0.05, na.rm = TRUE))
message("  DRGs with p.adj < 0.05: ", sum(drg_table$p.adjust < 0.05, na.rm = TRUE))

top_200 <- head(drg_table, 200)
fwrite(top_200, DRG_CSV)
message("  Top 200 DRGs saved to: ", DRG_CSV)

partners <- c("PDZD11", "HLCS", "BTD", "SLC5A7", "SLC19A2",
              "SLC22A12", "SLC26A4", "SLC5A3", "SLC23A1", "DPH2")
partner_hits <- drg_table[drg_table$gene %in% partners, ]
if (nrow(partner_hits) > 0) partner_hits <- partner_hits[order(-abs(distance))]
fwrite(partner_hits, PARTNER_CSV)
message("  Partner hits: ", nrow(partner_hits), " / ", length(partners))
if (nrow(partner_hits) > 0) {
  print(partner_hits[, .(gene, distance, p.value, p.adjust)])
}

# -- 4. Enrichment, Figures, Report -------------------------------------------
message("(4/4) Pathway enrichment and figures...")

drg_500 <- head(drg_table$gene, 500)

ego <- tryCatch({
  enrichGO(gene = drg_500, OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL", ont = "BP",
           pAdjustMethod = "BH", qvalueCutoff = 0.05)
}, error = function(e) { message("  GO failed: ", e$message); NULL })

ekegg <- tryCatch({
  enrichKEGG(gene = drg_500, organism = "hsa",
             pAdjustMethod = "BH", qvalueCutoff = 0.05)
}, error = function(e) { message("  KEGG failed: ", e$message); NULL })

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  fwrite(as.data.frame(ego), GO_CSV)
  message("  GO terms: ", nrow(as.data.frame(ego)))
}
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  fwrite(as.data.frame(ekegg), KEGG_CSV)
  message("  KEGG pathways: ", nrow(as.data.frame(ekegg)))
}

# Figures
message("  Generating figures...")

drg_table[, sig_label := fifelse(p.adjust < 0.05, "FDR < 0.05",
                                  fifelse(p.value < 0.05, "p < 0.05", "Not significant"))]

p_volcano <- ggplot(drg_table, aes(x = distance, y = -log10(p.value),
                                    color = sig_label)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("FDR < 0.05" = "#DC0000",
                                "p < 0.05" = "#4DBBD5",
                                "Not significant" = "gray60")) +
  theme_minimal(base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "SMVT Virtual KO -- Differentially Regulated Genes",
       subtitle = paste0("scTenifoldKnk on GSE178341 (", ncol(counts_matrix), " cells)"),
       x = "Distance (KO vs Control)", y = "-log10(p-value)", color = "Significance") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURE_DIR, "scTenifoldKnk_DRG_volcano.png"),
       p_volcano, width = 8, height = 6, dpi = 150)

top20 <- head(drg_table, 20)
top20[, gene := factor(gene, levels = rev(gene))]
p_top20 <- ggplot(top20, aes(x = distance, y = gene, fill = -log10(p.value))) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 11) +
  labs(title = "Top 20 DRGs -- SMVT Virtual KO",
       x = "Distance", y = NULL, fill = "-log10(p)")
ggsave(file.path(FIGURE_DIR, "scTenifoldKnk_top20_DRGs.png"),
       p_top20, width = 8, height = 6, dpi = 150)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  tryCatch({
    p_go <- dotplot(ego, showCategory = 20) + ggtitle("GO BP -- SMVT KO")
    ggsave(file.path(FIGURE_DIR, "scTenifoldKnk_enrichment_dotplot.png"),
           p_go, width = 10, height = 8, dpi = 150)
  }, error = function(e) message("  GO dotplot failed: ", e$message))
}

tryCatch({
  png(file.path(FIGURE_DIR, "scTenifoldKnk_network.png"),
      width = 10, height = 8, units = "in", res = 150)
  scTenifoldKnk::plotKO(ko_result)
  dev.off()
  message("  Network plot saved")
}, error = function(e) message("  Network plot failed: ", e$message))

message("\n=== Analysis Complete ===")
message("Result: ", RESULT_FILE)
message("DRGs: ", DRG_CSV)
message("Figures: ", FIGURE_DIR)
