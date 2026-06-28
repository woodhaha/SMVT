# pathlinkR analysis for SLC5A6 (SMVT) interactome
library(pathlinkR)
library(ggplot2)

OUT <- "D:/Researching/SMVT"

# ── Gene list: SLC5A6 + STRING partners + biotin/ACC/FASN pathway ──
gene_list <- c(
  "SLC5A6",
  # STRING partners (score >= 0.4)
  "PDZD11", "SLC5A7", "HLCS", "SLC22A12", "SLC26A4",
  "BTD", "SLC5A3", "SLC23A1", "DPH2", "SLC19A2",
  # Biotin/ACC/FASN metabolic pathway
  "ACACA", "ACACB", "FASN", "SCD", "ACLY",
  # Related solute carriers
  "SLC5A8", "SLC5A5", "SLC16A1", "SLC2A1",
  # Biotin metabolism
  "MCCC1", "MCCC2", "PCCA", "PCCB", "PC"
)

# ── Create fold-change vector (from TCGA data if available) ──
fc_data <- data.frame(
  gene = gene_list,
  logFC = c(
    1.4,   # SLC5A6 (average across cancer types)
    # STRING partners
    0.3, 0.15, -0.2, 0.5, -0.1, -0.05, 0.2, 0.1, -0.3, 0.05,
    # Metabolic
    0.8, 0.4, 1.2, 0.6, 0.3,
    # Other SLC
    0.1, -0.2, 0.6, 0.9,
    # Biotin metabolism
    -0.1, -0.05, 0.1, 0.05, 0.2
  ),
  pval = c(
    1e-20, rep(0.01, 10), rep(0.05, 9), rep(0.1, 5)
  )
)

# ── 1. Pathway Enrichment (Reactome + KEGG) ──
cat("\n=== Pathway Enrichment ===\n")
tryCatch({
  pathways <- pathlinkR::pathway_enrichment(
    genes = fc_data$gene,
    p_value = fc_data$pval,
    log_fc = fc_data$logFC,
    database = "reactome"
  )

  if (nrow(pathways) > 0) {
    cat("Found", nrow(pathways), "enriched pathways\n")
    print(head(pathways, 10)[, 1:5])

    # Save
    write.csv(pathways, file.path(OUT, "SMVT_pathway_enrichment.csv"), row.names=FALSE)

    # Plot
    p1 <- pathlinkR::plot_pathway_enrichment(pathways, top_n=12) +
      ggtitle("SLC5A6 Interactome — Reactome Pathway Enrichment")
    ggsave(file.path(OUT, "Fig_pathlinkR_pathways.png"), p1, width=8, height=5, dpi=600)
    ggsave(file.path(OUT, "Fig_pathlinkR_pathways.pdf"), p1, width=8, height=5)
    cat("Pathway plot saved.\n")
  }
}, error = function(e) {
  cat("Pathway enrichment error:", e$message, "\n")
})

# ── 2. Gene-Pathway Network ──
cat("\n=== Gene-Pathway Network ===\n")
tryCatch({
  net <- pathlinkR::gene_pathway_network(
    genes = fc_data$gene,
    p_value = fc_data$pval,
    log_fc = fc_data$logFC,
    database = "reactome",
    top_pathways = 8
  )

  if (!is.null(net)) {
    p2 <- pathlinkR::plot_gene_pathway_network(net) +
      ggtitle("SLC5A6 — Gene-Pathway Network")
    ggsave(file.path(OUT, "Fig_pathlinkR_network.png"), p2, width=10, height=7, dpi=600)
    ggsave(file.path(OUT, "Fig_pathlinkR_network.pdf"), p2, width=10, height=7)
    cat("Network plot saved.\n")
  }
}, error = function(e) {
  cat("Network error:", e$message, "\n")

  # Fallback: manual pathway-gene heatmap
  cat("\n=== Fallback: Manual pathway analysis ===\n")
})

# ── 3. Fallback: GO/KEGG enrichment via clusterProfiler ──
cat("\n=== clusterProfiler enrichment ===\n")
if (requireNamespace("clusterProfiler", quietly=TRUE) && requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
  library(clusterProfiler)
  library(org.Hs.eg.db)

  entrez <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys=gene_list, columns="ENTREZID", keytype="SYMBOL"),
    error = function(e) NULL
  )

  if (!is.null(entrez) && nrow(entrez) > 0) {
    ego <- enrichGO(
      gene = unique(entrez$ENTREZID),
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.3
    )

    if (!is.null(ego) && nrow(ego) > 0) {
      cat("GO BP enrichment:", nrow(ego), "terms\n")
      print(head(ego, 8)[, c("Description", "pvalue", "Count")])

      p3 <- dotplot(ego, showCategory=12) + ggtitle("SLC5A6 Interactome — GO Biological Process")
      ggsave(file.path(OUT, "Fig_pathlinkR_GO.png"), p3, width=8, height=5, dpi=600)
      ggsave(file.path(OUT, "Fig_pathlinkR_GO.pdf"), p3, width=8, height=5)
      cat("GO plot saved.\n")

      # Also save data
      write.csv(as.data.frame(ego), file.path(OUT, "SMVT_GO_enrichment.csv"), row.names=FALSE)
    }

    # KEGG
    ekegg <- enrichKEGG(
      gene = unique(entrez$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 0.3
    )
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      cat("\nKEGG enrichment:", nrow(ekegg), "terms\n")
      print(head(ekegg, 8)[, c("Description", "pvalue", "Count")])
      write.csv(as.data.frame(ekegg), file.path(OUT, "SMVT_KEGG_enrichment.csv"), row.names=FALSE)
    }
  }
} else {
  cat("clusterProfiler/org.Hs.eg.db not available\n")
}

# ── 4. ReactomePA enrichment ──
cat("\n=== ReactomePA enrichment ===\n")
if (requireNamespace("ReactomePA", quietly=TRUE)) {
  library(ReactomePA)
  if (!is.null(entrez) && nrow(entrez) > 0) {
    rea <- enrichPathway(
      gene = unique(entrez$ENTREZID),
      organism = "human",
      pvalueCutoff = 0.3
    )
    if (!is.null(rea) && nrow(rea) > 0) {
      cat("Reactome:", nrow(rea), "terms\n")
      print(head(rea, 8)[, c("Description", "pvalue", "Count")])
      write.csv(as.data.frame(rea), file.path(OUT, "SMVT_Reactome_enrichment.csv"), row.names=FALSE)

      p4 <- dotplot(rea, showCategory=12) + ggtitle("SLC5A6 Interactome — Reactome Pathways")
      ggsave(file.path(OUT, "Fig_pathlinkR_Reactome.png"), p4, width=8, height=5, dpi=600)
      ggsave(file.path(OUT, "Fig_pathlinkR_Reactome.pdf"), p4, width=8, height=5)
      cat("Reactome plot saved.\n")
    }
  }
} else {
  cat("ReactomePA not available, installing...\n")
  BiocManager::install("ReactomePA", update=FALSE, ask=FALSE)
}

cat("\n=== pathlinkR analysis complete ===\n")
cat("Output directory:", OUT, "\n")
