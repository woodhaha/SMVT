#!/usr/bin/env Rscript
# TCGA Survival Analysis — using TCGAbiolinks
# This pulls REAL clinical + expression data from GDC

if (!require("TCGAbiolinks")) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("TCGAbiolinks", ask = FALSE)
    library(TCGAbiolinks)
}
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")

library(TCGAbiolinks)
library(survival)
library(survminer)

setwd("D:/Researching/SMVT")
output_dir <- "03_Analysis/outputs"

# Cancer types to query
projects <- c(
    "TCGA-LUAD", "TCGA-LUSC", "TCGA-BLCA", "TCGA-COAD",
    "TCGA-LIHC", "TCGA-STAD", "TCGA-ESCA", "TCGA-BRCA",
    "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-PRAD",
    "TCGA-THCA", "TCGA-UCEC"
)

results <- list()

for (proj in projects) {
    cat(sprintf("\n[%s]\n", proj))

    # Query clinical data
    query_clin <- GDCquery(
        project = proj,
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR XML"
    )

    tryCatch({
        GDCdownload(query_clin, method = "client")
        clin <- GDCprepare_clinic(query_clin, clinical.info = "patient")

        if (is.null(clin) || nrow(clin) == 0) {
            cat("  No clinical data\n")
            next
        }

        # Get survival columns
        required <- c("patient_id", "days_to_death", "days_to_last_followup",
                      "vital_status", "age_at_initial_pathologic_diagnosis",
                      "gender", "tumor_stage")
        available <- intersect(required, names(clin))
        cat(sprintf("  Clinical fields: %d (%d patients)\n", length(available), nrow(clin)))

        # Query expression data for SLC5A6
        query_expr <- GDCquery(
            project = proj,
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification",
            workflow.type = "STAR - Counts",
            sample.type = "Primary Tumor"
        )

        # Only try to download if there are samples
        if (nrow(query_expr$results[[1]]) > 0) {
            # For efficiency, just get a subset
            n_samples <- min(20, nrow(query_expr$results[[1]]))
            query_expr$results[[1]] <- query_expr$results[[1]][1:n_samples, ]

            GDCdownload(query_expr, method = "client")
            expr <- GDCprepare(query_expr)

            if (!is.null(expr)) {
                # Check if SLC5A6 exists
                gene_name <- "SLC5A6"
                if (gene_name %in% rownames(expr)) {
                    cat(sprintf("  SLC5A6 expression: available\n"))
                } else {
                    cat("  SLC5A6 not found in expression matrix\n")
                }
            }
        } else {
            cat("  No expression files\n")
        }

    }, error = function(e) {
        cat(sprintf("  Error: %s\n", e$message))
    })
}

cat("\nDone.\n")
