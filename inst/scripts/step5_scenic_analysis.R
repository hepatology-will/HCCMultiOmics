#!/usr/bin/env Rscript
# =============================================================================
# Step 5: SCENIC Downstream Analysis
# Description: Import SCENIC results, perform TF enrichment, mediator screening, and visualization
# Note: Uses all malignant hepatocytes, not main trajectory cells
# Usage:
#   Rscript inst/scripts/step5_scenic_analysis.R --target EZH2 --downstream SLC7A11
# =============================================================================

library(HCCMultiOmics)
library(tidyverse)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)

get_standard_theme <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black", face = "plain"),
      axis.title = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, color = "black", face = "plain")
    )
}

# ============== Argument Parsing ==============
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("========== Step 5: SCENIC Downstream Analysis ==========\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript inst/scripts/step5_scenic_analysis.R --target EZH2 --downstream SLC7A11\n")
  cat("\n")
  cat("Parameters:\n")
  cat("  --target    Target gene (e.g., EZH2)\n")
  cat("  --downstream Downstream gene (e.g., SLC7A11)\n")
  cat("\n")
  quit(save = "no")
}

parse_args <- function() {
  opts <- list(target = "EZH2", downstream = "SLC7A11")

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--target", "-t")) {
      i <- i + 1
      if (i <= length(args)) opts$target <- args[i]
    }
    if (args[i] %in% c("--downstream", "-d")) {
      i <- i + 1
      if (i <= length(args)) opts$downstream <- args[i]
    }
    i <- i + 1
  }
  opts
}

args <- parse_args()
TARGET_GENE <- args$target
DOWNSTREAM_GENE <- args$downstream

OUTPUT_DIR <- "hcc_output/step5"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("========== Step 5: SCENIC Downstream Analysis ==========")
message("")

# ============== Load Malignant Hepatocyte Data ==============
message(">>> Step 5.1: Loading malignant hepatocyte data...")

# Prefer step3 output of all malignant hepatocytes (not main trajectory)
MALIGNANT_DATA_PATH <- "../../../hep_mal_reprocessed.rda"

if (!file.exists(MALIGNANT_DATA_PATH)) {
  # If not saved, try loading previous intermediate results
  alt_paths <- c(
    "../../../HCC_docker/hep_mal_reprocessed.rda",
    "inst/extdata/hep_mal_reprocessed.rda",
    "data/hep_mal_reprocessed.rda",
    "hcc_output/step3/hep_mal_reprocessed.rda"
  )
  for (p in alt_paths) {
    if (file.exists(p)) {
      MALIGNANT_DATA_PATH <- p
      break
    }
  }
}

if (!file.exists(MALIGNANT_DATA_PATH)) {
  stop("Malignant hepatocyte data not found. Please run step3 or ensure hep_mal_reprocessed.rda exists")
}

if (grepl("\\.rda$", MALIGNANT_DATA_PATH, ignore.case = TRUE)) {
  env <- new.env()
  load(MALIGNANT_DATA_PATH, envir = env)
  hep_malignant <- env$hep_mal_reprocessed
} else {
  hep_malignant <- readRDS(MALIGNANT_DATA_PATH)
}

message(sprintf(">>> Malignant hepatocyte count: %d", ncol(hep_malignant)))

# ============== Import SCENIC AUC Results ==============
message(">>> Step 5.2: Importing SCENIC AUC matrix...")

# First try step4 output
AUC_FILE <- "hcc_output/step4/AUC_matrix.csv"
if (!file.exists(AUC_FILE)) {
  # Fallback to existing results
  alt_auc_paths <- c(
    "../../../hep_mal_p1_AUC.csv",
    "inst/extdata/hep_mal_p1_AUC.csv"
  )
  for (p in alt_auc_paths) {
    if (file.exists(p)) {
      AUC_FILE <- p
      break
    }
  }
}

if (!file.exists(AUC_FILE)) {
  stop("SCENIC AUC file not found. Please run step4 first to generate AUC_matrix.csv")
}

message(">>> Loading SCENIC AUC matrix from: ", AUC_FILE)

hep_malignant <- import_scenic_results(
  seurat_obj = hep_malignant,
  auc_file = AUC_FILE
)

message(">>> SCENIC AUC imported successfully")
message(sprintf(">>> AUC Assay contains %d transcription factors", nrow(hep_malignant[["AUC"]]@data)))

# ============== TF Enrichment Analysis ==============
message(">>> Step 5.3: TF enrichment analysis...")

res_enrich <- run_tf_enrichment_analysis(hep_malignant, target_gene = TARGET_GENE)

pdf(file.path(OUTPUT_DIR, "tf_enrichment_lollipop.pdf"), width = 10, height = 8)
print(res_enrich$lollipop)
dev.off()

pdf(file.path(OUTPUT_DIR, "tf_enrichment_heatmap.pdf"), width = 12, height = 8)
print(res_enrich$heatmap)
dev.off()

# ============== Screen Mediator TFs ==============
message(">>> Step 5.4: Screening mediator transcription factors...")

screen_res <- screen_mediator_tfs(
  seurat_obj = hep_malignant,
  target_gene = TARGET_GENE,
  downstream_gene = DOWNSTREAM_GENE
)

pdf(file.path(OUTPUT_DIR, "mediator_screening.pdf"), width = 10, height = 10)
print(screen_res$plot)
dev.off()

write.csv(screen_res$results, file.path(OUTPUT_DIR, "mediator_tfs_results.csv"), row.names = FALSE)

top_tf <- screen_res$results$TF[1]
message(sprintf(">>> Best candidate TF: %s", top_tf))

# ============== Mediation Analysis ==============
message(">>> Step 5.5: Mediation analysis...")

med_res <- run_mediation_analysis(
  seurat_obj = hep_malignant,
  target_gene = TARGET_GENE,
  mediator_tf = top_tf,
  downstream_gene = DOWNSTREAM_GENE,
  n_sims = 1000
)

pdf(file.path(OUTPUT_DIR, "mediation_analysis.pdf"), width = 8, height = 6)
print(med_res$plot)
dev.off()

# ============== Visualization Validation ==============
message(">>> Step 5.6: Visualization validation...")

# Plot A: Target (RNA) vs TF (AUC)
p1 <- plot_earth_density(
  seurat_obj = hep_malignant,
  x_feature = TARGET_GENE,
  y_feature = top_tf,
  x_assay = "RNA",
  y_assay = "AUC"
) + ggplot2::labs(title = paste0("Step 1: ", TARGET_GENE, " activates ", top_tf))

# Plot B: TF (AUC) vs Downstream (RNA)
p2 <- plot_earth_density(
  seurat_obj = hep_malignant,
  x_feature = top_tf,
  y_feature = DOWNSTREAM_GENE,
  x_assay = "AUC",
  y_assay = "RNA"
) + ggplot2::labs(title = paste0("Step 2: ", top_tf, " activates ", DOWNSTREAM_GENE))

pdf(file.path(OUTPUT_DIR, "validation_panel.pdf"), width = 14, height = 6)
print(p1 + p2)
dev.off()

message("")
message("==========================================================")
message(">>> Step 5 Complete!")
message(sprintf(">>> Results saved to: %s", OUTPUT_DIR))
message("")
message("Output files:")
message("  - tf_enrichment_lollipop.pdf : TF enrichment lollipop plot")
message("  - tf_enrichment_heatmap.pdf  : TF enrichment heatmap")
message("  - mediator_screening.pdf     : Mediator TF screening quadrant plot")
message("  - mediator_tfs_results.csv  : Mediator TF results table")
message("  - mediation_analysis.pdf    : Mediation analysis plot")
message("  - validation_panel.pdf      : Validation panel")
message("")
message("==========================================================")
