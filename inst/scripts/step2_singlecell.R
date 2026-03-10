#!/usr/bin/env Rscript
# =============================================================================
# Step 2: Single-cell Analysis - Gene Expression Visualization
# Description: Select target gene and visualize expression in single-cell data
# Usage:
#   Rscript inst/scripts/step2_singlecell.R --gene EZH2
#   Rscript inst/scripts/step2_singlecell.R --gene EZH2 --data auto
#   Rscript inst/scripts/step2_singlecell.R --help
# =============================================================================

library(HCCMultiOmics)

# ============== Argument Parsing ==============
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("========== Step 2: Single-cell Gene Expression ==========\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript inst/scripts/step2_singlecell.R --gene EZH2           # Use local single-cell data\n")
  cat("  Rscript inst/scripts/step2_singlecell.R --gene EZH2 --data auto  # Auto-download single-cell data\n")
  cat("\n")
  cat("Description:\n")
  cat("  - Shows target gene expression across cell types\n")
  cat("  - Review plots in hcc_output/step2/\n")
  cat("  - Decide to use type 1 (mean expression) or type 2 (positive cell ratio) for downstream analysis\n")
  cat("\n")
  quit(save = "no")
}

parse_args <- function() {
  opts <- list(gene = NULL, data = "local")

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--gene", "-g")) {
      i <- i + 1
      if (i <= length(args)) opts$gene <- args[i]
    } else if (args[i] %in% c("--data", "-d")) {
      i <- i + 1
      if (i <= length(args)) opts$data <- args[i]
    }
    i <- i + 1
  }

  if (is.null(opts$gene)) {
    stop("Please specify --gene. Run --help for usage.")
  }
  opts
}

args <- parse_args()

OUTPUT_DIR <- "hcc_output/step2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("========== Step 2: Single-cell Gene Expression ==========")
message("")

TARGET_GENE <- args$gene
message(sprintf(">>> Target gene: %s", TARGET_GENE))
message("")

# ============== 1. Diagnostic Score Bar Plot ==============
message(">>> Step 2.1: Generating diagnostic score bar plot...")

res_diag_file <- "hcc_output/step1/diag_result.rds"
if (file.exists(res_diag_file)) {
  tryCatch({
    res_diag <- readRDS(res_diag_file)
    
    pdf(file.path(OUTPUT_DIR, sprintf("%s_diag_score.pdf", TARGET_GENE)), width = 8, height = 6)
    plot_gene_diag_score(res_diag, gene = TARGET_GENE)
    dev.off()
    message(sprintf("    Saved: %s_diag_score.pdf", TARGET_GENE))
  }, error = function(e) {
    message(sprintf("    Warning: Failed to generate diagnostic score plot - %s", e$message))
  })
} else {
  message(sprintf("    Warning: Diagnostic result file not found: %s", res_diag_file))
}

# ============== 2. Survival Curve ==============
message(">>> Step 2.2: Generating survival curve...")

# Load TCGA data
data_dir <- file.path(getwd(), "data")
if (file.exists(file.path(data_dir, "TCGA_Expr_Mat.rda"))) {
  load(file.path(data_dir, "TCGA_Expr_Mat.rda"), envir = environment())
  load(file.path(data_dir, "TCGA_Clin_Data.rda"), envir = environment())
}

pdf(file.path(OUTPUT_DIR, sprintf("%s_survival.pdf", TARGET_GENE)), width = 8, height = 6)
plot_gene_survival(gene = TARGET_GENE)
dev.off()
message(sprintf("    Saved: %s_survival.pdf", TARGET_GENE))

message("")

# ============== 3. Load Single-cell Data ==============
message(">>> Step 2.3: Loading single-cell data...")

SC_DATA_PATH <- "HCC_sc_data.rda"

download_sc_data <- function() {
  zenodo_id <- "18642056"
  zenodo_url <- sprintf("https://zenodo.org/records/%s/files/HCC_sc_data.rda", zenodo_id)
  message(sprintf(">>> Downloading from Zenodo (Record: %s)", zenodo_id))
  tryCatch({
    download.file(zenodo_url, SC_DATA_PATH, method = "wget", quiet = FALSE)
    message(">>> Download complete!")
  }, error = function(e) {
    stop(sprintf("Download failed: %s", e$message))
  })
}

if (!file.exists(SC_DATA_PATH)) {
  if (args$data == "auto") {
    download_sc_data()
  } else {
    alt_paths <- c("inst/scripts/HCC_sc_data.rda", "data/HCC_sc_data.rda")
    for (p in alt_paths) {
      if (file.exists(p)) { SC_DATA_PATH <- p; break }
    }
  }
}

if (!file.exists(SC_DATA_PATH)) {
  stop("Single-cell data file HCC_sc_data.rda not found")
}

message(sprintf(">>> Using single-cell data: %s", SC_DATA_PATH))

if (grepl("\\.rda$", SC_DATA_PATH, ignore.case = TRUE)) {
  env <- new.env()
  load(SC_DATA_PATH, envir = env)
  HCC_sc_data <- env$HCC_sc_data
} else {
  HCC_sc_data <- readRDS(SC_DATA_PATH)
}

message(sprintf(">>> Single-cell data: %d cells, %d features", ncol(HCC_sc_data), nrow(HCC_sc_data)))
message("")

# ============== 4. Single-cell Landscape ==============
message(">>> Step 2.4: Generating single-cell landscape...")

pdf(file.path(OUTPUT_DIR, "sc_landscape.pdf"), width = 12, height = 10)
plot_sc_landscape(HCC_sc_data)
dev.off()

# ============== 5. Gene Expression Visualization ==============
if (!TARGET_GENE %in% rownames(HCC_sc_data)) {
  stop(sprintf("Gene %s not found in single-cell data", TARGET_GENE))
}

message(sprintf(">>> Step 2.5: Visualizing %s expression...", TARGET_GENE))

# FeaturePlot + Violin
pdf(file.path(OUTPUT_DIR, sprintf("%s_landscape.pdf", TARGET_GENE)), width = 14, height = 6)
p_landscape <- plot_sc_gene(HCC_sc_data, gene = TARGET_GENE)
print(p_landscape)
dev.off()
message(sprintf("    Saved: %s_landscape.pdf", TARGET_GENE))

# Type 1: Mean Expression - Auto get highest expression cell type
res_type1 <- plot_sc_gene_stat(HCC_sc_data, gene = TARGET_GENE, plot_type = 1)
pdf(file.path(OUTPUT_DIR, sprintf("%s_stat_type1.pdf", TARGET_GENE)), width = 10, height = 6)
print(res_type1$plot)
dev.off()
message(sprintf("    Saved: %s_stat_type1.pdf (mean expression)", TARGET_GENE))
message(sprintf("    Highest expression cell type (type1): %s", res_type1$top_cell))

# Type 2: Positive Cell Ratio - Auto get highest expression cell type
res_type2 <- plot_sc_gene_stat(HCC_sc_data, gene = TARGET_GENE, plot_type = 2)
pdf(file.path(OUTPUT_DIR, sprintf("%s_stat_type2.pdf", TARGET_GENE)), width = 10, height = 6)
print(res_type2$plot)
dev.off()
message(sprintf("    Saved: %s_stat_type2.pdf (positive cell ratio)", TARGET_GENE))
message(sprintf("    Highest expression cell type (type2): %s", res_type2$top_cell))

message("")

# ============== 6. Save Results for Step3 ==============
# Save target gene
saveRDS(TARGET_GENE, file.path(OUTPUT_DIR, "target_gene.rds"))

# Save top expression cell types for both types
top_cells <- list(
  type1 = res_type1$top_cell,
  type2 = res_type2$top_cell
)
saveRDS(top_cells, file.path(OUTPUT_DIR, "top_cells.rds"))

message("")

message("==========================================================")
message(">>> Step 2 Complete!")
message(sprintf(">>> Results saved to: %s", OUTPUT_DIR))
message("")
message("Output files:")
message("  - sc_landscape.pdf        : Single-cell landscape")
message(sprintf("  - %s_diag_score.pdf    : Diagnostic score", TARGET_GENE))
message(sprintf("  - %s_survival.pdf    : Survival curve", TARGET_GENE))
message(sprintf("  - %s_landscape.pdf   : Gene FeaturePlot+Violin", TARGET_GENE))
message(sprintf("  - %s_stat_type1.pdf  : Mean expression statistics", TARGET_GENE))
message(sprintf("  - %s_stat_type2.pdf  : Positive cell ratio statistics", TARGET_GENE))
message("")
message(">>> Please select type 1 or type 2 based on results to run step3:")
message("   Rscript inst/scripts/step3_mechanism.R --type 1")
message("   Rscript inst/scripts/step3_mechanism.R --type 2")
message("==========================================================")
