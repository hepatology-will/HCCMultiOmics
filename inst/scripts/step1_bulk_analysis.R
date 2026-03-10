#!/usr/bin/env Rscript
# =============================================================================
# Step 1: Bulk RNA-seq Machine Learning Analysis
# Description: Select gene set and run machine learning modeling
# Usage:
#   Rscript inst/scripts/step1_bulk_analysis.R --builtin 1
#   Rscript inst/scripts/step1_bulk_analysis.R --custom genes.csv
#   Rscript inst/scripts/step1_bulk_analysis.R --builtin 1 --seed 123
#   Rscript inst/scripts/step1_bulk_analysis.R --help
# =============================================================================

library(HCCMultiOmics)

# ============== Argument Parsing ==============
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("========== Step 1: Bulk RNA-seq Analysis ==========\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript inst/scripts/step1_bulk_analysis.R --builtin 1        # Use built-in gene set\n")
  cat("  Rscript inst/scripts/step1_bulk_analysis.R --custom genes.csv # Custom gene set file\n")
  cat("  Rscript inst/scripts/step1_bulk_analysis.R --seed 123          # Set random seed\n")
  cat("\n")
  cat("Built-in Gene Sets:\n")
  cat("  1. Ferroptosis_FerrDb\n")
  cat("  2. Cuproptosis_FerrDb\n")
  cat("  3. Disulfidptosis\n")
  cat("  4. Autosis\n")
  cat("  5. immunogonic_cell_death\n")
  cat("  6. mitotic_death\n")
  cat("  7. parthanatos\n")
  cat("\n")
  cat("Supported Custom Formats:\n")
  cat("  - CSV/TXT: with columns all_genes, drivers, suppressors, syml\n")
  cat("  - RDS/RData: saved GeneSet object\n")
  cat("  - Plain text: one gene per line\n")
  cat("\n")
  quit(save = "no")
}

parse_args <- function() {
  opts <- list(builtin_index = 1, custom_file = NULL, seed = 123, mode = NULL)

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--builtin", "-b")) {
      i <- i + 1
      if (i <= length(args)) {
        opts$mode <- "builtin"
        opts$builtin_index <- as.integer(args[i])
      }
    } else if (args[i] %in% c("--custom", "-c")) {
      i <- i + 1
      if (i <= length(args)) {
        opts$mode <- "custom"
        opts$custom_file <- args[i]
      }
    } else if (args[i] %in% c("--seed", "-s")) {
      i <- i + 1
      if (i <= length(args)) {
        opts$seed <- as.integer(args[i])
      }
    }
    i <- i + 1
  }

  if (is.null(opts$mode)) {
    stop("Please specify --builtin or --custom. Run --help for usage.")
  }
  opts
}

args <- parse_args()

OUTPUT_DIR <- "hcc_output/step1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("========== Step 1: Bulk RNA-seq Analysis ==========")
message("")

# ============== Load Gene Set ==============
available_sets <- c(
  "Ferroptosis_FerrDb",
  "Cuproptosis_FerrDb",
  "Disulfidptosis",
  "Autosis",
  "immunogonic_cell_death",
  "mitotic_death",
  "parthanatos"
)

load_builtin_set <- function(index) {
  if (index < 1 || index > length(available_sets)) {
    stop(sprintf("Invalid gene set index: %d", index))
  }

  gene_set_name <- available_sets[index]
  data(list = gene_set_name, envir = environment())
  gene_set <- get(gene_set_name)

  message(sprintf(">>> Using built-in gene set: %s", gene_set_name))
  list(gene_set = gene_set, name = gene_set_name)
}

load_custom_set <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }

  ext <- tolower(path.expand(file_path))

  if (endsWith(ext, ".rds") || endsWith(ext, ".rda") || endsWith(ext, ".rdata")) {
    e <- new.env()
    load(file_path, envir = e)
    if (!exists("gene_set", envir = e)) {
      stop("'gene_set' object not found in file")
    }
    gene_set <- e$gene_set
    message(sprintf(">>> Loaded custom gene set: %s", basename(file_path)))
  } else {
    df <- tryCatch(
      readr::read_csv(file_path, show_col_types = FALSE),
      error = function(e) read.table(file_path, header = TRUE, sep = ",", fill = TRUE)
    )

    if (all(c("all_genes", "drivers", "suppressors", "syml") %in% colnames(df))) {
      gene_set <- create_gene_set(
        name = "Custom",
        all_genes = df$all_genes,
        drivers = df$drivers,
        suppressors = df$suppressors,
        syml = df$syml
      )
    } else if (all(c("all_genes", "drivers", "suppressors") %in% colnames(df))) {
      gene_set <- create_gene_set(
        name = "Custom",
        all_genes = df$all_genes,
        drivers = df$drivers,
        suppressors = df$suppressors
      )
    } else if ("all_genes" %in% colnames(df)) {
      gene_set <- create_gene_set(
        name = "Custom",
        all_genes = df$all_genes
      )
    } else {
      genes <- df[[1]]
      gene_set <- create_gene_set(
        name = "Custom",
        all_genes = genes
      )
    }
    message(sprintf(">>> Created custom gene set: %s", basename(file_path)))
  }

  list(gene_set = gene_set, name = basename(file_path))
}

if (args$mode == "builtin") {
  gene_info <- load_builtin_set(args$builtin_index)
} else {
  gene_info <- load_custom_set(args$custom_file)
}

gene_set <- gene_info$gene_set

message("")
message(sprintf(">>> Gene set: %s", gene_info$name))
message(sprintf(">>> Total genes: %d", length(gene_set$all_genes)))
if (!is.null(gene_set$drivers)) {
  message(sprintf(">>> Driver genes: %d", length(gene_set$drivers)))
}
if (!is.null(gene_set$suppressors)) {
  message(sprintf(">>> Suppressor genes: %d", length(gene_set$suppressors)))
}
if (!is.null(gene_set$syml)) {
  message(sprintf(">>> Syml genes: %d", length(gene_set$syml)))
}
message("")

print(gene_set)
message("")

ML_SEED <- args$seed
message(sprintf(">>> Random seed: %d", ML_SEED))
message("")
message(">>> Starting analysis...")

# ============== Run Analysis ==============
res <- run_intersection(gene_set)
message(sprintf(">>> Candidate gene count: %d", length(res$candidates)))

# Output PDF
pdf(file.path(OUTPUT_DIR, "venn_diagram.pdf"), width = 8, height = 8)
plot_venn(res)
dev.off()

message(">>> Starting machine learning modeling (time-consuming)...")
res_ml <- run_prognostic_model(res$candidates, seed = ML_SEED)

pdf(file.path(OUTPUT_DIR, "lasso_cv.pdf"), width = 8, height = 6)
plot_lasso_cv(res_ml)
dev.off()

pdf(file.path(OUTPUT_DIR, "lasso_coef.pdf"), width = 8, height = 6)
plot_lasso_coef(res_ml)
dev.off()

pdf(file.path(OUTPUT_DIR, "rsf_process.pdf"), width = 8, height = 6)
plot_rsf_process(res_ml)
dev.off()

pdf(file.path(OUTPUT_DIR, "rsf_vimp.pdf"), width = 8, height = 6)
plot_rsf_vimp(res_ml)
dev.off()

pdf(file.path(OUTPUT_DIR, "xgb_importance.pdf"), width = 8, height = 6)
plot_xgb_imp(res_ml)
dev.off()

pdf(file.path(OUTPUT_DIR, "model_venn.pdf"), width = 8, height = 8)
plot_model_venn(res_ml, top_n = 5)
dev.off()

res_diag <- run_diagnostic_model(res$candidates)

saveRDS(res_diag, file.path(OUTPUT_DIR, "diag_result.rds"))

pdf(file.path(OUTPUT_DIR, "diag_heatmap.pdf"), width = 10, height = 8)
plot_diag_heatmap(res_diag, 10)
dev.off()

# Save candidate gene list
write.csv(
  data.frame(gene = res$candidates),
  file.path(OUTPUT_DIR, "candidate_genes.csv"),
  row.names = FALSE
)

# Save gene set info for step3
gene_set_info <- list(name = gene_info$name)
saveRDS(gene_set_info, file.path(OUTPUT_DIR, "gene_set_info.rds"))

message("")
message("==========================================================")
message(">>> Step 1 Complete!")
message(sprintf(">>> Results saved to: %s", OUTPUT_DIR))
message("")
message(sprintf(">>> Candidate gene count: %d", length(res$candidates)))
message("")
message(">>> Next steps:")
message("   1. Review plots in hcc_output/step1/")
message("   2. Select a target gene")
message("   3. Run: Rscript inst/scripts/step2_singlecell_analysis.R --gene GENE_NAME")
message("==========================================================")
