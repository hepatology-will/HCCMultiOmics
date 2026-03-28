#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Spatial Transcriptomics Validation
# Description: Validate TF-target gene regulatory relationships using spatial transcriptomics
# Usage:
#   Rscript step6_spatial_validation.R --target EZH2 --tfs TF1,TF2
# =============================================================================

if (!require("HCCMultiOmics", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages("devtools")
  }
  devtools::install_github("hepatology-will/HCCMultiOmics")
}
library(HCCMultiOmics)

# ============== Argument Parsing ==============
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("========== Step 6: Spatial Transcriptomics Validation ==========\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript step6_spatial_validation.R --target TARGET_GENE --tfs TF1,TF2\n")
  cat("  Rscript step6_spatial_validation.R --target EZH2 --tfs YY1,EP300\n")
  cat("\n")
  cat("Parameters:\n")
  cat("  --target   Target gene from Step 2 (required)\n")
  cat("  --tfs      Comma-separated TF names from Step 5 (required)\n")
  cat("  --data     Spatial data: auto (download) or local path (optional)\n")
  cat("\n")
  cat("Note:\n")
  cat("  - Driver/Suppressor genes are loaded from Step 1 gene set\n")
  cat("  - Spatial data will be downloaded from Zenodo if not available locally\n")
  cat("  - Place downloaded data in: inst/scripts/sap_data/st_data.rda\n")
  cat("\n")
  quit(save = "no")
}

parse_args <- function() {
  opts <- list(target = NULL, tfs = NULL, data = "local")

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--target", "-t")) {
      i <- i + 1
      if (i <= length(args)) opts$target <- args[i]
    } else if (args[i] %in% c("--tfs", "-f")) {
      i <- i + 1
      if (i <= length(args)) opts$tfs <- strsplit(args[i], ",")[[1]]
    } else if (args[i] %in% c("--data", "-d")) {
      i <- i + 1
      if (i <= length(args)) opts$data <- args[i]
    }
    i <- i + 1
  }

  if (is.null(opts$target)) {
    stop("Please specify --target. Run --help for usage.")
  }
  if (is.null(opts$tfs)) {
    stop("Please specify --tfs. Run --help for usage.")
  }
  opts
}

args <- parse_args()

OUTPUT_DIR <- "hcc_output/step6"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAP_DATA_DIR <- "sap_data"
ZENODO_ID <- "19070970"

message("========== Step 6: Spatial Transcriptomics Validation ==========")
message("")

# ============== Load Target Gene ==============
TARGET_GENE <- args$target
message(sprintf(">>> Target gene: %s", TARGET_GENE))

# ============== Load TFs from Step 5 ==============
TF_LIST <- args$tfs
message(sprintf(">>> TFs to validate: %s", paste(TF_LIST, collapse = ", ")))

# ============== Load Gene Set (Driver/Suppressor) ==============
gene_set_info_file <- "hcc_output/step1/gene_set_info.rds"
if (file.exists(gene_set_info_file)) {
  gene_set_info <- readRDS(gene_set_info_file)
  gene_set_name <- gene_set_info$name
} else {
  gene_set_name <- NULL
}

driver_genes <- NULL
suppressor_genes <- NULL

if (!is.null(gene_set_name)) {
  gene_set_file <- paste0(gene_set_name, ".rda")
  data_dir <- file.path("..", "..", "data")

  if (file.exists(file.path(data_dir, gene_set_file))) {
    load(file.path(data_dir, gene_set_file), envir = environment())
    if (exists(gene_set_name)) {
      gene_set_obj <- get(gene_set_name)
      driver_genes <- gene_set_obj$drivers
      suppressor_genes <- gene_set_obj$suppressors
      message(sprintf(">>> Driver genes: %d", length(driver_genes)))
      message(sprintf(">>> Suppressor genes: %d", length(suppressor_genes)))
    }
  }
}

# ============== Load Spatial Data ==============
message(">>> Loading spatial transcriptomics data...")

ST_DATA_PATH <- file.path(SAP_DATA_DIR, "st_data.rda")

download_st_data <- function() {
  zenodo_url <- sprintf("https://zenodo.org/records/%s/files/st_data.rda", ZENODO_ID)
  message(sprintf(">>> Downloading from Zenodo (Record: %s)", ZENODO_ID))
  dir.create(SAP_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    download.file(zenodo_url, ST_DATA_PATH, method = "wget", quiet = FALSE)
    message(">>> Download complete!")
  }, error = function(e) {
    stop(sprintf("Download failed: %s", e$message))
  })
}

if (!file.exists(ST_DATA_PATH)) {
  if (args$data == "auto") {
    download_st_data()
  } else {
    alt_paths <- c(
      "inst/scripts/sap_data/st_data.rda",
      "sap_data/st_data.rda",
      "data/st_data.rda"
    )
    for (p in alt_paths) {
      if (file.exists(p)) {
        ST_DATA_PATH <- p
        break
      }
    }
  }
}

if (!file.exists(ST_DATA_PATH)) {
  stop("Spatial data file st_data.rda not found. Use --data auto to download or place file in inst/scripts/sap_data/")
}

message(sprintf(">>> Using spatial data: %s", ST_DATA_PATH))

if (grepl("\\.rda$", ST_DATA_PATH, ignore.case = TRUE)) {
  env <- new.env()
  load(ST_DATA_PATH, envir = env)
  st_data <- env$st_data
} else {
  st_data <- readRDS(ST_DATA_PATH)
}

message(sprintf(">>> Spatial data loaded: %d spots", ncol(st_data)))
message("")

# ============== Run Spatial Validation ==============
message(">>> Running spatial TF validation analysis...")

res <- run_spatial_validation(
  seurat_obj = st_data,
  target_gene = TARGET_GENE,
  tf_list = TF_LIST,
  driver_genes = driver_genes,
  suppressor_genes = suppressor_genes
)

# ============== Generate Plots ==============
message(">>> Generating spatial distribution plots...")

plot_spatial_distribution(
  plot_data = res$plot_data,
  output_dir = OUTPUT_DIR,
  target_gene = TARGET_GENE,
  tf_name = TF_LIST[1]
)

message(">>> Generating correlation hexbin plots...")

plot_correlation_hexbin(
  plot_data = res$plot_data,
  output_dir = OUTPUT_DIR
)

# ============== Save Results ==============
saveRDS(res, file.path(OUTPUT_DIR, "spatial_validation_result.rds"))

message("")
message("==========================================================")
message(">>> Step 6 Complete!")
message(sprintf(">>> Results saved to: %s", OUTPUT_DIR))
message("")
message("Output files:")
message("  - spatial_distribution.pdf : Spatial distribution of target, TF, ferro scores")
message("  - correlation_hexbin.pdf   : Hexbin correlation plots")
message("  - spatial_validation_result.rds : Full results")
message("")
message(">>> Workflow complete!")
message("==========================================================")
