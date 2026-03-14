#!/usr/bin/env Rscript
# =============================================================================
# Step 3: Mechanism Analysis - Workflow 1 (Hepatocytes)
# Description: User selects type, automatically inherits gene and data from step2
# Usage:
#   Rscript inst/scripts/step3_mechanism.R --type 1
#   Rscript inst/scripts/step3_mechanism.R --type 2
# =============================================================================

if (!require("HCCMultiOmics", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages("devtools")
  }
  devtools::install_github("hepatology-will/HCCMultiOmics")
}
if (!require("HiClimR", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages("devtools")
  }
  devtools::install_github("andymarty/HiClimR")
}
library(HCCMultiOmics)
library(tidyverse)
library(patchwork)
library(Seurat)

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
  cat("========== Step 3: Mechanism Analysis ==========\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript inst/scripts/step3_mechanism.R --type 1    # Use mean expression\n")
  cat("  Rscript inst/scripts/step3_mechanism.R --type 2    # Use positive cell ratio\n")
  cat("\n")
  quit(save = "no")
}

parse_args <- function() {
  opts <- list(type = 1)

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--type", "-t")) {
      i <- i + 1
      if (i <= length(args)) opts$type <- as.integer(args[i])
    }
    i <- i + 1
  }

  if (opts$type != 1 && opts$type != 2) {
    stop("type must be 1 or 2")
  }
  opts
}

args <- parse_args()

OUTPUT_DIR <- "hcc_output/step3"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("========== Step 3: Mechanism Analysis ==========")
message("")

# ============== Inherit info from step2 ==============
target_gene_file <- "hcc_output/step2/target_gene.rds"
if (!file.exists(target_gene_file)) {
  candidates_file <- "hcc_output/step1/candidate_genes.csv"
  if (!file.exists(candidates_file)) {
    stop("Candidate gene file not found. Please run step1 and step2 first.")
  }
  candidates <- read.csv(candidates_file, stringsAsFactors = FALSE)$gene
  TARGET_GENE <- candidates[1]
  message(">>> Warning: Target gene from step2 not found, using first candidate gene")
} else {
  TARGET_GENE <- readRDS(target_gene_file)
}
message(sprintf(">>> Target gene: %s", TARGET_GENE))

# Read top_cells info from step2
top_cells_file <- "hcc_output/step2/top_cells.rds"
if (file.exists(top_cells_file)) {
  top_cells <- readRDS(top_cells_file)
  inherited_type <- args$type
  top_cell_type <- if (inherited_type == 1) top_cells$type1 else top_cells$type2
  message(sprintf(">>> Using type%d highest expression cell type: %s", inherited_type, top_cell_type))
} else {
  message(">>> Warning: top_cells.rds not found")
}

# Load gene set (contains driver, suppressor, syml)
gene_set_info_file <- "hcc_output/step1/gene_set_info.rds"
if (file.exists(gene_set_info_file)) {
  gene_set_info <- readRDS(gene_set_info_file)
  gene_set_name <- gene_set_info$name
} else {
  gene_set_name <- NULL
}

gene_set_file <- paste0(gene_set_name, ".rda")
data_dir <- file.path("..", "..", "data")
if (file.exists(file.path(data_dir, gene_set_file))) {
  load(file.path(data_dir, gene_set_file), envir = environment())
  if (exists(gene_set_name)) {
    gene_set_obj <- get(gene_set_name)
    driver_genes <- gene_set_obj$drivers
    suppressor_genes <- gene_set_obj$suppressors
    syml_genes <- gene_set_obj$syml
    all_genes <- gene_set_obj$all_genes
    message(sprintf(">>> Driver gene count: %d", length(driver_genes)))
    message(sprintf(">>> Suppressor gene count: %d", length(suppressor_genes)))
    message(sprintf(">>> Syml gene count: %d", length(syml_genes)))
  } else {
    driver_genes <- NULL
    suppressor_genes <- NULL
    syml_genes <- NULL
    all_genes <- NULL
  }
} else {
  driver_genes <- NULL
  suppressor_genes <- NULL
  syml_genes <- NULL
  all_genes <- NULL
}

# Load single-cell data
SC_DATA_PATH <- "HCC_sc_data.rda"
if (!file.exists(SC_DATA_PATH)) {
  alt_paths <- c(
    "inst/scripts/HCC_sc_data.rda",
    "hcc_output/step2/HCC_sc_data.rda",
    "data/HCC_sc_data.rda"
  )
  for (p in alt_paths) {
    if (file.exists(p)) {
      SC_DATA_PATH <- p
      break
    }
  }
}

if (!file.exists(SC_DATA_PATH)) {
  stop("Single-cell data file not found. Please run step2 first.")
}

if (grepl("\\.rda$", SC_DATA_PATH, ignore.case = TRUE)) {
  env <- new.env()
  load(SC_DATA_PATH, envir = env)
  HCC_sc_data <- env$HCC_sc_data
} else {
  HCC_sc_data <- readRDS(SC_DATA_PATH)
}

message(sprintf(">>> Single-cell data: %d cells", ncol(HCC_sc_data)))
message("")

# ============== Workflow 1: Hepatocyte Mechanism Analysis ==============
if (grepl("Hepatocyte|hepatic|liver", top_cell_type, ignore.case = TRUE)) {
  message(">>> Target gene primarily expressed in hepatocytes, starting hepatocyte mechanism analysis workflow (Workflow 1)...")

  # Step 1: Extract hepatocytes
  message(">>> Step 1: Extracting hepatocytes...")
  hep_cells <- WhichCells(HCC_sc_data, expression = SingleR_Labels == top_cell_type)
  hep_obj <- subset(HCC_sc_data, cells = hep_cells)
  message(sprintf(">>> Hepatocyte count: %d", ncol(hep_obj)))

  # Step 2: Visualization by site
  message(">>> Step 2: Visualizing by site...")
  GROUP_COL <- "site"

  if (GROUP_COL %in% colnames(hep_obj@meta.data)) {
    res_site <- plot_sc_gene_by_site(
      hep_obj,
      gene = TARGET_GENE,
      cell_type = top_cell_type,
      group_col = GROUP_COL,
      plot_type = args$type
    )

    pdf(file.path(OUTPUT_DIR, sprintf("%s_by_site.pdf", TARGET_GENE)), width = 12, height = 5)
    print(res_site$umap + res_site$barplot)
    dev.off()
  }

  # Step 3: Identify malignant cells
  message(">>> Step 3: Identifying malignant cells (CopyKAT)...")
  hep_obj <- identify_malignant_cells(hep_obj, cell_type = top_cell_type)

  pdf(file.path(OUTPUT_DIR, "malignant_stats.pdf"), width = 10, height = 6)
  p_stats <- plot_copykat_stats(hep_obj, group_col = "site")
  print(p_stats)
  dev.off()

  pdf(file.path(OUTPUT_DIR, "malignant_umap.pdf"), width = 10, height = 8)
  p_umap <- plot_copykat_umap(hep_obj, group_col = "site")
  print(p_umap)
  dev.off()

  # Step 4: Extract malignant hepatocytes
  message(">>> Step 4: Extracting malignant hepatocytes...")
  malignant_cells <- WhichCells(hep_obj, expression = cnv_status == "aneuploid")
  hep_malignant <- subset(hep_obj, cells = malignant_cells)
  message(sprintf(">>> Malignant hepatocyte count: %d", ncol(hep_malignant)))

  # Step 5: Re-clustering
  message(">>> Step 5: Re-clustering malignant hepatocytes...")
  hep_malignant <- reprocess_malignant_cells(hep_malignant, batch_col = "patient", theta = 3, resolution = 0.5)

  pdf(file.path(OUTPUT_DIR, "malignant_clusters.pdf"), width = 10, height = 8)
  p_clusters <- plot_malignant_clusters(hep_malignant)
  print(p_clusters)
  dev.off()

  # Step 6: CytoTRACE stemness analysis
  message(">>> Step 6: CytoTRACE stemness analysis...")
  res_stemness <- run_stemness_analysis(hep_malignant)
  hep_malignant <- res_stemness$obj

  pdf(file.path(OUTPUT_DIR, "cyto_score.pdf"), width = 10, height = 8)
  print(res_stemness$plot_umap)
  dev.off()

  pdf(file.path(OUTPUT_DIR, "cyto_violin.pdf"), width = 10, height = 6)
  print(res_stemness$plot_vln)
  dev.off()

  # Step 7: Select cluster with minimum stemness score as starting point
  message(">>> Step 7: Selecting cluster with minimum stemness score as starting point...")
  cluster_cyto <- aggregate(hep_malignant$cyto_score, by = list(hep_malignant$seurat_clusters), FUN = mean)
  start_cluster <- cluster_cyto$Group.1[which.min(cluster_cyto$x)]
  message(sprintf(">>> Cluster with lowest stemness score: %s", start_cluster))

  # Step 8: Slingshot pseudotime analysis
  message(">>> Step 8: Slingshot pseudotime analysis...")
  res_trajectory <- run_trajectory_analysis(hep_malignant, start_cluster = start_cluster)

  pdf(file.path(OUTPUT_DIR, "trajectory.pdf"), width = 12, height = 10)
  print(res_trajectory$plot_all_curves)
  dev.off()

  # Step 9: Extract main trajectory cells
  message(">>> Step 9: Extracting main trajectory cells...")
  res_subset <- subset_trajectory(hep_malignant, res_trajectory$sds, curve_id = 1)
  hep_trajectory <- res_subset$obj
  message(sprintf(">>> Main trajectory cell count: %d", ncol(hep_trajectory)))

  pdf(file.path(OUTPUT_DIR, "pseudotime.pdf"), width = 10, height = 8)
  print(res_subset$plot_final)
  dev.off()

  # Step 10-11: Calculate Driver/Suppressor scores and visualize correlation with target gene
  message(">>> Step 10-11: Driver/Suppressor scores vs target gene correlation...")

  # Calculate Driver score
  if (!is.null(driver_genes) && length(driver_genes) > 0) {
    hep_trajectory <- calc_ferro_score(hep_trajectory, gene_list = driver_genes, name = "DriverScore")

    pdf(file.path(OUTPUT_DIR, "driver_cor.pdf"), width = 8, height = 6)
    p_driver_cor <- plot_mechanism_cor(hep_trajectory, x_feature = TARGET_GENE, y_feature = "DriverScore")
    print(p_driver_cor)
    dev.off()

    # Plot Driver trend
    pdf(file.path(OUTPUT_DIR, "driver_trend.pdf"), width = 10, height = 6)
    p_driver_trend <- plot_pseudotime_trend_bins(
      seurat_obj = hep_trajectory,
      feature = "DriverScore",
      pseudotime_col = "pseudotime",
      color = "#3A8E6A",
      title = "Ferroptosis Sensitivity"
    )
    print(p_driver_trend)
    dev.off()
  }

  # Calculate Suppressor score
  if (!is.null(suppressor_genes) && length(suppressor_genes) > 0) {
    hep_trajectory <- calc_ferro_score(hep_trajectory, gene_list = suppressor_genes, name = "SuppressorScore")

    pdf(file.path(OUTPUT_DIR, "suppressor_cor.pdf"), width = 8, height = 6)
    p_suppressor_cor <- plot_mechanism_cor(hep_trajectory, x_feature = TARGET_GENE, y_feature = "SuppressorScore")
    print(p_suppressor_cor)
    dev.off()

    # Plot Suppressor trend
    pdf(file.path(OUTPUT_DIR, "suppressor_trend.pdf"), width = 10, height = 6)
    p_suppressor_trend <- plot_pseudotime_trend_bins(
      seurat_obj = hep_trajectory,
      feature = "SuppressorScore",
      pseudotime_col = "pseudotime",
      color = "#C98648",
      title = "Ferroptosis Suppression"
    )
    print(p_suppressor_trend)
    dev.off()
  }

  # Step 12: Screen genes with highest correlation to syml genes
  message(">>> Step 12: Screening genes with highest correlation to syml genes...")

  if (!is.null(syml_genes) && length(syml_genes) > 0) {
    screen_res <- screen_ferro_correlations(
      seurat_obj = hep_trajectory,
      target_gene = TARGET_GENE,
      db_genes = syml_genes,
      top_n = 25
    )

    pdf(file.path(OUTPUT_DIR, "syml_screening.pdf"), width = 10, height = 8)
    print(screen_res$plot)
    dev.off()

    top_syml_gene <- screen_res$top_gene
    message(sprintf(">>> Gene with highest syml correlation: %s", top_syml_gene))

    # Step 13: Plot density
    message(">>> Step 13: Plotting density...")
    pdf(file.path(OUTPUT_DIR, sprintf("%s_%s_density.pdf", TARGET_GENE, top_syml_gene)), width = 8, height = 6)
    p_density <- plot_nature_density(hep_trajectory, x_gene = TARGET_GENE, y_gene = top_syml_gene)
    print(p_density)
    dev.off()

    # Step 14: tradeSeq analysis
    message(">>> Step 14: tradeSeq analysis...")

    if (!is.null(all_genes) && length(all_genes) > 0) {
      valid_all <- intersect(all_genes, rownames(hep_trajectory))
      tradeseq_genes <- c(TARGET_GENE, valid_all)
      message(sprintf(">>> Using %d genes for tradeSeq analysis", length(tradeseq_genes)))

      res_tradeseq <- run_ferro_tradeseq(
        seurat_obj = hep_trajectory,
        gene_list = tradeseq_genes,
        pseudotime_col = "pseudotime",
        n_cores = 1
      )
      sce_ts <- res_tradeseq$sce

      # Heatmap
      tryCatch({
        pdf(file.path(OUTPUT_DIR, "tradeseq_heatmap.pdf"), width = 10, height = 8)
        p_heat <- plot_tradeseq_heatmap(sce_obj = sce_ts, sig_genes = valid_all, key_genes = c(TARGET_GENE, top_syml_gene))
        print(p_heat)
        dev.off()
      }, error = function(e) {
        message(">>> tradeSeq heatmap skipped: ", e$message)
      })

      # Trends
      pdf(file.path(OUTPUT_DIR, "tradeseq_trends.pdf"), width = 12, height = 8)
      p_trends <- plot_tradeseq_trends(sce_obj = sce_ts, plot_genes = c(TARGET_GENE, top_syml_gene))
      print(p_trends)
      dev.off()
    }
  }

  # Step 15: Target gene trend curve
  message(">>> Step 15: Target gene trend curve...")
  pdf(file.path(OUTPUT_DIR, sprintf("%s_pseudotime_trend.pdf", TARGET_GENE)), width = 10, height = 6)
  p_final_trend <- plot_pseudotime_trend_bins(
    seurat_obj = hep_trajectory,
    feature = TARGET_GENE,
    pseudotime_col = "pseudotime",
    color = "#C98648",
    title = paste0("Epigenetic Activation (", TARGET_GENE, ")")
  )
  print(p_final_trend)
  dev.off()

  # Step 16: Export SCENIC input format
  # Use ALL malignant hepatocytes (not just trajectory cells) for SCENIC analysis
  message(">>> Step 16: Exporting SCENIC input format...")
  message(">>> Using ALL malignant hepatocytes (n = ", ncol(hep_malignant), ")")
  export_scenic_data(
    hep_malignant,
    output_dir = OUTPUT_DIR,
    file_name = "scenic_input.csv",
    target_gene = TARGET_GENE
  )
  
  # Save malignant hepatocytes for downstream analysis (step5)
  message(">>> Saving malignant hepatocytes for downstream analysis...")
  hep_mal_reprocessed <- hep_malignant
  save(hep_mal_reprocessed, file = file.path(OUTPUT_DIR, "hep_mal_reprocessed.rda"))
  message(">>> Saved to: ", file.path(OUTPUT_DIR, "hep_mal_reprocessed.rda"))

} else {
  # ============== Workflow 2: Microenvironment Analysis ==============
  message(sprintf(">>> Target gene primarily expressed in %s, starting microenvironment analysis workflow (Workflow 2)...", top_cell_type))

  # Step 1: Ligand-receptor interaction analysis
  message(">>> Step 2-1: Ligand-receptor interaction analysis...")
  res_ligand <- auto_analyze_ligand(
    HCC_sc_data,
    gene = TARGET_GENE,
    group_col = "SingleR_Labels",
    db_use = "human"
  )

  pdf(file.path(OUTPUT_DIR, "ligand_bubble.pdf"), width = 10, height = 8)
  print(res_ligand$bubble)
  dev.off()

  # Step 2: Stromal cell function analysis
  message(">>> Step 2-2: Stromal cell function analysis...")
  res_stromal <- auto_stromal_function(
    HCC_sc_data,
    gene = TARGET_GENE,
    cell_type = top_cell_type,
    group_col = "SingleR_Labels",
    species = "Homo sapiens",
    top_n = 10
  )

  pdf(file.path(OUTPUT_DIR, "stromal_function.pdf"), width = 10, height = 8)
  print(res_stromal$plot)
  dev.off()

  # Step 3: Immune microenvironment analysis
  message(">>> Step 2-3: Immune microenvironment analysis...")
  res_immune <- auto_immune_landscape(
    seurat_obj = HCC_sc_data,
    gene = TARGET_GENE,
    stromal_class = top_cell_type,
    group_col = "SingleR_Labels",
    sample_col = "orig.ident",
    plot_type = "all"
  )

  pdf(file.path(OUTPUT_DIR, "immune_landscape.pdf"), width = 12, height = 10)
  print(res_immune$landscape)
  dev.off()

  pdf(file.path(OUTPUT_DIR, "immune_scatter_grid.pdf"), width = 14, height = 10)
  print(res_immune$scatter_grid)
  dev.off()
}

message("")
message("==========================================================")
message(">>> Step 3 Complete!")
message(sprintf(">>> Results saved to: %s", OUTPUT_DIR))
message("")
message("==========================================================")
