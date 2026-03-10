#' @importFrom ggplot2 ggplot aes stat_density_2d scale_fill_gradient geom_point geom_density_2d geom_smooth annotate labs theme_classic theme element_blank element_line element_text after_stat geom_line scale_x_continuous
#' @importFrom stats cor.test sd
NULL

# ==============================================================================
# Part 1: Calculate Functional Scores (AddModuleScore Wrapper)
# ==============================================================================

#' Calculate Ferroptosis (or Custom) Module Scores
#'
#' Calculates module scores for a given gene set (e.g., Ferroptosis Drivers).
#'
#' @param seurat_obj A Seurat object.
#' @param gene_list A vector of gene symbols (e.g., ferroptosis drivers).
#' @param name The name for the metadata column (e.g., "FerroScore").
#' @return A Seurat object with the calculated score added to metadata.
#' @export
calc_ferro_score <- function(seurat_obj, gene_list, name = "FerroScore") {

  # 1. Check if genes exist in matrix
  valid_genes <- intersect(gene_list, rownames(seurat_obj))

  if(length(valid_genes) == 0) {
    stop("Error: None of the genes in the list are found in the Seurat object.")
  }

  message(paste(">>> Calculating score using", length(valid_genes), "valid genes..."))

  # 2. Run AddModuleScore
  seurat_obj <- Seurat::AddModuleScore(
    object = seurat_obj,
    features = list(valid_genes),
    name = name
  )

  # 3. Fix column name (remove '1' that Seurat automatically adds)
  raw_col_name <- paste0(name, "1")
  if (raw_col_name %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[name]] <- seurat_obj@meta.data[[raw_col_name]]
    seurat_obj@meta.data[[raw_col_name]] <- NULL
  }

  message(paste(">>> Score calculated. Added metadata column:", name))
  return(seurat_obj)
}

# ==============================================================================
# Part 2: Correlation Analysis Plotting (Nature-style density scatter plot)
# ==============================================================================

#' Plot Mechanism Correlation (Density + Regression)
#'
#' Visualizes the correlation between a Target Gene (e.g., EZH2) and a Functional Score
#' (e.g., Ferroptosis) with a high-quality "Nature-style" density plot.
#'
#' @param seurat_obj A Seurat object with necessary metadata/expression.
#' @param x_feature Name of the gene/metadata for X axis (e.g., "EZH2").
#' @param y_feature Name of the gene/metadata for Y axis (e.g., "FerroScore").
#' @param filter_zeros Boolean. Whether to filter out cells with 0 expression in X (default: TRUE).
#' @return A ggplot object.
#' @export
plot_mechanism_cor <- function(seurat_obj, x_feature = "EZH2", y_feature = "FerroScore", filter_zeros = TRUE) {

  # 1. Data extraction
  df_cor <- Seurat::FetchData(seurat_obj, vars = c(x_feature, y_feature))
  colnames(df_cor) <- c("X_Val", "Y_Val")

  # 2. Filter zero values
  if (filter_zeros) {
    df_cor <- df_cor %>% dplyr::filter(X_Val > 0)
  }

  # 3. Calculate correlation (Spearman)
  stat_res <- stats::cor.test(df_cor$X_Val, df_cor$Y_Val, method = "spearman")

  rho_val <- round(unname(stat_res$estimate), 2)
  p_val_str <- ifelse(stat_res$p.value < 0.001, " < 0.001", paste0(" = ", signif(stat_res$p.value, 2)))
  label_plotmath <- paste0("Spearman~rho == ", rho_val, "~~~italic(P)", p_val_str)

  # 4. Define theme and color scheme
  line_color <- "#8EADD9"
  fill_color <- "#AFCCD9"
  density_low <- "#F9F9F9"
  density_high <- "#D4E4FC"

  theme_nature_elegant <- ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 1.0, color = "#333333"),
      axis.text = ggplot2::element_text(size = 14, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  # 5. Plot
  p <- ggplot2::ggplot(df_cor, ggplot2::aes(x = X_Val, y = Y_Val)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)), geom = "polygon", alpha = 0.4) +
    ggplot2::scale_fill_gradient(low = density_low, high = density_high) +
    ggplot2::geom_point(color = "grey60", size = 1.2, alpha = 0.2, stroke = 0) +
    ggplot2::geom_density_2d(color = "grey50", linewidth = 0.4, bins = 5) +
    ggplot2::geom_smooth(method = "lm", color = line_color, fill = fill_color, linewidth = 1.8) +
    ggplot2::annotate("text", x = max(df_cor$X_Val), y = max(df_cor$Y_Val),
                      label = label_plotmath, parse = TRUE,
                      size = 5.5, fontface = "bold", hjust = 1, vjust = 1) +
    ggplot2::labs(x = paste(x_feature, "Expression"), y = paste(y_feature, "Score")) +
    theme_nature_elegant

  return(p)
}
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr ntile group_by summarise mutate
NULL

#' @importFrom ggplot2 ggplot aes geom_line geom_point annotate scale_x_continuous labs theme_classic theme element_text element_line element_blank
#' @importFrom dplyr mutate ntile group_by summarise
#' @importFrom stats cor.test
NULL

#' Plot Trend along Pseudotime Bins (Nature Style)
#'
#' Modified to plot a SINGLE feature with Stage I-IV x-axis and correlation stats.
#'
#' @param seurat_obj A Seurat object.
#' @param feature The feature to plot (Gene name or Metadata column).
#' @param pseudotime_col The metadata column name for pseudotime.
#' @param bins Number of bins (default: 4 for Stage I-IV).
#' @param color The color for the line and points (hex code).
#' @param title Custom title for the plot.
#' @return A ggplot object.
#' @export
plot_pseudotime_trend_bins <- function(seurat_obj, feature, pseudotime_col,
                                       bins = 4,
                                       color = "#333333",
                                       title = NULL) {

  # 1. Check if feature exists (compatible with Gene and Metadata)
  if (feature %in% rownames(seurat_obj)) {
    type <- "gene"
  } else if (feature %in% colnames(seurat_obj@meta.data)) {
    type <- "meta"
  } else {
    stop(paste("Feature", feature, "not found in RNA assay or Metadata."))
  }

  # 2. Extract data
  df <- Seurat::FetchData(seurat_obj, vars = c(pseudotime_col, feature))
  colnames(df) <- c("Time", "Value")

  # Filter invalid values
  df <- df[df$Value > 0 & !is.na(df$Time), ]

  # 3. Calculate correlation (Spearman) for legend
  cor_res <- stats::cor.test(df$Time, df$Value, method = "spearman")
  rho_val <- round(cor_res$estimate, 2)
  p_val <- cor_res$p.value
  p_text <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", round(p_val, 3)))
  stats_label <- paste0("Rho = ", rho_val, ", ", p_text)

  # 4. Bin data (Calculating Mean per Bin)
  df_summary <- df %>%
    dplyr::mutate(Bin = dplyr::ntile(Time, bins)) %>%
    dplyr::group_by(Bin) %>%
    dplyr::summarise(Mean_Val = mean(Value), .groups = "drop")

  # 5. Generate X-axis labels (Stage I, II, III, IV...)
  roman_labels <- as.character(utils::as.roman(1:bins))
  stage_labels <- paste("Stage", roman_labels)

  # 6. Plot (Strictly replicate reference style)
  p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = Bin, y = Mean_Val)) +
    # Line
    ggplot2::geom_line(color = color, linewidth = 1.2) +
    # Hollow circles: fill="white" (white background), color=color (colored border), stroke=2 (thick border)
    ggplot2::geom_point(size = 4.5, shape = 21, fill = "white", color = color, stroke = 2) +
    # Stats label (adaptive position)
    ggplot2::annotate("text", x = 1.2, y = max(df_summary$Mean_Val),
                      label = stats_label,
                      color = color, size = 5, fontface = "bold", hjust = 0, vjust = -0.5) +
    # Axis adjustment
    ggplot2::scale_x_continuous(breaks = 1:bins, labels = stage_labels) +
    ggplot2::labs(
      x = "Pseudotime Progression",
      y = ifelse(type == "gene", paste(feature, "Expression"), paste(feature, "Score")),
      title = ifelse(is.null(title), feature, title)
    ) +
    # Theme adjustment
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      axis.line = ggplot2::element_line(linewidth = 0.8),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}

#' Batch Screen Ferroptosis Gene Correlations (Corrected: Properly extract Top 1)
#'
#' @param seurat_obj Seurat object
#' @param target_gene Target gene (default "EZH2")
#' @param db_genes Ferroptosis gene list
#' @param top_n How many to display (default 25)
#' @return List: plot (ggplot), results (all results), top_gene (Top 1 gene name)
#' @export
screen_ferro_correlations <- function(seurat_obj, target_gene = "EZH2", db_genes, top_n = 25) {

  # 1. Preprocessing
  valid_genes <- intersect(unique(db_genes), rownames(seurat_obj))
  valid_genes <- setdiff(valid_genes, target_gene)

  if(length(valid_genes) == 0) stop("No valid genes found.")
  message(paste(">>> Screening correlations for", length(valid_genes), "genes..."))

  # 2. Batch calculation
  res_list <- list()
  vec_target <- Seurat::FetchData(seurat_obj, vars = target_gene)[, 1]

  # Batch extract data
  mat_genes <- Seurat::FetchData(seurat_obj, vars = valid_genes)

  for (g in valid_genes) {
    vec_g <- mat_genes[[g]]
    # Double non-zero filtering
    mask <- (vec_target > 0) & (vec_g > 0)

    if (sum(mask) > 30) {
      stat <- stats::cor.test(vec_target[mask], vec_g[mask], method = "spearman")
      res_list[[g]] <- data.frame(
        Gene = g,
        Rho = stat$estimate,
        Pval = stat$p.value,
        Correlation_Type = ifelse(stat$estimate > 0, "Positive", "Negative")
      )
    }
  }

  if(length(res_list) == 0) stop("No significant correlations found.")
  df_res <- do.call(rbind, res_list)
  rownames(df_res) <- NULL

  # 3. Sort and extract Top 1 (Key fix)
  # Sort by absolute value descending -> first row is the strongest
  df_plot <- df_res %>%
    dplyr::arrange(dplyr::desc(abs(Rho))) %>%
    utils::head(top_n)

  # [Fix point] Take the first row directly, this is the one with highest correlation!
  best_gene <- as.character(df_plot$Gene[1])

  # Set factor levels: rev() is to make the top rank appear at the top of the plot after coord_flip
  df_plot <- df_plot %>%
    dplyr::mutate(Gene = factor(Gene, levels = rev(Gene)))

  # P-value labels
  df_plot$P_label <- ifelse(df_plot$Pval < 0.001, "***",
                            ifelse(df_plot$Pval < 0.01, "**", "*"))

  # 4. Plot
  colors_type <- c("Positive" = "#BC3C29", "Negative" = "#0072B5")

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Gene, y = Rho)) +
    ggplot2::geom_segment(ggplot2::aes(x = Gene, xend = Gene, y = 0, yend = Rho, color = Correlation_Type), linewidth = 1.2) +
    # Highlight Top 1
    ggplot2::geom_point(ggplot2::aes(color = Correlation_Type, size = (Gene == best_gene)), show.legend = FALSE) +
    ggplot2::scale_size_manual(values = c("FALSE" = 4, "TRUE" = 6)) +
    ggplot2::geom_text(ggplot2::aes(label = P_label, y = Rho + 0.02 * sign(Rho)), size = 5, vjust = 0.7) +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = colors_type) +
    ggplot2::labs(y = paste("Spearman Correlation with", target_gene), x = NULL,
                  title = paste("Top Candidates Correlated with", target_gene)) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linetype = "dashed"),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )

  return(list(plot = p, results = df_res, top_gene = best_gene))
}
# ==============================================================================
# Function 2: Nature-style Density Plot (Validation Function)
# ==============================================================================

#' Plot Nature-style Density Scatter Plot (Validate Correlation)
#'
#' @param seurat_obj Seurat object
#' @param x_gene X-axis gene (usually EZH2)
#' @param y_gene Y-axis gene (usually the screened Top 1)
#' @return ggplot object
#' @export
plot_nature_density <- function(seurat_obj, x_gene, y_gene) {

  # 1. Extract data
  df <- Seurat::FetchData(seurat_obj, vars = c(x_gene, y_gene))
  colnames(df) <- c("X_Val", "Y_Val")

  # 2. Double filter (only show cells expressing both)
  df_nonzero <- df[df$X_Val > 0 & df$Y_Val > 0, ]

  # 3. Calculate correlation annotation
  stat <- stats::cor.test(df_nonzero$X_Val, df_nonzero$Y_Val, method = "spearman")
  rho_val <- round(stat$estimate, 2)
  p_val_str <- ifelse(stat$p.value < 0.001, " < 0.001", paste0(" = ", signif(stat$p.value, 2)))
  label_plotmath <- paste0("Spearman~rho == ", rho_val, "~~~italic(P)", p_val_str)

  # 4. Define strict color scheme (replicate your requirement)
  line_color  <- "#0072B5"   # Blue line
  fill_color  <- "#AFCCD9"   # Light blue fill
  point_color <- "#888888"   # Gray scatter points
  density_low <- "#FFFFFF"
  density_high <- "#E5E7E9"  # Silver gray background

  # 5. Plot
  p <- ggplot2::ggplot(df_nonzero, ggplot2::aes(x = X_Val, y = Y_Val)) +
    # Density layer
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)), geom = "polygon", alpha = 0.5) +
    ggplot2::scale_fill_gradient(low = density_low, high = density_high) +
    # Scatter layer
    ggplot2::geom_point(color = point_color, size = 1.8, alpha = 0.4, stroke = 0) +
    # Fit layer
    ggplot2::geom_smooth(method = "lm", color = line_color, fill = fill_color, linewidth = 2.0) +
    # Annotation layer
    ggplot2::annotate("text", x = max(df_nonzero$X_Val), y = max(df_nonzero$Y_Val),
                      label = label_plotmath, parse = TRUE,
                      size = 6, fontface = "bold", hjust = 1, vjust = 1) +
    # Label layer
    ggplot2::labs(x = paste(x_gene, "Expression"),
                  y = bquote(bold(.(y_gene) ~ "Expression"))) +
    # Theme layer
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 1.0, color = "#333333"),
      axis.text = ggplot2::element_text(size = 14, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  return(p)
}



#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by ungroup arrange desc
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual theme_classic labs theme element_blank element_text element_line
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom stats p.adjust

NULL

# ==============================================================================
# Part 1: tradeSeq Model Fitting (Memory Optimization + Fix set.seed Error)
# ==============================================================================

#' Run tradeSeq Analysis for Ferroptosis Genes
#'
#' Converts Seurat object to SingleCellExperiment (subsetted to specific genes to save RAM),
#' fits GAM models using tradeSeq.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_list A vector of gene symbols (e.g., Ferroptosis genes).
#' @param pseudotime_col The column name in metadata representing pseudotime.
#' @param nknots Number of knots for GAM fitting (default: 6).
#' @param n_cores Number of cores for parallel processing (default: 1).
#' @return A list containing:
#' \item{sce}{The fitted SingleCellExperiment object.}
#' \item{asso_res}{A data frame of association test results, sorted by p-value.}
#' @export
run_ferro_tradeseq <- function(seurat_obj, gene_list, pseudotime_col = "pseudotime", nknots = 6, n_cores = 1) {

  .check_bioc_packages(c("tradeSeq", "SingleCellExperiment", "SummarizedExperiment"))
  .check_cran_packages("BiocParallel")

  # 1. Check pseudotime column
  if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", pseudotime_col, "not found in metadata."))
  }

  # 2. Gene filtering: only keep genes present in Seurat
  # This step is critical: filter before extracting large matrix to avoid 2.8GB memory overflow
  valid_genes <- intersect(unique(gene_list), rownames(seurat_obj))
  if (length(valid_genes) < 2) stop("Too few valid genes found in the matrix.")

  message(paste(">>> Extracting data for", length(valid_genes), "genes (Memory Optimized)..."))

  # 3. Only extract target genes sparse matrix and convert to dense matrix
  # Seurat V5 compatible syntax
  counts_sub <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[valid_genes, ]
  counts_sub <- as.matrix(counts_sub)

  # 4. Build SCE object
  sce_ts <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_sub)
  )

  # 5. Add metadata (use SummarizedExperiment to avoid namespace errors)
  # Ensure cell order is consistent
  meta_data <- seurat_obj@meta.data[colnames(counts_sub), ]

  SummarizedExperiment::colData(sce_ts)$pseudotime <- meta_data[[pseudotime_col]]
  SummarizedExperiment::colData(sce_ts)$cluster <- meta_data$seurat_clusters

  # 6. Filter out cells without pseudotime
  valid_cells <- !is.na(SummarizedExperiment::colData(sce_ts)$pseudotime)
  sce_ts <- sce_ts[, valid_cells]

  message(paste(">>> Analyzing", length(valid_genes), "genes across", ncol(sce_ts), "cells..."))

  # 7. Build tradeSeq weight matrix
  pt <- SummarizedExperiment::colData(sce_ts)$pseudotime
  ptime_mat <- matrix(pt, ncol = 1)
  colnames(ptime_mat) <- "lin1"

  cell_w <- matrix(1, nrow = ncol(sce_ts), ncol = 1)
  colnames(cell_w) <- "lin1"

  # 8. Run fitGAM
  # set.seed is from base package, call directly
  set.seed(123)
  message(paste(">>> Fitting GAM with nknots =", nknots, "using", n_cores, "cores..."))

  sce_ts <- tryCatch({
    tradeSeq::fitGAM(
      counts      = SingleCellExperiment::counts(sce_ts),
      pseudotime  = ptime_mat,
      cellWeights = cell_w,
      nknots      = nknots,
      verbose     = TRUE,
      parallel    = (n_cores > 1),
      BPPARAM     = if(n_cores > 1) BiocParallel::MulticoreParam(workers = n_cores) else NULL
    )
  }, error = function(e) {
    stop("fitGAM failed: ", e$message)
  })

  # 9. Association test
  message(">>> Running associationTest...")
  asso_res <- tradeSeq::associationTest(sce_ts)
  asso_res$gene <- rownames(asso_res)
  asso_res <- asso_res[order(asso_res$pvalue), ]

  return(list(sce = sce_ts, asso_res = asso_res))
}

# ==============================================================================
# Part 2: Heatmap Visualization (No changes needed, ensure Imports exist)
# ==============================================================================

#' Plot tradeSeq Heatmap
#' @param sce_obj The fitted SingleCellExperiment object.
#' @param sig_genes Significant genes to plot.
#' @param key_genes Specific genes to label (e.g., EZH2).
#' @param n_points Smoothing points (default: 100).
#' @export
plot_tradeseq_heatmap <- function(sce_obj, sig_genes, key_genes = c("EZH2", "SLC7A11"), n_points = 100) {

  valid_sig <- intersect(sig_genes, rownames(sce_obj))
  if(length(valid_sig) == 0) stop("No valid significant genes found.")

  message(">>> Clustering expression patterns...")

  # clusterExpressionPatterns uses SummarizedExperiment here
  clus_res <- tradeSeq::clusterExpressionPatterns(
    sce_obj, nPoints = n_points, genes = valid_sig,
    reduceMethod = "PCA", verbose = FALSE
  )

  heatmap_mat <- clus_res$yhatScaled

  # Handle row labels: only show key_genes
  gene_labels <- ifelse(rownames(heatmap_mat) %in% key_genes, rownames(heatmap_mat), "")

  pheatmap::pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    labels_row = gene_labels,
    color = grDevices::colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
    fontsize_row = 10,
    border_color = NA,
    main = "Gene Expression Dynamics along Pseudotime"
  )
}

# ==============================================================================
# Part 3: Trend Plot (No changes)
# ==============================================================================

#' Plot Smooth Trends
#' @param sce_obj The fitted SCE object.
#' @param plot_genes Genes to plot.
#' @param colors Named vector of colors.
#' @export
plot_tradeseq_trends <- function(sce_obj, plot_genes, colors = NULL) {

  valid_genes <- intersect(plot_genes, rownames(sce_obj))
  if(length(valid_genes) == 0) stop("Genes not found in the fitted model.")

  message(">>> Predicting smooth curves...")
  yhat_smooth <- tradeSeq::predictSmooth(sce_obj, gene = valid_genes, nPoints = 100, tidy = TRUE)

  # Z-score standardization
  yhat_visual <- yhat_smooth %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(scaled_expr = scale(yhat)) %>%
    dplyr::ungroup()

  # Color logic
  if (is.null(colors)) {
    base_cols <- c("EZH2" = "#D58D4E", "SLC7A11" = "#86C495")
    extra_genes <- setdiff(valid_genes, names(base_cols))
    if(length(extra_genes)>0) {
      extra_cols <- grDevices::rainbow(length(extra_genes))
      names(extra_cols) <- extra_genes
      colors <- c(base_cols, extra_cols)
    } else { colors <- base_cols }
  }

  ggplot2::ggplot(yhat_visual, ggplot2::aes(x = time, y = scaled_expr, color = gene)) +
    ggplot2::geom_line(linewidth = 2, alpha = 0.9) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(x = "Pseudotime", y = "Standardized Expression") +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 14, face = "bold"),
      axis.line = ggplot2::element_line(linewidth = 1.0),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 14)
    )
}
