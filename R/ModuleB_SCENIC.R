
#' @importFrom data.table fwrite fread
#' @importFrom Matrix t
#' @importFrom stats cor.test lm cor quantile sd
#' @importFrom dplyr filter arrange desc mutate group_by summarise select everything %>% across all_of n
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_vline geom_hline theme_bw labs geom_segment coord_flip scale_color_manual scale_fill_gradient2 scale_size scale_fill_manual annotate geom_density aes_string stat_density_2d scale_fill_gradientn geom_smooth theme_classic theme element_blank element_line element_text margin geom_rect alpha guides element_rect guide_legend scale_color_viridis_c
#' @importFrom ggrepel geom_text_repel

#' @importFrom tidyr pivot_longer
NULL

# ==============================================================================
# Part 1: Data I/O (Export/Import)
# ==============================================================================

#' Export Count Matrix for SCENIC
#' @export
export_scenic_data <- function(seurat_obj, output_dir = ".", file_name = "scenic_input.csv", target_gene = NULL) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_path <- file.path(output_dir, file_name)

  if (!is.null(target_gene)) {
    if (!target_gene %in% rownames(seurat_obj)) stop(paste("Error:", target_gene, "not found."))
    expr_vals <- tryCatch(
      Seurat::GetAssayData(seurat_obj, layer = "data")[target_gene, ],
      error = function(e) Seurat::GetAssayData(seurat_obj, layer = "data")[target_gene, ]
    )
    pos_cells <- names(expr_vals)[expr_vals > 0]
    message(paste(">>> Filtering cells: Keeping", length(pos_cells), "cells with", target_gene, "> 0."))
    seurat_obj <- subset(seurat_obj, cells = pos_cells)
  }

  message(">>> Extracting count matrix...")
  counts_mat <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  )
  message(">>> Transposing matrix (Cells x Genes)...")
  counts_t <- Matrix::t(counts_mat)
  counts_df <- as.data.frame(as.matrix(counts_t))
  counts_df <- cbind(CellID = rownames(counts_df), counts_df)

  message(paste(">>> Saving matrix to:", out_path))
  data.table::fwrite(counts_df, file = out_path, row.names = FALSE)
  message("Export complete!")
}

#' Import SCENIC AUC Matrix
#' @export
import_scenic_results <- function(seurat_obj, auc_file) {
  if (!file.exists(auc_file)) stop(paste("Error: AUC file not found at", auc_file))

  message(">>> Reading SCENIC AUC matrix...")
  auc_df <- data.table::fread(auc_file, data.table = FALSE)
  rownames(auc_df) <- auc_df[, 1]
  auc_df <- auc_df[, -1]

  colnames(auc_df) <- gsub("\\(\\+\\)", "", colnames(auc_df))

  auc_mat <- as.matrix(t(auc_df))

  # [Key Step] Take intersection (Step 2 logic)
  common_cells <- intersect(colnames(seurat_obj), colnames(auc_mat))

  if (length(common_cells) == 0) stop("Error: No matching cell barcodes.")

  message(paste(">>> Intersecting cells. Retained:", length(common_cells)))

  # Align matrix and object
  auc_mat <- auc_mat[, common_cells, drop = FALSE]
  seurat_obj <- subset(seurat_obj, cells = common_cells)

  # Add Assay
  seurat_obj[["AUC"]] <- Seurat::CreateAssayObject(data = auc_mat)
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = "AUC")

  Seurat::DefaultAssay(seurat_obj) <- "AUC"
  return(seurat_obj)
}

# ==============================================================================
# Part 2: Mediator TF Screening (Screening) - Corrected cell filtering logic
# ==============================================================================

#' Screen for Mediator Transcription Factors
#'
#' @param seurat_obj Seurat object with "AUC" assay.
#' @param target_gene The upstream driver gene (e.g., "EZH2").
#' @param downstream_gene The downstream effector gene (e.g., "SLC7A11").
#' @return A list containing the results data frame (`res_df`) and a screening plot (`plot`).
#' @export
screen_mediator_tfs <- function(seurat_obj, target_gene = "EZH2", downstream_gene = "SLC7A11") {

  message(paste(">>> Screening TFs for axis:", target_gene, "-> TF ->", downstream_gene))

  # --- 1. Data extraction ---
  get_data_layer <- function(obj, assay, feature) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ],
      error = function(e) Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ]
    )
  }

  vec_ezh2 <- get_data_layer(seurat_obj, "RNA", target_gene)
  vec_slc  <- get_data_layer(seurat_obj, "RNA", downstream_gene)

  # Extract AUC matrix
  auc_mat <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "AUC", layer = "data"),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = "AUC", layer = "data")
  )

  # --- 2. Filter valid cells (Step 1 + Step 3 Combined Logic) ---
  # Local code Step 1: Filter EZH2 > 0
  cells_target_pos <- names(vec_ezh2[vec_ezh2 > 0])

  # Local code Step 3: Filter SLC7A11 > 0
  cells_down_pos <- names(vec_slc[vec_slc > 0])

  # [Core Fix] Take intersection! Ensure cells are double-positive to replicate your local Step 3 operation on top of Step 1
  valid_cells <- intersect(cells_target_pos, cells_down_pos)

  # Fallback strategy
  if(length(valid_cells) < 50) {
    warning("Double-positive cells < 50. Using all cells instead (Comparison might differ).")
    valid_cells <- colnames(seurat_obj)
  }

  message(paste0(">>> Valid cells (", target_gene, "+ & ", downstream_gene, "+): ", length(valid_cells)))

  # --- 3. Loop through TFs to calculate ---
  tf_list <- rownames(auc_mat)
  tf_list <- tf_list[tf_list != target_gene] # Remove EZH2 itself

  res_df <- data.frame(TF = tf_list,
                       Cor_EZH2 = NA, P_EZH2 = NA,
                       Cor_SLC = NA, P_SLC = NA,
                       Score = NA)

  # Pre-extract data
  vec_ezh2_val <- vec_ezh2[valid_cells]
  vec_slc_val  <- vec_slc[valid_cells]
  auc_mat_val  <- auc_mat[, valid_cells]

  message(">>> Looping through TFs...")
  for(i in 1:length(tf_list)) {
    tf <- tf_list[i]
    vec_tf <- auc_mat_val[tf, ]

    # Path A: EZH2 -> TF
    res1 <- cor.test(vec_ezh2_val, vec_tf, method = "spearman")

    # Path B: TF -> SLC7A11
    res2 <- cor.test(vec_tf, vec_slc_val, method = "spearman")

    res_df$Cor_EZH2[i] <- as.numeric(res1$estimate)
    res_df$P_EZH2[i]   <- res1$p.value
    res_df$Cor_SLC[i]  <- as.numeric(res2$estimate)
    res_df$P_SLC[i]    <- res2$p.value

    res_df$Score[i] <- abs(as.numeric(res1$estimate) * as.numeric(res2$estimate))
  }

  # --- 4. Filter significant results ---
  res_sig <- res_df %>%
    dplyr::filter(P_EZH2 < 0.05 & P_SLC < 0.05) %>%
    dplyr::arrange(dplyr::desc(Score))

  message("Top 10 Candidates:")
  print(head(res_sig, 10))

  # --- 5. Plot ---
  p <- plot_quadrant_screening_v2(res_sig, target_gene, downstream_gene)

  return(list(results = res_sig, plot = p))
}

# --- Internal plotting function ---
plot_quadrant_screening_v2 <- function(res_sig, target, down) {

  # Top 5 Label
  top_n_show <- 5
  res_plot <- res_sig %>%
    dplyr::arrange(dplyr::desc(Score)) %>%
    dplyr::mutate(Label = ifelse(dplyr::row_number() <= top_n_show, TF, NA))

  # Color scheme
  color_low  <- "#76C6C3"
  color_mid  <- "white"
  color_high <- "#E698B9"
  mid_point <- mean(range(res_plot$Score, na.rm = TRUE))

  p <- ggplot2::ggplot(res_plot, ggplot2::aes(x = Cor_EZH2, y = Cor_SLC)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", size = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", size = 0.8) +

    # Shape 21 points
    ggplot2::geom_point(ggplot2::aes(fill = Score, size = Score), shape = 21, color = "grey30", stroke = 0.2, alpha = 0.9) +

    # Text Repel
    ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 4.5, fontface = "bold.italic",
                             color = "black", box.padding = 0.6, max.overlaps = 50,
                             segment.color = "grey50", seed = 123) +

    # Colors
    ggplot2::scale_fill_gradient2(low = color_low, mid = color_mid, high = color_high,
                                  midpoint = mid_point, name = "Combined Score") +
    ggplot2::scale_size(range = c(2, 7), guide = "none") +

    # Labs & Theme
    ggplot2::labs(x = paste("Correlation with", target, "(RNA)"),
                  y = paste("Correlation with", down, "(RNA)")) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = c(0.85, 0.2),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9)
    )

  return(p)
}

# ==============================================================================
# Part 3: Mediation Analysis & Part 4: Validation Plotting (Unchanged)
# ==============================================================================

#' Run Mediation Analysis with Bootstrap
#' @export
run_mediation_analysis <- function(seurat_obj, target_gene, mediator_tf, downstream_gene, n_sims = 500) {
  .check_cran_packages("mediation")

  message(paste(">>> Running Mediation:", target_gene, "->", mediator_tf, "->", downstream_gene))

  get_data_layer <- function(obj, assay, feature) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ],
      error = function(e) Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ]
    )
  }

  vec_x <- get_data_layer(seurat_obj, "RNA", target_gene)
  vec_y <- get_data_layer(seurat_obj, "RNA", downstream_gene)
  vec_m <- get_data_layer(seurat_obj, "AUC", mediator_tf)

  # Apply double filtering logic
  cells_x <- names(vec_x[vec_x > 0])
  cells_y <- names(vec_y[vec_y > 0])
  valid_cells <- intersect(cells_x, cells_y)

  if(length(valid_cells) < 50) valid_cells <- colnames(seurat_obj)

  df_med <- data.frame(
    X = vec_x[valid_cells],
    Y = vec_y[valid_cells],
    M = vec_m[valid_cells]
  )

  message("   Fitting regression models...")
  model_m <- lm(M ~ X, data = df_med)
  model_y <- lm(Y ~ X + M, data = df_med)

  message(paste("   Bootstrapping (", n_sims, "sims )..."))
  set.seed(123)
  med_out <- mediation::mediate(model_m, model_y, treat = "X", mediator = "M", boot = TRUE, sims = n_sims)

  p <- plot_mediation_density(med_out)
  message(paste("Mediation complete. Prop. Mediated:", round(med_out$n0 * 100, 1), "%"))
  return(list(model = med_out, plot = p))
}

plot_mediation_density <- function(med_out) {
  sim_data <- data.frame(
    Indirect = as.numeric(med_out$d0.sims),
    Direct   = as.numeric(med_out$z0.sims)
  ) %>% tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Effect_Type", values_to = "Estimate")

  mean_vals <- sim_data %>% dplyr::group_by(Effect_Type) %>% dplyr::summarise(Mean_Val = mean(Estimate))
  prop_mediated <- round(med_out$n0 * 100, 1)

  color_indirect <- "#C89F9C"
  color_direct   <- "#8FA39D"

  p <- ggplot2::ggplot(sim_data, ggplot2::aes(x = Estimate, fill = Effect_Type)) +
    ggplot2::geom_density(alpha = 0.7, color = NA) +
    ggplot2::geom_vline(data = mean_vals, ggplot2::aes(xintercept = Mean_Val, color = Effect_Type), linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("Indirect" = color_indirect, "Direct" = color_direct)) +
    ggplot2::scale_color_manual(values = c("Indirect" = "#A87F7C", "Direct" = "#6F837D")) +
    ggplot2::annotate("text", x = max(sim_data$Estimate) * 0.9, y = Inf,
                      label = paste0("Mediated Proportion\n", prop_mediated, "%"),
                      hjust = 1, vjust = 1.5, size = 5, fontface = "bold.italic", color = color_indirect) +
    ggplot2::labs(x = "Effect Size Estimate (Bootstrap)", y = "Density") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = c(0.15, 0.85), legend.title = ggplot2::element_blank())
  return(p)
}

#' Plot Earth-Tone Density Correlation
#' @export
plot_earth_density <- function(seurat_obj, x_feature, y_feature, x_assay = "RNA", y_assay = "AUC") {
  get_data_layer <- function(obj, assay, feature) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ],
      error = function(e) Seurat::GetAssayData(obj, assay = assay, layer = "data")[feature, ]
    )
  }
  x_val <- get_data_layer(seurat_obj, x_assay, x_feature)
  y_val <- get_data_layer(seurat_obj, y_assay, y_feature)
  df_plot <- data.frame(X = x_val, Y = y_val)
  df_plot <- df_plot %>% dplyr::filter(X > 0 & Y > 0)
  res <- stats::cor.test(df_plot$X, df_plot$Y, method = "spearman")
  rho_val <- round(res$estimate, 3)
  p_text <- ifelse(res$p.value < 0.001, "P < 0.001", paste0("P = ", signif(res$p.value, 2)))
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = X, y = Y)) +
    ggplot2::stat_density_2d(geom = "polygon", ggplot2::aes(fill = ..level..), contour = TRUE) +
    ggplot2::scale_fill_gradientn(colors = c("white", "#469D76", "#D58D4E"), name = "Density") +
    ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed", size = 0.8, se = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("Spearman R = ", rho_val, "\n", p_text),
                      hjust = -0.1, vjust = 1.2, size = 5, fontface = "bold") +
    ggplot2::labs(x = x_feature, y = y_feature) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = "right")
  return(p)
}
# ==============================================================================
# Part 5: Linkage Analysis (Lollipop + Heatmap) - Append to file end
# ==============================================================================

#' Run TF Enrichment Analysis (Lollipop + Heatmap Linkage)
#'
#' @param seurat_obj Seurat object (Raw object, function will filter target_gene > 0 internally).
#' @param target_gene Gene to correlate with (default "EZH2").
#' @return A list containing `lollipop_plot` and `heatmap_obj`.
#' @export
run_tf_enrichment_analysis <- function(seurat_obj, target_gene = "EZH2") {

  message(paste0(">>> Starting Linked Analysis for ", target_gene, "..."))

  # --- 0. Simulate local environment: only extract cells with Target > 0 ---
  # This ensures results match exactly what you get locally with hep_ezh2_pos
  vec_target_all <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")[target_gene, ],
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")[target_gene, ]
  )

  pos_cells <- names(vec_target_all[vec_target_all > 0])
  message(paste(">>> Internally subsetting:", length(pos_cells), "cells with", target_gene, "> 0."))

  if(length(pos_cells) == 0) stop("No positive cells found for target gene.")

  # --- 1. Data preparation (based on pos_cells) ---
  # mat_auc: rows=TF, columns=Cell
  mat_auc <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "AUC", layer = "data")[, pos_cells],
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = "AUC", layer = "data")[, pos_cells]
  )

  # vec_target: vector
  vec_target <- vec_target_all[pos_cells]

  # Remove Target itself
  if(target_gene %in% rownames(mat_auc)) {
    mat_auc <- mat_auc[rownames(mat_auc) != target_gene, ]
  }

  # --- 2. Fast correlation calculation (Matrix operations - replicate local code) ---
  message(">>> Calculating correlations...")
  cor_val <- stats::cor(t(as.matrix(mat_auc)), vec_target, method = "spearman")
  cor_df <- data.frame(TF = rownames(cor_val), R = as.numeric(cor_val))

  # --- 3. Filter Top (replicate local code) ---
  # Note: must use utils::head to prevent dplyr conflict
  top_pos <- cor_df %>% dplyr::arrange(dplyr::desc(R)) %>% utils::head(10) %>% dplyr::mutate(Type = "Positive")
  top_neg <- cor_df %>% dplyr::arrange(R) %>% utils::head(5) %>% dplyr::mutate(Type = "Negative")

  plot_data <- rbind(top_pos, top_neg)
  # Lock factor order
  plot_data$TF <- factor(plot_data$TF, levels = plot_data$TF[order(plot_data$R)])

  selected_tfs <- as.character(plot_data$TF)
  message(paste(">>> Selected TFs:", paste(selected_tfs, collapse = ", ")))

  # --- 4. Draw advanced lollipop chart (replicate local ggplot code) ---
  colors_lolly <- c("Positive" = "#B69FD6", "Negative" = "#C8B6A5")

  p_lolly <- ggplot2::ggplot(plot_data, ggplot2::aes(x = TF, y = R)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = TF, xend = TF, y = 0, yend = R), color = "grey70", size = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = Type), size = 5) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = colors_lolly) +
    ggplot2::labs(x = "", y = paste("Spearman Correlation with", target_gene)) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )

  # --- 5. Draw heatmap (replicate local Heatmap code) ---
  # Use selected_tfs and pos_cells data from above

  .check_cran_packages(c("ComplexHeatmap", "circlize"))

  valid_tfs <- intersect(selected_tfs, rownames(mat_auc))
  sub_auc <- t(as.matrix(mat_auc[valid_tfs, ]))

  df_smooth <- as.data.frame(sub_auc)
  df_smooth$Target_Raw <- vec_target

  # Binning (replicate local logic: quantile)
  # No need to filter >0 here since pos_cells is already >0
  df_smooth$Bin <- tryCatch({
    cut(df_smooth$Target_Raw,
        breaks = stats::quantile(df_smooth$Target_Raw, probs = seq(0, 1, 0.05)),
        include.lowest = TRUE)
  }, error = function(e) {
    # Fallback: if too few data points for quantile, use equal bins
    cut(df_smooth$Target_Raw, breaks = 20)
  })

  mat_heatmap <- df_smooth %>%
    dplyr::group_by(Bin) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(valid_tfs), \(x) mean(x, na.rm = TRUE))) %>%
    dplyr::select(-Bin) %>%
    as.matrix() %>%
    t()

  mat_heatmap_scaled <- t(scale(t(mat_heatmap)))

  # Color scheme (replicate local)
  col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#0072B5", "white", "#D58D4E"))
  anno_col = circlize::colorRamp2(c(1, 20), c("#F0F0F0", "#BC3C29"))

  ha_top = ComplexHeatmap::HeatmapAnnotation(
    Trend = 1:20,
    col = list(Trend = anno_col),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.3, "cm"),
    border = FALSE
  )

  hm_obj <- ComplexHeatmap::Heatmap(mat_heatmap_scaled,
                                    name = "TF Activity\n(Z-score)",
                                    col = col_fun,
                                    cluster_columns = FALSE,
                                    cluster_rows = TRUE,
                                    show_column_names = FALSE,
                                    row_names_gp = grid::gpar(fontsize = 12, fontface = "italic"),
                                    top_annotation = ha_top,
                                    rect_gp = grid::gpar(col = "white", lwd = 1),
                                    row_dend_width = grid::unit(1.5, "cm"),
                                    border = FALSE
  )

  message(">>> Analysis Complete! Returning plots list.")
  return(list(lollipop = p_lolly, heatmap = hm_obj, top_tfs = selected_tfs))
}
