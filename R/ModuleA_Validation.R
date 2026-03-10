#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual scale_y_continuous theme_bw labs theme element_blank element_text
#' @importFrom dplyr select mutate everything filter
#' @importFrom tidyr pivot_longer
#' @importFrom survival survdiff survfit Surv
#' @importFrom survminer ggsurvplot
#' @importFrom stats median pchisq
NULL

# ==============================================================================
# Part 1: Single Gene Diagnostic Value Visualization (Figure I)
# ==============================================================================

#' Plot Single Gene Diagnostic Importance (Bar Chart)
#'
#' Visualizes the importance score of a specific gene across RF, SVM, GBM and Average.
#'
#' @param diag_result Result object from `run_diagnostic_model`.
#' @param gene The gene symbol to visualize (e.g., "EZH2").
#' @return A ggplot object.
#' @export
plot_gene_diag_score <- function(diag_result, gene = "EZH2") {

  # 1. Extract data
  final_imp <- diag_result$importance

  # Check if gene exists
  if(!gene %in% final_imp$Gene) {
    stop(paste("Gene", gene, "not found in the diagnostic model results."))
  }

  # 2. Extract and transform data (wide -> long)
  # Extract this gene's row
  gene_data <- final_imp[final_imp$Gene == gene, ]
  rank_val <- which(final_imp$Gene == gene)
  total_genes <- nrow(final_imp)

  # Use tidyr to reshape data for ggplot color mapping
  plot_dat <- gene_data %>%
    dplyr::select(-Gene) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "Model_Type",
      values_to = "Score"
    ) %>%
    # Remove _Score suffix
    dplyr::mutate(Model_Type = gsub("_Score", "", Model_Type)) %>%
    # Rename Average to "Average" (in case previous step didn't handle it)
    dplyr::mutate(Model_Type = ifelse(Model_Type == "Average", "Average", Model_Type))

  # 3. Set plot order: single models first, Average last
  model_levels <- c(setdiff(plot_dat$Model_Type, "Average"), "Average")
  plot_dat$Model_Type <- factor(plot_dat$Model_Type, levels = model_levels)

  # 4. Define colors (replicating your color scheme)
  # Average uses red (#E64B35), others use cool colors
  my_colors <- c("RF" = "#4DBBD5", "SVM" = "#00A087", "GBM" = "#3C5488", "Average" = "#E64B35")

  # 5. Plot
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = Model_Type, y = Score, fill = Model_Type)) +
    ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    # Value labels
    ggplot2::geom_text(ggplot2::aes(label = round(Score, 1)), vjust = -0.5, size = 5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = my_colors) +
    # Expand Y-axis slightly to prevent labels being cut off
    ggplot2::scale_y_continuous(limits = c(0, 115), expand = c(0,0)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = paste0("Diagnostic Importance: ", gene),
      subtitle = paste0("Overall Rank: ", rank_val, " / ", total_genes),
      x = "",
      y = "Importance Score (0-100)"
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray30", face = "italic"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(face = "bold")
    )

  return(p)
}


# ==============================================================================
# Part 2: Single Gene Survival Analysis (Figure J)
# ==============================================================================

#' Plot Single Gene Survival Analysis (KM Curve)
#'
#' Performs KM analysis for a single gene using internal TCGA data.
#' Splits samples by median expression (High vs Low).
#'
#' @param gene The gene symbol to analyze.
#' @param expr_mat Expression matrix (optional, defaults to internal data).
#' @param clin_data Clinical data (optional, defaults to internal data).
#' @return A ggsurvplot object.
#' @export
plot_gene_survival <- function(gene = "EZH2", expr_mat = NULL, clin_data = NULL) {

  # 1. Auto-load data
  if(is.null(expr_mat)) expr_mat <- get("TCGA_Expr_Mat")
  if(is.null(clin_data)) clin_data <- get("TCGA_Clin_Data")

  # 2. Check gene
  if(!gene %in% rownames(expr_mat)) {
    stop(paste("Gene", gene, "not found in expression matrix."))
  }

  # 3. Data integration and cleaning
  # Extract Tumor samples (Group == "Tumor")
  tumor_clin <- clin_data[clin_data$Group == "Tumor" & !is.na(clin_data$OS.time) & !is.na(clin_data$OS), ]
  common_samples <- intersect(tumor_clin$SampleID, colnames(expr_mat))

  # Extract expression
  gene_expr <- as.numeric(expr_mat[gene, common_samples])

  # Build survival analysis dataframe
  surv_df <- tumor_clin[match(common_samples, tumor_clin$SampleID), ]
  surv_df$Expression <- gene_expr

  # Time unit conversion: Days -> Months (divide by ~30.4 or 30, using 30 for readability)
  surv_df$Time_Months <- surv_df$OS.time

  # 4. Grouping (High/Low by Median)
  cut_point <- stats::median(surv_df$Expression)
  surv_df$Group <- ifelse(surv_df$Expression > cut_point, "High", "Low")
  # Ensure Low is baseline (usually blue), High is comparison level (usually red/yellow)
  surv_df$Group <- factor(surv_df$Group, levels = c("Low", "High"))

  # 5. Calculate P-value (Log-Rank Test)
  diff <- survival::survdiff(survival::Surv(Time_Months, OS) ~ Group, data = surv_df)
  pval <- 1 - stats::pchisq(diff$chisq, length(diff$n) - 1)

  # Format P-value text
  if(pval < 0.001) {
    p_label <- "P < 0.001"
  } else {
    p_label <- paste0("P = ", sprintf("%.3f", pval))
  }

  # 6. Fit survival curve
  fit <- survival::survfit(survival::Surv(Time_Months, OS) ~ Group, data = surv_df)

  # 7. Plot (using survminer to replicate your style)
  p <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    pval = p_label,        # Display P-value
    pval.size = 6,
    conf.int = TRUE,       # Display confidence interval (shadow)
    risk.table = TRUE,     # Display risk table
    risk.table.col = "strata",
    palette = c("#2E9FDF", "#E7B800"), # Classic blue-yellow color scheme (JCO style variant)
    legend.labs = c("Low", "High"),
    legend.title = gene,
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    xlim = c(0, 120),      # Limit X-axis to 10 years (120 months)
    break.time.by = 20,    # Tick every 20 months
    ggtheme = ggplot2::theme_classic(base_size = 14),

    # Risk table and censoring plot details
    ncensor.plot = TRUE,   # Display small vertical lines for censored points
    ncensor.plot.height = 0.2,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE # Hide Y-axis text in risk table for cleaner look
  )

  return(p)
}

#p_bar <- plot_gene_diag_score(res_diag, gene = "EZH2")
#plot_gene_survival(gene = "EZH2")
#print(p_bar)

# 3. Plot Figure J: Survival curve (auto-calls built-in TCGA data)
# This step plots directly, no print needed
#plot_gene_survival(gene = "EZH2")
