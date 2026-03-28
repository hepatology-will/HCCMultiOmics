#' @importFrom Seurat GetAssayData AddModuleScore FetchData GetTissueCoordinates
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn theme_void theme
#' @importFrom ggplot2 element_text element_blank element_rect unit
#' @importFrom dplyr filter select arrange
#' @importFrom viridis magma plasma viridis
#' @importFrom patchwork plot_layout
#' @importFrom ggpubr stat_cor
NULL

# ==============================================================================
# Module C: Spatial Transcriptomics Validation
# ==============================================================================

# Helper function to check and install required packages
.check_spatial_packages <- function() {
  pkgs_needed <- c("dorothea", "decoupleR")
  
  for (pkg in pkgs_needed) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(sprintf(">>> Installing %s...", pkg))
      if (!require("BiocManager", quietly = TRUE)) {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
        install.packages("BiocManager")
      }
      BiocManager::install(c("dorothea", "decoupleR"), ask = FALSE)
    }
  }
  
  library(dorothea, character.only = TRUE)
  library(decoupleR, character.only = TRUE)
}

#' Run Spatial TF Validation Analysis
#'
#' This function validates transcription factor activity from Step 5 using spatial
#' transcriptomics data. It calculates TF activity scores using DoRothEA database
#' and generates spatial distribution plots and correlation plots.
#'
#' @param seurat_obj A Seurat object containing spatial transcriptomics data
#' @param target_gene Target gene name (used as upstream gene, e.g., from Step 2)
#' @param tf_list Vector of transcription factors to validate (from Step 5)
#' @param driver_genes Vector of driver genes (for ferroptosis suppressor score)
#' @param suppressor_genes Vector of suppressor genes (for ferroptosis suppressor score)
#' @return A list containing modified Seurat object and validation results
#' @export
run_spatial_validation <- function(seurat_obj, target_gene, tf_list = NULL,
                                   driver_genes = NULL, suppressor_genes = NULL) {

  # Check and install required packages
  .check_spatial_packages()

  message(">>> Running Spatial TF Validation...")

  # 1. Calculate TF activity using DoRothEA
  if (!is.null(tf_list) && length(tf_list) > 0) {
    message(">>> Calculating TF activity scores...")

    # Load DoRothEA database
    data(dorothea_hs, package = "dorothea")

    # Filter high confidence network (A, B grades)
    net <- dorothea_hs %>%
      dplyr::filter(confidence %in% c("A", "B")) %>%
      dplyr::select(source = tf, target = target, mor = mor)

    # Prepare expression matrix
    mat <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "SCT", slot = "data"))

    # Run decoupleR (ulm algorithm)
    acts <- decoupleR::run_ulm(mat = mat, network = net,
                                 .source = 'source', .target = 'target', .mor = 'mor')

    # Add TF activity scores to Seurat object for each TF
    for (tf in tf_list) {
      tf_res <- acts %>%
        dplyr::filter(source == tf) %>%
        dplyr::filter(statistic == "ulm")

      tf_scores <- tf_res$score
      names(tf_scores) <- tf_res$condition

      col_name <- paste0(tf, "_TF_Activity")
      seurat_obj[[col_name]] <- tf_scores[colnames(seurat_obj)]
    }
  }

  # 2. Calculate target gene module score (instead of PRC2)
  # User's target gene serves as the upstream "EZH2 complex"
  message(sprintf(">>> Calculating module score for target gene: %s", target_gene))
  target_gene_list <- list(c(target_gene))
  names(target_gene_list) <- "Target_Gene"
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = target_gene_list, name = "TargetGene")

  # 3. Calculate ferroptosis suppressor score using user-provided genes
  # Combine driver and suppressor genes as ferroptosis-related genes
  ferro_genes <- c()
  if (!is.null(driver_genes) && length(driver_genes) > 0) {
    ferro_genes <- c(ferro_genes, driver_genes)
  }
  if (!is.null(suppressor_genes) && length(suppressor_genes) > 0) {
    ferro_genes <- c(ferro_genes, suppressor_genes)
  }

  if (length(ferro_genes) > 0) {
    message(sprintf(">>> Calculating ferroptosis suppressor score with %d genes...", length(ferro_genes)))
    ferro_list <- list(ferro_genes)
    names(ferro_list) <- "Ferro_Suppressor"
    seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = ferro_list, name = "Ferro_Suppressor")
  }

  # 4. Prepare plot data
  message(">>> Preparing plot data...")

  vars_to_fetch <- c("TargetGene1")

  if (!is.null(tf_list) && length(tf_list) > 0) {
    tf_col <- paste0(tf_list[1], "_TF_Activity")
    vars_to_fetch <- c(vars_to_fetch, tf_col)
  }

  if (length(ferro_genes) > 0) {
    vars_to_fetch <- c(vars_to_fetch, "Ferro_Suppressor1")
  }

  plot_data <- Seurat::FetchData(seurat_obj, vars = vars_to_fetch)

  if (!is.null(tf_list) && length(tf_list) > 0) {
    tf_col <- paste0(tf_list[1], "_TF_Activity")
    colnames(plot_data) <- c("Upstream_Target", "Mid_TF_Activity", "Down_Ferro_Score")
  } else {
    colnames(plot_data) <- c("Upstream_Target", "Down_Ferro_Score")
  }

  # Get spatial coordinates
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  if (!("imagecol" %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("imagerow", "imagecol")
  }
  plot_final <- cbind(coords, plot_data)

  result <- list(
    seurat_obj = seurat_obj,
    plot_data = plot_final,
    tf_activity = if (!is.null(tf_list)) acts else NULL
  )

  return(result)
}

#' Plot Spatial Distribution
#'
#' Generates spatial distribution plots for target gene, TF activity, and ferroptosis score
#'
#' @param plot_data Data frame with spatial coordinates and scores
#' @param output_dir Output directory for plots
#' @param target_gene Target gene name
#' @param tf_name TF name used (first TF if multiple)
#' @export
plot_spatial_distribution <- function(plot_data, output_dir, target_gene, tf_name = NULL) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Define color schemes (Ink-style high contrast)
  col_ink_blue <- colorRampPalette(c("#F0F4F8", "#B0C4DE", "#607D8B", "#2C3E50", "#1A2530"))(100)
  col_ink_red <- colorRampPalette(c("#F9F2F2", "#D8BFD8", "#A0522D", "#800000", "#3E0000"))(100)
  col_ink_green <- colorRampPalette(c("#F0F5F0", "#8FBC8F", "#556B2F", "#006400", "#003300"))(100)

  plot_spatial_ink <- function(data, col_name, custom_colors, title) {
    data <- data %>% dplyr::arrange(!!dplyr::sym(col_name))

    p <- ggplot2::ggplot(data, ggplot2::aes(x = imagecol, y = -imagerow, color = !!dplyr::sym(col_name))) +
      ggplot2::geom_point(size = 2.8, alpha = 1) +
      ggplot2::scale_color_gradientn(colors = custom_colors) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 0, 0, 0),
        legend.position = "right",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8, color = "grey30"),
        legend.key.height = ggplot2::unit(0.8, "cm"),
        legend.key.width = ggplot2::unit(0.2, "cm")
      ) +
      ggplot2::coord_fixed()

    return(p)
  }

  # Plot target gene (upstream)
  p1 <- plot_spatial_ink(plot_data, "Upstream_Target", col_ink_blue,
                          paste0("Upstream: ", target_gene))

  plots_list <- list(p1)

  # Plot TF activity (mid) if available
  if ("Mid_TF_Activity" %in% colnames(plot_data)) {
    p2 <- plot_spatial_ink(plot_data, "Mid_TF_Activity", col_ink_red,
                            paste0("Mediator: ", tf_name, " Activity"))
    plots_list <- c(plots_list, list(p2))
  }

  # Plot ferroptosis suppressor score (downstream)
  if ("Down_Ferro_Score" %in% colnames(plot_data)) {
    p3 <- plot_spatial_ink(plot_data, "Down_Ferro_Score", col_ink_green,
                            "Effector: Ferroptosis Suppressor")
    plots_list <- c(plots_list, list(p3))
  }

  # Combine and save
  if (length(plots_list) > 1) {
    combined_plot <- plots_list[[1]]
    for (i in 2:length(plots_list)) {
      combined_plot <- combined_plot + plots_list[[i]]
    }
  } else {
    combined_plot <- plots_list[[1]]
  }

  pdf_file <- file.path(output_dir, "spatial_distribution.pdf")
  pdf(pdf_file, width = 6 * length(plots_list), height = 6)
  print(combined_plot)
  dev.off()

  message(sprintf(">>> Saved: %s", pdf_file))

  return(invisible(pdf_file))
}

#' Plot Correlation Hexbin
#'
#' Generates hexbin correlation plots between upstream, mediator, and downstream
#'
#' @param plot_data Data frame with scores
#' @param output_dir Output directory for plots
#' @export
plot_correlation_hexbin <- function(plot_data, output_dir) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Define color schemes
  col_ink_blue <- colorRampPalette(c("#F0F4F8", "#B0C4DE", "#607D8B", "#2C3E50", "#1A2530"))(100)
  col_ink_red <- colorRampPalette(c("#F9F2F2", "#D8BFD8", "#A0522D", "#800000", "#3E0000"))(100)
  col_ink_green <- colorRampPalette(c("#F0F5F0", "#8FBC8F", "#556B2F", "#006400", "#003300"))(100)

  plot_hex_ink <- function(data, x_var, y_var, color_palette, x_lab = NULL, y_lab = NULL) {

    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!dplyr::sym(x_var), y = !!dplyr::sym(y_var))) +
      ggplot2::geom_hex(bins = 55, color = "white", size = 0.1) +
      ggplot2::scale_fill_gradientn(colors = color_palette) +
      ggplot2::geom_smooth(method = "lm", color = "#2F2F2F", size = 0.8, linetype = "solid", se = FALSE) +
      ggpubr::stat_cor(method = "pearson", size = 4, color = "#2F2F2F",
                       label.x.npc = "left", label.y.npc = "top") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(color = "#D0D0D0", size = 0.5),
        axis.line.y = ggplot2::element_line(color = "#D0D0D0", size = 0.5),
        axis.text = ggplot2::element_text(color = "grey40", size = 10),
        axis.title = ggplot2::element_text(color = "grey30", size = 11, face = "plain"),
        legend.position = "none"
      ) +
      ggplot2::labs(x = x_lab, y = y_lab)

    return(p)
  }

  plots_list <- list()

  # Get column names
  cols <- colnames(plot_data)
  cols <- cols[cols %in% c("Upstream_Target", "Mid_TF_Activity", "Down_Ferro_Score")]

  if (length(cols) >= 2) {
    # Upstream vs Mid
    if ("Upstream_Target" %in% cols && "Mid_TF_Activity" %in% cols) {
      c1 <- plot_hex_ink(plot_data, "Upstream_Target", "Mid_TF_Activity",
                         col_ink_blue,
                         x_lab = paste0("Target Gene (", cols[1], ")"), y_lab = "TF Activity")
      plots_list <- c(plots_list, list(c1))
    }

    # Mid vs Down
    if ("Mid_TF_Activity" %in% cols && "Down_Ferro_Score" %in% cols) {
      c2 <- plot_hex_ink(plot_data, "Mid_TF_Activity", "Down_Ferro_Score",
                         col_ink_red,
                         x_lab = "TF Activity", y_lab = "Ferro Suppressor")
      plots_list <- c(plots_list, list(c2))
    }

    # Upstream vs Down
    if ("Upstream_Target" %in% cols && "Down_Ferro_Score" %in% cols) {
      c3 <- plot_hex_ink(plot_data, "Upstream_Target", "Down_Ferro_Score",
                         col_ink_green,
                         x_lab = paste0("Target Gene (", cols[1], ")"), y_lab = "Ferro Suppressor")
      plots_list <- c(plots_list, list(c3))
    }
  }

  if (length(plots_list) > 0) {
    combined_plot <- plots_list[[1]]
    for (i in 2:length(plots_list)) {
      combined_plot <- combined_plot + plots_list[[i]]
    }

    pdf_file <- file.path(output_dir, "correlation_hexbin.pdf")
    pdf(pdf_file, width = 4 * length(plots_list), height = 4)
    print(combined_plot)
    dev.off()

    message(sprintf(">>> Saved: %s", pdf_file))
  }

  return(invisible(NULL))
}
