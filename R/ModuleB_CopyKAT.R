#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual scale_y_continuous labs theme_classic theme element_blank element_line element_text geom_col position_fill position_stack ggtitle
#' @importFrom dplyr group_by summarise ungroup filter %>% mutate n
#' @importFrom scales percent
#' @importFrom patchwork wrap_plots
#' @importFrom stats reorder
NULL

# ==============================================================================
# Helper 1: Ultra-low saturation Morandi color scheme (reuse previous definition, for Site coloring)
# ==============================================================================
get_morandi_colors <- function(n) {
  muted_morandi_12 <- c(
    "#8F9EAB", "#C59993", "#9FAFA3", "#D6C6A8",
    "#A796B0", "#9C9C9C", "#7D8C8E", "#B0A49D",
    "#A3B6C2", "#D1B8B4", "#B4BDA6", "#C8C3BC"
  )
  if (n > length(muted_morandi_12)) colorRampPalette(muted_morandi_12)(n) else muted_morandi_12[1:n]
}

# ==============================================================================
# Helper 2: CopyKAT semantic color scheme (CNV specific)
# ==============================================================================
get_copykat_colors <- function() {
  c(
    "aneuploid"   = "#B48283",  # Muted red-brown (Malignant) - echoing FeaturePlot high values
    "diploid"     = "#8F9EAB",  # Foggy blue-gray (Normal) - echoing main theme
    "not.defined" = "#E0E0E0"   # Very light gray (Noise) - blending into background
  )
}

get_standard_theme <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black", face = "plain"),
      axis.title = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12)
    )
}

# ==============================================================================
# Part 1: Malignant cell identification and extraction
# ==============================================================================

#' Identify Malignant Hepatocytes using CopyKAT Data
#' @export
identify_malignant_cells <- function(seurat_obj, pred_data = NULL, cell_type = "Hepatocytes") {

  # 1. Load prediction data
  if (is.null(pred_data)) {
    # Assume copykat_data is available in the environment
    pred_data <- get("copykat_data")
  }

  # 2. Check intersection
  common_cells <- intersect(colnames(seurat_obj), rownames(pred_data))
  if (length(common_cells) == 0) {
    stop("Error: No common cell barcodes found between Seurat object and CopyKAT predictions.")
  }
  message(paste("Matched CopyKAT predictions for", length(common_cells), "cells."))

  # 3. Map metadata
  cnv_vec <- pred_data[colnames(seurat_obj), "copykat_pred"]
  names(cnv_vec) <- colnames(seurat_obj)

  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = cnv_vec, col.name = "cnv_status")

  # 4. Extract hepatocytes
  if (!"SingleR_Labels" %in% colnames(seurat_obj@meta.data)) {
    stop("Error: 'SingleR_Labels' not found in metadata.")
  }

  message(paste("Subsetting for", cell_type, "..."))
  obj_hep <- subset(seurat_obj, subset = SingleR_Labels == cell_type)

  return(obj_hep)
}

# ==============================================================================
# Part 2: Statistical plotting (stacked bar chart)
# ==============================================================================

#' Plot CopyKAT CNV Statistics (Barplots)
#' @export
plot_copykat_stats <- function(hep_obj, group_col = "site") {

  if (!"cnv_status" %in% colnames(hep_obj@meta.data)) stop("Error: 'cnv_status' not found.")

  df_plot <- data.frame(
    site = hep_obj@meta.data[[group_col]],
    cnv_status = hep_obj@meta.data[["cnv_status"]]
  )

  df_clean <- df_plot %>%
    dplyr::filter(!is.na(cnv_status)) %>%
    dplyr::group_by(site, cnv_status) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  # Ensure Aneuploid is in the most prominent position (order in stacked chart, does not affect color mapping)
  df_clean$cnv_status <- factor(df_clean$cnv_status, levels = c("not.defined", "diploid", "aneuploid"))

  colors_copykat <- get_copykat_colors()

  # --- Figure A: Percentage stacked chart ---
  p_percent <- ggplot2::ggplot(df_clean, ggplot2::aes(x = site, y = n, fill = cnv_status)) +
    ggplot2::geom_bar(stat = "identity", position = "fill", color = NA, width = 0.7) +
    ggplot2::scale_fill_manual(values = colors_copykat) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Proportion of Cells (%)") +
    get_standard_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black", vjust = 0.5))

  # --- Figure B: Absolute count stacked chart ---
  p_count <- ggplot2::ggplot(df_clean, ggplot2::aes(x = site, y = n, fill = cnv_status)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", color = NA, width = 0.7) +
    ggplot2::scale_fill_manual(values = colors_copykat) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, sum(df_clean$n) * 1.05)) +
    ggplot2::labs(x = NULL, y = "Number of Hepatocytes") +
    get_standard_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black", vjust = 0.5))

  return(list(percent_plot = p_percent, count_plot = p_count))
}

# ==============================================================================
# Part 3: UMAP comparison plot (CNV vs Site)
# ==============================================================================

#' Plot CopyKAT UMAP Comparison
#' @export
plot_copykat_umap <- function(hep_obj, group_col = "site") {

  # 1. CNV UMAP (malignant vs normal)
  colors_copykat <- get_copykat_colors()

  p1 <- Seurat::DimPlot(hep_obj, reduction = "umap", group.by = "cnv_status",
                        cols = colors_copykat, pt.size = 0.8, order = TRUE) +
    ggplot2::ggtitle("CopyKAT Classification") +
    get_standard_theme() +
    ggplot2::theme(legend.position = "bottom")

  # 2. Site UMAP (tissue origin)
  # Key modification: without specifying cols, default rainbow colors would be used which is jarring.
  # We call get_morandi_colors to generate matching cool-gray color scheme.
  n_sites <- length(unique(hep_obj@meta.data[[group_col]]))
  site_cols <- get_morandi_colors(n_sites)

  p2 <- Seurat::DimPlot(hep_obj, reduction = "umap", group.by = group_col,
                        cols = site_cols, pt.size = 0.8) +
    ggplot2::ggtitle("Tissue Site") +
    get_standard_theme() +
    ggplot2::theme(legend.position = "bottom")

  # 3. Composite
  p_final <- patchwork::wrap_plots(p2, p1, ncol = 2)
  return(p_final)
}
