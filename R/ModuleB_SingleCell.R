#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual scale_y_continuous labs theme_classic theme element_blank element_line element_text geom_col ggtitle
#' @importFrom dplyr group_by summarise %>%
#' @importFrom patchwork wrap_plots
#' @importFrom stats reorder
NULL

# ==============================================================================
# Helper: Define ultra-low saturation Morandi color scheme (Ultra-Low Saturation Morandi)
# ==============================================================================
get_morandi_colors <- function(n) {
  # This set of colors all have heavy gray tones, visually very unified and restrained
  muted_morandi_12 <- c(
    "#8F9EAB",  # Foggy blue-gray (main color: subdued)
    "#C59993",  # Dried rose (contrast: extremely soft)
    "#9FAFA3",  # Sage green (cool-toned gray-green)
    "#D6C6A8",  # Linen gray-yellow (only warm tone, but not bright)
    "#A796B0",  # Lavender gray (elegant purple tone)
    "#9C9C9C",  # Neutral warm gray (transition)
    "#7D8C8E",  # Deep sea algae gray (increase contrast)
    "#B0A49D",  # Brownish gray
    "#A3B6C2",  # Light steel blue
    "#D1B8B4",  # Light lotus pink
    "#B4BDA6",  # Light moss green
    "#C8C3BC"   # Even close to white gray
  )

  if (n > length(muted_morandi_12)) {
    return(colorRampPalette(muted_morandi_12)(n))
  } else {
    return(muted_morandi_12[1:n])
  }
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
# Part 1: Single cell atlas display (Landscape)
# ==============================================================================

#' Plot Single Cell Landscape (UMAP)
#' @export
plot_sc_landscape <- function(seurat_obj) {
  if (!"umap" %in% names(seurat_obj@reductions)) stop("Error: No 'umap' reduction.")
  if (!"SingleR_Labels" %in% colnames(seurat_obj@meta.data)) stop("Error: No 'SingleR_Labels'.")

  n_groups <- length(unique(seurat_obj$SingleR_Labels))
  plot_cols <- get_morandi_colors(n_groups)

  p <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_Labels",
                       cols = plot_cols, pt.size = 0.8, label = TRUE, label.size = 5,
                       label.box = FALSE, repel = TRUE) +
    ggplot2::ggtitle(NULL) +
    get_standard_theme() +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 1, color = "grey30"),
                   axis.text = ggplot2::element_text(color = "grey30"),
                   axis.title = ggplot2::element_text(color = "grey30"))
  return(p)
}

# ==============================================================================
# Part 2: Core gene expression localization (Feature + Violin)
# ==============================================================================

#' Plot Target Gene Expression
#' @note FeaturePlot colors changed to: very light gray (#F0F0F0) -> low saturation red-brown (#B48283)
#' @export
plot_sc_gene <- function(seurat_obj, gene = "EZH2", color_low = "#F0F0F0", color_high = "#B48283") {
  if (!gene %in% rownames(seurat_obj)) stop(paste("Gene", gene, "not found."))

  # 1. FeaturePlot: Colors more muted, remove harsh red
  p1 <- Seurat::FeaturePlot(seurat_obj, features = gene, cols = c(color_low, color_high), order = TRUE) +
    ggplot2::labs(title = paste0(gene, " Distribution")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # 2. VlnPlot: Also use new low saturation color scheme
  n_groups <- length(unique(seurat_obj$SingleR_Labels))
  vln_cols <- get_morandi_colors(n_groups)

  p2 <- Seurat::VlnPlot(seurat_obj, features = gene, group.by = "SingleR_Labels",
                        pt.size = 0, cols = vln_cols) +
    ggplot2::labs(title = paste0(gene, " Expression"), x = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ==============================================================================
# Part 3: Gene cell type statistics (Barplot)
# ==============================================================================

#' Plot Gene Statistics per Cell Type
#' @export
plot_sc_gene_stat <- function(seurat_obj, gene = "EZH2", plot_type = 2) {
  if (!gene %in% rownames(seurat_obj)) stop(paste("Gene", gene, "not found."))

  n_groups <- length(unique(seurat_obj$SingleR_Labels))
  plot_cols <- get_morandi_colors(n_groups)

  expr_data <- Seurat::GetAssayData(seurat_obj, layer = "data")[gene, ]
  df_calc <- data.frame(celltype = seurat_obj$SingleR_Labels, expression = expr_data)

  if (plot_type == 2) {
    stat_df <- df_calc %>% dplyr::group_by(celltype) %>% dplyr::summarise(Value = mean(expression > 0) * 100)
    y_label <- paste0(gene, "+ Cells Proportion (%)")
  } else {
    stat_df <- df_calc %>% dplyr::group_by(celltype) %>% dplyr::summarise(Value = mean(expression))
    y_label <- paste0("Average Expression of ", gene)
  }

  # Find cell type with highest expression
  top_cell <- as.character(stat_df$celltype[which.max(stat_df$Value)])

  p <- ggplot2::ggplot(stat_df, ggplot2::aes(x = stats::reorder(celltype, -Value), y = Value, fill = celltype)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = plot_cols) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(stat_df$Value) * 1.1)) +
    ggplot2::labs(x = NULL, y = y_label) +
    get_standard_theme() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))

  return(list(plot = p, top_cell = top_cell))
}
# ==============================================================================
# Part 4: Site distribution and gene statistics for specific cell type (Subset Analysis)
# ==============================================================================

#' Plot Gene Statistics by Site in Specific Cell Type
#' @param plot_type Integer. 1 for Average Expression, 2 for Proportion (%). Default is 2.
#' @export
plot_sc_gene_by_site <- function(seurat_obj, gene = "EZH2", cell_type = "Hepatocytes",
                                 group_col = "site", plot_type = 2) {

  # 1. Validate input
  if (!gene %in% rownames(seurat_obj)) stop(paste("Gene", gene, "not found."))
  if (!group_col %in% colnames(seurat_obj@meta.data)) stop(paste("Column", group_col, "not found in meta.data."))
  if (!cell_type %in% seurat_obj$SingleR_Labels) stop(paste("Cell type", cell_type, "not found."))

  # 2. Extract subset of specific cell type
  message(paste("Subsetting", cell_type, "..."))
  obj_sub <- subset(seurat_obj, subset = SingleR_Labels == cell_type)

  # 3. Draw UMAP (left plot, show distribution of this cell type)
  sites <- unique(obj_sub@meta.data[[group_col]])
  plot_cols <- get_morandi_colors(length(sites))

  p_umap <- Seurat::DimPlot(obj_sub, reduction = "umap", group.by = group_col,
                            cols = plot_cols, pt.size = 1.0, label = FALSE, repel = TRUE) +
    ggplot2::ggtitle(paste0(cell_type, " Distribution")) +
    get_standard_theme() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # 4. Prepare calculation data
  expr_vals <- Seurat::GetAssayData(obj_sub, layer = "data")[gene, ]
  group_vals <- obj_sub@meta.data[[group_col]]
  df_calc <- data.frame(group = group_vals, expression = expr_vals)

  # 5. Calculate statistics based on plot_type (consistent with Part 3 logic)
  if (plot_type == 2) {
    # Mode 2: Calculate expression proportion
    stat_df <- df_calc %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(Value = mean(expression > 0) * 100)

    y_label <- paste0(gene, "+ ", cell_type, " Proportion (%)")

  } else {
    # Mode 1: Calculate average expression
    stat_df <- df_calc %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(Value = mean(expression))

    y_label <- paste0("Average Expression of ", gene, " in ", cell_type)
  }

  # 6. Draw bar plot (right plot)
  p_bar <- ggplot2::ggplot(stat_df, ggplot2::aes(x = group, y = Value, fill = group)) +
    ggplot2::geom_bar(stat = "identity", color = NA, width = 0.6) +
    ggplot2::scale_fill_manual(values = plot_cols) +
    # Dynamically adjust Y-axis upper limit, leave 15% space to prevent label overflow
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(stat_df$Value) * 1.15)) +
    ggplot2::labs(x = NULL, y = y_label) +
    get_standard_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(vjust = 0.5, color = "black")
    )

  return(list(umap = p_umap, barplot = p_bar))
}
