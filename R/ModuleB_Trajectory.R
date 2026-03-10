#' @importFrom ggplot2 ggplot aes geom_boxplot geom_path ggtitle scale_color_gradientn labs theme_classic theme element_blank element_line element_text arrow unit annotate guides guide_legend
#' @importFrom patchwork wrap_plots
#' @importFrom stats reorder


NULL

# ==============================================================================
# Helper: Standard theme and Morandi color scheme
# ==============================================================================

# 1. Discrete colors (for Cluster)
get_morandi_colors <- function(n) {
  muted_morandi_12 <- c(
    "#8F9EAB", "#C59993", "#9FAFA3", "#D6C6A8",
    "#A796B0", "#9C9C9C", "#7D8C8E", "#B0A49D",
    "#A3B6C2", "#D1B8B4", "#B4BDA6", "#C8C3BC"
  )
  if (n > length(muted_morandi_12)) colorRampPalette(muted_morandi_12)(n) else muted_morandi_12[1:n]
}

# 2. Continuous gradient colors (for CytoTRACE / Pseudotime)
# Logic: Low (cool gray-blue) -> Middle (almost white) -> High (warm red-brown)
get_morandi_gradient <- function() {
  c("#8F9EAB", "#F0F0F0", "#B48283")
}

get_standard_theme <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black", face = "plain"),
      axis.title = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10, color = "black", face = "plain")
    )
}

# ==============================================================================
# Part 0: Malignant cell reclustering (Harmony + UMAP)
# ==============================================================================

#' Reprocess Malignant Cells (Harmony Integration)
#' @export
reprocess_malignant_cells <- function(seurat_obj, batch_col = "patient", theta = 3, resolution = 0.5) {
  .check_cran_packages("harmony")

  message(">>> Running Normalization and PCA...")
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)

  message(paste(">>> Running Harmony (theta =", theta, ")..."))
  if (!requireNamespace("harmony", quietly = TRUE)) stop("Package 'harmony' is required.")

  seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = batch_col, plot_convergence = FALSE, theta = theta)

  message(">>> Running UMAP and Clustering on Harmony reduction...")
  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = "harmony", dims = 1:10)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:10)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution)

  return(seurat_obj)
}

#' Plot Malignant Subclusters
#' @export
plot_malignant_clusters <- function(seurat_obj) {

  n_clus <- length(unique(seurat_obj$seurat_clusters))
  plot_cols <- get_morandi_colors(n_clus) # Use new color scheme

  p_base <- Seurat::DimPlot(seurat_obj, reduction = "umap", cols = plot_cols,
                            pt.size = 1.2, label = FALSE, repel = TRUE) +
    ggplot2::ggtitle(NULL) +
    get_standard_theme() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))

  p_final <- Seurat::LabelClusters(plot = p_base, id = "ident", color = "black",
                                   size = 5, box = FALSE, fontface = "plain", repel = TRUE)

  return(p_final)
}

# ==============================================================================
# Part 1: Stemness assessment (CytoTRACE)
# ==============================================================================

#' Run CytoTRACE Stemness Analysis
#' @export
run_stemness_analysis <- function(seurat_obj, ncores = 1) {

  if (!requireNamespace("CytoTRACE", quietly = TRUE)) stop("Package 'CytoTRACE' is required.")

  message(">>> Extracting count matrix...")
  mat <- tryCatch(
    as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")),
    error = function(e) as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))
  )

  message(">>> Running CytoTRACE (enableFast=TRUE)...")
  results <- CytoTRACE::CytoTRACE(mat, ncores = ncores, subsamplesize = 1000, enableFast = TRUE)

  seurat_obj$cyto_score <- results$CytoTRACE[colnames(seurat_obj)]

  # 1. FeaturePlot: Use Morandi cool-warm gradient
  colors_gradient <- get_morandi_gradient()

  p1 <- Seurat::FeaturePlot(seurat_obj, features = "cyto_score", reduction = "umap", pt.size = 1.0) +
    ggplot2::scale_color_gradientn(colors = colors_gradient) +
    ggplot2::ggtitle("CytoTRACE Score (Stemness)") +
    get_standard_theme() +
    ggplot2::labs(color = "Stemness")

  # 2. Violin Plot: Use discrete Morandi colors
  n_clus <- length(unique(seurat_obj$seurat_clusters))
  plot_cols <- get_morandi_colors(n_clus)

  p2 <- Seurat::VlnPlot(seurat_obj, features = "cyto_score", group.by = "seurat_clusters",
                        cols = plot_cols, pt.size = 0) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, size = 0.3) +
    ggplot2::ggtitle("Stemness per Cluster") +
    get_standard_theme() +
    ggplot2::labs(x = "Cluster", y = "CytoTRACE Score") +
    ggplot2::theme(legend.position = "none")

  return(list(obj = seurat_obj, plot_umap = p1, plot_vln = p2))
}

# ==============================================================================
# Part 2: Pseudotime trajectory analysis (Slingshot)
# ==============================================================================

#' Run Slingshot Trajectory Analysis
#' @export
run_trajectory_analysis <- function(seurat_obj, start_cluster = "1") {
  .check_bioc_packages(c("slingshot", "SingleCellExperiment"))
  .check_cran_packages("harmony")

  message(paste(">>> Running Slingshot with start.clus =", start_cluster))

  emb <- Seurat::Embeddings(seurat_obj, "umap")
  clu <- seurat_obj$seurat_clusters

  message(">>> Preparing SingleCellExperiment object...")
  counts_mat <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = data.frame(cluster = clu)
  )
  SingleCellExperiment::reducedDims(sce)$UMAP <- emb

  message(">>> Running slingshot...")
  sce <- slingshot::slingshot(sce, clusterLabels = "cluster", reducedDim = "UMAP", start.clus = start_cluster)
  sds <- slingshot::SlingshotDataSet(sce)

  # Background uses CytoTRACE or Cluster colors; keeping CytoTRACE gradient as background is informative
  colors_gradient <- get_morandi_gradient()

  p_base <- Seurat::FeaturePlot(seurat_obj, features = "cyto_score", reduction = "umap", pt.size = 1.0) +
    ggplot2::scale_color_gradientn(colors = colors_gradient) +
    ggplot2::ggtitle(NULL) +
    get_standard_theme() +
    ggplot2::labs(color = "Stemness")

  num_curves <- length(sds@curves)
  p_all <- p_base

  # Trajectory lines remain black - academic standard and clearest on low-saturation background
  for(i in 1:num_curves){
    curve_coords <- sds@curves[[i]]$s
    curve_df <- data.frame(UMAP_1 = curve_coords[,1], UMAP_2 = curve_coords[,2])
    p_all <- p_all +
      ggplot2::geom_path(data = curve_df, ggplot2::aes(x = UMAP_1, y = UMAP_2),
                         color = "black", size = 1.0,
                         arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                         inherit.aes = FALSE) +
      ggplot2::annotate("text", x = tail(curve_df$UMAP_1, 1), y = tail(curve_df$UMAP_2, 1),
                        label = paste0("C", i), size = 6, fontface = "bold")
  }

  message(paste(">>> Found", num_curves, "trajectories. Check plot to select one."))
  return(list(sds = sds, plot_all_curves = p_all))
}

# ==============================================================================
# Part 3: Extract specific trajectory (Subset Curve)
# ==============================================================================

#' Subset Cells on a Specific Trajectory
#' @export
subset_trajectory <- function(seurat_obj, sds, curve_id = 1) {

  message(paste(">>> Subsetting Curve", curve_id, "..."))

  weights <- slingshot::slingCurveWeights(sds)
  pt_matrix <- slingshot::slingPseudotime(sds)

  if (curve_id > ncol(weights)) stop("Error: curve_id exceeds number of curves.")

  cells_on_curve <- rownames(weights)[weights[, curve_id] > 0 & !is.na(weights[, curve_id])]
  message(paste(">>> Cells on curve:", length(cells_on_curve)))

  obj_sub <- subset(seurat_obj, cells = cells_on_curve)
  obj_sub$pseudotime <- pt_matrix[colnames(obj_sub), curve_id]
  obj_sub <- obj_sub[, order(obj_sub$pseudotime)]

  # Pseudotime is also a continuous variable, use same Morandi gradient
  colors_gradient <- get_morandi_gradient()

  curve_coords <- sds@curves[[curve_id]]$s
  curve_df <- data.frame(UMAP_1 = curve_coords[,1], UMAP_2 = curve_coords[,2])

  p_final <- Seurat::FeaturePlot(obj_sub, features = "pseudotime", reduction = "umap", pt.size = 1.0) +
    ggplot2::scale_color_gradientn(colors = colors_gradient, na.value = "grey90") +
    ggplot2::ggtitle(NULL) +
    get_standard_theme() +
    ggplot2::labs(color = "Pseudotime") +
    ggplot2::geom_path(data = curve_df, ggplot2::aes(x = UMAP_1, y = UMAP_2),
                       color = "black", size = 1.2,
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
                       inherit.aes = FALSE)

  return(list(obj = obj_sub, plot_final = p_final))
}
