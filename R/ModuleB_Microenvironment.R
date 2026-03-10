# ==============================================================================
# Module B: Automated Microenvironment Analysis (Final Corrected)
# ==============================================================================
# 1. Auto-Mapping: Gene -> Pathway.
# 2. Hybrid Plotting: Draws Chord plot immediately; Returns ggplot objects.
# 3. Style: Morandi Aesthetic.
# ==============================================================================

#' @import ggplot2
#' @import dplyr
#' @importFrom methods is
#' @importFrom graphics par title
NULL

# ==============================================================================
# Main Function: One-Click Automation
# ==============================================================================

#' Auto-Analyze Ligand-Receptor Interaction
#'
#' @description
#' The Master Function.
#' 1. Identifies the pathway for the input gene.
#' 2. Draws the Network Structure (Chord Diagram) IMMEDIATELY to the plot pane.
#' 3. Returns a list containing the Bubble Plot and Centrality Plot (ggplot2).
#'
#' @param object Seurat or CellChat object.
#' @param gene Character. Target gene (e.g., "ANGPT2").
#' @param group_col Character. Metadata column. Default "SingleR_Labels".
#' @param db_use Character. "human" or "mouse".
#'
#' @return A list with: $bubble (ggplot), $centrality (ggplot), $object (CellChat).
#' @export
auto_analyze_ligand <- function(object, gene = "ANGPT2", group_col = "SingleR_Labels", db_use = "human") {

  # --- 1. Calculation & Mapping ---
  cellchat <- run_cellchat_core(object, group_col, db_use)

  message(paste0("[Auto] Mapping gene '", gene, "' to pathway..."))
  pathway <- get_pathway_from_gene(cellchat, gene)

  if (is.null(pathway)) stop(paste0("Gene '", gene, "' not found in CellChat DB."))
  if (!pathway %in% cellchat@netP$pathways) stop(paste0("Pathway '", pathway, "' is not significant in this dataset."))

  message(paste0("[Auto] Target Pathway: ", pathway))

  # --- 2. Compute Centrality (Essential for Bar Plot) ---
  # Important: The slot.name must be 'netP' for pathway level analysis
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  # --- 3. EXECUTE PLOTS ---

  # [Action 1] Draw Base Plot (Chord) IMMEDIATELY
  # Base plots cannot be stored in a list easily without weird tricks.
  # So we draw it as a "side effect" so the user sees it instantly.
  message("[Vis] Drawing Network Structure (Chord Diagram) to plot pane...")
  tryCatch({
    graphics::par(mfrow = c(1, 1), xpd = TRUE) # Reset par settings
    CellChat::netVisual_aggregate(
      cellchat,
      signaling = pathway,
      layout = "circle",
      pt.title = 20,
      title.space = 0.2,
      remove.isolate = TRUE
    )
    graphics::title(paste0(pathway, " Network Structure"), line = -1)
  }, error = function(e) message("[Warn] Chord plot failed: ", e$message))

  # [Action 2] Generate GGPLOTs (to return)
  message("[Vis] Generating Bubble and Centrality plots...")
  p_bubble <- plot_morandi_bubble(cellchat, pathway)
  p_bar <- plot_morandi_centrality(cellchat, pathway)

  # --- 4. Return ---
  message("[Done] Chord diagram is on screen. Bubble and Bar plots are in the returned list.")
  return(list(
    bubble = p_bubble,      # Bubble plot (ggplot)
    centrality = p_bar,    # Bar plot (ggplot)
    object = cellchat       # Object (for downstream analysis)
  ))
}

# ==============================================================================
# Helper Functions (Internal)
# ==============================================================================

# Internal: Core Calculation Logic
run_cellchat_core <- function(object, group_col, db_use) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required but not installed.\n",
         "Please install it manually from GitHub:\n",
         "  devtools::install_github('sqjin/CellChat')")
  }

  if (methods::is(object, "CellChat")) return(object)

  message("[Core] Initializing CellChat (this may take time)...")
  data_input <- Seurat::GetAssayData(object, assay = "RNA", layer = "data")
  cc <- CellChat::createCellChat(data_input, meta = object@meta.data, group.by = group_col)
  cc@DB <- if(db_use == "human") CellChat::CellChatDB.human else CellChat::CellChatDB.mouse

  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "triMean")
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  return(cc)
}

# Internal: Map Gene -> Pathway
get_pathway_from_gene <- function(cc, gene) {
  DB <- cc@DB$interaction
  hit <- DB[which(DB$ligand == gene | DB$gene_name == gene), ]
  if (nrow(hit) == 0) hit <- DB[grep(gene, DB$interaction_name)[1], ]
  if (nrow(hit) == 0 || is.na(hit$pathway_name[1])) return(NULL)
  return(unique(hit$pathway_name)[1])
}

# Internal: Bubble Plot (Morandi Style)
plot_morandi_bubble <- function(cc, signaling) {
  df <- tryCatch({
    CellChat::netVisual_bubble(cc, signaling = signaling, remove.isolate = TRUE, return.data = TRUE)$communication
  }, error = function(e) NULL)

  if (is.null(df)) return(NULL)

  df <- df %>%
    dplyr::filter(prob > 0) %>%
    dplyr::mutate(log_p = -log10(pval + 1e-5), pair = paste0(source, " -> ", target))

  # Morandi Blue-Red Palette
  cols <- c("#5B7A8C", "#F5F5F5", "#C85C5C")

  ggplot2::ggplot(df, ggplot2::aes(x = pair, y = interaction_name_2)) +
    ggplot2::geom_point(ggplot2::aes(size = log_p, color = prob)) +
    ggplot2::scale_color_gradientn(colors = cols, name = "Prob.") +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "-log10(P)") +
    ggplot2::labs(title = paste0(signaling, ": Interactions"), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      panel.border = ggplot2::element_rect(color = "grey80", fill = NA)
    )
}

# Internal: Centrality Bar Plot (Morandi Style)
plot_morandi_centrality <- function(cc, signaling) {
  # Extract centrality scores safely
  cent <- cc@netP$centr[[signaling]]
  if (is.null(cent)) return(NULL)

  df_send <- data.frame(Cell = names(cent$outdeg), Score = cent$outdeg, Role = "Sender")
  df_recv <- data.frame(Cell = names(cent$indeg), Score = cent$indeg, Role = "Receiver")
  df_all <- rbind(df_send, df_recv)

  # Morandi Palette: Red for Sender, Blue for Receiver
  cols <- c("Sender" = "#C85C5C", "Receiver" = "#5B7A8C")

  ggplot2::ggplot(df_all, ggplot2::aes(x = stats::reorder(Cell, Score), y = Score, fill = Role)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Role, scales = "free_x") +
    ggplot2::labs(title = paste0(signaling, ": Network Roles"), x = NULL, y = "Centrality Score") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "none")
}
# ==============================================================================
# ==============================================================================
# MODULE C: Stromal Functional State Decoding (Universal Mechanism)
# ==============================================================================
# Focus: Correlating Target Gene with Global Hallmark Functions.
# Logic: Gene High -> Function Score High = Mechanism confirmed.
# Style: Morandi Aesthetic.
# ==============================================================================

#' @importFrom stats cor.test reorder
NULL

# ==============================================================================
# Helper: Get All Hallmark Pathways (Comprehensive)
# ==============================================================================
get_all_hallmarks <- function(species = "Homo sapiens") {
  .check_cran_packages("msigdbr")

  # Fetch ALL Hallmark pathways from msigdbr
  # Ensure 'msigdbr' is installed: install.packages("msigdbr")
  m_df <- msigdbr::msigdbr(species = species, category = "H")
  path_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
  return(path_list)
}

# ==============================================================================
# Main Function: Auto-Analyze Functional State
# ==============================================================================

#' Analyze Stromal Functional Mechanism (Comprehensive)
#'
#' @description
#' Decodes the mechanism of a stromal gene by correlating it with ALL Hallmark pathways.
#' Automatically highlights the top positively and negatively correlated functions.
#'
#' @param seurat_obj A Seurat object.
#' @param gene Character. The target gene (e.g., "ANGPT2").
#' @param cell_type Character. The specific stromal cell type (e.g., "Endothelial_cells").
#' @param group_col Character. Metadata column for cell types.
#' @param species Character. "Homo sapiens" or "Mus musculus".
#' @param top_n Integer. How many top positive/negative pathways to show? Default 10.
#'
#' @return A list containing the Lollipop Plot (ggplot) and Correlation Data.
#' @export
auto_stromal_function <- function(seurat_obj, gene, cell_type, group_col = "SingleR_Labels",
                                  species = "Homo sapiens", top_n = 10) {

  message(paste0("[Module C] Decoding functional state for ", gene, " in ", cell_type, "..."))

  # --- 1. Subset Data ---
  # Validate cell type existence
  if (!cell_type %in% seurat_obj@meta.data[[group_col]]) {
    stop(paste0("Cell type '", cell_type, "' not found in column '", group_col, "'."))
  }

  target_obj <- subset(seurat_obj, subset = !!sym(group_col) == cell_type)

  if (!gene %in% rownames(target_obj)) stop(paste0("Gene '", gene, "' not found in dataset."))

  # --- 2. Calculate Functional Scores (All Hallmarks) ---
  message("[Module C] Scoring ALL Hallmark pathways (this may take a moment)...")
  pathways <- get_all_hallmarks(species)

  # Clean names (remove HALLMARK_ prefix for cleaner plots)
  names(pathways) <- gsub("HALLMARK_", "", names(pathways))

  # AddModuleScore
  # Seurat adds scores as "HallmarkScore_1", "HallmarkScore_2"... we need to track names
  target_obj <- Seurat::AddModuleScore(target_obj, features = pathways, name = "HallmarkScore_")

  # --- 3. Correlation Analysis ---
  gene_expr <- Seurat::GetAssayData(target_obj, layer = "data")[gene, ]

  results <- data.frame(Pathway = names(pathways), R = NA, P = NA)

  # Loop through all pathways to calculate correlation
  for (i in seq_along(pathways)) {
    score_col <- paste0("HallmarkScore_", i)
    path_score <- target_obj@meta.data[[score_col]]

    # Spearman correlation (Robust to outliers)
    tryCatch({
      test <- stats::cor.test(gene_expr, path_score, method = "spearman")
      results$R[i] <- test$estimate
      results$P[i] <- test$p.value
    }, error = function(e) {
      results$R[i] <- 0
      results$P[i] <- 1
    })
  }

  # --- 4. Filter Top Pathways (Smart Selection) ---
  # Select Top N Positive and Top N Negative for visualization
  top_pos <- results %>% dplyr::arrange(dplyr::desc(R)) %>% head(top_n)
  top_neg <- results %>% dplyr::arrange(R) %>% head(top_n)

  plot_data <- rbind(top_pos, top_neg) %>%
    dplyr::mutate(LogP = -log10(P + 1e-10)) %>%
    dplyr::arrange(R) # Sort for plotting order

  # --- 5. Morandi Visualization (Lollipop Plot) ---
  message("[Vis] Generating Functional Association Plot...")

  # Morandi Colors: Red for Pos, Blue for Neg
  cols <- c("#B69FD6", "#C8B6A5")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(Pathway, R), y = R, color = R > 0)) +
    ggplot2::geom_segment(ggplot2::aes(xend = Pathway, yend = 0), size = 1.2) +
    ggplot2::geom_point(ggplot2::aes(size = LogP), alpha = 1) +

    ggplot2::scale_color_manual(values = cols, labels = c("Negative", "Positive"), name = "Correlation") +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "-log10(P)") +

    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +

    ggplot2::labs(
      title = paste0("Mechanism: ", gene, " in ", cell_type),
      subtitle = paste0("Top ", top_n*2, " Associated Hallmark Pathways"),
      x = NULL, y = "Spearman Correlation (R)"
    ) +

    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "grey80", fill = NA),
      axis.text.y = ggplot2::element_text(size = 10, color = "#333333")
    )

  return(list(plot = p, data = results))
}
# ==============================================================================
# MODULE D: Immune Landscape Analysis (The "Gatekeeper" Panorama)
# ==============================================================================
# Focus: Global correlation between a Stromal Gene and ALL Immune Cell types.
# Output 1: Landscape Lollipop Plot (Overview of Recruitment/Exclusion).
# Output 2: Faceted Scatter Plot (Detailed trends).
# ==============================================================================

#' @import ggplot2
#' @import dplyr
#' @importFrom stats cor.test reorder
#' @importFrom reshape2 melt
NULL

#' Auto-Analyze Immune Landscape
#'
#' @description
#' Calculates the correlation between a stromal gene and the abundance of ALL other cell types.
#' Generates a "Landscape" plot showing who is recruited (Positive) vs excluded (Negative).
#'
#' @param seurat_obj A Seurat object.
#' @param gene Character. The stromal gene (e.g., "ANGPT2").
#' @param stromal_class Character. The source cell type (e.g., "Endothelial_cells").
#' @param group_col Character. Metadata column for cell types.
#' @param sample_col Character. Metadata column for Sample IDs.
#' @param plot_type Character. "all" (return list of plots), "summary" (lollipop only), "scatter" (grid only).
#'
#' @return A list containing the Landscape Plot, Scatter Grid, and Data.
#' @export
auto_immune_landscape <- function(seurat_obj, gene, stromal_class,
                                  group_col = "SingleR_Labels", sample_col = "orig.ident",
                                  plot_type = "all") {

  message(paste0("[Module D] Mapping Immune Landscape for ", gene, " in ", stromal_class, "..."))

  # --- 1. Data Prep: Calculate Proportions for ALL cell types ---
  meta <- seurat_obj@meta.data

  # Total cells per sample (Force numeric)
  total_counts <- as.numeric(table(meta[[sample_col]]))
  names(total_counts) <- names(table(meta[[sample_col]]))

  # Identify all cell types except the stromal one
  all_cells <- unique(meta[[group_col]])
  target_cells <- setdiff(all_cells, stromal_class)

  # Create a matrix to store proportions
  # Rows = Samples, Cols = Cell Types
  prop_mat <- matrix(NA, nrow = length(total_counts), ncol = length(target_cells))
  rownames(prop_mat) <- names(total_counts)
  colnames(prop_mat) <- target_cells

  for (cell in target_cells) {
    counts_raw <- table(meta[meta[[group_col]] == cell, sample_col])
    # Match samples and fill 0s
    counts_vec <- setNames(rep(0, length(total_counts)), names(total_counts))
    counts_vec[names(counts_raw)] <- as.numeric(counts_raw)

    # Calculate %
    prop_mat[, cell] <- (counts_vec / total_counts) * 100
  }

  # --- 2. Data Prep: Calculate Stromal Gene Expression ---
  stromal_obj <- subset(seurat_obj, subset = !!sym(group_col) == stromal_class)

  # Handle Seurat v5 AverageExpression return type
  avg_res <- Seurat::AverageExpression(stromal_obj, features = gene, group.by = sample_col, assay = "RNA")
  avg_mat <- if(is.list(avg_res) && "RNA" %in% names(avg_res)) avg_res$RNA else avg_res

  # Align samples
  common_samples <- intersect(rownames(prop_mat), colnames(avg_mat))
  if (length(common_samples) < 3) stop("Not enough common samples for correlation.")

  expr_vec <- as.numeric(avg_mat[gene, common_samples])
  prop_mat <- prop_mat[common_samples, ]

  # --- 3. Correlation Loop ---
  results <- data.frame(CellType = colnames(prop_mat), R = NA, P = NA)

  for (i in 1:nrow(results)) {
    ct <- results$CellType[i]
    test <- stats::cor.test(expr_vec, prop_mat[, ct], method = "spearman")
    results$R[i] <- test$estimate
    results$P[i] <- test$p.value
  }

  results <- results %>%
    dplyr::mutate(
      LogP = -log10(P + 1e-10),
      Direction = ifelse(R > 0, "Recruitment", "Exclusion"),
      Significance = ifelse(P < 0.05, "Significant", "NS")
    ) %>%
    dplyr::arrange(R) # Sort by R for plotting

  # --- 4. Plot A: The Landscape (Lollipop) ---
  # This is the "Publication Quality" Summary

  cols <- c("Exclusion" = "#86C495", "Recruitment" = "#F2B37B")

  p_landscape <- ggplot2::ggplot(results, ggplot2::aes(x = stats::reorder(CellType, R), y = R, color = Direction)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
    ggplot2::geom_segment(ggplot2::aes(xend = CellType, yend = 0), size = 1) +
    # Use fill for the point to make it look like a bubble
    ggplot2::geom_point(ggplot2::aes(size = LogP), alpha = 1) +

    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "-log10(P)") +
    ggplot2::coord_flip() +

    ggplot2::labs(
      title = paste0("Immune Landscape: ", gene),
      subtitle = paste0("Stromal-Immune Crosstalk in ", stromal_class),
      x = NULL, y = "Spearman Correlation (R)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right"
    )

  # --- 5. Plot B: The Detail Grid (Faceted Scatter) ---
  # We melt the data to make it compatible with facet_wrap

  # Construct a long dataframe for plotting
  plot_df <- data.frame(Sample = common_samples, GeneExpr = expr_vec)
  plot_df <- cbind(plot_df, prop_mat)
  long_df <- reshape2::melt(plot_df, id.vars = c("Sample", "GeneExpr"),
                            variable.name = "CellType", value.name = "Proportion")

  # Merge with results to get P-values for labels
  long_df <- merge(long_df, results[, c("CellType", "R", "P", "Direction")], by = "CellType")

  # Optional: Filter to show only cell types with |R| > 0.1 or P < 0.1 to avoid clutter?
  # For now, let's keep all but order them by R
  long_df$CellType <- factor(long_df$CellType, levels = rev(results$CellType))

  # Create the label text for each facet
  long_df$Label <- paste0("R=", round(long_df$R, 2), "\np=", formatC(long_df$P, digits = 1, format = "g"))

  p_grid <- ggplot2::ggplot(long_df, ggplot2::aes(x = GeneExpr, y = Proportion)) +
    ggplot2::geom_smooth(method = "lm", ggplot2::aes(color = Direction, fill = Direction), alpha = 0.1) +
    ggplot2::geom_point(ggplot2::aes(color = Direction), size = 1.5, alpha = 0.7) +

    # Facet Wrap with Free Scales is key for "Not Ugly"
    ggplot2::facet_wrap(~CellType, scales = "free_y", ncol = 4) +

    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +

    # Add stats to each plot
    ggplot2::geom_text(
      data = unique(long_df[, c("CellType", "Label", "GeneExpr", "Proportion")]),
      ggplot2::aes(label = Label),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
      size = 3, fontface = "italic", color = "black", inherit.aes = FALSE
    ) +

    ggplot2::labs(
      x = paste0(gene, " Expression"),
      y = "Cell Proportion (%)",
      title = "Correlations by Cell Type"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none" # Color is self-explanatory
    )

  return(list(landscape = p_landscape, scatter_grid = p_grid, data = results))
}
