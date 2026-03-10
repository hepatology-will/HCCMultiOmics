#' @importFrom glmnet cv.glmnet
#' @importFrom randomForestSRC tune.rfsrc rfsrc
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix xgb.importance
#' @importFrom survival Surv
#' @importFrom stats na.omit predict
#' @importFrom grid grid.draw grid.newpage grid.text gpar
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point geom_vline annotate labs theme_classic theme element_blank element_line element_text geom_col coord_flip scale_fill_manual scale_y_continuous geom_line geom_hline geom_segment scale_x_discrete
#' @importFrom dplyr mutate filter arrange pull desc
NULL
#' Perform Intersection Analysis for Target Discovery
#'
#' This function intersects the user-provided gene set with internal TCGA-LIHC DEGs
#' and GeneCards HCC-related genes to prioritize candidates.
#'
#' @param gene_set_obj An object of class `HCCGeneSet`.
#' @param degs_list Vector of strings. Background DEGs (defaults to internal data).
#' @param hcc_list Vector of strings. Background HCC genes (defaults to internal data).
#'
#' @return A list containing:
#' \item{candidates}{Vector of overlapping gene symbols (The Center intersection)}
#' \item{venn_lists}{A named list of the three gene sets used for plotting}
#' \item{counts}{A list of counts for each region}
#' @export
run_intersection <- function(gene_set_obj,
                             degs_list = NULL,
                             hcc_list = NULL) {

  # 1. Check input
  if (!inherits(gene_set_obj, "HCCGeneSet")) {
    stop("Input must be an HCCGeneSet object.")
  }

  # 2. Load background data (use built-in data if not provided by user)
  # Note: In package development, built-in data is usually available as LazyData
  if (is.null(degs_list)) degs_list <- get("TCGA_LIHC_DEGs")
  if (is.null(hcc_list)) hcc_list <- get("HCC_GeneCards_List")

  user_genes <- gene_set_obj$all_genes

  # 3. Calculate intersection genes (Center)
  # Intersection: DEGs & HCC & User
  candidates <- intersect(intersect(degs_list, hcc_list), user_genes)

  # 4. Prepare output
  results <- list(
    candidates = candidates,
    venn_lists = list(
      Deg = degs_list,
      Hcc = hcc_list,
      User = user_genes,
      User_Name = gene_set_obj$name # Save name for plotting
    )
  )

  class(results) <- "HCCDiscoveryResult"
  return(results)
}

#' Plot Discovery Venn Diagram
#'
#' Draws a publication-ready Venn diagram matching the specific style of the HCC study.
#' Colors: DEGs (Blue), HCC (Gold), User/Ferroptosis (Green).
#'
#' @param discovery_result The result object from `run_intersection()`.
#' @return Draws a plot on the current device.
#' @export
plot_venn <- function(discovery_result) {
  .check_cran_packages("VennDiagram")

  if (!inherits(discovery_result, "HCCDiscoveryResult")) {
    stop("Input must be the result from run_intersection().")
  }

  # --- 1. Extract data ---
  lists <- discovery_result$venn_lists

  # Set A: DEGs
  set_deg <- lists$Deg
  # Set B: HCC
  set_hcc <- lists$Hcc
  # Set C: User (e.g., Ferroptosis)
  set_user <- lists$User
  user_name <- if (!is.null(lists$User_Name)) lists$User_Name else "User"

  # --- 2. Dynamically calculate values for Venn diagram ---
  # draw.triple.venn requires area sizes and intersection sizes for each region

  area1 <- length(set_deg)
  area2 <- length(set_hcc)
  area3 <- length(set_user)

  n12 <- length(intersect(set_deg, set_hcc))
  n23 <- length(intersect(set_hcc, set_user))
  n13 <- length(intersect(set_deg, set_user))
  n123 <- length(intersect(intersect(set_deg, set_hcc), set_user))

  # --- 3. Define colors (Strictly following your snippet) ---
  # Deg (Blue), Hcc (Gold), User (Green)
  my_colors <- c("#4575B4", "#D9863E", "#1B9E77")

  # --- 4. Plot ---
  grid::grid.newpage()

  venn_plot <- VennDiagram::draw.triple.venn(
    # --- Data ---
    area1 = area1,
    area2 = area2,
    area3 = area3,
    n12   = n12,
    n23   = n23,
    n13   = n13,
    n123  = n123,

    # --- Labels (Category names) ---
    category = c("Deg", "Hcc", user_name),

    # --- Fill colors (Key step) ---
    fill = my_colors,
    alpha = c(0.5, 0.5, 0.5), # Maintain transparency

    # --- Border settings ---
    col = my_colors, # Border color matches fill
    lwd = 2,         # Border width
    lty = "solid",

    # --- Number style ---
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",

    # --- Category name style ---
    cat.cex = 1.8,
    cat.fontface = "bold",
    cat.col = my_colors, # Category names color auto-match

    # --- Layout adjustment (Strictly following your snippet) ---
    # Ensure Deg is top-left, Hcc is top-right, User is bottom
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.05, 0.05, 0.05),

    # --- Scale control ---
    # Turn off euler.d to prevent circles from disappearing due to magnitude differences
    euler.d = FALSE,
    scaled = FALSE,

    # Prevent writing log files
    ind = FALSE
  )

grid::grid.draw(venn_plot)
}
