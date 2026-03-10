#' @importFrom pheatmap pheatmap
#' @importFrom dplyr arrange desc
#' @importFrom stats na.omit
NULL

# ==============================================================================
# Part 1: Core Computation Module (Diagnostic Modeling)
# ==============================================================================

#' Run Multi-Model Diagnostic Feature Selection (RF + SVM + GBM)
#'
#' Builds three diagnostic classifiers (Random Forest, SVM-Radial, GBM) to distinguish
#' Tumor vs Normal samples. Returns aggregated feature importance.
#'
#' @param candidates Vector of gene symbols (e.g., from Prognostic Model).
#' @param expr_mat Expression matrix (Genes x Samples).
#' @param clin_data Clinical data frame (SampleID, Group).
#' @param seed Random seed for reproducibility.
#' @return An object of class `HCCDiagnosticResult`.
#' @export
run_diagnostic_model <- function(candidates, expr_mat = NULL, clin_data = NULL, seed = 123) {

  .check_cran_packages(c("caret", "gbm", "kernlab"))

  # Force load gbm package to prevent caret from failing to find relative.influence
  library(gbm)

  # --- 1. Data Preparation ---
  if(is.null(expr_mat)) expr_mat <- get("TCGA_Expr_Mat")
  if(is.null(clin_data)) clin_data <- get("TCGA_Clin_Data")

  # Extract Tumor and Normal samples
  valid_samples <- intersect(clin_data$SampleID, colnames(expr_mat))
  valid_genes <- intersect(candidates, rownames(expr_mat))

  if(length(valid_genes) < 2) stop("Not enough candidate genes found.")
  message(paste("Starting Diagnostic Workflow | Samples:", length(valid_samples), "| Genes:", length(valid_genes)))

  # Build wide matrix (Samples x Genes)
  dat_expr <- t(expr_mat[valid_genes, valid_samples])
  dat_clin <- clin_data[match(valid_samples, clin_data$SampleID), "Group"]

  # Combine as data.frame
  model_dat <- data.frame(dat_expr, cato = dat_clin)

  # --- 2. Strict Data Cleaning ---
  # 2.1 Remove genes with >30% zeros (consistent with prognostic model)
  gene_cols <- setdiff(colnames(model_dat), "cato")
  zero_fraction <- colMeans(model_dat[, gene_cols] == 0, na.rm = TRUE)
  genes_to_keep <- names(zero_fraction[zero_fraction <= 0.3])

  if(length(gene_cols) - length(genes_to_keep) > 0) {
    message(paste(">>> Preprocessing: Dropped", length(gene_cols) - length(genes_to_keep), "genes (>30% zeros)."))
  }
  model_dat <- model_dat[, c(genes_to_keep, "cato")]

  # 2.2 Factorize labels (Normal/Tumor)
  model_dat$cato <- factor(model_dat$cato, levels = c("Normal", "Tumor"))

  # 2.3 Clean gene names
  colnames(model_dat) <- make.names(colnames(model_dat))

  message(paste(">>> Data Ready. Final Features:", length(genes_to_keep)))

  # --- 3. Set Training Parameters (10-fold CV, 3 repeats) ---
  fitControl <- caret::trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = TRUE,
    verboseIter = FALSE
  )

  # --- 4. Train Three Major Models in Batch ---

  # Model A: Random Forest
  message(">>> Training Random Forest...")
  set.seed(seed)
  fit_rf <- caret::train(cato ~ ., data = model_dat,
                         method = "rf",
                         metric = "ROC",
                         trControl = fitControl,
                         importance = TRUE)

  # Model B: SVM (Radial)
  message(">>> Training SVM (Radial)...")
  set.seed(seed)
  fit_svm <- caret::train(cato ~ ., data = model_dat,
                          method = "svmRadial",
                          metric = "ROC",
                          trControl = fitControl)

  # Model C: GBM
  message(">>> Training GBM...")
  set.seed(seed)
  gbmGrid <- expand.grid(
    interaction.depth = c(1, 3),
    n.trees = c(100, 300),
    shrinkage = c(0.1, 0.01),
    n.minobsinnode = 10
  )
  fit_gbm <- caret::train(cato ~ ., data = model_dat,
                          method = "gbm",
                          metric = "ROC",
                          trControl = fitControl,
                          tuneGrid = gbmGrid,
                          verbose = FALSE)

  # --- 5. Extract Importance ---
  message(">>> Aggregating Feature Importance...")

  imp_list <- list()

  # RF Importance
  tmp_rf <- caret::varImp(fit_rf, scale = TRUE)$importance
  # Take mean in case there are multiple columns
  imp_list[["RF"]] <- data.frame(Gene = rownames(tmp_rf), RF_Score = rowMeans(tmp_rf))

  # SVM Importance
  tmp_svm <- caret::varImp(fit_svm, scale = TRUE)$importance
  imp_list[["SVM"]] <- data.frame(Gene = rownames(tmp_svm), SVM_Score = rowMeans(tmp_svm))

  # GBM Importance
  # library(gbm is now loaded, this should definitely work
  tmp_gbm <- caret::varImp(fit_gbm, scale = TRUE)$importance
  imp_list[["GBM"]] <- data.frame(Gene = rownames(tmp_gbm), GBM_Score = rowMeans(tmp_gbm))

  # Merge All
  final_imp <- Reduce(function(x, y) merge(x, y, by = "Gene"), imp_list)
  final_imp$Average_Score <- rowMeans(final_imp[, -1])
  final_imp <- final_imp %>% dplyr::arrange(dplyr::desc(Average_Score))

  res <- list(
    models = list(rf = fit_rf, svm = fit_svm, gbm = fit_gbm),
    importance = final_imp
  )
  class(res) <- "HCCDiagnosticResult"
  return(res)
}

# ==============================================================================
# Part 2: Visualization Module (Heatmap)
# ==============================================================================

#' Plot Diagnostic Feature Importance Heatmap
#'
#' Visualizes the importance scores of top features across RF, SVM, and GBM models.
#'
#' @param res Object of class `HCCDiagnosticResult`.
#' @param top_n Number of top genes to display (default: 20).
#' @export
plot_diag_heatmap <- function(res, top_n = 20) {

  # 1. Prepare data
  final_imp <- res$importance
  heatmap_dat <- final_imp[, c("RF_Score", "SVM_Score", "GBM_Score")]

  # 2. Beautify column names
  colnames(heatmap_dat) <- c("Random Forest", "SVM", "GBM")
  rownames(heatmap_dat) <- final_imp$Gene

  # 3. Get Top N
  top_dat <- head(heatmap_dat, top_n)

  # 4. Plot
  pheatmap::pheatmap(
    top_dat,
    cluster_rows = FALSE, # Already sorted, no clustering
    cluster_cols = FALSE, # Keep model order
    display_numbers = TRUE,
    number_format = "%.1f",
    color = colorRampPalette(c("white", "firebrick3"))(50),
    main = paste0("Diagnostic Feature Importance (Top ", top_n, ")"),
    fontsize = 12,
    cellwidth = 40,
    cellheight = 20,
    angle_col = 45
  )
}
#plot_diag_heatmap(res_diag,10)
