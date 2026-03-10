#' @importFrom glmnet cv.glmnet
#' @importFrom randomForestSRC tune.rfsrc rfsrc impute
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix xgb.importance
#' @importFrom survival Surv
#' @importFrom stats na.omit predict coef reorder
#' @importFrom grid grid.draw grid.newpage grid.text gpar
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point geom_vline annotate labs theme_classic theme_bw theme element_blank element_line element_text geom_col coord_flip scale_fill_manual scale_y_continuous geom_line geom_hline geom_segment scale_x_discrete
#' @importFrom dplyr mutate filter arrange pull desc
NULL

# ==============================================================================
# Part 1: Core calculation module
# ==============================================================================

#' Run Ultimate Rigorous Prognostic Screening
#'
#' @param candidates Vector of gene symbols.
#' @param expr_mat Expression matrix.
#' @param clin_data Clinical data frame.
#' @param seed Random seed for reproducibility (Default: 123). Change this if result varies!
#' @return An object of class `HCCPrognosticResult`.
#' @export
run_prognostic_model <- function(candidates, expr_mat = NULL, clin_data = NULL, seed = 123) {

  # --- 1. Data preparation ---
  if(is.null(expr_mat)) expr_mat <- get("TCGA_Expr_Mat")
  if(is.null(clin_data)) clin_data <- get("TCGA_Clin_Data")

  # Filter Tumor samples
  tumor_clin <- subset(clin_data, Group == "Tumor" & !is.na(OS.time) & !is.na(OS))
  valid_samples <- intersect(tumor_clin$SampleID, colnames(expr_mat))
  valid_genes <- intersect(candidates, rownames(expr_mat))

  if(length(valid_genes) < 2) stop("Not enough candidate genes found.")
  message(paste("Starting ML Workflow | Seed:", seed, "| Samples:", length(valid_samples), "| Initial Genes:", length(valid_genes)))

  # Build wide matrix
  dat_expr <- t(expr_mat[valid_genes, valid_samples])
  dat_clin <- tumor_clin[match(valid_samples, tumor_clin$SampleID), c("OS.time", "OS")]
  dat <- data.frame(dat_expr, OS = dat_clin$OS, OS.time = dat_clin$OS.time)
  colnames(dat) <- make.names(colnames(dat))

  # ============================================================================
  # Part 0: Data preprocessing
  # ============================================================================

  # 1. Remove low expression
  gene_cols <- setdiff(colnames(dat), c("OS", "OS.time"))
  zero_fraction <- colMeans(dat[, gene_cols] == 0, na.rm = TRUE)
  genes_to_keep <- names(zero_fraction[zero_fraction <= 0.3]) # 30% threshold
  dat <- dat[, c(genes_to_keep, "OS", "OS.time")]
  gene_cols <- genes_to_keep

  # 2. Random Forest imputation (affected by Seed)
  if(any(is.na(dat))) {
    message(">>> Preprocessing: Missing values detected. Running RF Imputation...")
    set.seed(seed) # Set seed
    dat <- randomForestSRC::impute(Surv(OS.time, OS) ~ ., data = dat, ntree = 100, nimpute = 1, splitrule = "random", do.trace = FALSE)
    message(">>> Imputation complete.")
  }

  x_full <- as.matrix(dat[, gene_cols])
  y_full <- survival::Surv(dat$OS.time, dat$OS)

  message(paste(">>> Data ready. Final Gene Count:", length(gene_cols)))

  # ============================================================================
  # Part 1: LASSO Cox
  # ============================================================================
  message(">>> Running LASSO Cox...")
  set.seed(seed)
  cv_lasso <- glmnet::cv.glmnet(x_full, y_full, family = "cox", type.measure = "C", nfolds = 10)
  coef_lasso <- stats::coef(cv_lasso, s = "lambda.min")
  lasso_genes <- rownames(coef_lasso)[which(as.numeric(coef_lasso) != 0)]

  df_lasso_plot <- data.frame(
    log_lambda = log(cv_lasso$lambda),
    c_index    = cv_lasso$cvm,
    cv_up      = cv_lasso$cvup,
    cv_lo      = cv_lasso$cvlo,
    n_genes    = cv_lasso$nzero
  )

  # ============================================================================
  # Part 2: RSF
  # ============================================================================
  message(">>> Running RSF Tuning...")
  set.seed(seed)
  tune_rsf <- randomForestSRC::tune.rfsrc(Surv(OS.time, OS) ~ ., data = dat,
                                          ntree = 500, stepFactor = 2, improve = 1e-3,
                                          trace = FALSE, doBest = FALSE)
  best_mtry <- tune_rsf$optimal["mtry"]
  best_nodesize <- tune_rsf$optimal["nodesize"]

  message(">>> Building Final RSF Model...")
  set.seed(seed)
  rf_final <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = dat,
                                     ntree = 1000, mtry = best_mtry, nodesize = best_nodesize,
                                     importance = "perm", splitrule = "logrank", block.size = 1)

  vimp_scores <- rf_final$importance
  rsf_genes <- names(sort(vimp_scores[vimp_scores > 0], decreasing = TRUE))
  if(length(rsf_genes) > 30) rsf_genes <- rsf_genes[1:30]

  err_vals <- rf_final$err.rate
  if(is.matrix(err_vals)) err_vals <- err_vals[, 1]
  df_rsf_err <- data.frame(Trees = 1:length(err_vals), Error = as.numeric(err_vals))

  # ============================================================================
  # Part 3: XGBoost (Robust + nthread control)
  # ============================================================================
  message(">>> Running XGBoost CV...")
  y_xgb <- dat$OS.time * ifelse(dat$OS == 1, 1, -1)
  dall <- xgboost::xgb.DMatrix(data = x_full, label = y_xgb)

  # Key modification: add nthread = 1 to ensure maximum reproducibility
  params <- list(
    booster = "gbtree",
    objective = "survival:cox",
    eta = 0.01,
    gamma = 0.1,
    max_depth = 5,
    min_child_weight = 1,
    subsample = 0.7,
    colsample_bytree = 0.7,
    nthread = 1 # Force single thread, reduce concurrent randomness
  )

  set.seed(seed) # Ensure XGBoost internal random numbers are consistent
  max_rounds <- 2000
  cv_xgb <- xgboost::xgb.cv(
    params = params, data = dall, nrounds = max_rounds, nfold = 5,
    early_stopping_rounds = 50, verbose = 0, metrics = "cox-nloglik"
  )

  best_nrounds <- cv_xgb$best_iteration
  if (is.null(best_nrounds) || length(best_nrounds) == 0) best_nrounds <- max_rounds

  message(paste(">>> Training Final XGBoost (nrounds:", best_nrounds, ")..."))
  final_xgb <- xgboost::xgb.train(params = params, data = dall, nrounds = best_nrounds, evals = list(train = dall), verbose = 0)

  imp_matrix <- xgboost::xgb.importance(feature_names = colnames(x_full), model = final_xgb)
  xgb_genes <- imp_matrix$Feature[1:min(20, nrow(imp_matrix))]

  # ============================================================================
  # Part 4: Package results
  # ============================================================================
  res <- list(
    genes = list(lasso = lasso_genes, rsf = rsf_genes, xgb = xgb_genes),
    models = list(lasso = cv_lasso, rsf = rf_final, xgb = final_xgb),
    plot_data = list(
      lasso_cv = df_lasso_plot,
      lasso_lambda_min = cv_lasso$lambda.min,
      lasso_lambda_1se = cv_lasso$lambda.1se,
      rsf_err = df_rsf_err,
      rsf_vimp = rf_final$importance,
      xgb_imp = imp_matrix
    )
  )
  class(res) <- "HCCPrognosticResult"
  return(res)
}

# ==============================================================================
# Part 2: Visualization functions (unchanged)
# ==============================================================================

#' Plot LASSO Cross-Validation Curve
#' @param res Object of class HCCPrognosticResult
#' @export
plot_lasso_cv <- function(res) {
  df <- res$plot_data$lasso_cv
  l_min <- log(res$plot_data$lasso_lambda_min)
  l_1se <- log(res$plot_data$lasso_lambda_1se)
  n_genes_min <- df$n_genes[which.min(abs(df$log_lambda - l_min))]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log_lambda, y = c_index)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = cv_lo, ymax = cv_up), color = "grey60", width = 0.1, alpha = 0.6) +
    ggplot2::geom_point(color = "#BC3C29", size = 2.5) +
    ggplot2::geom_vline(xintercept = l_min, linetype = "dashed", color = "black", size = 0.8) +
    ggplot2::geom_vline(xintercept = l_1se, linetype = "dotted", color = "grey40", size = 0.8) +
    ggplot2::annotate("text", x = l_min, y = min(df$cv_lo) * 1.01, label = paste0("Optimal: ", n_genes_min, " genes"), color = "#BC3C29", fontface = "bold", vjust = 1, hjust = 0.5) +
    ggplot2::labs(x = expression(Log(lambda)), y = "C-index (10-fold CV)") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"), axis.text = ggplot2::element_text(color = "black", size = 12), axis.title = ggplot2::element_text(face = "bold", size = 14))
  return(p)
}

#' Plot LASSO Coefficients
#' @param res Object of class HCCPrognosticResult
#' @export
plot_lasso_coef <- function(res) {
  tmp_coeffs <- stats::coef(res$models$lasso, s = "lambda.min")
  lasso_df <- data.frame(Gene = rownames(tmp_coeffs), Coef = as.numeric(tmp_coeffs))
  lasso_df <- lasso_df[lasso_df$Coef != 0 & lasso_df$Gene != "(Intercept)", ]
  lasso_df$Association <- ifelse(lasso_df$Coef > 0, "Positive", "Negative")
  colors_lasso <- c("Positive" = "#D58D4E", "Negative" = "#0072B5")

  p <- ggplot2::ggplot(lasso_df, ggplot2::aes(x = stats::reorder(Gene, Coef), y = Coef, fill = Association)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.8) +
    ggplot2::geom_col(width = 0.7, color = NA) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors_lasso) +
    ggplot2::labs(x = "", y = "LASSO Coefficient") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), axis.line = ggplot2::element_line(size = 0.8, color = "black"), axis.text = ggplot2::element_text(color = "black", size = 12), axis.text.y = ggplot2::element_text(face = "bold.italic"), axis.title = ggplot2::element_text(face = "bold", size = 14), legend.position = "top", legend.title = ggplot2::element_blank())
  return(p)
}

#' Plot RSF Convergence
#' @param res Object of class HCCPrognosticResult
#' @export
plot_rsf_process <- function(res) {
  df <- res$plot_data$rsf_err
  final_err <- tail(df$Error, 1)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Trees, y = Error)) +
    ggplot2::geom_line(color = "#469D76", linewidth = 1) +
    ggplot2::geom_hline(yintercept = final_err, linetype = "dashed", color = "#BC3C29", linewidth = 0.8) +
    ggplot2::annotate("text", x = max(df$Trees) * 0.6, y = max(df$Error), label = paste0("Final Error: ", round(final_err, 3)), color = "black", fontface = "bold", size = 4) +
    ggplot2::labs(x = "Number of Trees", y = "OOB Error Rate") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"), axis.text = ggplot2::element_text(color = "black", size = 12), axis.title = ggplot2::element_text(face = "bold", size = 14))
  return(p)
}

#' Plot RSF Variable Importance (Lollipop)
#' @param res Object of class HCCPrognosticResult
#' @export
plot_rsf_vimp <- function(res) {
  vimp_all <- res$plot_data$rsf_vimp
  top_genes <- res$genes$rsf

  df <- data.frame(Gene = names(vimp_all), VIMP = vimp_all)
  df <- df[df$Gene %in% top_genes, ]
  df <- df[order(df$VIMP), ]
  if(nrow(df) > 20) df <- tail(df, 20)
  df$Gene <- factor(df$Gene, levels = df$Gene)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Gene, y = VIMP)) +
    ggplot2::geom_segment(ggplot2::aes(x = Gene, xend = Gene, y = 0, yend = VIMP), color = "grey70", size = 1) +
    ggplot2::geom_point(color = "#0072B5", fill = "#0072B5", size = 4, shape = 21, stroke = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Variable Importance (Permutation)") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), axis.line = ggplot2::element_line(linewidth = 0.8, color = "black"), panel.grid.major.x = ggplot2::element_line(color = "grey90", linetype = "dashed"))
  return(p)
}

#' Plot XGBoost Importance
#' @param res Object of class HCCPrognosticResult
#' @export
plot_xgb_imp <- function(res) {
  xgb_plot_df <- res$plot_data$xgb_imp
  xgb_plot_df <- xgb_plot_df %>% dplyr::arrange(dplyr::desc(Gain)) %>% head(20)

  p <- ggplot2::ggplot(xgb_plot_df, ggplot2::aes(x = stats::reorder(Feature, Gain), y = Gain)) +
    ggplot2::geom_col(fill = "#D58D4E", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "Gain (Contribution to Model)") +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  return(p)
}

#' Plot Top Feature Intersection (Venn)
#' @param res Object of class HCCPrognosticResult
#' @param top_n Number of top features to select from each model for intersection.
#' @export
plot_model_venn <- function(res, top_n = 5) {
  .check_cran_packages("VennDiagram")

  lasso_cf <- stats::coef(res$models$lasso, s = "lambda.min")
  lasso_df <- data.frame(Gene = rownames(lasso_cf), Val = abs(as.numeric(lasso_cf)))
  lasso_top <- lasso_df[lasso_df$Val > 0 & lasso_df$Gene != "(Intercept)", ]
  lasso_top <- lasso_top[order(lasso_top$Val, decreasing = TRUE), ]
  genes_lasso <- head(lasso_top$Gene, top_n)

  rsf_imp <- res$plot_data$rsf_vimp
  genes_rsf <- names(sort(rsf_imp, decreasing = TRUE))[1:top_n]

  xgb_imp <- res$plot_data$xgb_imp
  genes_xgb <- head(xgb_imp$Feature, top_n)

  list_genes <- list(LASSO = genes_lasso, RSF = genes_rsf, XGBoost = genes_xgb)
  common_genes <- Reduce(intersect, list_genes)

  grid::grid.newpage()
  venn.plot <- VennDiagram::venn.diagram(
    x = list_genes, filename = NULL,
    col = "transparent", fill = c("#E64B35", "#00A087", "#3C5488"), alpha = 0.50,
    label.col = "white", cex = 1.5, fontfamily = "sans", fontface = "bold",
    cat.default.pos = "text", cat.col = c("#E64B35", "#00A087", "#3C5488"),
    cat.cex = 1.5, cat.fontfamily = "sans", cat.fontface = "bold",
    cat.pos = c(-27, 27, 180), cat.dist = c(0.055, 0.055, 0.055),
    margin = 0.1
  )
  grid::grid.draw(venn.plot)

  if(length(common_genes) > 0) {
    label_txt <- paste(head(common_genes, 3), collapse = "\n")
    if(length(common_genes) > 3) label_txt <- paste0(label_txt, "\n...")
    grid::grid.text(label = label_txt, x = 0.5, y = 0.5, gp = grid::gpar(col = "grey90", fontsize = 12, fontface = "bold"))
  }
}
#res_ml=res_123
#p1 <- plot_lasso_cv(res_ml)
#p1
#p2 <- plot_lasso_coef(res_ml)
#print(p2)
#p3 <- plot_rsf_process(res_ml)
# print(p3)
# p4 <- plot_rsf_vimp(res_ml)
# print(p4)
# p3
#p5 <- plot_xgb_imp(res_ml)
#p5
# plot_model_venn(res_ml, top_n = 5)
