# ==============================================================================
# MODULE: Prediction
# Inputs:  test_expr_mat  (genes x test_samples numeric matrix; normalised,
#                          batch-corrected in the same way as training data)
#          train_obj      (list returned by train_ecmexplainer())
#          avail_features (data.frame; same used for training)
# Outputs: list(
#            final_report  = data.frame: Sample_ID, Predicted_Phenotype, posterior %
#            val_pcs       = data.frame of PC scores for test samples (for LDA plot)
#            val_lda_coords = matrix of LD scores for test samples
#          )
# ==============================================================================

#' Project new samples into the trained ECMExplainer geometric space and predict.
#'
#' Automatically dispatches to \code{predict_ktsp()} when the model bundle was
#' trained with the binary k-TSP method (\code{train_obj$encoding == "sign"}).
#'
#' @param test_expr_mat         genes x test_samples numeric matrix; log2(CPM+1)
#' @param train_obj             list returned by \code{train_ecmexplainer()} or
#'                              loaded from a bundle .rds file
#' @param avail_features        data.frame; columns Panel, Gene_A, Gene_B
#' @param precomputed_ratio_mat optional features x test_samples matrix; skips
#'                              ratio computation (continuous mode only)
#' @return named list with $final_report, $val_pcs, $val_lda_coords
predict_ecmexplainer <- function(test_expr_mat,
                                 train_obj,
                                 avail_features,
                                 precomputed_ratio_mat = NULL) {

  # ---- Dispatch: binary k-TSP (sign-encoded) model --------------------------
  if (identical(train_obj$encoding, "sign")) {
    res <- predict_ktsp(test_expr_mat, train_obj)
    return(res)
  }

  ad_rows  <- train_obj$ad_rows
  sig_rows <- train_obj$sig_rows
  pca_ad   <- train_obj$pca_ad
  pca_sig  <- train_obj$pca_sig
  n_pcs_ad  <- train_obj$n_pcs_ad
  n_pcs_sig <- train_obj$n_pcs_sig

  # ---- Compute ratios for test samples (or use pre-computed matrix) ----
  if (!is.null(precomputed_ratio_mat)) {
    # Gene-mode: caller pre-computed the feature matrix (features x test_samples)
    val_ratio <- precomputed_ratio_mat
  } else {
    # Standard ratio mode
    n_avail <- nrow(avail_features)
    ratio_names <- paste0(avail_features$Gene_A, " / ", avail_features$Gene_B)
    val_ratio <- matrix(
      NA_real_,
      nrow     = n_avail,
      ncol     = ncol(test_expr_mat),
      dimnames = list(ratio_names, colnames(test_expr_mat))
    )
    for (k in seq_len(n_avail)) {
      val_ratio[k, ] <- test_expr_mat[avail_features$Gene_A[k], ] -
                        test_expr_mat[avail_features$Gene_B[k], ]
    }
  }

  mat_ad_val  <- t(val_ratio[ad_rows,  , drop = FALSE])
  mat_sig_val <- if (length(sig_rows) > 0)
    t(val_ratio[sig_rows, , drop = FALSE])
  else
    matrix(0, nrow = ncol(test_expr_mat), ncol = 0)

  # ---- Project into training PCA space (no data leakage) ----
  project_into_pca <- function(val_mat, pca_obj, n_pcs) {
    if (is.null(pca_obj) || ncol(val_mat) == 0) return(NULL)
    # Use only the genes that entered prcomp (those with var > 0 during training)
    train_genes <- names(pca_obj$center)
    val_subset  <- val_mat[, train_genes, drop = FALSE]
    scaled_val  <- scale(val_subset,
                         center = pca_obj$center,
                         scale  = pca_obj$scale)
    scaled_val %*% pca_obj$rotation[, seq_len(n_pcs), drop = FALSE]
  }

  proj_ad  <- project_into_pca(mat_ad_val,  pca_ad,  n_pcs_ad)
  proj_sig <- project_into_pca(mat_sig_val, pca_sig, n_pcs_sig)

  # ---- Assemble val_pcs data.frame matching training column layout ----
  if (!is.null(proj_ad)) {
    proj_ad_df <- as.data.frame(proj_ad)
    colnames(proj_ad_df) <- paste0("AD_PC", seq_len(n_pcs_ad))
  } else {
    stop("[Predict] AD PCA projection failed.")
  }

  if (!is.null(proj_sig) && ncol(proj_sig) > 0) {
    proj_sig_df <- as.data.frame(proj_sig)
    colnames(proj_sig_df) <- paste0("SIG_PC", seq_len(n_pcs_sig))
    val_pcs <- cbind(proj_ad_df, proj_sig_df)
  } else {
    val_pcs <- proj_ad_df
  }

  rownames(val_pcs) <- colnames(test_expr_mat)

  # ---- LDA prediction ----
  preds          <- predict(train_obj$lda_model, val_pcs)
  val_lda_coords <- preds$x
  post_probs     <- as.data.frame(round(preds$posterior * 100, 2))

  final_report <- data.frame(
    Sample_ID           = colnames(test_expr_mat),
    Predicted_Phenotype = as.character(preds$class),
    stringsAsFactors    = FALSE
  )
  final_report <- cbind(final_report, post_probs)

  list(
    final_report   = final_report,
    val_pcs        = val_pcs,
    val_lda_coords = val_lda_coords
  )
}
