# ==============================================================================
# MODULE: Training (PCA + LDA) on Available Ratios
# Inputs:  ratio_mat      (ratios x training_samples numeric matrix)
#          avail_features (data.frame with column Panel, Gene_A, Gene_B)
#          phenotype_vec  (factor; one entry per training sample)
#          metadata_df    (optional data.frame; used for ratio-level batch correction)
#          batch_col      (optional character; column in metadata_df for batch variable)
#          protect_cols   (optional character vector; biological columns to protect)
# Outputs: list(
#            lda_model           = trained lda object,
#            pca_ad              = prcomp object for AD panel,
#            pca_sig             = prcomp object for SIG panel,
#            n_pcs_ad            = integer; PCs retained for AD,
#            n_pcs_sig           = integer; PCs retained for SIG,
#            train_pcs           = data.frame of PC scores for training samples,
#            lda_coords          = matrix of LD scores for training samples,
#            lda_var             = numeric vector % variance per LD axis,
#            global_pca          = prcomp on all ratios (for global overview plot),
#            pca_var             = numeric vector % variance for global PCA,
#            phenotype_vec       = factor used for training (aligned),
#            loocv_acc           = numeric; leave-one-out CV accuracy (%),
#            batch_correction_level = character; "RATIO" or "NONE",
#            mat_ad_train_corr   = samples x AD-ratios corrected matrix,
#            mat_sig_train_corr  = samples x SIG-ratios corrected matrix (or NULL)
#          )
# ==============================================================================

#' Retrain the PCA + LDA ECMExplainer model using available ratios.
#'
#' When \code{metadata_df} and \code{batch_col} are supplied, ratio-level batch
#' correction is applied via \code{apply_ratio_batch_correction()} (V19 pipeline)
#' before PCA. LOOCV accuracy is always computed.
#'
#' @param ratio_mat      ratios x training_samples numeric matrix
#' @param avail_features data.frame; must have column Panel ("AD" or "SIG")
#' @param phenotype_vec  factor with one level per phenotype class
#' @param metadata_df    optional data.frame aligned to columns of ratio_mat;
#'                       required for ratio-level batch correction
#' @param batch_col      optional character; column in metadata_df to use as batch
#'                       (default "Study_ID")
#' @param protect_cols   optional character vector; columns in metadata_df to protect
#'                       as biological covariates during batch correction
#' @return named list with model objects and coordinates (see header)
train_ecmexplainer <- function(ratio_mat,
                               avail_features,
                               phenotype_vec,
                               metadata_df   = NULL,
                               batch_col     = "Study_ID",
                               protect_cols  = NULL) {

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for LDA training.")

  phenotype_vec <- as.factor(phenotype_vec)

  # ---- Split ratio matrix by panel type ----
  panel_col <- avail_features$Panel
  # Handle both "A/D" (from Back Up CSV) and "AD" / "SIG" naming conventions
  is_ad  <- panel_col %in% c("A/D", "AD")
  is_sig <- panel_col %in% c("SIG", "Signalling", "SIGNALLING")

  ad_rows  <- rownames(ratio_mat)[is_ad]
  sig_rows <- rownames(ratio_mat)[is_sig]

  # Graceful fallback: if one panel is empty, lump all ratios into ad
  if (length(ad_rows) == 0 && length(sig_rows) > 0) {
    ad_rows  <- sig_rows
    sig_rows <- character(0)
    message("[Train] No AD ratios found; treating all ratios as AD panel.")
  }

  mat_ad_train  <- t(ratio_mat[ad_rows,  , drop = FALSE])   # samples x AD-ratios
  mat_sig_train <- if (length(sig_rows) > 0)
    t(ratio_mat[sig_rows, , drop = FALSE])
  else
    matrix(0, nrow = nrow(mat_ad_train), ncol = 0)

  # ---- Ratio-level batch correction (V19 pipeline) ----
  # Applied to the ratio matrix before PCA if metadata + batch_col are supplied.
  # Uses limma::removeBatchEffect to remove study-level technical variation
  # while protecting biological phenotype signal.
  batch_correction_level <- "NONE"
  mat_ad_train_corr  <- mat_ad_train
  mat_sig_train_corr <- mat_sig_train

  if (!is.null(metadata_df) && batch_col %in% colnames(metadata_df) &&
      nlevels(factor(metadata_df[[batch_col]])) > 1) {

    ratio_mat_raw <- t(rbind(
      t(mat_ad_train),
      if (ncol(mat_sig_train) > 0) t(mat_sig_train) else matrix(0, 0, nrow(mat_ad_train))
    ))

    bc_result <- tryCatch(
      apply_ratio_batch_correction(
        ratio_mat    = ratio_mat,
        metadata_df  = metadata_df,
        batch_col    = batch_col,
        protect_cols = protect_cols
      ),
      error = function(e) {
        message("[Train] Ratio batch correction failed: ", e$message, " — proceeding without.")
        list(corrected_mat = ratio_mat, batch_used = NA)
      }
    )

    corrected_ratio_mat <- bc_result$corrected_mat
    mat_ad_train_corr   <- t(corrected_ratio_mat[ad_rows, , drop = FALSE])
    mat_sig_train_corr  <- if (length(sig_rows) > 0)
      t(corrected_ratio_mat[sig_rows, , drop = FALSE])
    else
      matrix(0, nrow = nrow(mat_ad_train_corr), ncol = 0)
    batch_correction_level <- "RATIO"
    message("[Train] Ratio-level batch correction applied (", batch_col, ").")
  }

  # ---- Global PCA (all ratios, for overview plot) — use corrected data ----
  ratio_for_pca <- if (batch_correction_level == "RATIO") corrected_ratio_mat else ratio_mat
  global_pca <- prcomp(t(ratio_for_pca), center = TRUE, scale. = TRUE)
  pca_var    <- round(100 * global_pca$sdev^2 / sum(global_pca$sdev^2), 1)

  # ---- Dual PCA (AD + SIG panels separately, corrected matrices) ----
  run_pca <- function(mat) {
    var_mask <- apply(mat, 2, var, na.rm = TRUE) > 0
    if (sum(var_mask) == 0) {
      return(prcomp(mat, center = FALSE, scale. = FALSE))
    }
    prcomp(mat[, var_mask, drop = FALSE], center = TRUE, scale. = TRUE)
  }

  pca_ad  <- run_pca(mat_ad_train_corr)
  pca_sig <- if (ncol(mat_sig_train_corr) > 0) run_pca(mat_sig_train_corr) else NULL

  # Number of PCs explaining >= 95% variance
  n_pcs <- function(pca_obj) {
    if (is.null(pca_obj) || length(pca_obj$sdev) == 0) return(1L)
    cum_var <- cumsum(pca_obj$sdev^2 / sum(pca_obj$sdev^2))
    max(1L, which(cum_var >= 0.95)[1], na.rm = TRUE)
  }

  n_pcs_ad  <- n_pcs(pca_ad)
  n_pcs_sig <- n_pcs(pca_sig)

  # ---- Build combined PC data.frame for LDA ----
  pc_ad_df <- as.data.frame(pca_ad$x[, seq_len(n_pcs_ad), drop = FALSE])
  colnames(pc_ad_df) <- paste0("AD_PC", seq_len(n_pcs_ad))

  if (!is.null(pca_sig) && n_pcs_sig > 0 && ncol(mat_sig_train_corr) > 0) {
    pc_sig_df <- as.data.frame(pca_sig$x[, seq_len(n_pcs_sig), drop = FALSE])
    colnames(pc_sig_df) <- paste0("SIG_PC", seq_len(n_pcs_sig))
    train_pcs <- cbind(pc_ad_df, pc_sig_df)
  } else {
    train_pcs <- pc_ad_df
  }

  # ---- Train LDA ----
  lda_model  <- MASS::lda(phenotype_vec ~ ., data = train_pcs)
  lda_out    <- predict(lda_model, train_pcs)
  lda_coords <- lda_out$x

  n_ld    <- ncol(lda_coords)
  lda_var <- if (n_ld > 0)
    round(100 * lda_model$svd^2 / sum(lda_model$svd^2), 1)
  else
    numeric(0)

  # ---- LOOCV accuracy ----
  loocv_acc <- tryCatch({
    lda_cv <- suppressWarnings(
      MASS::lda(train_pcs, grouping = phenotype_vec, CV = TRUE)
    )
    round(mean(lda_cv$class == phenotype_vec) * 100, 1)
  }, error = function(e) {
    message("[Train] LOOCV failed: ", e$message)
    NA_real_
  })

  list(
    lda_model              = lda_model,
    pca_ad                 = pca_ad,
    pca_sig                = pca_sig,
    n_pcs_ad               = n_pcs_ad,
    n_pcs_sig              = n_pcs_sig,
    train_pcs              = train_pcs,
    lda_coords             = lda_coords,
    lda_var                = lda_var,
    global_pca             = global_pca,
    pca_var                = pca_var,
    phenotype_vec          = phenotype_vec,
    ad_rows                = ad_rows,
    sig_rows               = sig_rows,
    loocv_acc              = loocv_acc,
    batch_correction_level = batch_correction_level,
    mat_ad_train_corr      = mat_ad_train_corr,
    mat_sig_train_corr     = mat_sig_train_corr
  )
}
