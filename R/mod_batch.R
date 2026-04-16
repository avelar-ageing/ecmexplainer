# ==============================================================================
# MODULE: Batch Correction
#
# Two functions:
#
#  apply_ratio_batch_correction()   [PRIMARY — V19 pipeline]
#    Operates on a RATIO × samples matrix using limma::removeBatchEffect.
#    This is the method used during training in the V19 pipeline: batch
#    correction is deferred until ratio computation, not applied at the
#    gene level.  No sva dependency required.
#
#  apply_batch_correction()         [LEGACY / optional user data]
#    Operates on a GENE × samples matrix using sva::ComBat.
#    Retained for cases where a user explicitly requests gene-level
#    correction on their own uploaded data (optional step in Test Wizard).
# ==============================================================================

#' Apply ComBat batch correction to a logCPM matrix.
#'
#' @param log_cpm_mat  genes x samples numeric matrix
#' @param metadata_df  data.frame; rows must align with columns of log_cpm_mat
#' @param batch_col    column name in metadata_df to use as batch variable
#' @param protect_cols character vector of column names to include in the
#'                     biological model (model.matrix covariates). NULL = intercept only.
#' @return list with $corrected_mat and $batch_used
apply_batch_correction <- function(log_cpm_mat,
                                   metadata_df,
                                   batch_col    = "Study_ID",
                                   protect_cols = NULL) {

  if (!requireNamespace("sva", quietly = TRUE))
    stop("Package 'sva' is required for batch correction.")

  stopifnot(is.matrix(log_cpm_mat))
  stopifnot(batch_col %in% colnames(metadata_df))
  stopifnot(nrow(metadata_df) == ncol(log_cpm_mat))

  batch <- factor(metadata_df[[batch_col]])

  # Need at least 2 batches to run ComBat
  if (nlevels(batch) < 2) {
    message("[BatchCorrection] Only one batch level found – skipping ComBat.")
    return(list(corrected_mat = log_cpm_mat, batch_used = batch_col))
  }

  # Build covariate model matrix
  if (!is.null(protect_cols) && length(protect_cols) > 0) {
    valid_protect <- intersect(protect_cols, colnames(metadata_df))
    if (length(valid_protect) > 0) {
      formula_str <- paste("~", paste(valid_protect, collapse = " + "))
      mod <- tryCatch(
        model.matrix(as.formula(formula_str), data = metadata_df),
        error = function(e) {
          message("[BatchCorrection] Could not build model matrix: ", e$message,
                  " – using intercept only.")
          model.matrix(~1, data = metadata_df)
        }
      )
    } else {
      mod <- model.matrix(~1, data = metadata_df)
    }
  } else {
    mod <- model.matrix(~1, data = metadata_df)
  }

  # Remove zero-variance genes before ComBat
  gene_vars  <- apply(log_cpm_mat, 1, var, na.rm = TRUE)
  keep_genes <- !is.na(gene_vars) & gene_vars > 0
  if (sum(keep_genes) == 0)
    stop("[BatchCorrection] No genes with non-zero variance found.")

  mat_filt <- log_cpm_mat[keep_genes, , drop = FALSE]

  corrected_filt <- tryCatch(
    sva::ComBat(
      dat         = mat_filt,
      batch       = batch,
      mod         = mod,
      par.prior   = TRUE,
      prior.plots = FALSE
    ),
    error = function(e) {
      warning("[BatchCorrection] ComBat failed: ", e$message,
              " – returning uncorrected matrix.")
      mat_filt
    }
  )

  # Reconstruct full matrix (zero-variance genes pass through uncorrected)
  out_mat <- log_cpm_mat
  out_mat[keep_genes, ] <- corrected_filt

  list(corrected_mat = out_mat, batch_used = batch_col)
}


# ==============================================================================
# V19 PRIMARY: Ratio-level batch correction via limma::removeBatchEffect
# ==============================================================================

#' Apply ratio-level batch correction using limma::removeBatchEffect.
#'
#' This mirrors the V19 training pipeline where correction is applied to the
#' computed ratio matrix (not the raw gene expression), preserving the
#' biological signal more faithfully.
#'
#' @param ratio_mat    ratios × samples numeric matrix
#' @param metadata_df  data.frame; rows must align with columns of ratio_mat
#' @param batch_col    column name in metadata_df to use as batch variable (default "Study_ID")
#' @param protect_cols character vector of column names to include in the
#'                     biological model (model.matrix covariates). NULL = intercept only.
#' @return list with $corrected_mat (ratios × samples) and $batch_used
apply_ratio_batch_correction <- function(ratio_mat,
                                         metadata_df,
                                         batch_col    = "Study_ID",
                                         protect_cols = NULL) {

  if (!requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' is required for ratio-level batch correction.")

  stopifnot(is.matrix(ratio_mat))
  stopifnot(batch_col %in% colnames(metadata_df))
  stopifnot(nrow(metadata_df) == ncol(ratio_mat))

  batch <- factor(metadata_df[[batch_col]])

  # Need at least 2 batches to correct
  if (nlevels(batch) < 2) {
    message("[RatioBatch] Only one batch level found – skipping removeBatchEffect.")
    return(list(corrected_mat = ratio_mat, batch_used = batch_col))
  }

  # Build covariate design matrix
  if (!is.null(protect_cols) && length(protect_cols) > 0) {
    valid_protect <- intersect(protect_cols, colnames(metadata_df))
    if (length(valid_protect) > 0) {
      formula_str <- paste("~", paste(valid_protect, collapse = " + "))
      design_mat  <- tryCatch(
        model.matrix(as.formula(formula_str), data = metadata_df),
        error = function(e) {
          message("[RatioBatch] Could not build model matrix: ", e$message,
                  " – using intercept only.")
          model.matrix(~1, data = metadata_df)
        }
      )
    } else {
      design_mat <- model.matrix(~1, data = metadata_df)
    }
  } else {
    design_mat <- model.matrix(~1, data = metadata_df)
  }

  corrected <- tryCatch(
    limma::removeBatchEffect(ratio_mat, batch = batch, design = design_mat),
    error = function(e) {
      warning("[RatioBatch] removeBatchEffect failed: ", e$message,
              " – returning uncorrected ratio matrix.")
      ratio_mat
    }
  )

  list(corrected_mat = corrected, batch_used = batch_col)
}

