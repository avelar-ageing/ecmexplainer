# ==============================================================================
# MODULE: Ratio Detection & Computation
# Inputs:  expr_mat      (genes x samples numeric matrix)
#          features_df   (data.frame with columns: Panel, Gene_A, Gene_B)
# Outputs: list(
#            avail_features = subset of features_df with both genes present,
#            ratio_mat      = ratios x samples matrix  (Gene_A expr - Gene_B expr),
#            n_total        = integer total features in features_df,
#            n_avail        = integer available features found in expr_mat
#          )
# Notes:   Ratio is computed as log-difference (logCPM(A) - logCPM(B)),
#          which is equivalent to log2 fold-change in expression space.
# ==============================================================================

#' Detect available ECM ratios and compute the ratio matrix.
#'
#' @param expr_mat    genes x samples numeric matrix; rownames = gene symbols
#' @param features_df data.frame with columns Panel, Gene_A, Gene_B
#' @return named list: $avail_features, $ratio_mat, $n_total, $n_avail
compute_ratios <- function(expr_mat, features_df) {

  stopifnot(is.matrix(expr_mat))
  stopifnot(all(c("Gene_A", "Gene_B") %in% colnames(features_df)))

  avail <- features_df[
    features_df$Gene_A %in% rownames(expr_mat) &
    features_df$Gene_B %in% rownames(expr_mat), ,
    drop = FALSE
  ]

  n_avail <- nrow(avail)

  if (n_avail == 0) {
    stop(paste0(
      "[Ratios] No ECMExplainer features could be found in the supplied expression matrix.\n",
      "Ensure gene symbols match those used during training (human HGNC symbols)."
    ))
  }

  ratio_names <- paste0(avail$Gene_A, " / ", avail$Gene_B)
  ratio_mat   <- matrix(
    NA_real_,
    nrow     = n_avail,
    ncol     = ncol(expr_mat),
    dimnames = list(ratio_names, colnames(expr_mat))
  )

  for (k in seq_len(n_avail)) {
    ratio_mat[k, ] <- expr_mat[avail$Gene_A[k], ] - expr_mat[avail$Gene_B[k], ]
  }

  list(
    avail_features = avail,
    ratio_mat      = ratio_mat,
    n_total        = nrow(features_df),
    n_avail        = n_avail
  )
}
