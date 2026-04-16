# ==============================================================================
# MODULE: Normalization
# RNA-seq only (V19 pipeline — microarray integration removed).
#
# Inputs:  counts_matrix   (genes x samples, raw integer counts OR CPM)
# Outputs: list(
#            log_cpm = numeric matrix genes x samples (log2(CPM + 1))
#          )
# ==============================================================================

#' Scale raw base-pair coverage counts to read-count equivalents (recount3-style).
#'
#' recount3 stores gene counts as the sum of base-level coverage across each gene body.
#' These values are much larger than conventional read counts and must be scaled before
#' CPM normalisation. This function divides each sample's counts by (colSums / 1e6),
#' converting coverage-based integers into read-count-scale equivalents.
#'
#' @param counts_matrix  genes × samples integer matrix of raw coverage counts
#' @return               matrix of the same dimensions with scaled integer values
transform_counts <- function(counts_matrix) {
  counts_matrix <- as.matrix(counts_matrix)
  scale_factors <- colSums(counts_matrix) / 1e6
  scale_factors[scale_factors == 0] <- 1   # guard against empty samples
  round(sweep(counts_matrix, 2, scale_factors, "/"))
}


#' Normalise a raw counts matrix to log2(CPM + 1).
#'
#' All input is treated as RNA-seq. Pass \code{skip_rnaseq_cpm = TRUE} when
#' the matrix is already in CPM units (e.g. pre-normalised uploads).
#'
#' @param counts_matrix  genes x samples integer/numeric matrix; rownames = gene symbols
#' @param skip_rnaseq_cpm logical (default FALSE). When TRUE, the data is treated as
#'                        already in CPM units — library-size normalisation is skipped
#'                        and only log2(x + 1) is applied.
#' @return named list with $log_cpm (matrix)
normalize_counts <- function(counts_matrix,
                             skip_rnaseq_cpm = FALSE) {

  stopifnot(is.matrix(counts_matrix) || is.data.frame(counts_matrix))
  counts_matrix <- as.matrix(counts_matrix)

  if (isTRUE(skip_rnaseq_cpm)) {
    # Data already in CPM — just apply log2(CPM + 1)
    log_cpm <- log2(counts_matrix + 1)
  } else {
    # Standard: library-size normalisation -> log2(CPM + 1)
    lib_sizes <- colSums(counts_matrix)
    lib_sizes[lib_sizes == 0] <- 1    # guard against zero-library samples
    cpm_mat <- sweep(counts_matrix, 2, lib_sizes / 1e6, "/")
    log_cpm <- log2(cpm_mat + 1)
  }

  list(log_cpm = log_cpm)
}
