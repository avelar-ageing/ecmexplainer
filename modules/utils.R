# ==============================================================================
# utils.R – Shared constants, colour palette, and small helpers
# Inputs : none (sourced at app start)
# Outputs: ECM_COLORS, add_alpha(), safe_prcomp(), pcs_for_95pct(), %||%
# ==============================================================================

# ── Null-coalescing operator ──────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Publication colour palette for phenotypes ─────────────────────────────────
ECM_COLORS <- c(
  Apoptosis      = "#E41A1C",
  Fibrosis       = "#377EB8",
  Inflammation   = "#4DAF4A",
  OxidativeStress = "#984EA3",
  Proliferating  = "#FF7F00",
  Quiescence     = "#A65628",
  Senescence     = "#333333"
)

# Add alpha channel to a named colour vector
add_alpha <- function(colors, alpha = 0.85) {
  r <- col2rgb(colors, alpha = TRUE)
  r[4, ] <- round(r[4, ] * alpha)
  grDevices::rgb(r[1, ], r[2, ], r[3, ], r[4, ], maxColorValue = 255)
}

# Number of PCs explaining >= pct variance (default 95 %)
pcs_for_95pct <- function(pca_obj, pct = 0.95) {
  cumvar <- cumsum(pca_obj$sdev^2 / sum(pca_obj$sdev^2))
  max(1L, which(cumvar >= pct)[1L], na.rm = TRUE)
}

# Safe prcomp: remove zero-variance columns before calling prcomp
safe_prcomp <- function(mat, ...) {
  keep <- apply(mat, 2, var) > 0
  if (sum(keep) < 1) stop("No variable features found after variance filter.")
  list(
    pca    = prcomp(mat[, keep, drop = FALSE], ...),
    keep   = keep
  )
}

# Project new samples into an existing PCA space (no data leakage)
project_pca <- function(new_mat, pca_obj, keep_mask, n_pcs) {
  new_sub <- new_mat[, keep_mask, drop = FALSE]
  scale(new_sub,
        center = pca_obj$center,
        scale  = pca_obj$scale) %*% pca_obj$rotation[, seq_len(n_pcs), drop = FALSE]
}
