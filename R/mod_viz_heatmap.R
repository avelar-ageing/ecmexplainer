# ==============================================================================
# mod_viz_heatmap.R  —  Sign matrix heatmap + Distance heatmap
# ==============================================================================
#
# Implements two visualisations from the k-TSP manuscript:
#
#  1. plot_sign_heatmap()
#     Sign matrix {-1, 0, +1} as a ggplot2 tile heatmap.
#     blue = -1 (Gene_A < Gene_B), cream/white = 0, red = +1 (Gene_A > Gene_B)
#     Rows = ratios (Ratio_ID), columns = samples, annotated by phenotype.
#
#  2. plot_distance_heatmap()
#     Per-sample Euclidean distance to each training phenotype centroid in
#     LDA space.  Columns = phenotypes, rows = samples (or conditions).
#     Useful for visualising which centroid each sample is closest to.
#
# Both functions return ggplot objects; the caller is responsible for rendering.
# ==============================================================================

# Colour constants matching the manuscript
.SIGN_COLORS <- c("-1" = "#3182BD", "0" = "#F7F4F0", "+1" = "#E31A1C")

ECM_PHENO_COLORS <- c(
  Fibrosis      = "#377EB8",
  Hypoxia       = "#E41A1C",
  Inflammation  = "#4DAF4A",
  Proliferating = "#FF7F00",
  Quiescence    = "#A65628",
  Senescence    = "#333333"
)

# Fallback colour palette for phenotypes not in the default list
.pheno_palette <- function(levels_vec) {
  base_pal <- c("#377EB8","#E41A1C","#4DAF4A","#FF7F00","#A65628","#333333",
                "#984EA3","#F781BF","#999999","#A6CEE3")
  clrs <- ECM_PHENO_COLORS[levels_vec]
  missing <- is.na(clrs)
  if (any(missing)) {
    extra <- base_pal[seq_len(sum(missing))]
    clrs[missing] <- extra
  }
  clrs
}


# ==============================================================================
# plot_sign_heatmap()
# ==============================================================================
#' Render the sign matrix as a ggplot2 tile heatmap.
#'
#' @param sign_mat   [ratios × samples] integer matrix; values ∈ {-1, 0, +1}.
#'                   rownames = Ratio_ID; colnames = Sample_ID.
#' @param pheno_vec  named character vector (names = Sample_ID → phenotype label),
#'                   OR a factor/character vector aligned to columns of sign_mat.
#' @param title      character; plot title.
#' @param cluster_rows  logical; hierarchically cluster ratios? (default TRUE)
#' @param cluster_cols  logical; hierarchically cluster samples within phenotype?
#'                      (default FALSE — columns are ordered by phenotype)
#' @return ggplot object
plot_sign_heatmap <- function(sign_mat,
                              pheno_vec     = NULL,
                              title         = "ECMExplainer Sign Matrix",
                              cluster_rows  = TRUE,
                              cluster_cols  = FALSE) {

  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr",    quietly = TRUE)) stop("tidyr required")

  # ---- Cluster rows (ratios) if requested ----
  row_order <- seq_len(nrow(sign_mat))
  if (cluster_rows && nrow(sign_mat) > 2) {
    tryCatch({
      hc <- hclust(dist(sign_mat, method = "manhattan"), method = "ward.D2")
      row_order <- hc$order
    }, error = function(e) NULL)
  }
  sign_mat <- sign_mat[row_order, , drop = FALSE]

  # ---- Align phenotype vector to columns ----
  if (!is.null(pheno_vec)) {
    if (!is.null(names(pheno_vec))) {
      pheno_aligned <- pheno_vec[colnames(sign_mat)]
    } else {
      pheno_aligned <- pheno_vec
    }
    pheno_aligned[is.na(pheno_aligned)] <- "Unknown"
  } else {
    pheno_aligned <- rep("Sample", ncol(sign_mat))
  }

  # ---- Order columns by phenotype (then optionally cluster within) ----
  col_order <- order(pheno_aligned)
  if (cluster_cols) {
    pheno_levels <- unique(pheno_aligned[col_order])
    col_order_clust <- unlist(lapply(pheno_levels, function(p) {
      idx <- which(pheno_aligned == p)
      if (length(idx) <= 2) return(idx)
      sub_mat <- sign_mat[, idx, drop = FALSE]
      tryCatch({
        hc2 <- hclust(dist(t(sub_mat), method = "manhattan"), method = "ward.D2")
        idx[hc2$order]
      }, error = function(e) idx)
    }))
    col_order <- col_order_clust
  }
  sign_mat      <- sign_mat[, col_order, drop = FALSE]
  pheno_aligned <- pheno_aligned[col_order]

  # ---- Melt to long format ----
  df_long <- data.frame(
    Ratio     = rep(factor(rownames(sign_mat), levels = rev(rownames(sign_mat))),
                    times = ncol(sign_mat)),
    Sample_ID = rep(colnames(sign_mat), each = nrow(sign_mat)),
    Sign      = factor(as.vector(sign_mat), levels = c(-1L, 0L, 1L),
                       labels = c("-1", "0", "+1")),
    Phenotype = rep(pheno_aligned, each = nrow(sign_mat)),
    stringsAsFactors = FALSE
  )
  # Preserve column order in the factor
  df_long$Sample_ID <- factor(df_long$Sample_ID,
                               levels = unique(df_long$Sample_ID))

  pheno_levels <- unique(pheno_aligned)
  pheno_clrs   <- .pheno_palette(pheno_levels)

  # ---- Build ggplot ----
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample_ID, y = Ratio, fill = Sign)) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_manual(
      values = .SIGN_COLORS,
      name   = "sign(A \u2212 B)",
      labels = c("\u22121 (A < B)", "0 (tied)", "+1 (A > B)")
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(size = 7),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold", size = 13)
    )

  # Add phenotype colour bar as a strip on top (using facet hack)
  if (length(unique(pheno_aligned)) > 1) {
    df_long$Phenotype <- factor(df_long$Phenotype, levels = pheno_levels)
    p <- p +
      ggplot2::facet_grid(. ~ Phenotype, scales = "free_x", space = "free_x",
                          switch = "x") +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = NA, colour = NA),
        strip.text.x     = ggplot2::element_text(size = 9, angle = 45, hjust = 1),
        panel.spacing    = ggplot2::unit(0.5, "pt")
      )
  }

  p
}


# ==============================================================================
# plot_distance_heatmap()
# ==============================================================================
#' Render per-sample Euclidean distance to each training centroid.
#'
#' @param val_lda_coords  [n_samples × n_LDs] matrix of LDA scores for test data.
#' @param train_obj       model bundle (list with $lda_model and $train_phenotype).
#' @param pheno_vec       named/aligned character vector of test sample phenotypes
#'                        (used for row annotation; can be NULL).
#' @param condition_vec   optional character vector (length = n_samples);
#'                        if provided, rows are aggregated by condition (median dist).
#' @param title           character; plot title.
#' @return ggplot object
plot_distance_heatmap <- function(val_lda_coords,
                                  train_obj,
                                  pheno_vec     = NULL,
                                  condition_vec = NULL,
                                  title         = "Distance to Training Centroids (LDA space)") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  # ---- Compute training centroids in LDA space ----
  lda_scores_train <- predict(train_obj$lda_model,
                              newdata = as.data.frame(
                                if (!is.null(train_obj$mat_ad_train_corr)) {
                                  # Re-project training data
                                  cbind(
                                    predict(train_obj$pca_ad)[, seq_len(train_obj$n_pcs_ad), drop = FALSE],
                                    predict(train_obj$pca_sig)[, seq_len(train_obj$n_pcs_sig), drop = FALSE]
                                  )
                                } else {
                                  matrix(0, nrow = train_obj$n_training, ncol = 1)
                                }
                              ))$x

  pheno_levels <- train_obj$phenotypes %||% levels(train_obj$train_phenotype)
  train_pheno  <- as.character(train_obj$train_phenotype)

  # Use stored LDA scores if available (avoids re-projecting)
  if (!is.null(train_obj$lda_coords_train)) {
    lda_scores_train <- train_obj$lda_coords_train
  }

  centroids <- do.call(rbind, lapply(pheno_levels, function(p) {
    idx <- which(train_pheno == p)
    if (length(idx) == 0) return(rep(NA_real_, ncol(val_lda_coords)))
    colMeans(lda_scores_train[idx, seq_len(ncol(val_lda_coords)), drop = FALSE])
  }))
  rownames(centroids) <- pheno_levels

  # ---- Distance from each sample to each centroid ----
  n_s    <- nrow(val_lda_coords)
  n_p    <- length(pheno_levels)
  n_dim  <- min(ncol(val_lda_coords), ncol(centroids), na.rm = TRUE)
  dist_mat <- matrix(NA_real_, nrow = n_s, ncol = n_p,
                     dimnames = list(rownames(val_lda_coords), pheno_levels))
  for (p in pheno_levels) {
    if (all(is.na(centroids[p, ]))) next
    diff_sq <- sweep(val_lda_coords[, seq_len(n_dim), drop = FALSE],
                     2,
                     centroids[p, seq_len(n_dim)],
                     "-")^2
    dist_mat[, p] <- sqrt(rowSums(diff_sq))
  }

  # ---- Optionally aggregate by condition (median per condition) ----
  row_labels <- rownames(val_lda_coords)
  row_pheno  <- if (!is.null(pheno_vec)) {
    if (!is.null(names(pheno_vec))) pheno_vec[rownames(val_lda_coords)] else pheno_vec
  } else rep(NA_character_, n_s)

  if (!is.null(condition_vec)) {
    conds   <- unique(condition_vec)
    agg_mat <- matrix(NA_real_, nrow = length(conds), ncol = n_p,
                      dimnames = list(conds, pheno_levels))
    agg_phe <- character(length(conds))
    for (i in seq_along(conds)) {
      rows_i         <- which(condition_vec == conds[i])
      agg_mat[i, ]   <- apply(dist_mat[rows_i, , drop = FALSE], 2, median, na.rm = TRUE)
      agg_phe[i]     <- names(sort(table(row_pheno[rows_i]), decreasing = TRUE))[1]
    }
    dist_mat   <- agg_mat
    row_labels <- conds
    row_pheno  <- agg_phe
  }

  # ---- Melt to long format ----
  df_long <- data.frame(
    Sample    = rep(row_labels, times = n_p),
    Phenotype = rep(pheno_levels, each = nrow(dist_mat)),
    Distance  = as.vector(dist_mat),
    stringsAsFactors = FALSE
  )
  df_long$Sample    <- factor(df_long$Sample,    levels = row_labels)
  df_long$Phenotype <- factor(df_long$Phenotype, levels = pheno_levels)

  # Scale each column independently so colours reflect within-phenotype variation
  df_long$Dist_scaled <- ave(df_long$Distance, df_long$Phenotype, FUN = function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) return(rep(0.5, length(x)))
    (x - r[1]) / diff(r)
  })

  pheno_clrs <- .pheno_palette(pheno_levels)

  p <- ggplot2::ggplot(df_long,
         ggplot2::aes(x = Phenotype, y = Sample, fill = Dist_scaled)) +
    ggplot2::geom_tile(colour = "white", size = 0.3) +
    ggplot2::scale_fill_gradient(low = "#D73027", high = "#4575B4",
                                 na.value = "grey90",
                                 name = "Distance\n(scaled)") +
    ggplot2::labs(title = title, x = "Training phenotype centroid", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 40, hjust = 1, size = 10),
      axis.text.y     = ggplot2::element_text(size = 8),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      plot.title      = ggplot2::element_text(face = "bold", size = 13)
    )

  # Annotate each cell with rounded distance value
  df_long$Label <- ifelse(is.na(df_long$Distance), "",
                          sprintf("%.1f", df_long$Distance))
  p <- p +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 2.8,
                       colour = "black", data = df_long)
  p
}
