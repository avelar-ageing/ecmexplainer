# Suppress R CMD CHECK notes for ggplot2 NSE column references
utils::globalVariables(c(
  "Phenotype", "Sample_ID", "Set", "Probability", ".data",
  "ColourBy", "Label"
))

# ==============================================================================
# MODULE: Plots
# All plot functions return a ggplot2 object (or plotly where specified).
# Standardised colour palette matches the publication figures.
#
# Exported functions:
#   ecm_colors()             - named colour palette
#   plot_pca_training()      - ggplot: global PCA of training data
#   plot_lda_training()      - ggplot: LDA of training data
#   plot_pca_combined()      - ggplot: PCA training + test samples overlaid
#   plot_lda_combined()      - ggplot: LDA training + test samples overlaid
#   plot_prediction_probs()  - ggplot: stacked bar of posterior probabilities
# ==============================================================================

# ---- Colour palette ----

# ---- Colour palette ----

# Maximum-distinctness palette: Okabe-Ito (colorblind-safe, 8 colours) extended
# with Kelly high-contrast colours (positions 9-22).  Colours at positions
# 1-22 are chosen to maximise perceptual separation across the full set.
# The 7 built-in ECM phenotype names are pinned to fixed positions so they
# always render with the same colour regardless of which other classes appear.
.ECM_PALETTE_SEQ <- c(
  # ---- Okabe-Ito (positions 1-8, colorblind-safe) ----
  "#E69F00",  # 1  orange
  "#56B4E9",  # 2  sky blue
  "#009E73",  # 3  bluish green
  "#F0E442",  # 4  yellow
  "#0072B2",  # 5  blue
  "#D55E00",  # 6  vermillion
  "#CC79A7",  # 7  reddish purple
  "#000000",  # 8  black
  # ---- Kelly extensions (positions 9-22, maximum contrast) ----
  "#F2F3F4",  # 9  white (border visible)
  "#222222",  # 10 dark grey  (skip if too close to black, see below)
  "#F3C300",  # 11 vivid yellow
  "#875692",  # 12 strong purple
  "#F38400",  # 13 vivid orange
  "#A1CAF1",  # 14 very light blue
  "#BE0032",  # 15 vivid red
  "#C2B280",  # 16 grayish yellow
  "#848482",  # 17 medium grey
  "#008856",  # 18 vivid green
  "#E68FAC",  # 19 strong purplish pink
  "#0067A5",  # 20 strong blue
  "#F99379",  # 21 strong yellowish pink
  "#604E97"   # 22 strong violet
)

# Fixed colour assignments for the 7 built-in ECM phenotype classes.
# These are deliberately spread across the palette for maximum separation.
.ECM_NAMED <- c(
  Apoptosis       = "#D55E00",   # vermillion  (pos 6)
  Fibrosis        = "#0072B2",   # blue        (pos 5)
  Inflammation    = "#009E73",   # bluish green(pos 3)
  OxidativeStress = "#875692",   # strong purple(pos 12)
  Proliferating   = "#E69F00",   # orange      (pos 1)
  Quiescence      = "#56B4E9",   # sky blue    (pos 2)
  Senescence      = "#BE0032"    # vivid red   (pos 15)
)

#' Return the standard ECMExplainER phenotype colour palette.
#' @return named character vector of hex colours
ecm_colors <- function() .ECM_NAMED

# ---- Internal helpers ----
.ensure_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required for plotting.")
}

# Return a palette that covers every level in `levels_vec`.
# Named ECM classes always get their fixed colour; any unknown levels are
# assigned colours from the maximally-distinct sequence, skipping slots
# already used by named classes, so adjacent colours are never similar.
.pal_for_levels <- function(levels_vec) {
  lvls   <- as.character(unique(levels_vec))
  result <- setNames(rep(NA_character_, length(lvls)), lvls)

  # 1. Assign fixed colours to known ECM classes
  known <- intersect(lvls, names(.ECM_NAMED))
  result[known] <- .ECM_NAMED[known]

  # 2. Assign palette slots (in order) to unknown classes,
  #    skipping any slots whose colour is already used by a named class
  used_cols <- .ECM_NAMED[known]
  pool      <- .ECM_PALETTE_SEQ[!.ECM_PALETTE_SEQ %in% used_cols]
  unknown   <- lvls[is.na(result)]
  if (length(unknown) > 0) {
    n_assign <- min(length(unknown), length(pool))
    result[unknown[seq_len(n_assign)]] <- pool[seq_len(n_assign)]
    # Fallback: if we somehow exceeded the pool, repeat from the start
    if (length(unknown) > n_assign) {
      extras <- unknown[(n_assign + 1):length(unknown)]
      result[extras] <- .ECM_PALETTE_SEQ[((seq_along(extras) - 1) %% length(.ECM_PALETTE_SEQ)) + 1]
    }
  }
  result
}

# ---- Training PCA ----

#' Plot global PCA of training samples coloured by Phenotype.
#'
#' @param train_obj  list from train_ecmexplainer()
#' @param metadata   data.frame with at least columns Sample_ID and Phenotype
#' @param pc_x       integer; which PC to use as x-axis (default 1)
#' @param pc_y       integer; which PC to use as y-axis (default 2)
#' @return ggplot object
plot_pca_training <- function(train_obj, metadata, pc_x = 1L, pc_y = 2L) {
  .ensure_ggplot()

  gp    <- train_obj$global_pca
  pvar  <- train_obj$pca_var
  df    <- as.data.frame(gp$x)
  df$Sample_ID <- rownames(df)
  df    <- merge(df, metadata[, c("Sample_ID", "Phenotype")], by = "Sample_ID")
  pal   <- .pal_for_levels(unique(df$Phenotype))

  x_lab <- paste0("PC", pc_x, " (", pvar[pc_x], "%)")
  y_lab <- paste0("PC", pc_y, " (", pvar[pc_y], "%)")

  ggplot2::ggplot(df, ggplot2::aes(
    x      = .data[[paste0("PC", pc_x)]],
    y      = .data[[paste0("PC", pc_y)]],
    fill   = Phenotype,
    label  = Sample_ID
  )) +
    ggplot2::geom_point(shape = 21, colour = "black", size = 3.5, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
    ggplot2::labs(title = "Training Global PCA", x = x_lab, y = y_lab) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

# ---- Training LDA ----

#' Plot LDA of training samples coloured by Phenotype.
#'
#' @param train_obj  list from train_ecmexplainer()
#' @param metadata   data.frame with at least columns Sample_ID and Phenotype
#' @param ld_x       integer; which LD to use as x-axis (default 1)
#' @param ld_y       integer; which LD to use as y-axis (default 2)
#' @return ggplot object
plot_lda_training <- function(train_obj, metadata, ld_x = 1L, ld_y = 2L) {
  .ensure_ggplot()

  lda_c <- as.data.frame(train_obj$lda_coords)
  lv    <- train_obj$lda_var

  # Guard: if only 1 LD axis
  if (ncol(lda_c) < 2) {
    lda_c$LD2 <- 0
    lv <- c(lv, 0)
  }

  lda_c$Sample_ID <- rownames(lda_c)
  lda_c <- merge(lda_c, metadata[, c("Sample_ID", "Phenotype")], by = "Sample_ID")
  pal   <- .pal_for_levels(unique(lda_c$Phenotype))

  x_lab <- paste0("LD", ld_x, " (", lv[ld_x], "%)")
  y_lab <- paste0("LD", ld_y, " (", lv[ld_y], "%)")

  ggplot2::ggplot(lda_c, ggplot2::aes(
    x     = .data[[paste0("LD", ld_x)]],
    y     = .data[[paste0("LD", ld_y)]],
    fill  = Phenotype,
    label = Sample_ID
  )) +
    ggplot2::geom_point(shape = 21, colour = "black", size = 3.5, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
    ggplot2::labs(title = "Training LDA", x = x_lab, y = y_lab) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

# ---- Combined PCA (training + test) ----

#' Overlay test samples onto training global PCA.
#'
#' @param train_obj   list from train_ecmexplainer()
#' @param pred_obj    list from predict_ecmexplainer()
#' @param metadata    data.frame with Sample_ID + Phenotype for training
#' @param pc_x        integer x PC
#' @param pc_y        integer y PC
#' @return ggplot object
plot_pca_combined <- function(train_obj, pred_obj, metadata, pc_x = 1L, pc_y = 2L) {
  .ensure_ggplot()
  pal <- ecm_colors()

  gp   <- train_obj$global_pca
  pvar <- train_obj$pca_var

  # Training coordinates
  train_df <- as.data.frame(gp$x[, c(paste0("PC", pc_x), paste0("PC", pc_y)), drop = FALSE])
  train_df$Sample_ID <- rownames(train_df)
  train_df <- merge(train_df, metadata[, c("Sample_ID", "Phenotype")], by = "Sample_ID")
  train_df$Set <- "Training"

  # Test coordinates: project val_pcs back through global PCA is complex; use prediction LDA
  # For PCA overlay we compute test ratios projection directly from the global PCA object
  # using the AD ratio matrix stored in pred_obj if available.
  # Simpler: annotate with predicted phenotype
  test_df <- pred_obj$final_report[, c("Sample_ID", "Predicted_Phenotype")]
  colnames(test_df)[2] <- "Phenotype"

  # We do not have global PCA coords for test samples directly without the ratio_mat.
  # Return only training PCA with a note; the combined plot is in LDA space.
  # If val_pcs is available, use it.
  # This function is best shown in LDA; see plot_lda_combined().
  message("[Plots] PCA combined overlay not yet available – use plot_lda_combined().")
  plot_pca_training(train_obj, metadata, pc_x, pc_y)
}

# ---- Combined LDA (training + test) ----

#' Overlay test samples onto training LDA coordinates.
#'
#' @param train_obj   list from train_ecmexplainer()
#' @param pred_obj    list from predict_ecmexplainer()
#' @param train_meta  data.frame with Sample_ID + Phenotype for training
#' @param ld_x        integer x LD
#' @param ld_y        integer y LD
#' @return ggplot object
plot_lda_combined <- function(train_obj, pred_obj, train_meta, ld_x = 1L, ld_y = 2L,
                              highlight_ids = character(0)) {
  .ensure_ggplot()

  lv <- train_obj$lda_var
  ld_cols_train <- paste0("LD", seq_len(ncol(train_obj$lda_coords)))

  n_test <- if (!is.null(pred_obj$val_lda_coords)) nrow(pred_obj$val_lda_coords) else 0L

  # Build test df — may be zero rows when all samples are filtered out
  if (n_test > 0) {
    ld_cols_test <- paste0("LD", seq_len(ncol(pred_obj$val_lda_coords)))
    test_df <- as.data.frame(pred_obj$val_lda_coords)
    colnames(test_df) <- ld_cols_test[seq_len(ncol(test_df))]
    test_df$Sample_ID <- pred_obj$final_report$Sample_ID
    test_df$Phenotype <- pred_obj$final_report$Predicted_Phenotype
    test_df$Set <- "Test (Predicted)"
  } else {
    ld_cols_test <- ld_cols_train   # use train cols as reference
    test_df <- data.frame(matrix(ncol = length(ld_cols_test) + 3, nrow = 0))
    colnames(test_df) <- c(ld_cols_test, "Sample_ID", "Phenotype", "Set")
  }

  # Build training df — only if train_meta has rows
  if (!is.null(train_meta) && nrow(train_meta) > 0) {
    train_df <- as.data.frame(train_obj$lda_coords)
    colnames(train_df) <- ld_cols_train
    train_df$Sample_ID <- rownames(train_df)
    # Use phenotype_vec directly — it is guaranteed row-aligned with lda_coords by construction
    train_df$Phenotype <- as.character(train_obj$phenotype_vec)
    train_df$Set <- "Training"
    shared_ld <- intersect(colnames(train_df), colnames(test_df))
    combined <- rbind(
      train_df[, c(shared_ld, "Sample_ID", "Phenotype", "Set")],
      test_df[,  c(shared_ld, "Sample_ID", "Phenotype", "Set")]
    )
  } else {
    combined <- test_df[, c(intersect(ld_cols_test, ld_cols_train[seq_along(ld_cols_test)]),
                            "Sample_ID", "Phenotype", "Set"), drop = FALSE]
    # ensure LD1/LD2 exist
    if (!"LD1" %in% colnames(combined)) combined$LD1 <- 0
    if (!"LD2" %in% colnames(combined)) combined$LD2 <- 0
  }

  if (ncol(as.data.frame(train_obj$lda_coords)) < 2) {
    combined$LD2 <- 0
    lv <- c(lv, 0)
  }

  x_col <- paste0("LD", ld_x)
  y_col <- paste0("LD", ld_y)
  x_lab <- paste0(x_col, " (", lv[ld_x], "%)")
  y_lab <- paste0(y_col, " (", lv[ld_y], "%)")

  # Build palette to cover every phenotype class present in this combined dataset
  pal <- .pal_for_levels(unique(combined$Phenotype))

  # Highlighted points overlay
  hl_df <- combined[combined$Sample_ID %in% highlight_ids, , drop = FALSE]

  p <- ggplot2::ggplot(combined, ggplot2::aes(
    x     = .data[[x_col]],
    y     = .data[[y_col]],
    fill  = Phenotype,
    shape = Set,
    label = Sample_ID
  )) +
    ggplot2::geom_point(colour = "black", size = 3.5, alpha = 0.88) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
    ggplot2::scale_shape_manual(values = c("Training" = 21, "Test (Predicted)" = 24)) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4)),
      shape = ggplot2::guide_legend(override.aes = list(fill  = "grey50", size = 4))
    ) +
    ggplot2::labs(title = "LDA: Training + Test Samples", x = x_lab, y = y_lab,
                  fill = "Phenotype", shape = "Dataset") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "right")

  # Draw highlight ring around selected samples
  if (nrow(hl_df) > 0) {
    p <- p + ggplot2::geom_point(
      data        = hl_df,
      mapping     = ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
      inherit.aes = FALSE,
      colour      = "#f39c12",
      fill        = NA,
      shape       = 21,
      size        = 6.5,
      stroke      = 2
    )
  }

  p
}

# ---- Prediction probability bar chart ----

#' Stacked/grouped bar chart of posterior probabilities for test samples.
#'
#' @param pred_obj  list from predict_ecmexplainer()
#' @return ggplot object
plot_prediction_probs <- function(pred_obj,
                                  plot_df    = NULL,
                                  facet_cols = NULL) {
  .ensure_ggplot()

  # Use pre-filtered/annotated df if provided, otherwise use full report
  report <- if (!is.null(plot_df)) plot_df else pred_obj$final_report

  prob_cols <- setdiff(
    colnames(report),
    c("Sample_ID", "Predicted_Phenotype", facet_cols)
  )

  long_df <- tidyr::pivot_longer(
    report,
    cols      = tidyr::all_of(prob_cols),
    names_to  = "Phenotype",
    values_to = "Probability"
  )

  # Build palette from the actual phenotype levels in this dataset
  pal <- .pal_for_levels(unique(long_df$Phenotype))

  # Order samples by their top predicted phenotype so bars cluster nicely
  sample_order <- report$Sample_ID[
    order(report$Predicted_Phenotype, report$Sample_ID)
  ]
  long_df$Sample_ID <- factor(long_df$Sample_ID, levels = sample_order)

  n_samples  <- length(unique(long_df$Sample_ID))
  x_text_size <- ifelse(n_samples > 80, 5,
                  ifelse(n_samples > 40, 7,
                  ifelse(n_samples > 20, 9, 11)))

  p <- ggplot2::ggplot(long_df, ggplot2::aes(
    x    = Sample_ID,
    y    = Probability,
    fill = Phenotype
  )) +
    ggplot2::geom_col(position = "stack", colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::labs(
      title = "ECMExplainER Phenotype Probabilities",
      x     = "Sample",
      y     = "Posterior Probability (%)",
      fill  = "Phenotype"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                              size  = x_text_size),
      legend.position = "right"
    )

  # Faceting
  if (!is.null(facet_cols) && length(facet_cols) > 0) {
    # Keep only facet cols that actually exist in long_df
    valid_facets <- intersect(facet_cols, colnames(long_df))
    if (length(valid_facets) > 0) {
      facet_formula <- stats::as.formula(
        paste("~", paste(valid_facets, collapse = " + "))
      )
      p <- p + ggplot2::facet_wrap(
        facet_formula,
        scales = "free_x",
        ncol   = 2
      ) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
          strip.text       = ggplot2::element_text(face = "bold", size = 10),
          panel.spacing    = ggplot2::unit(0.8, "lines")
        )
    }
  }

  p
}

# ---- Training PCA coloured by arbitrary metadata column ----

#' Plot global training PCA coloured by any metadata column.
#'
#' @param train_obj   list from train_ecmexplainer()
#' @param metadata    data.frame with at least Sample_ID and colour_by column
#' @param colour_by   character; column name in metadata to use for colouring
#' @param pc_x        integer x PC
#' @param pc_y        integer y PC
#' @return ggplot object
plot_pca_training_meta <- function(train_obj, metadata, colour_by,
                                   pc_x = 1L, pc_y = 2L) {
  .ensure_ggplot()

  # Accept either integer (1) or string ("PC1") for pc_x / pc_y
  pc_x_col <- if (is.character(pc_x) && grepl("^PC", pc_x)) pc_x else paste0("PC", pc_x)
  pc_y_col <- if (is.character(pc_y) && grepl("^PC", pc_y)) pc_y else paste0("PC", pc_y)
  # Extract the numeric index for variance look-up
  pc_x_idx <- as.integer(sub("PC", "", pc_x_col))
  pc_y_idx <- as.integer(sub("PC", "", pc_y_col))

  gp   <- train_obj$global_pca
  pvar <- train_obj$pca_var

  df <- as.data.frame(gp$x)
  df$Sample_ID <- rownames(df)

  if (!colour_by %in% colnames(metadata))
    stop(sprintf("[Plots] Column '%s' not found in training metadata.", colour_by))

  df <- merge(df, metadata[, c("Sample_ID", colour_by)], by = "Sample_ID")
  colnames(df)[colnames(df) == colour_by] <- "ColourBy"

  x_lab <- paste0(pc_x_col, " (", pvar[pc_x_idx], "%)")
  y_lab <- paste0(pc_y_col, " (", pvar[pc_y_idx], "%)")

  # Use phenotype palette if colouring by Phenotype, otherwise a discrete palette
  if (colour_by == "Phenotype") {
    colour_scale <- ggplot2::scale_fill_manual(
      values   = .pal_for_levels(unique(df$ColourBy)),
      na.value = "grey70"
    )
  } else {
    colour_scale <- ggplot2::scale_fill_brewer(palette = "Set1", na.value = "grey70")
  }

  ggplot2::ggplot(df, ggplot2::aes(
    x     = .data[[pc_x_col]],
    y     = .data[[pc_y_col]],
    fill  = ColourBy,
    label = Sample_ID
  )) +
    ggplot2::geom_point(shape = 21, colour = "black", size = 3.5, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    colour_scale +
    ggplot2::labs(
      title = paste0("Training PCA – coloured by: ", colour_by),
      x     = x_lab,
      y     = y_lab,
      fill  = colour_by
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

# ---- Ratio PCA of user samples (optionally with training overlay) ----

#' PCA on the ECM ratio matrix of user samples, coloured by any metadata column.
#' Optionally overlays training samples (coloured by Phenotype).
#'
#' When \code{train_obj} is supplied the user samples are projected onto the
#' training \code{global_pca} space (which was fitted on batch-corrected
#' training ratios) rather than computing a new PCA from user data alone.
#' This ensures the PCA axes are on a batch-corrected reference frame.
#'
#' @param ratio_res       list from compute_ratios() on user data
#' @param user_meta       data.frame with Sample_ID and at least one other column
#' @param colour_by       character; column name in user_meta for colouring
#' @param train_ratio_mat optional ratios x training_samples matrix to overlay
#' @param train_meta      optional training metadata data.frame (needs Sample_ID + Phenotype)
#' @param show_training   logical; whether to show training samples
#' @param pc_x            integer x PC
#' @param pc_y            integer y PC
#' @param train_obj       optional list from train_ecmexplainer(); when supplied, user
#'                        samples are projected onto the pre-computed (batch-corrected)
#'                        \code{global_pca} rather than fitting a fresh PCA
#' @return ggplot object
plot_ratio_pca_user <- function(ratio_res,
                                user_meta,
                                colour_by,
                                train_ratio_mat = NULL,
                                train_meta      = NULL,
                                show_training   = FALSE,
                                pc_x = 1L,
                                pc_y = 2L,
                                train_obj       = NULL) {
  .ensure_ggplot()

  user_ratio_mat <- ratio_res$ratio_mat   # ratios x user_samples

  if (!colour_by %in% colnames(user_meta))
    stop(sprintf("[Plots] Column '%s' not found in user metadata.", colour_by))

  # ---- Determine PCA space ----
  # When train_obj is available, project user samples onto the training global PCA
  # (which was fitted on batch-corrected training ratios). This keeps the PCA
  # reference frame consistent with the corrected training space.
  # Otherwise fall back to a fresh prcomp() on user + optional training data.

  if (!is.null(train_obj) && !is.null(train_obj$global_pca)) {
    gp      <- train_obj$global_pca
    pca_var <- round(100 * gp$sdev^2 / sum(gp$sdev^2), 1)

    # Rows shared between user ratios and the training PCA features
    shared_ratios <- intersect(rownames(user_ratio_mat), rownames(gp$rotation))
    if (length(shared_ratios) == 0)
      stop("[Plots] No shared ratios between user data and training PCA loadings.")

    # Project user samples: stats::predict() handles centering/scaling automatically
    user_proj <- stats::predict(gp, newdata = t(user_ratio_mat[shared_ratios, , drop = FALSE]))
    coords    <- as.data.frame(user_proj)
    coords$Sample_ID <- rownames(coords)
    coords$Origin    <- "User"

    # Training samples — already stored in gp$x
    if (show_training) {
      train_coords <- as.data.frame(gp$x)
      train_coords$Sample_ID <- rownames(train_coords)
      train_coords$Origin    <- "Training"
      coords <- rbind(
        coords[, intersect(colnames(coords), colnames(train_coords)), drop = FALSE],
        train_coords[, intersect(colnames(coords), colnames(train_coords)), drop = FALSE]
      )
    }

    pc_x_col <- paste0("PC", pc_x)
    pc_y_col <- paste0("PC", pc_y)

    # Annotate user samples
    user_ann <- user_meta[, c("Sample_ID", colour_by), drop = FALSE]
    colnames(user_ann)[2] <- "ColourBy"
    user_ann$Set <- "User"
    user_df <- merge(
      coords[coords$Origin == "User", , drop = FALSE],
      user_ann, by = "Sample_ID"
    )

    if (show_training && !is.null(train_meta)) {
      tm_ids <- intersect(train_coords$Sample_ID, train_meta$Sample_ID)
      train_ann <- train_meta[train_meta$Sample_ID %in% tm_ids,
                              c("Sample_ID", "Phenotype"), drop = FALSE]
      colnames(train_ann)[2] <- "ColourBy"
      train_ann$Set <- "Training"
      t_df <- merge(
        train_coords[train_coords$Sample_ID %in% tm_ids, , drop = FALSE],
        train_ann, by = "Sample_ID"
      )
      plot_df <- rbind(user_df[, intersect(colnames(user_df), colnames(t_df))],
                       t_df[,    intersect(colnames(user_df), colnames(t_df))])
    } else {
      plot_df <- user_df
    }

  } else {
    # ---- Fallback: fresh prcomp on user (+ optionally training) data ----
    # Build combined matrix for PCA
    if (show_training && !is.null(train_ratio_mat)) {
      shared_ratios <- intersect(rownames(user_ratio_mat), rownames(train_ratio_mat))
      if (length(shared_ratios) == 0)
        stop("[Plots] No shared ratios between user and training data for PCA.")
      combined_mat <- cbind(
        train_ratio_mat[shared_ratios, , drop = FALSE],
        user_ratio_mat[shared_ratios,  , drop = FALSE]
      )
    } else {
      combined_mat <- user_ratio_mat
    }

    # ---- PCA ----
    var_mask    <- apply(combined_mat, 1, var, na.rm = TRUE) > 0
    pca_obj     <- prcomp(t(combined_mat[var_mask, , drop = FALSE]),
                          center = TRUE, scale. = TRUE)
    pca_var     <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 1)

    coords      <- as.data.frame(pca_obj$x)
    coords$Sample_ID <- rownames(coords)

    pc_x_col <- paste0("PC", pc_x)
    pc_y_col <- paste0("PC", pc_y)

    # Annotate user samples
    user_ann <- user_meta[, c("Sample_ID", colour_by), drop = FALSE]
    colnames(user_ann)[2] <- "ColourBy"
    user_ann$Set <- "User"
    user_df <- merge(
      coords[coords$Sample_ID %in% user_meta$Sample_ID, ],
      user_ann, by = "Sample_ID"
    )

    # Annotate training samples (if shown)
    if (show_training && !is.null(train_meta)) {
      train_ids <- colnames(train_ratio_mat)
      train_ann <- train_meta[train_meta$Sample_ID %in% train_ids,
                              c("Sample_ID", "Phenotype"), drop = FALSE]
      colnames(train_ann)[2] <- "ColourBy"
      train_ann$Set <- "Training"
      train_df <- merge(
        coords[coords$Sample_ID %in% train_ids, ],
        train_ann, by = "Sample_ID"
      )
      plot_df <- rbind(user_df, train_df)
    } else {
      plot_df <- user_df
    }
  }  # end if/else train_obj branch

  x_lab <- paste0(pc_x_col, " (", pca_var[pc_x], "%)")
  y_lab <- paste0(pc_y_col, " (", pca_var[pc_y], "%)")

  # Shape: user = circle, training = triangle
  shape_vals <- c("User" = 21, "Training" = 24)

  # Colour scale: phenotype palette if colouring by Phenotype, else brewer
  if (colour_by == "Phenotype") {
    colour_scale <- ggplot2::scale_fill_manual(
      values   = .pal_for_levels(unique(plot_df$ColourBy)),
      na.value = "grey70"
    )
  } else {
    n_levels <- length(unique(plot_df$ColourBy))
    if (n_levels <= 12) {
      colour_scale <- ggplot2::scale_fill_brewer(palette = "Set1", na.value = "grey70")
    } else {
      colour_scale <- ggplot2::scale_fill_viridis_d(na.value = "grey70")
    }
  }

  title_str <- paste0("Ratio PCA – user samples, coloured by: ", colour_by)
  if (show_training) title_str <- paste0(title_str, " (+ training overlay)")

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x     = .data[[pc_x_col]],
    y     = .data[[pc_y_col]],
    fill  = ColourBy,
    shape = Set,
    label = Sample_ID
  )) +
    ggplot2::geom_point(colour = "black", size = 3.5, alpha = 0.88) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    colour_scale +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4),
                                    title = colour_by),
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey50", size = 4))
    ) +
    ggplot2::labs(title = title_str, x = x_lab, y = y_lab, shape = "Dataset") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "right")

  if (show_training) {
    p <- p + ggplot2::scale_shape_manual(values = shape_vals)
  } else {
    p <- p + ggplot2::scale_shape_manual(values = c("User" = 21), guide = "none")
  }

  p
}
