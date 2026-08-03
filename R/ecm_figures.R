# ==============================================================================
# ecm_figures.R — manuscript-style figures.
#
# Palette and geometry follow the published figures so app output is directly
# comparable to the paper. Pure functions: they take data, return a ggplot.
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2) })

ECMX_COLORS <- c(Fibrosis = "#377EB8", Hypoxia = "#E41A1C", Inflammation = "#4DAF4A",
                 Proliferating = "#FF7F00", Quiescence = "#A65628", Senescence = "#333333")

#' Colours for an arbitrary phenotype vocabulary: keep the manuscript's colours
#' for the six reference states, generate distinct ones for user-added labels.
ecmx_palette <- function(levels) {
  known <- intersect(levels, names(ECMX_COLORS))
  extra <- setdiff(levels, names(ECMX_COLORS))
  pal <- ECMX_COLORS[known]
  if (length(extra)) {
    hues <- grDevices::hcl(h = seq(15, 375, length.out = length(extra) + 1)[-1],
                           l = 55, c = 90)
    pal <- c(pal, setNames(hues, extra))
  }
  pal[levels]
}

ecmx_theme <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))
}

#' Figure 1-style LDA space: training samples with 85% ellipses, centroid stars,
#' and (optionally) projected new samples overlaid in the same frozen space.
ecmx_fig_lda <- function(train_coords, train_pheno, centroids,
                         new_coords = NULL, new_label = NULL,
                         dims = c(1, 2), title = "ECM phenotype space") {
  stopifnot(ncol(train_coords) >= max(dims))
  df <- data.frame(D1 = train_coords[, dims[1]], D2 = train_coords[, dims[2]],
                   Phenotype = factor(train_pheno))
  pal <- ecmx_palette(levels(df$Phenotype))
  cen <- data.frame(D1 = centroids[, dims[1]], D2 = centroids[, dims[2]],
                    Phenotype = rownames(centroids))
  ok <- names(which(table(df$Phenotype) >= 4))

  p <- ggplot(df, aes(D1, D2, colour = Phenotype)) +
    stat_ellipse(data = df[df$Phenotype %in% ok, ], aes(fill = Phenotype),
                 geom = "polygon", type = "norm", level = 0.85,
                 alpha = 0.12, colour = NA, show.legend = FALSE) +
    geom_point(size = 2.4, alpha = 0.75) +
    geom_point(data = cen, shape = 8, size = 4.5, stroke = 1.2, colour = "black") +
    scale_colour_manual(values = pal, drop = FALSE) +
    scale_fill_manual(values = pal, guide = "none", drop = FALSE) +
    ecmx_theme() +
    labs(title = title, x = paste0("LD", dims[1]), y = paste0("LD", dims[2]))

  if (!is.null(new_coords) && nrow(new_coords)) {
    nd <- data.frame(D1 = new_coords[, dims[1]], D2 = new_coords[, dims[2]])
    nd$Group <- if (is.null(new_label)) "Your samples" else new_label
    p <- p +
      geom_point(data = nd, aes(D1, D2, shape = Group), inherit.aes = FALSE,
                 size = 3.1, stroke = 1.1, colour = "black", fill = "white") +
      scale_shape_manual(values = rep(c(21, 24, 22, 23, 25), length.out =
                                        length(unique(nd$Group))), name = NULL)
  }
  p
}

#' Per-sample posterior probabilities as stacked bars, grouped by condition
#' (Methods 5.5).
ecmx_fig_posteriors <- function(posterior, sample_ids, group = NULL) {
  df <- data.frame(
    SampleID  = rep(sample_ids, times = ncol(posterior)),
    Phenotype = rep(colnames(posterior), each = nrow(posterior)),
    P         = as.vector(posterior),
    Group     = if (is.null(group)) "All" else rep(group, times = ncol(posterior)),
    stringsAsFactors = FALSE
  )
  df$Phenotype <- factor(df$Phenotype)
  ggplot(df, aes(SampleID, P, fill = Phenotype)) +
    geom_col(width = 0.9) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = ecmx_palette(levels(df$Phenotype))) +
    scale_y_continuous(expand = c(0, 0)) +
    ecmx_theme() +
    labs(title = "Posterior probability per sample", x = NULL, y = "Posterior") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          panel.spacing = unit(2, "pt"))
}

#' Six-column distance heatmap: median distance per condition to each training
#' centroid (Methods 5.6). Proximity summary only — direction is not encoded.
ecmx_fig_distance_heatmap <- function(distances, group) {
  agg <- aggregate(distances, by = list(Group = group), FUN = stats::median)
  m <- as.matrix(agg[, -1, drop = FALSE]); rownames(m) <- agg$Group
  df <- data.frame(
    Group     = rep(rownames(m), times = ncol(m)),
    Phenotype = rep(colnames(m), each  = nrow(m)),
    Distance  = as.vector(m), stringsAsFactors = FALSE
  )
  ggplot(df, aes(Phenotype, Group, fill = Distance)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", Distance)), size = 3.1,
              colour = ifelse(df$Distance > stats::median(df$Distance), "white", "black")) +
    scale_fill_gradientn(colours = c("#a50026", "#f7f7f7", "#313695")) +
    ecmx_theme() +
    labs(title = "Median distance to training centroids",
         x = NULL, y = NULL, fill = "Euclidean\ndistance") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

#' Protractor: angle between an observed trajectory and each phenotype
#' direction, with bootstrap confidence fans (Methods 5.6).
ecmx_fig_protractor <- function(angles, title = "Trajectory orientation") {
  a <- angles[!is.na(angles$Angle), , drop = FALSE]
  if (!nrow(a)) return(ggplot() + ecmx_theme() + labs(title = "Not enough samples"))
  a$r <- a$Distance / max(a$Distance, na.rm = TRUE)
  pal <- ecmx_palette(a$Phenotype)
  p <- ggplot(a) +
    geom_segment(aes(x = Angle, xend = Angle, y = 0, yend = r, colour = Phenotype),
                 linewidth = 1.4, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
    geom_point(aes(Angle, r, colour = Phenotype), size = 2.6) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, 30)) +
    coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    ecmx_theme() +
    labs(title = title, x = "Angle to trajectory (degrees)", y = "Relative distance")
  if (all(is.finite(a$Lower))) {
    p <- p + geom_rect(aes(xmin = Lower, xmax = Upper, ymin = 0, ymax = r,
                           fill = Phenotype), alpha = 0.15, colour = NA) +
      scale_fill_manual(values = pal, guide = "none")
  }
  p
}

#' Two-anchor trajectory plot: distance to centroid A against distance to
#' centroid B, with medoid points joined by arrows in the given order
#' (Methods 5.7). The gold triangles mark where the two training centroids sit.
ecmx_fig_trajectory <- function(distances, group, anchor_a, anchor_b,
                                order_levels = NULL, centroid_sep = NULL) {
  df <- data.frame(A = distances[, anchor_a], B = distances[, anchor_b],
                   Group = as.character(group), stringsAsFactors = FALSE)
  if (is.null(order_levels)) order_levels <- unique(df$Group)
  med <- do.call(rbind, lapply(order_levels, function(g) {
    sub <- df[df$Group == g, c("A", "B"), drop = FALSE]
    if (!nrow(sub)) return(NULL)
    i <- ecmx_medoid(as.matrix(sub))
    data.frame(Group = g, A = sub$A[i], B = sub$B[i], n = nrow(sub))
  }))
  med$Group <- factor(med$Group, levels = order_levels)
  med <- med[order(med$Group), ]

  p <- ggplot(df, aes(A, B)) +
    geom_point(colour = "grey72", size = 1.9, alpha = 0.6) +
    geom_path(data = med, aes(A, B), linewidth = 1, colour = "grey25",
              arrow = arrow(length = unit(8, "pt"), type = "closed")) +
    geom_point(data = med, aes(A, B, colour = Group), size = 4.2) +
    ggrepel::geom_text_repel(data = med, aes(A, B, label = paste0(Group, " (n=", n, ")")),
                             size = 3.2, seed = 1) +
    ecmx_theme() +
    labs(title = sprintf("Trajectory in the %s-%s plane", anchor_a, anchor_b),
         x = paste("Distance to", anchor_a), y = paste("Distance to", anchor_b))
  if (!is.null(centroid_sep)) {
    gold <- data.frame(A = c(0, centroid_sep), B = c(centroid_sep, 0))
    p <- p + geom_point(data = gold, aes(A, B), shape = 24, size = 3.6,
                        fill = "#DAA520", colour = "black", inherit.aes = FALSE)
  }
  p
}

#' ECM Aging Score distributions per group (Methods 5.9). Higher = older.
ecmx_fig_aging_score <- function(scores, group) {
  df <- data.frame(Score = scores, Group = factor(group))
  df$Group <- stats::reorder(df$Group, df$Score, FUN = stats::median)
  ggplot(df, aes(Group, Score, fill = Group)) +
    geom_violin(alpha = 0.35, colour = NA, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.09, size = 1.5, alpha = 0.7) +
    ecmx_theme() +
    labs(title = "ECM Aging Score", x = NULL,
         y = expression(S^sign~"(higher = older)")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1))
}

#' Mean sign per phenotype across the panel (the Figure 1 heatmap), rendered as
#' a ggplot so it composes with the rest of the app.
ecmx_fig_panel_heatmap <- function(sign_matrix, phenotype, features, max_rows = 60) {
  z <- t(apply(sign_matrix, 1, function(x) tapply(x, phenotype, mean)))
  sel <- features$Selecting_Phenotype[match(rownames(z), features$Ratio_ID)]
  ord <- order(factor(sel, levels = sort(unique(sel))))
  z <- z[ord, , drop = FALSE]; sel <- sel[ord]
  if (nrow(z) > max_rows) {                 # keep it legible in a browser
    keep <- round(seq(1, nrow(z), length.out = max_rows))
    z <- z[keep, , drop = FALSE]; sel <- sel[keep]
  }
  df <- data.frame(
    Ratio     = factor(rep(rownames(z), times = ncol(z)), levels = rev(rownames(z))),
    Phenotype = rep(colnames(z), each = nrow(z)),
    MeanSign  = as.vector(z), stringsAsFactors = FALSE
  )
  ggplot(df, aes(Phenotype, Ratio, fill = MeanSign)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#313695", "white", "#a50026"), limits = c(-1, 1)) +
    ecmx_theme(11) +
    labs(title = sprintf("Mean sign per phenotype (%d of %d ratios shown)",
                         nrow(z), nrow(sign_matrix)),
         x = NULL, y = NULL, fill = "Mean sign") +
    theme(axis.text.y = element_text(size = 5),
          axis.text.x = element_text(angle = 30, hjust = 1))
}
