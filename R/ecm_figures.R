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

#' Add to any figure drawn on a phone-width screen: legends move under the plot
#' (side by side, each one column) and text shrinks, so the panel keeps the width.
ecmx_compact <- function() {
  theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "vertical",
        legend.justification = "left", legend.box.just = "top",
        text = element_text(size = 10), plot.title = element_text(size = 12, face = "bold"),
        plot.title.position = "plot", plot.subtitle = element_text(size = 9.5),
        legend.title = element_text(size = 9.5, face = "bold"), legend.text = element_text(size = 8.5),
        legend.key.size = unit(11, "pt"), plot.margin = margin(4, 6, 4, 4))
}

ecmx_theme <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))
}

#' Facet strips that stay readable with many or long group names. Labels are
#' vertical and wrapped, because a facet is only as wide as its samples and a
#' horizontal label does not fit a small group. Above 20 groups no label can fit,
#' so strips are left blank and a caption says so. Add AFTER ecmx_theme(): a
#' complete theme added later would reset these settings.
ecmx_facet_strips <- function(n_groups) {
  many <- n_groups > 20
  list(facet_grid(~ Group, scales = "free_x", space = "free_x",
                  labeller = labeller(Group = label_wrap_gen(width = 22))),
       theme(strip.text = if (many) element_blank() else element_text(angle = 90, size = 9, hjust = 0),
             strip.background = if (many) element_blank() else element_rect(colour = "black"),
             strip.clip = "off", panel.spacing = unit(if (many) 0.5 else 2, "pt")),
       if (many) labs(caption = sprintf(
         "%d groups: too many to label. Group by a column with fewer values, or filter samples.", n_groups)))
}

#' Figure 1-style LDA space: training samples with 85% ellipses, centroid stars,
#' and (optionally) projected new samples overlaid in the same frozen space.
#'
#' The panel is drawn with a fixed 1:1 aspect ratio. That is not only cosmetic —
#' the LDA discriminant space is sphered, so a unit along LD1 means the same thing
#' as a unit along LD2, and letting the browser stretch the axes independently
#' misrepresents the distances the classifier actually uses.
#'
#' Projected samples are coloured by `new_label`, not merely shaped by it: with
#' more than about five groups, recycled point shapes are indistinguishable and
#' changing the grouping appears to do nothing.
ecmx_fig_lda <- function(train_coords, train_pheno, centroids,
                         new_coords = NULL, new_label = NULL, group_title = NULL,
                         dims = c(1, 2), title = "ECM phenotype space",
                         fixed_aspect = TRUE, show_groups = FALSE, point_scale = 1) {
  stopifnot(ncol(train_coords) >= max(dims))
  df <- data.frame(D1 = train_coords[, dims[1]], D2 = train_coords[, dims[2]],
                   Phenotype = factor(train_pheno))
  pal <- ecmx_palette(levels(df$Phenotype))
  cen <- data.frame(D1 = centroids[, dims[1]], D2 = centroids[, dims[2]],
                    Phenotype = rownames(centroids))
  ok <- names(which(table(df$Phenotype) >= 4))

  # The ellipses carry phenotype identity through their outline, not a fill
  # mapping. That keeps the `fill` aesthetic free for the projected-sample legend
  # below — a second scale_fill_manual would otherwise replace this one and
  # silently strip the ellipse colours.
  p <- ggplot(df, aes(D1, D2, colour = Phenotype)) +
    stat_ellipse(data = df[df$Phenotype %in% ok, ], geom = "polygon",
                 type = "norm", level = 0.85, fill = "grey70", alpha = 0.10,
                 linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.4 * point_scale, alpha = 0.7) +
    geom_point(data = cen, shape = 8, size = 4.5 * point_scale, stroke = 1.2 * point_scale, colour = "black") +
    scale_colour_manual(values = pal, drop = FALSE, name = "Training phenotype") +
    ecmx_theme() +
    labs(title = title, x = paste0("LD", dims[1]), y = paste0("LD", dims[2]))

  if (!is.null(new_coords) && nrow(new_coords)) {
    nd <- data.frame(D1 = new_coords[, dims[1]], D2 = new_coords[, dims[2]],
                     Group = if (is.null(new_label)) "Your samples"
                             else as.character(new_label))
    nd$Group <- factor(nd$Group)
    # A deliberately darker palette, so the overlay never reads as a training
    # phenotype even when the grouping column *is* Phenotype.
    gcols <- grDevices::hcl(
      h = seq(20, 380, length.out = nlevels(nd$Group) + 1)[-1], l = 42, c = 105)

    # Convex hull per group, drawn *under* the points. Without this, changing the
    # grouping only recolours the same markers — which on a cohort where the
    # grouping columns are confounded (studies that are each mostly one
    # phenotype) is indistinguishable from a legend change. The hulls make the
    # partition itself visible.
    if (show_groups) {
      hulls <- do.call(rbind, lapply(split(nd, nd$Group), function(s) {
        if (nrow(s) < 3) return(NULL)
        s[grDevices::chull(s$D1, s$D2), , drop = FALSE]
      }))
      if (!is.null(hulls) && nrow(hulls))
        p <- p + geom_polygon(data = hulls, aes(D1, D2, fill = Group),
                              inherit.aes = FALSE, alpha = 0.16, colour = NA)
    }

    p <- p +
      geom_point(data = nd, aes(D1, D2, fill = Group), inherit.aes = FALSE,
                 shape = 21, size = 3.4 * point_scale, stroke = 0.9 * point_scale, colour = "black") +
      scale_fill_manual(values = gcols,
                        name = if (is.null(group_title)) "Your samples" else group_title,
                        guide = guide_legend(override.aes = list(shape = 21, size = 3.6, alpha = 1)))
  }
  if (fixed_aspect) p <- p + coord_fixed(ratio = 1)
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
  # Within each group, order samples by their most likely phenotype and its
  # posterior, so blocks of the same call sit together.
  top <- colnames(posterior)[max.col(posterior, ties.method = "first")]
  ord <- order(if (is.null(group)) rep("All", length(sample_ids)) else group, top, -apply(posterior, 1, max))
  df$SampleID <- factor(df$SampleID, levels = unique(sample_ids[ord]))
  n <- length(sample_ids); show_names <- n <= 50
  ggplot(df, aes(SampleID, P, fill = Phenotype)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = ecmx_palette(levels(df$Phenotype))) +
    scale_y_continuous(expand = c(0, 0)) +
    ecmx_theme() +
    labs(title = "Posterior probability per sample", x = NULL, y = "Posterior",
         caption = if (!show_names) sprintf("%d samples; names omitted for legibility (see the Results table).", n)) +
    ecmx_facet_strips(length(unique(df$Group))) +
    theme(axis.text.x = if (show_names) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6) else element_blank(),
          axis.ticks.x = if (show_names) element_line() else element_blank(),
          legend.position = "bottom")
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
                        fill = "#DAA520", colour = "black", inherit.aes = FALSE) +
      labs(caption = sprintf("Gold triangles = exact positions of the %s and %s training centroids. Grey = individual samples; large points = condition medoids.",
                             anchor_a, anchor_b))
  }
  p
}

#' ECM Aging Score distributions per group (Methods 5.9). Higher = older.
#' @param reference optional data frame (Bin, Median, Q25, Q75): the GSE113957
#'   age bands, drawn as dashed medians with the Young interquartile range shaded.
ecmx_fig_aging_score <- function(scores, group, reference = NULL) {
  df <- data.frame(Score = scores, Group = factor(group))
  df$Group <- stats::reorder(df$Group, df$Score, FUN = stats::median)
  p <- ggplot(df, aes(Group, Score, fill = Group))
  if (!is.null(reference)) {
    young <- reference[grepl("^Young", reference$Bin), ]
    if (nrow(young)) p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = young$Q25,
                                       ymax = young$Q75, fill = "#3b4cca", alpha = 0.08)
    p <- p + geom_hline(data = reference, aes(yintercept = Median), linetype = "dashed",
                        colour = "grey45", linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = Inf, y = reference$Median, label = reference$Bin, hjust = -0.05,
               size = 3, colour = "grey30") +
      coord_cartesian(clip = "off") +
      labs(caption = "Dashed lines: median score of healthy dermal fibroblasts per age group (GSE113957); shaded: Young interquartile range.")
  }
  p +
    geom_violin(alpha = 0.35, colour = NA, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.09, size = 1.5, alpha = 0.7) +
    ecmx_theme() +
    labs(title = "ECM Aging Score", x = NULL,
         y = expression(S^sign~"(higher = older)")) +
    theme(legend.position = "none",
          plot.margin = margin(5.5, if (is.null(reference)) 5.5 else 105, 5.5, 5.5),
          axis.text.x = element_text(angle = 25, hjust = 1))
}

#' Mean sign per phenotype across the panel (the Figure 1 heatmap), rendered as
#' a ggplot so it composes with the rest of the app.
#' @param selecting one or more phenotypes: show only the ratios selected *for*
#'   those phenotypes. Filtering by phenotype shows all of that phenotype's
#'   ratios, which is the point — the whole 180-row panel cannot be labelled
#'   legibly in a browser, so an unfiltered view is necessarily thinned.
#' @param max_rows thinning cap for the unfiltered view; set NULL for no cap.
ecmx_fig_panel_heatmap <- function(sign_matrix, phenotype, features,
                                   selecting = NULL, max_rows = 60) {
  z <- t(apply(sign_matrix, 1, function(x) tapply(x, phenotype, mean)))
  sel <- features$Selecting_Phenotype[match(rownames(z), features$Ratio_ID)]

  if (!is.null(selecting) && length(selecting)) {
    keep <- which(sel %in% selecting)
    if (!length(keep))
      return(ggplot() + ecmx_theme() +
               labs(title = sprintf("No ratios select for %s",
                                    paste(selecting, collapse = ", "))))
    z <- z[keep, , drop = FALSE]; sel <- sel[keep]
    max_rows <- NULL                        # a single phenotype fits; show it all
  }

  ord <- order(factor(sel, levels = sort(unique(sel))))
  z <- z[ord, , drop = FALSE]; sel <- sel[ord]

  thinned <- FALSE
  if (!is.null(max_rows) && nrow(z) > max_rows) {
    keep <- round(seq(1, nrow(z), length.out = max_rows))
    z <- z[keep, , drop = FALSE]; sel <- sel[keep]
    thinned <- TRUE
  }

  df <- data.frame(
    Ratio     = factor(rep(rownames(z), times = ncol(z)), levels = rev(rownames(z))),
    Phenotype = rep(colnames(z), each = nrow(z)),
    MeanSign  = as.vector(z), stringsAsFactors = FALSE
  )
  subtitle <- if (thinned)
    sprintf("%d of %d ratios shown (evenly thinned for legibility)", nrow(z), nrow(sign_matrix))
  else
    sprintf("all %d ratios%s", nrow(z),
            if (is.null(selecting)) "" else paste0(" selecting for ", paste(selecting, collapse = ", ")))

  ggplot(df, aes(Phenotype, Ratio, fill = MeanSign)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#313695", "white", "#a50026"), limits = c(-1, 1)) +
    ecmx_theme(11) +
    labs(title = "Mean sign per phenotype", subtitle = subtitle,
         x = NULL, y = NULL, fill = "Mean sign") +
    theme(axis.text.y = element_text(size = if (nrow(z) > 40) 7 else 10),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey35"))
}

# ---- panel comparison and validation figures (Figure 7, SI Table 6) ----------

# App palette: ECMExplainer pink, controls blue and grey.
ECMX_PANEL_COLORS <- c("ADSIG-180 (ECMExplainer)" = "#d6336c", "POOL-180" = "#8c96a0",
                       "ALLGENE-180" = "#3b4cca")
ecmx_panel_palette <- function(labels) {
  known <- intersect(labels, names(ECMX_PANEL_COLORS))
  extra <- setdiff(labels, names(ECMX_PANEL_COLORS))
  pal <- ECMX_PANEL_COLORS[known]
  if (length(extra))
    pal <- c(pal, setNames(grDevices::hcl(seq(150, 330, length.out = length(extra) + 1)[-1],
                                          l = 50, c = 80), extra))
  pal[labels]
}

#' Figure 7B: recall per truth phenotype for each panel, as a radar chart.
#' @param recall data frame with Panel, Phenotype, Recall_pct
#' Drawn in Cartesian coordinates so the edges between spokes are straight.
ecmx_fig_recall_radar <- function(recall) {
  recall$Panel <- factor(recall$Panel, levels = unique(recall$Panel))
  ph <- sort(unique(recall$Phenotype)); k <- length(ph)
  ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k
  xy <- function(r, i) data.frame(x = r * cos(ang[i]), y = r * sin(ang[i]))
  recall$i <- match(recall$Phenotype, ph)
  pts <- cbind(recall, xy(recall$Recall_pct, recall$i))
  pts <- pts[order(pts$Panel, pts$i), ]
  rings <- do.call(rbind, lapply(c(25, 50, 75, 100), function(r)
    cbind(xy(r, c(seq_len(k), 1)), r = r)))
  spokes <- data.frame(x = 0, y = 0, xend = 100 * cos(ang), yend = 100 * sin(ang))
  lab <- cbind(xy(112, seq_len(k)), Phenotype = ph)
  # Anchor side labels outward so they never sit on the 100% point.
  lab$hjust <- ifelse(cos(ang) > 0.3, 0, ifelse(cos(ang) < -0.3, 1, 0.5))
  ggplot() +
    geom_path(data = rings, aes(x, y, group = r), colour = "grey85") +
    geom_segment(data = spokes, aes(x, y, xend = xend, yend = yend), colour = "grey85") +
    annotate("text", x = 3, y = c(25, 50, 75, 100), label = c("25", "50", "75", "100%"),
             hjust = 0, size = 2.8, colour = "grey50") +
    geom_polygon(data = pts, aes(x, y, colour = Panel, fill = Panel, group = Panel),
                 alpha = 0.08, linewidth = 0.9) +
    geom_point(data = pts, aes(x, y, colour = Panel), size = 2.4) +
    geom_text(data = lab, aes(x, y, label = Phenotype, hjust = hjust), size = 3.6) +
    scale_colour_manual(values = ecmx_panel_palette(levels(recall$Panel)), name = NULL) +
    scale_fill_manual(values = ecmx_panel_palette(levels(recall$Panel)), name = NULL) +
    coord_fixed(xlim = c(-165, 165), ylim = c(-135, 135), clip = "off") +
    theme_void(base_size = 12) +
    labs(title = "Recall per phenotype") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
}

#' Accuracy per panel with Wilson 95% intervals, micro and balanced side by side.
#' @param summary data frame with Panel, Micro_accuracy_pct, Micro_CI95_low,
#'   Micro_CI95_high, Balanced_accuracy_pct, N_scored
ecmx_fig_accuracy <- function(summary) {
  summary$Panel <- factor(summary$Panel, levels = rev(unique(summary$Panel)))
  ggplot(summary, aes(y = Panel)) +
    geom_errorbar(aes(xmin = Micro_CI95_low, xmax = Micro_CI95_high), width = 0.18,
                  orientation = "y", colour = "grey45") +
    geom_point(aes(x = Micro_accuracy_pct, colour = "Micro accuracy (95% CI)"), size = 3.4) +
    geom_point(aes(x = Balanced_accuracy_pct, colour = "Balanced accuracy"), size = 3.4, shape = 18) +
    scale_colour_manual(values = c("Micro accuracy (95% CI)" = "#3b4cca",
                                   "Balanced accuracy" = "#d6336c"), name = NULL) +
    scale_x_continuous(limits = c(0, 100)) +
    ecmx_theme(12) +
    labs(title = sprintf("Accuracy on %d labelled samples", max(summary$N_scored)),
         x = "%", y = NULL) +
    theme(legend.position = "bottom")
}

#' Truth x prediction counts, one facet per panel (SI Table 6 sheet F).
#' @param conf data frame with Panel, Truth, Predicted, Freq
ecmx_fig_confusion <- function(conf) {
  conf$Panel <- factor(conf$Panel, levels = unique(conf$Panel))
  has_truth <- tapply(conf$Freq, conf$Truth, sum) > 0       # drop label rows nobody carries
  conf <- conf[conf$Truth %in% names(has_truth)[has_truth], , drop = FALSE]
  conf <- do.call(rbind, lapply(split(conf, list(conf$Panel, conf$Truth), drop = TRUE), function(d) {
    d$Share <- if (sum(d$Freq)) d$Freq / sum(d$Freq) else 0; d }))
  ggplot(conf, aes(Predicted, Truth, fill = Share)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 3.3,
              colour = ifelse(conf$Share > 0.55, "white", "black")) +
    facet_wrap(~ Panel) +
    scale_fill_gradient(low = "#fbf3f7", high = "#a61e4d", limits = c(0, 1),
                        name = "Share of\ntruth row") +
    ecmx_theme(11) +
    labs(title = "Confusion matrices", x = "Predicted", y = "Known phenotype") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

#' Accuracy within each dataset for each panel (SI Table 6 sheet D).
#' @param per_ds data frame with Panel, Dataset, N, Accuracy_pct
ecmx_fig_dataset_accuracy <- function(per_ds) {
  per_ds$Panel <- factor(per_ds$Panel, levels = unique(per_ds$Panel))
  per_ds$Label <- sprintf("%s (n=%d)", per_ds$Dataset, per_ds$N)
  per_ds$Label <- factor(per_ds$Label, levels = rev(sort(unique(per_ds$Label))))
  ggplot(per_ds, aes(Panel, Label, fill = Accuracy_pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", Accuracy_pct)), size = 3.2,
              colour = ifelse(per_ds$Accuracy_pct > 60, "white", "black")) +
    scale_fill_gradient(low = "#fbf3f7", high = "#a61e4d", limits = c(0, 100), name = "Accuracy\n(%)") +
    ecmx_theme(11) +
    labs(title = "Accuracy per dataset", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

#' Predicted phenotype for every sample under every panel, so disagreements
#' between panels are visible sample by sample.
#' @param preds data frame with SampleID, Group, Panel, Predicted
ecmx_fig_panel_agreement <- function(preds) {
  preds$Panel <- factor(preds$Panel, levels = rev(unique(preds$Panel)))
  ph <- sort(unique(preds$Predicted))
  ggplot(preds, aes(SampleID, Panel, fill = Predicted)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = ecmx_palette(ph), name = "Predicted") +
    ecmx_theme(11) +
    labs(title = "Prediction per sample and panel", x = NULL, y = NULL) +
    ecmx_facet_strips(length(unique(preds$Group))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")
}
