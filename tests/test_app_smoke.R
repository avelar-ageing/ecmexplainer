# ==============================================================================
# test_app_smoke.R — exercises everything the app does, without a browser.
#
#   Rscript tests/test_app_smoke.R
#
# Renders each figure to a temporary PNG so ggplot build errors surface here
# rather than as a blank panel in the UI.
# ==============================================================================
suppressPackageStartupMessages({ library(MASS); library(ggplot2) })
source("R/ecm_core.R")
source("R/ecm_figures.R")

fails <- 0L
step <- function(label, expr) {
  res <- tryCatch({ force(expr); "PASS" },
                  error = function(e) paste0("FAIL: ", conditionMessage(e)))
  cat(sprintf("%-52s %s\n", label, res))
  if (!identical(res, "PASS")) fails <<- fails + 1L
  invisible(res)
}
render <- function(p) { f <- tempfile(fileext = ".png")
  suppressWarnings(ggsave(f, p, width = 8, height = 5, dpi = 72)); unlink(f); TRUE }

MODEL <- readRDS("data/TrainedModel_ADSIG_180.rds")
CEN   <- ecmx_centroids(MODEL)
TRC   <- local({
  pcs <- as.data.frame(MODEL$pca$x[, seq_len(MODEL$n_pcs), drop = FALSE])
  colnames(pcs) <- paste0("PC", seq_len(MODEL$n_pcs))
  predict(MODEL$lda, pcs)$x
})

cat("\n== data loading ==\n")
step("read gzipped demo counts", { m <<- ecmx_read_expression("data/demo_counts.csv.gz"); stopifnot(ncol(m) == 24) })
step("demo metadata",            { md <<- read.csv("data/demo_metadata.csv", stringsAsFactors = FALSE) })
step("counts detected as unlogged", stopifnot(!ecmx_looks_logged(m)))
step("log2-CPM transform",       { lg <<- ecmx_logcpm(m); stopifnot(ecmx_looks_logged(lg)) })

cat("\n== projection ==\n")
step("project through frozen model", { pr <<- ecmx_project(MODEL, lg) })
step("full ratio coverage",      stopifnot(pr$coverage == 1))
step("one call per sample",      stopifnot(length(pr$class) == 24))
step("posteriors sum to 1",      stopifnot(all(abs(rowSums(pr$posterior) - 1) < 1e-8)))
step("distances to 6 centroids", { d <<- ecmx_distances(pr$coords, CEN); stopifnot(dim(d) == c(24, 6)) })

cat("\n== robustness ==\n")
step("unknown genes are ignored", {
  x <- lg; rownames(x)[1:50] <- paste0("NOTAGENE", 1:50)
  p <- ecmx_project(MODEL, x); stopifnot(length(p$class) == 24)
})
step("lower-case symbols still match", {
  x <- lg; rownames(x) <- tolower(rownames(x))
  p <- ecmx_project(MODEL, x); stopifnot(p$coverage == 1)
})
step("half the panel missing is survivable", {
  keep <- setdiff(rownames(lg), ecmx_parse_ratios(ecmx_model_ratios(MODEL))$Gene_A)
  p <- ecmx_project(MODEL, lg[keep, , drop = FALSE]); stopifnot(p$coverage == 0)
})
step("single sample projects", {
  p <- ecmx_project(MODEL, lg[, 1, drop = FALSE]); stopifnot(length(p$class) == 1)
})

cat("\n== figures ==\n")
g <- md$Phenotype[match(colnames(lg), md$SampleID)]
step("LDA phenotype space", render(ecmx_fig_lda(TRC, MODEL$train_phenotype, CEN, pr$coords, g)))
step("posterior bars",      render(ecmx_fig_posteriors(pr$posterior, colnames(lg), g)))
step("distance heatmap",    render(ecmx_fig_distance_heatmap(d, g)))
step("panel heatmap (all, thinned)", render(ecmx_fig_panel_heatmap(
  MODEL$sign_matrix_train, MODEL$train_phenotype, MODEL$features)))
step("panel heatmap (one phenotype, unthinned)", {
  f <- MODEL$features
  ph <- names(sort(table(f$Selecting_Phenotype), decreasing = TRUE))[1]
  p <- ecmx_fig_panel_heatmap(MODEL$sign_matrix_train, MODEL$train_phenotype, f,
                              selecting = ph)
  # every ratio for that phenotype must be present, not a thinned subset
  n_rows <- length(unique(ggplot2::ggplot_build(p)$data[[1]]$y))
  stopifnot(n_rows == sum(f$Selecting_Phenotype == ph))
  render(p)
})
step("panel heatmap (unknown phenotype degrades gracefully)",
     render(ecmx_fig_panel_heatmap(MODEL$sign_matrix_train, MODEL$train_phenotype,
                                   MODEL$features, selecting = "NoSuchState")))
step("LDA overlay differs when grouping changes", {
  a <- ecmx_fig_lda(TRC, MODEL$train_phenotype, CEN, pr$coords, g)
  b <- ecmx_fig_lda(TRC, MODEL$train_phenotype, CEN, pr$coords,
                    md$Study[match(colnames(lg), md$SampleID)])
  fa <- ggplot2::ggplot_build(a)$data; fb <- ggplot2::ggplot_build(b)$data
  # the overlay layer is the last one; its fill must actually change
  stopifnot(!identical(fa[[length(fa)]]$fill, fb[[length(fb)]]$fill))
  render(a); render(b)
})
step("LDA panel is drawn 1:1 (not stretched)", {
  # ggplot2 3.x gives coord_fixed() the class CoordFixed; 4.x returns a
  # CoordCartesian carrying a `ratio` field. Assert the field, not the class.
  co <- ggplot2::ggplot_build(
    ecmx_fig_lda(TRC, MODEL$train_phenotype, CEN, pr$coords, g))$layout$coord
  stopifnot(isTRUE(all.equal(co$ratio, 1)))
  co2 <- ggplot2::ggplot_build(
    ecmx_fig_lda(TRC, MODEL$train_phenotype, CEN, pr$coords, g,
                 fixed_aspect = FALSE))$layout$coord
  stopifnot(is.null(co2$ratio))          # opt-out still works
})
step("trajectory plot",     render(ecmx_fig_trajectory(d, g, "Proliferating", "Senescence",
  centroid_sep = sqrt(sum((CEN["Proliferating", ] - CEN["Senescence", ])^2)))))
step("protractor + bootstrap", {
  a <- ecmx_trajectory_angles(pr$coords[g == "Proliferating", , drop = FALSE],
                              pr$coords[g == "Senescence", , drop = FALSE],
                              CEN, n_boot = 200)
  stopifnot(nrow(a) == 6, all(a$Angle >= 0 & a$Angle <= 180))
  render(ecmx_fig_protractor(a))
})

cat("\n== aging score ==\n")
step("fit betas and score", {
  set.seed(1); ages <- seq(20, 80, length.out = ncol(lg))
  b <- ecmx_fit_age_betas(lg, ages, ecmx_model_ratios(MODEL), fdr = 0.5)
  stopifnot(!is.null(b), nrow(b) > 0)
  s <- ecmx_aging_score(lg, b); stopifnot(length(s) == ncol(lg), all(is.finite(s)))
  render(ecmx_fig_aging_score(s, g))
})

cat("\n== retraining with user-defined labels ==\n")
step("train on arbitrary 3-class labels", {
  cnt <- readRDS("data/salmon_counts_Combined_FINAL.rds")
  mt  <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
  mt  <- mt[match(colnames(cnt), mt$SampleID), ]
  # Merge the six reference states into three user-named ones: arbitrary labels,
  # but real biological signal, which is the realistic user scenario.
  mt$Phenotype <- c(Proliferating = "Growing", Quiescence = "Growing",
                    Senescence = "Arrested", Fibrosis = "Arrested",
                    Inflammation = "Stressed", Hypoxia = "Stressed")[mt$Phenotype]
  gl  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
  f <- ecmx_train(cnt, mt, gl, design = "ADSIG", n_per_pheno = 5L)
  stopifnot(setequal(f$phenotypes, c("Growing", "Arrested", "Stressed")), f$panel_size == 15)
  p <- ecmx_project(f, lg); stopifnot(length(p$class) == 24)
})
step("eta2 ranking produces a different panel", {
  cnt <- readRDS("data/salmon_counts_Combined_FINAL.rds")
  mt  <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
  gl  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
  a <- ecmx_train(cnt, mt, gl, n_per_pheno = 6L, rank_by = "tsp_delta")
  b <- ecmx_train(cnt, mt, gl, n_per_pheno = 6L, rank_by = "eta2")
  stopifnot(!identical(sort(a$features$Ratio_ID), sort(b$features$Ratio_ID)))
  cat("     (tsp_delta vs eta2 overlap: ",
      length(intersect(a$features$Ratio_ID, b$features$Ratio_ID)), "/36)\n", sep = "")
})
step("min_per_module floor is honoured", {
  cnt <- readRDS("data/salmon_counts_Combined_FINAL.rds")
  mt  <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
  gl  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
  f <- ecmx_train(cnt, mt, gl, n_per_pheno = 6L, min_per_module = 5L)
  stopifnot(sum(f$features$Module == "AD") >= 5, sum(f$features$Module == "SIG") >= 5)
})

cat("\n== guard rails ==\n")
step("single-sample class is rejected", {
  cnt <- readRDS("data/salmon_counts_Combined_FINAL.rds")
  mt  <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
  mt  <- mt[match(colnames(cnt), mt$SampleID), ]; mt$Phenotype[1] <- "Loner"
  gl  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
  r <- tryCatch({ ecmx_train(cnt, mt, gl, n_per_pheno = 3L); "no error" },
                error = function(e) "errored")
  stopifnot(identical(r, "errored"))
})

cat("\n", strrep("=", 70), "\n", sep = "")
if (fails == 0L) cat("ALL SMOKE TESTS PASSED\n") else cat(fails, "FAILED\n")
quit(status = if (fails == 0L) 0L else 1L)
