# ==============================================================================
# test_app_server.R — drives the real Shiny server headlessly (shiny::testServer)
# through the paths a user takes: the validation cohorts, filters, accuracy,
# every figure, an upload, and retraining with a custom gene list.
#
#   Rscript tests/test_app_server.R        # ~1 min
# ==============================================================================
suppressPackageStartupMessages({ library(shiny); library(ggplot2) })
Sys.setenv(ECMX_WORKERS = "1")

fails <- 0L
check <- function(label, ok) {
  cat(sprintf("%-62s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) fails <<- fails + 1L
}
renders <- function(p) {
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({ suppressWarnings(ggsave(f, p, width = 9, height = 7, dpi = 60)); TRUE },
                 error = function(e) { message(conditionMessage(e)); FALSE })
  unlink(f); ok
}
upload <- function(path) data.frame(name = basename(path), size = file.size(path),
                                    type = "text/csv", datapath = path)

app <- shinyAppDir(".")
testServer(app, {
  cat("\n== Classify: manuscript validation cohorts ==\n")
  session$setInputs(src = "validation", panels = c("ADSIG180", "POOL180", "ALLGENE180"),
                    primary = "ADSIG180", group_col = "Condition",
                    truth_col = "Truth_Phenotype", ds_col = "Dataset", viz = "lda")
  check("327 samples scored by 3 panels",
        length(scored()) == 3 && ncol(expr_mat()) == 327)
  a <- acc_tables()$summary
  check(sprintf("ADSIG-180 balanced %.1f%% / micro %.1f%% (SI 72.7 / 83.6)",
                a$Balanced_accuracy_pct[1], a$Micro_accuracy_pct[1]),
        round(a$Balanced_accuracy_pct[1], 1) == 72.7 && round(a$Micro_accuracy_pct[1], 1) == 83.6)
  check("ALLGENE-180 micro 62.6% (SI)", round(a$Micro_accuracy_pct[3], 1) == 62.6)
  check("McNemar ADSIG vs ALLGENE p = 2.18e-06 (SI)",
        signif(acc_tables()$mcnemar$McNemar_exact_p[2], 3) == 2.18e-06)
  for (f in c("fig_lda", "fig_post", "fig_dist", "fig_heat", "fig_cmp",
              "fig_acc_radar", "fig_acc_conf", "fig_acc_ds"))
    check(paste("renders", f), renders(get(f)()))
  check("results table has every panel's columns",
        all(c("Pred_ADSIG180", "Mahal_POOL180", "OOD95_ALLGENE180") %in% names(results_tbl())))

  cat("\n== Filters ==\n")
  session$setInputs(cf_col1 = "Dataset", cf_val1 = c("GSE273882", "GSE302114"))
  check("filter to two datasets leaves 14 samples", ncol(expr_mat()) == 14)
  check("both datasets 100% accurate under ADSIG-180",
        acc_tables()$summary$Micro_accuracy_pct[1] == 100)
  session$setInputs(cf_col1 = "", cf_val1 = NULL)
  check("clearing the filter restores 327", ncol(expr_mat()) == 327)

  cat("\n== Trajectories ==\n")
  session$setInputs(tr_group_col = "Condition", anchor_a = "Proliferating",
                    anchor_b = "Fibrosis", prot_from = "Untreated", prot_to = "TGFb_Treated",
                    n_boot = 200)
  check("protractor angles computed", nrow(angles_tbl()) == 6)

  cat("\n== Upload ==\n")
  session$setInputs(src = "upload", expr_file = upload("data/examples/example_expression.csv"),
                    meta_file = upload("data/examples/example_metadata.csv"),
                    truth_col = "Truth_Phenotype", ds_col = "Dataset")
  check("example upload: 14 samples", ncol(expr_mat()) == 14)
  check("example upload reproduces the bundled predictions",
        acc_tables()$summary$Micro_accuracy_pct[1] == 100)
})

cat("\n== Retrain (synchronous core path) ==\n")
source("R/ecm_core.R")
cnt <- readRDS("data/salmon_counts_Combined_FINAL.rds")
md  <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
gl  <- data.frame(Gene = c("COL1A1", "COL1A2", "COL3A1", "FN1", "MMP1", "MMP2", "MMP3", "TIMP1",
                           "TIMP3", "SERPINE1", "CTGF", "CCN2", "IL6", "CXCL8", "VEGFA", "LOX",
                           "POSTN", "SPARC", "BGN", "DCN", "LUM", "TGFB1", "TGFBI", "THBS1",
                           "CDKN1A", "MKI67", "TOP2A", "CCNB1", "HMOX1", "NDRG1"))
fit <- run_training <- NULL
t0 <- Sys.time()
prep <- ecmx_prepare(cnt, md, gl, design = "POOL")
fit <- ecmx_train(design = "POOL", n_per_pheno = 3L, prep = prep)
check(sprintf("custom 30-gene POOL panel trains (%d ratios, LOOCV %.1f%%)",
              fit$panel_size, fit$loocv_acc), fit$panel_size >= 6)
lo <- ecmx_loso(prep, n_per_pheno = 3L)
check(sprintf("nested LOSO runs on it (%.1f%%)", lo$accuracy), is.finite(lo$accuracy))
sub <- md[md$Phenotype %in% c("Proliferating", "Senescence", "Quiescence"), ]
fit2 <- ecmx_train(cnt[, sub$SampleID], sub, design = "ALLGENE", n_per_pheno = 10L, top_n = 1000L)
check(sprintf("ALLGENE on a 3-class subset with custom labels (%d ratios)", fit2$panel_size),
      fit2$panel_size == 30 && identical(fit2$phenotypes, c("Proliferating", "Quiescence", "Senescence")))
p <- ecmx_project(fit2, readRDS("data/validation_expr.rds"))
check("a retrained model scores new samples", length(p$class) == 327)
cat(sprintf("retrain checks took %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n", strrep("=", 78), "\n", sep = "")
if (fails == 0L) cat("ALL CHECKS PASSED\n") else cat(fails, "CHECK(S) FAILED\n")
quit(status = if (fails == 0L) 0L else 1L)
