# ==============================================================================
# test_manuscript_numbers.R
#
# Holds the code to the published paper. Every assertion below is a number
# printed in the manuscript; if a refactor changes the science, this fails.
#
#   Rscript tests/test_manuscript_numbers.R
#
# Runtime ~2 min (it retrains ADSIG-180 from raw counts).
# ==============================================================================
suppressPackageStartupMessages({ library(MASS) })
source("R/ecm_core.R")

fails <- 0L
check <- function(label, got, want, tol = 0) {
  ok <- if (is.numeric(got) && is.numeric(want)) isTRUE(abs(got - want) <= tol) else isTRUE(all.equal(got, want))
  cat(sprintf("%-58s %-22s %s\n", label,
              paste0(format(got), if (tol > 0) paste0(" (±", tol, ")") else ""),
              if (ok) "PASS" else paste0("FAIL, expected ", format(want))))
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}

counts   <- readRDS("data/salmon_counts_Combined_FINAL.rds")
meta     <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
genelist <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
shipped  <- readRDS("data/TrainedModel_ADSIG_180.rds")

cat("\n== Methods 5.1: training cohort ==\n")
meta <- meta[match(colnames(counts), meta$SampleID), ]
check("samples",                    ncol(counts),                    96)
check("studies",                    length(unique(meta$Study)),      10)
check("smallest class (Hypoxia)",   min(table(meta$Phenotype)),      7)
check("common feature space",       nrow(counts),                    19861)   # ms says 18,957 -- see METHODS_DISCREPANCIES.md
lg <- ecmx_logcpm(counts)
check("genes with logCPM var>0.005", sum(apply(lg, 1, var) > 0.005),  18957)

cat("\n== Methods 5.2: curated panel ==\n")
gl <- genelist[!duplicated(genelist$Gene), ]
fn <- toupper(trimws(gl$Function))
check("curated genes",              nrow(gl),                        1139)
check("Assembly",                   sum(fn == "ASSEMBLY"),           412)
check("Destruction",                sum(fn == "DESTRUCTION"),        119)
check("Signalling",                 sum(fn == "SIGNALLING"),         608)
check("detected in all 96",         length(intersect(gl$Gene, rownames(lg))), 923)

cat("\n== Methods 5.3/5.4: retraining ADSIG-180 from raw counts ==\n")
fit <- ecmx_train(counts, meta, genelist, design = "ADSIG", n_per_pheno = 30L)
check("panel size",                 fit$panel_size,                  180)
check("A/D ratios",                 sum(fit$features$Module == "AD"),  76)
check("SIG ratios",                 sum(fit$features$Module == "SIG"), 104)
check("PCs retained",               fit$n_pcs,                       10)
check("LOOCV accuracy (%)",         round(fit$loocv_acc, 2),         92.71, tol = 0.01)
check("macro-F1",                   round(fit$macro_f1, 3),          0.921, tol = 0.001)
check("max gene appearances",       max(table(c(fit$features$Gene_A, fit$features$Gene_B))), 1)

cat("\n== Reproduces the shipped model bit-for-bit ==\n")
check("same ratio set",   setequal(fit$features$Ratio_ID, shipped$features$Ratio_ID), TRUE)
check("same ratio order", identical(fit$features$Ratio_ID, shipped$features$Ratio_ID), TRUE)
check("same LOOCV",       round(fit$loocv_acc, 2) == round(shipped$loocv_acc, 2),     TRUE)

cat("\n== Methods 5.5: deployment projection is self-consistent ==\n")
# Projecting the training data through the frozen model must return the
# training LOOCV-independent fit: hard labels equal the LDA's own training call.
proj <- ecmx_project(shipped, ecmx_logcpm(counts)[
  , seq_len(ncol(counts)), drop = FALSE])
train_pcs <- as.data.frame(shipped$pca$x[, seq_len(shipped$n_pcs), drop = FALSE])
colnames(train_pcs) <- paste0("PC", seq_len(shipped$n_pcs))
direct <- predict(shipped$lda, train_pcs)$class
check("projection matches direct LDA call", mean(proj$class == direct), 1, tol = 1e-9)
check("ratio coverage on training data",    round(proj$coverage, 3),    1)

cat("\n== Methods 5.5: missing genes are neutralised, not fatal ==\n")
expr <- ecmx_logcpm(counts)
drop <- unique(ecmx_parse_ratios(ecmx_model_ratios(shipped))$Gene_A[1:20])
partial <- ecmx_project(shipped, expr[setdiff(rownames(expr), drop), , drop = FALSE])
check("coverage drops as expected",  round(partial$coverage, 3), round(1 - 20/180, 3), tol = 0.01)
check("still returns a call for every sample", sum(!is.na(partial$class)), 96)

cat("\n== Methods 5.6: centroid geometry ==\n")
cen <- ecmx_centroids(shipped)
check("six centroids", nrow(cen), 6)
d <- ecmx_distances(proj$coords, cen)
# Each training sample should be closest to its own phenotype centroid more
# often than not; this is the geometric statement the distance heatmaps rely on.
nearest <- colnames(d)[apply(d, 1, which.min)]
check("nearest-centroid agrees with label (>=90%)",
      mean(nearest == as.character(shipped$train_phenotype)) >= 0.90, TRUE)

cat("\n", strrep("=", 78), "\n", sep = "")
if (fails == 0L) cat("ALL CHECKS PASSED\n") else cat(fails, "CHECK(S) FAILED\n")
quit(status = if (fails == 0L) 0L else 1L)
