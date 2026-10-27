# ==============================================================================
# test_validation.R
#
# Scores the 327 bundled validation samples (data/validation_expr.rds) through
# the app's projection code with all three bundled panels, and asserts that
# every prediction, nearest centroid and OOD flag, and the headline accuracies,
# equal Supplementary Table 6.
#
#   Rscript tests/test_validation.R        # ~5 s
# ==============================================================================
suppressPackageStartupMessages({ library(MASS); library(readxl) })
source("R/ecm_core.R")

fails <- 0L
check <- function(label, ok) {
  cat(sprintf("%-66s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) fails <<- fails + 1L
}

expr <- readRDS("data/validation_expr.rds")
meta <- read.csv("data/validation_metadata.csv", stringsAsFactors = FALSE)
si6  <- "supplementary/SI_Table6.xlsx"
A <- as.data.frame(read_excel(si6, "A_MasterWide_AllSamples", skip = 2))
O <- as.data.frame(read_excel(si6, "O_OOD_and_Confidence", skip = 2))
B <- as.data.frame(read_excel(si6, "B_PanelSummary", skip = 2))
C <- as.data.frame(read_excel(si6, "C_PerPhenotypeRecall", skip = 2))
stopifnot(identical(colnames(expr), A$SampleID), identical(meta$SampleID, A$SampleID))

panels <- c(ADSIG180 = "ADSIG", POOL180 = "POOL", ALLGENE180 = "ALLGENE")
for (tag in names(panels)) {
  m  <- readRDS(sprintf("data/TrainedModel_%s_180.rds", panels[[tag]]))
  p  <- ecmx_project(m, expr)
  cf <- ecmx_confidence(m, p)
  si <- A[[paste0("Pred_", tag)]]; in_si <- !is.na(si)
  cat(sprintf("\n== %s (%d samples in SI Table 6) ==\n", tag, sum(in_si)))
  check("predicted phenotype identical for every sample",
        all(as.character(p$class)[in_si] == si[in_si]))
  check("nearest centroid identical",
        all(cf$Closest[in_si] == O[[paste0("Closest_", tag)]][in_si]))
  check("max posterior equal to 4 dp",
        all(abs(round(cf$MaxPost, 4) - as.numeric(O[[paste0("MaxPost_", tag)]])) [in_si] < 1e-9))
  check("Mahalanobis distance equal to the table's rounding",
        all(abs(cf$Mahal - as.numeric(O[[paste0("Mahal_", tag)]]))[in_si] <= 0.0051))

  mt <- ecmx_metrics(as.character(p$class), meta$Truth_Phenotype, meta$Dataset)
  b  <- B[B$Column_in_MasterWide == tag, ]
  s  <- mt$summary
  check(sprintf("micro accuracy %.1f%% = SI %.1f%%", s$Micro_accuracy_pct, b$Micro_accuracy_pct),
        round(s$Micro_accuracy_pct, 1) == b$Micro_accuracy_pct)
  check(sprintf("balanced accuracy %.1f%% = SI %.1f%%", s$Balanced_accuracy_pct, b$Balanced_accuracy_pct),
        round(s$Balanced_accuracy_pct, 1) == b$Balanced_accuracy_pct)
  check(sprintf("Wilson CI %.1f-%.1f = SI", s$Micro_CI95_low, s$Micro_CI95_high),
        round(s$Micro_CI95_low, 1) == b$Micro_CI95_low && round(s$Micro_CI95_high, 1) == b$Micro_CI95_high)
  check("mean per-dataset accuracy = SI",
        round(s$Mean_perDataset_accuracy_pct, 1) == b$Mean_perDataset_accuracy_pct)
  check("stress sensitivity = SI", round(s$Stress_sensitivity_pct, 1) == b$Stress_sensitivity_pct)
  check("per-phenotype recall = SI sheet C (Figure 7B)",
        all(round(mt$recall$Recall_pct, 1) == C[[tag]][match(mt$recall$Phenotype, C$Phenotype)]))
}

cat("\n", strrep("=", 78), "\n", sep = "")
if (fails == 0L) cat("ALL CHECKS PASSED\n") else cat(fails, "CHECK(S) FAILED\n")
quit(status = if (fails == 0L) 0L else 1L)
