# ==============================================================================
# build_models.R — rebuild the three bundled panels from the training cohort.
#
#   Rscript scripts/build_models.R
#
# ADSIG-180 is the deployed ECMExplainer model; POOL-180 (all pairs within the
# curated matrisome list) and ALLGENE-180 (all pairs within the 2,000 most
# variable genes genome-wide) are the size-matched controls of Supplementary
# Tables 5 and 6. Each is checked against the supplementary tables before it is
# written; the script stops on any mismatch. ADSIG-180 is compared with the
# shipped model, which it must reproduce ratio for ratio, and is not rewritten.
# ==============================================================================
suppressPackageStartupMessages({ library(MASS); library(readxl) })
source("R/ecm_core.R")

counts   <- readRDS("data/salmon_counts_Combined_FINAL.rds")
meta     <- read.csv("data/salmon_metadata_Combined_FINAL.csv", stringsAsFactors = FALSE)
genelist <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
si <- function(file, sheet, skip = 0)
  as.data.frame(read_excel(file.path("supplementary", file), sheet, skip = skip))

expect <- list(
  ADSIG   = list(ratios = si("SI_Table3.xlsx", "A_Ratio_Panel")$Ratio,
                 loocv = 92.71, loso = 47.92),
  POOL    = list(ratios = si("SI_Table5.xlsx", "A_POOL180_Ratios", 2)$Ratio,
                 loocv = 97.92, loso = 41.67),
  ALLGENE = list(ratios = si("SI_Table5.xlsx", "I_ALLGENE180_Ratios", 2)$Ratio,
                 loocv = 100.00, loso = 42.71))

cv <- list()
for (d in names(expect)) {
  prep <- ecmx_prepare(counts, meta, genelist, design = d)
  fit  <- ecmx_train(design = d, n_per_pheno = 30L, prep = prep)
  lo   <- ecmx_loso(prep, n_per_pheno = 30L)
  e    <- expect[[d]]
  ok <- setequal(fit$features$Ratio_ID, e$ratios) &&
        round(fit$loocv_acc, 2) == e$loocv && round(lo$accuracy, 2) == e$loso
  cat(sprintf("%-8s %s  %d ratios | LOOCV %.2f%% | nested LOSO %.2f%% (restricted %.2f%%)\n",
              d, if (ok) "OK  " else "FAIL", fit$panel_size, fit$loocv_acc,
              lo$accuracy, lo$restricted))
  if (!ok) stop(d, "-180 does not reproduce the supplementary tables.")
  fit$loso <- lo[c("accuracy", "restricted", "balanced", "macro_f1", "per_class")]
  cv[[d]] <- data.frame(Panel = c(ADSIG = "ADSIG-180 (ECMExplainer)", POOL = "POOL-180",
                                  ALLGENE = "ALLGENE-180")[[d]],
                        Ratios = fit$panel_size, PCs = fit$n_pcs,
                        `LOOCV (%)` = round(fit$loocv_acc, 2),
                        `Nested LOSO (%)` = round(lo$accuracy, 2),
                        `Nested LOSO, phenotypes in >= 3 studies (%)` = round(lo$restricted, 2),
                        `Nested LOSO macro-F1` = round(lo$macro_f1, 3), check.names = FALSE)
  fit$label <- c(ADSIG = "ADSIG-180 (ECMExplainer)", POOL = "POOL-180",
                 ALLGENE = "ALLGENE-180")[[d]]
  if (d == "ADSIG") {
    shipped <- readRDS("data/TrainedModel_ADSIG_180.rds")
    stopifnot(identical(fit$features$Ratio_ID, shipped$features$Ratio_ID))
  } else {
    fit$trained_at <- NULL       # keep the file byte-stable across rebuilds
    saveRDS(fit, sprintf("data/TrainedModel_%s_180.rds", d))
  }
}
write.csv(do.call(rbind, cv), "data/panel_cross_validation.csv", row.names = FALSE)
cat("Models and data/panel_cross_validation.csv written.\n")
