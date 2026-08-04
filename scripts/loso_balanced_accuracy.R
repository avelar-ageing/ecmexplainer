# ==============================================================================
# loso_balanced_accuracy.R
#
# Recomputes the LOSO sweep metrics from the per-sample predictions, separating
# plain accuracy from balanced accuracy.
#
# Why this exists: 002c_nested_loso_sweep.R reports LOSO as
#   mean(loso_pred == phenotype, na.rm = TRUE)
# which is *plain* (unweighted) accuracy. Methods 5.4 calls the resulting 57.3%
# a "balanced accuracy". It is not: the training set is 42% Proliferating
# (40/96), and Proliferating is the best-recovered class, so the unweighted mean
# is pulled well above the mean of per-class recalls.
#
#   Rscript scripts/loso_balanced_accuracy.R [path/to/PanelSweep_perSample.csv]
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "Salmon_kTSP_FAST_NESTED/PanelSweep_perSample.csv"
if (!file.exists(path))
  stop("Per-sample sweep predictions not found at ", path,
       "\nRun scripts/002c_nested_loso_sweep.R first.")

ps <- read.csv(path, stringsAsFactors = FALSE)

summarise_one <- function(d) {
  d   <- d[!is.na(d$LOSO_pred), , drop = FALSE]
  hit <- d$LOSO_pred == d$True
  rec <- tapply(hit, d$True, mean)
  nst <- tapply(d$Study, d$True, function(x) length(unique(x)))
  data.frame(
    Design            = d$Design[1],
    N_target          = d$N_target[1],
    n_scored          = nrow(d),
    plain_accuracy    = round(mean(hit) * 100, 2),
    balanced_accuracy = round(mean(rec, na.rm = TRUE) * 100, 2),
    bal_acc_ge2_study = round(mean(rec[names(rec)[nst >= 2]], na.rm = TRUE) * 100, 2),
    bal_acc_ge3_study = round(mean(rec[names(rec)[nst >= 3]], na.rm = TRUE) * 100, 2)
  )
}

out <- do.call(rbind, lapply(
  split(ps, list(ps$Design, ps$N_target), drop = TRUE), summarise_one))
out <- out[order(out$Design, out$N_target), ]
rownames(out) <- NULL
print(out, row.names = FALSE)

cat("\nPer-class LOSO recall (%), ADSIG at the largest panel:\n")
d <- ps[ps$Design == "ADSIG" & ps$N_target == max(ps$N_target) & !is.na(ps$LOSO_pred), ]
print(round(tapply(d$LOSO_pred == d$True, d$True, mean) * 100, 1))
cat("\nSamples per class:\n"); print(table(d$True))
cat("\nSource studies per class:\n")
print(tapply(d$Study, d$True, function(x) length(unique(x))))

cat("\nNote: `plain_accuracy` is the column the sweep script writes as `LOSO`.\n",
    "`balanced_accuracy` is the mean of per-class recalls, which is what the\n",
    "term normally denotes.\n", sep = "")
