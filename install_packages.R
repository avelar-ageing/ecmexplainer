# Dependencies for ECMExplainer. CRAN only — no Bioconductor required, because
# the deployed method needs no batch correction (sign features are per-sample
# batch-invariant by construction).
pkgs <- c("shiny", "ggplot2", "MASS", "ggrepel", "commonmark")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All dependencies already present.")
}
invisible(lapply(pkgs, function(p)
  message(sprintf("  %-12s %s", p, as.character(packageVersion(p))))))
