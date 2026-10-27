# Dependencies for ECMExplainer. CRAN only — no Bioconductor is needed to run the
# app, because sign features need no batch correction.
#   shiny, ggplot2, MASS, ggrepel, commonmark   the app
#   future, promises                            retraining in the background
#   readxl                                      tests and scripts that read the SI tables
pkgs <- c("shiny", "ggplot2", "MASS", "ggrepel", "commonmark", "future", "promises", "readxl")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All dependencies already present.")
}
if (utils::packageVersion("shiny") < "1.8.1")
  stop("ECMExplainer needs shiny >= 1.8.1; run install.packages('shiny').")
invisible(lapply(pkgs, function(p)
  message(sprintf("  %-12s %s", p, as.character(packageVersion(p))))))
