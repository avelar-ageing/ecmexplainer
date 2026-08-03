# ==============================================================================
# install_packages.R
#
# Run this script ONCE before launching the ECMExplainer Shiny app to install
# all required R packages.
#
# ── RECOMMENDED (reproducible): use renv ──────────────────────────────────────
# A renv.lock file is included. To restore the exact package versions used
# during development, run from INSIDE the shiny_app/ directory:
#
#   install.packages("renv")   # only needed once
#   renv::restore()            # installs everything from renv.lock
#
# After that, launch the app with:
#   shiny::runApp(".")         # from inside shiny_app/
#   # OR
#   shiny::runApp("shiny_app") # from the project root
#
# ── ALTERNATIVE (install latest versions) ─────────────────────────────────────
# If you prefer to install the latest package versions (not pinned to lockfile),
# run this script directly:
#
# Usage (from an R console or RStudio):
#   source("shiny_app/install_packages.R")
#
# Or from a terminal:
#   Rscript shiny_app/install_packages.R
# ==============================================================================

cat("=======================================================\n")
cat("  ECMExplainer – Package Installer\n")
cat("=======================================================\n\n")

# ---- CRAN packages ----
cran_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "DT",
  "ggplot2",
  "ggrepel",
  "tidyr",
  "dplyr",
  "MASS",
  "jsonlite",
  "waiter"
)

cran_missing <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]

if (length(cran_missing) > 0) {
  cat("Installing missing CRAN packages:", paste(cran_missing, collapse = ", "), "\n\n")
  install.packages(cran_missing, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

# ---- Bioconductor packages ----
bioc_packages <- c("sva", "limma")

bioc_missing <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]

if (length(bioc_missing) > 0) {
  cat("\nInstalling missing Bioconductor packages:", paste(bioc_missing, collapse = ", "), "\n\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# ---- Final check ----
cat("\n--- Verification ---\n")
all_packages <- c(cran_packages, bioc_packages)
status <- sapply(all_packages, requireNamespace, quietly = TRUE)

for (pkg in all_packages) {
  cat(sprintf("  %-20s %s\n", pkg, if (status[pkg]) "[OK]" else "[MISSING - check errors above]"))
}

if (all(status)) {
  cat("\n✓ All packages installed. You can now run the app:\n")
  cat("    shiny::runApp(\"shiny_app\")\n\n")
} else {
  cat("\n✗ Some packages failed to install. Please check the error messages above.\n\n")
}
