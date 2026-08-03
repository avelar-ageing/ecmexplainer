# ==============================================================================
# app.R  –  ECMExplainer Shiny App Entry Point (golem-compatible)
#
# Run with:   shiny::runApp("shiny_app")   from the project root, OR
#             shiny::runApp()              from inside the shiny_app/ directory.
# ==============================================================================

# ---- Source all R modules ----
# Determine the directory containing this app.R
app_dir <- tryCatch(
  dirname(normalizePath(sys.frame(0)$ofile, mustWork = TRUE)),
  error = function(e) getwd()
)
r_dir <- file.path(app_dir, "R")
if (!dir.exists(r_dir)) r_dir <- file.path(getwd(), "R")

for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

# ---- Required packages ----
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(DT)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(MASS)
  library(sva)
  library(dplyr)
  library(waiter)
})

# ---- Upload size limit (set to 2 GB) ----
options(shiny.maxRequestSize = 2 * 1024^3)

# ---- Launch ----
shiny::shinyApp(ui = app_ui, server = app_server)
