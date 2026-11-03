# ==============================================================================
# ECMExplainer — Shiny application
#
# Positions fibroblast transcriptomes in the ECMExplainer functional matrisome
# phenotype space with the manuscript's three sign-ratio panels (ADSIG-180, the
# ECMExplainer model, and the POOL-180 and ALLGENE-180 controls), scores them
# against known labels, and retrains panels on new data or new gene lists.
#
# Run locally (no limits):        shiny::runApp()
# Run as a public instance:       ECMX_MODE=web shiny::runApp()
#
# Web mode caps upload size, sample count and retraining size so a shared
# instance stays responsive. Nothing scientific differs between the modes — the
# models, the projection and the figures are identical.
# ==============================================================================
suppressPackageStartupMessages({
  library(shiny); library(ggplot2); library(MASS)
})
source("R/ecm_core.R")
source("R/ecm_figures.R")

ECMX_MODE <- tolower(Sys.getenv("ECMX_MODE", "local"))
WEB <- identical(ECMX_MODE, "web")
LIM <- list(
  upload_mb     = as.numeric(Sys.getenv("ECMX_MAX_UPLOAD_MB", if (WEB) "50"  else "4000")),
  max_samples   = as.integer(Sys.getenv("ECMX_MAX_SAMPLES",   if (WEB) "500" else "100000")),
  max_train     = as.integer(Sys.getenv("ECMX_MAX_TRAIN",     if (WEB) "200" else "100000")),
  max_per_pheno = as.integer(Sys.getenv("ECMX_MAX_PER_PHENO", if (WEB) "60"  else "500"))
)
options(shiny.maxRequestSize = LIM$upload_mb * 1024^2)
REPO_URL <- "https://github.com/avelar-ageing/ecmexplainer"
APP_DIR  <- normalizePath(".")

# base R gained %||% in 4.4.0; define it if we are on something older.
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

# Retraining runs in a background R process when `future` is available, so one
# visitor's retrain never stalls the app for everyone else on a shared server.
# ECMX_WORKERS=0 runs retraining in-process (for environments without local sockets).
ASYNC <- requireNamespace("future", quietly = TRUE) && requireNamespace("promises", quietly = TRUE) &&
  as.integer(Sys.getenv("ECMX_WORKERS", "2")) > 0
if (ASYNC) future::plan(future::multisession,
                        workers = as.integer(Sys.getenv("ECMX_WORKERS", "2")))

# ---- bundled panels ------------------------------------------------------------
#' Attach what every view needs once, so switching panels costs nothing.
prepare_panel <- function(m, label, key) {
  m$label <- label; m$key <- key
  pcs <- as.data.frame(m$pca$x[, seq_len(m$n_pcs), drop = FALSE])
  colnames(pcs) <- paste0("PC", seq_len(m$n_pcs))
  m$train_coords <- predict(m$lda, pcs)$x
  m$centroids <- ecmx_centroids(m)
  f <- m$features
  if (is.null(f$Ratio_ID) && !is.null(f$Ratio)) f$Ratio_ID <- f$Ratio
  if (is.null(f$Module)   && !is.null(f$Category)) f$Module <- f$Category
  if (is.null(f$eta_sq)) f$eta_sq <- NA_real_
  m$features <- f
  m
}
PANELS <- list(
  ADSIG180   = prepare_panel(readRDS("data/TrainedModel_ADSIG_180.rds"),
                             "ADSIG-180 (ECMExplainer)", "ADSIG180"),
  POOL180    = prepare_panel(readRDS("data/TrainedModel_POOL_180.rds"), "POOL-180", "POOL180"),
  ALLGENE180 = prepare_panel(readRDS("data/TrainedModel_ALLGENE_180.rds"), "ALLGENE-180", "ALLGENE180"))
PANEL_NOTES <- c(
  ADSIG180   = "ECMExplainer: Assembly/Destruction and Signalling/Signalling ratios from the 1,139-gene curated matrisome list",
  POOL180    = "Control: any pair within the same curated list, without the functional pairing rule",
  ALLGENE180 = "Control: any pair within the 2,000 most variable genes genome-wide")

GENELIST  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
TRAIN_CNT <- "data/salmon_counts_Combined_FINAL.rds"
TRAIN_MET <- "data/salmon_metadata_Combined_FINAL.csv"
VAL_EXPR  <- readRDS("data/validation_expr.rds")
VAL_META  <- read.csv("data/validation_metadata.csv", stringsAsFactors = FALSE)
EXAMPLE_DS <- "SRP113329"     # "Load example": one validation cohort, proliferating vs RAS-induced senescence
CONTROL_KEYS <- c("POOL180", "ALLGENE180")
AGE_WEIGHTS <- read.csv("data/aging_weights_SI4.csv", stringsAsFactors = FALSE)
AGE_REF     <- read.csv("data/aging_reference_bins.csv", stringsAsFactors = FALSE)
# Set CITATION_URL once the preprint or paper has a DOI; until then no link is shown.
CITATION     <- paste("Sheekey E, Avelar RA, Fujiwara H. ECMExplainer: a functional matrisome framework",
                      "for decoding cellular stress states and aging trajectories. Manuscript in preparation, 2026.")
CITATION_URL <- ""
# Column names shown to people, not R identifiers.
pretty_col <- function(x) {
  map <- c(Truth_Phenotype = "Known phenotype", Expected_Phenotype = "Expected phenotype",
           Exclusion_Status = "Exclusion status")
  ifelse(x %in% names(map), map[x], gsub("_", " ", x))
}
col_choices <- function(cols) setNames(cols, pretty_col(cols))
SI_FILES  <- c("SI_Table1.xlsx" = "Curated matrisome gene list",
               "SI_Table2.xlsx" = "Training samples (n = 96)",
               "SI_Table3.xlsx" = "ADSIG-180 ratio panel, panel-size sweep, gene-level shifts",
               "SI_Table4.xlsx" = "Age-associated ratios and the ECM Aging Score",
               "SI_Table5.xlsx" = "POOL-180 and ALLGENE-180 control panels",
               "SI_Table6.xlsx" = "Validation across 18 cohorts, all panels")

# ---- the manuscript's validation, recomputed once at start-up ---------------------
score_panels <- function(panels, expr) {
  lapply(panels, function(m) {
    p <- ecmx_project(m, expr)
    list(panel = m, proj = p, conf = ecmx_confidence(m, p))
  })
}
validation_tables <- function(scored, truth, dataset) {
  ms <- lapply(scored, function(s) ecmx_metrics(as.character(s$proj$class), truth, dataset))
  lab <- vapply(scored, function(s) s$panel$label, character(1))
  summary <- do.call(rbind, Map(function(m, l) cbind(Panel = l, m$summary), ms, lab))
  recall  <- do.call(rbind, Map(function(m, l) cbind(Panel = l, m$recall), ms, lab))
  per_ds  <- do.call(rbind, Map(function(m, l) if (!is.null(m$per_dataset)) cbind(Panel = l, m$per_dataset), ms, lab))
  conf    <- do.call(rbind, Map(function(m, l) cbind(Panel = l, as.data.frame(m$confusion)), ms, lab))
  mcn <- NULL
  if (length(ms) > 1) {
    ref <- ms[[1]]$correct
    mcn <- do.call(rbind, lapply(seq_along(ms)[-1], function(i) {
      cmp <- ms[[i]]$correct; common <- intersect(names(ref), names(cmp))
      t <- ecmx_mcnemar(ref[common], cmp[common])
      data.frame(Comparison = paste(lab[1], "vs", lab[i]), N_paired = length(common),
                 Accuracy_reference_pct = 100 * mean(ref[common]),
                 Accuracy_comparator_pct = 100 * mean(cmp[common]),
                 Only_reference_correct = t$only_a, Only_comparator_correct = t$only_b,
                 McNemar_exact_p = signif(t$p, 3))
    }))
    mcn$BH_FDR <- signif(stats::p.adjust(mcn$McNemar_exact_p, "BH"), 3)
  }
  list(summary = summary, recall = recall, per_dataset = per_ds, confusion = conf, mcnemar = mcn)
}
MS_SCORED <- score_panels(PANELS, VAL_EXPR)
MS_TABLES <- validation_tables(MS_SCORED, VAL_META$Truth_Phenotype, VAL_META$Dataset)

# Training, at top level so a background worker receives it without the session.
run_training <- function(args) {
  prep <- do.call(ecmx_prepare, args$prep)
  fit <- do.call(ecmx_train, c(list(prep = prep), args$train))
  if (isTRUE(args$loso) && !is.null(prep$study) && length(unique(prep$study)) > 1)
    fit$loso <- do.call(ecmx_loso, c(list(prep = prep, n_per_pheno = args$train$n_per_pheno),
                                     args$train[c("eta_min", "rank_by", "min_per_module")]))
  fit
}
PANEL_CV <- read.csv("data/panel_cross_validation.csv", stringsAsFactors = FALSE, check.names = FALSE)

# ---- small helpers -------------------------------------------------------------------
notify_err <- function(e) showNotification(conditionMessage(e), type = "error", duration = 12)
info_box   <- function(...) div(class = "ecmx-note", ...)
warn_box   <- function(...) div(class = "ecmx-warn", ...)
#' Format numbers for display: whole numbers without decimals, p-values and
#' FDRs to 3 significant figures, everything else to `digits` places.
fmt_tbl <- function(df, digits = 1) {
  for (nm in names(df)) {
    x <- df[[nm]]
    if (!is.numeric(x)) next
    df[[nm]] <- if (grepl("(^|_)(p|P)$|_p$|FDR", nm)) formatC(x, digits = 3, format = "g")
                else if (all(is.na(x) | x == round(x))) formatC(x, format = "d", big.mark = ",")
                else formatC(x, digits = digits, format = "f")
    df[[nm]][is.na(x)] <- "–"
  }
  names(df) <- gsub("_", " ", names(df))
  df
}

#' Read a metadata CSV and make sure it has a SampleID column (or use the first).
read_meta <- function(path) {
  md <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  id <- which(tolower(names(md)) %in% c("sampleid", "sample_id", "sample"))[1]
  if (is.na(id)) id <- 1L
  names(md)[id] <- "SampleID"
  md$SampleID <- as.character(md$SampleID)
  md
}

#' Read a gene list: a CSV with a Gene column (and optionally Function), a
#' one-column file of symbols, or pasted text.
read_genelist <- function(path = NULL, text = NULL) {
  if (!is.null(text)) {
    g <- unique(trimws(unlist(strsplit(text, "[,;[:space:]]+"))))
    return(data.frame(Gene = g[nzchar(g)], stringsAsFactors = FALSE))
  }
  gl <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gcol <- which(tolower(names(gl)) %in% c("gene", "symbol", "gene_symbol", "hgnc_symbol"))[1]
  if (is.na(gcol)) {
    # No header we recognise: treat the whole first column, header included, as symbols.
    gl <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE)
    gcol <- 1L
  }
  names(gl)[gcol] <- "Gene"
  gl$Gene <- trimws(as.character(gl$Gene))
  gl[nzchar(gl$Gene), , drop = FALSE]
}

# Two filter slots per data set: a column and the values to keep.
filter_slots_ui <- function(prefix, md, n = 2) {
  if (is.null(md)) return(helpText("Add sample metadata to filter samples."))
  cols <- setdiff(names(md), "SampleID")
  tagList(lapply(seq_len(n), function(i) fluidRow(
    column(5, selectInput(paste0(prefix, "_col", i), if (i == 1) "Filter by" else NULL,
                          c("(no filter)" = "", cols))),
    column(7, uiOutput(paste0(prefix, "_val", i, "_ui"))))))
}
apply_filters <- function(input, prefix, md, n = 2) {
  if (is.null(md)) return(NULL)
  keep <- rep(TRUE, nrow(md))
  for (i in seq_len(n)) {
    col <- input[[paste0(prefix, "_col", i)]]; val <- input[[paste0(prefix, "_val", i)]]
    if (!is.null(col) && nzchar(col) && col %in% names(md) && length(val))
      keep <- keep & as.character(md[[col]]) %in% val
  }
  md$SampleID[keep]
}
filter_value_ui <- function(input, prefix, i, md) {
  col <- input[[paste0(prefix, "_col", i)]]
  if (is.null(md) || is.null(col) || !nzchar(col) || !col %in% names(md)) return(NULL)
  v <- sort(unique(as.character(md[[col]])))
  selectizeInput(paste0(prefix, "_val", i), if (i == 1) "keep" else NULL, v, selected = v,
                 multiple = TRUE, options = list(plugins = list("remove_button")))
}

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  title = "ECMExplainer",
  tags$head(
    tags$meta(name = "description", content = paste(
      "ECMExplainer: a functional matrisome framework for decoding cellular stress",
      "states and aging trajectories. Sheekey, Avelar & Fujiwara.")),
    tags$style(HTML("
    @font-face {
      font-family:'Space Grotesk'; font-style:normal; font-weight:300 700; font-display:swap;
      src:url('fonts/SpaceGrotesk-Variable.woff2') format('woff2');
    }
    :root {
      --ecmx-navy:#22224a; --ecmx-blue:#3b4cca; --ecmx-deep:#1f2a6b;
      --ecmx-pink:#d6336c; --ecmx-pink-dark:#a61e4d; --ecmx-pink-soft:#ffb3cf;
      --ecmx-light:var(--ecmx-pink); --ecmx-page:#fbf8fb; --ecmx-rule:#ecdfe8;
    }
    body { background:var(--ecmx-page); }
    .container-fluid { max-width:1680px; }

    .ecmx-head {
      background:linear-gradient(115deg, var(--ecmx-deep) 0%, var(--ecmx-blue) 55%, var(--ecmx-pink) 100%);
      color:#fff; padding:16px 22px 14px; margin:-15px -15px 18px -15px;
      border-bottom:4px solid var(--ecmx-pink);
    }
    .ecmx-head .wordmark {
      font-family:'Space Grotesk', 'Helvetica Neue', Arial, sans-serif;
      font-size:2.15em; font-weight:700; letter-spacing:-.02em; line-height:1;
    }
    .ecmx-head .wordmark .ecm { color:#fff; }
    .ecmx-head .wordmark .explainer { color:var(--ecmx-pink-soft); font-weight:500; }
    .ecmx-head .wordmark .dot { color:var(--ecmx-pink-soft); }
    .ecmx-head .tagline  { color:#eef0fb; font-size:13px; margin-left:12px; }
    .ecmx-head .authors  { color:#fff; font-size:13px; margin-top:7px; }
    .ecmx-head .affil    { color:#e3e5f5; font-size:11.5px; margin-top:3px; line-height:1.5; }
    .ecmx-head .cite     { opacity:.85; font-style:italic; }
    .ecmx-head sup { font-size:70%; }
    .ecmx-head a { color:#fff; text-decoration:underline; }

    .well {
      background:#fff; border:1px solid var(--ecmx-rule);
      border-top:3px solid var(--ecmx-light); border-radius:4px;
      box-shadow:0 1px 2px rgba(34,34,74,.06);
    }
    .well h4 { color:var(--ecmx-navy); font-weight:700; font-size:14px;
               letter-spacing:.02em; margin-top:14px; }
    .well h4:first-child { margin-top:0; }
    .well h5 { color:var(--ecmx-pink-dark); font-weight:700; font-size:13px; }
    .tab-content { background:#fff; border:1px solid var(--ecmx-rule);
                   border-top:none; border-radius:0 0 4px 4px; padding:6px 14px 14px; }
    .tab-content .tab-content { border:none; padding:6px 0 0; }
    .nav-tabs { border-bottom:1px solid var(--ecmx-rule); }
    .nav-tabs > li > a { color:#5a6b7d; font-weight:600; font-size:13px; }
    .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover, .nav-tabs > li.active > a:focus {
      color:var(--ecmx-pink-dark); border-top:2px solid var(--ecmx-pink);
    }
    .ecmx-note { background:#f4f5fd; border:1px solid #d3d7f3;
                 border-left:4px solid var(--ecmx-blue); padding:9px 12px;
                 font-size:12.5px; margin:10px 0; border-radius:3px; color:#33475b; }
    .ecmx-warn { background:#fff8e1; border:1px solid #f5dfa8;
                 border-left:4px solid #f39c12; padding:9px 12px;
                 font-size:12.5px; margin:10px 0; border-radius:3px; color:#6b5320; }
    .ecmx-ok   { background:#eef8f0; border:1px solid #bfe3c7;
                 border-left:4px solid #2e8b57; padding:9px 12px;
                 font-size:12.5px; margin:10px 0; border-radius:3px; color:#24513a; }
    .btn-primary { background:var(--ecmx-pink); border-color:var(--ecmx-pink-dark); }
    .btn-primary:hover, .btn-primary:focus { background:var(--ecmx-pink-dark); border-color:var(--ecmx-pink-dark); }
    .btn-default { border-color:#e2c6d4; color:var(--ecmx-navy); }
    .btn-default:hover { border-color:var(--ecmx-pink); color:var(--ecmx-pink-dark); }
    a { color:var(--ecmx-pink-dark); } a:hover, a:focus { color:var(--ecmx-blue); }
    .ecmx-head a { color:#fff; }
    .checkbox input[type=checkbox], .radio input[type=radio] { accent-color:var(--ecmx-pink); }
    label { color:#33475b; font-weight:600; font-size:12.5px; }
    .help-block { font-size:11.5px; color:#616e78; }
    .shiny-output-error-validation { color:#8a6d1f; font-size:13px; }
    table.data, .table { font-size:12px; }
    .table > thead > tr > th { color:var(--ecmx-navy); border-bottom:2px solid var(--ecmx-rule); }
    .ecmx-steps { font-size:12.5px; color:#33475b; padding-left:18px; }
    .ecmx-steps li { margin-bottom:3px; }
    .ecmx-scroll { overflow-x:auto; }
    .ecmx-phone-only { display:none; }
    @media (max-width: 767px) {
      .ecmx-head { padding:10px 14px 9px; margin-bottom:10px; }
      .ecmx-head .wordmark { font-size:1.75em; }
      .ecmx-head .tagline { display:block; margin:3px 0 0; font-size:12px; }
      .ecmx-head .authors { font-size:12px; margin-top:4px; }
      .ecmx-head .affil { display:none; }              /* affiliations and citation: About tab */
      /* one swipeable row of tabs instead of three wrapped rows */
      .nav-tabs { display:flex; flex-wrap:nowrap; overflow-x:auto; overflow-y:hidden;
                  -webkit-overflow-scrolling:touch; scrollbar-width:none; }
      .nav-tabs::-webkit-scrollbar { display:none; }
      .nav-tabs > li { float:none; flex:0 0 auto; }
      .nav-tabs > li > a { padding:8px 10px; white-space:nowrap; margin-right:0; }
      .tab-content { padding:6px 8px 10px; }
      /* results first; the controls follow */
      .tab-pane > .row { display:flex; flex-direction:column; }
      .tab-pane > .row > .col-sm-9 { order:-1; }
      .ecmx-phone-only { display:block; }
    }
  "))),

  div(class = "ecmx-head",
      # One string, so no whitespace separates the two halves of the name.
      div(span(class = "wordmark", HTML('<span class="ecm">ECM</span><span class="explainer">Explainer</span>')),
          span(class = "tagline",
               "A functional matrisome framework for decoding cellular stress states and aging trajectories")),
      div(class = "authors",
          HTML("Eleanor Sheekey<sup>1,3</sup>, Roberto A. Avelar<sup>2</sup>, Hironobu Fujiwara<sup>1</sup>")),
      div(class = "affil",
          HTML(paste0(
            "<sup>1</sup> RIKEN Center for Biosystems Dynamics Research, Kobe, Japan &nbsp;·&nbsp; ",
            "<sup>2</sup> Institute for Biostatistics and Informatics in Medicine and Ageing Research, ",
            "Rostock University Medical Center, Germany &nbsp;·&nbsp; ",
            "<sup>3</sup> Corresponding author &nbsp;·&nbsp; ",
            "<a href='", REPO_URL, "' target='_blank' rel='noopener'>Code and local install</a>"))),
      div(class = "affil cite",
          HTML(paste0("<b>Cite:</b> ", htmltools::htmlEscape(CITATION),
                      if (nzchar(CITATION_URL)) paste0(" <a href='", CITATION_URL,
                                                       "' target='_blank' rel='noopener'>", CITATION_URL, "</a>"))))),

  if (WEB) warn_box(
    strong("Shared public instance. "),
    sprintf("Uploads are capped at %g MB and %d samples; retraining at %d samples. ",
            LIM$upload_mb, LIM$max_samples, LIM$max_train),
    "Uploaded files are held in memory for your session only and are never written to disk. ",
    "For unrestricted use, ", a("run the app on your own computer", href = REPO_URL, target = "_blank"), "."),

  tabsetPanel(
    id = "tabs",

    # ---- Classify ----------------------------------------------------------------
    tabPanel(
      "Classify",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("1 · Expression data"),
          radioButtons("src", NULL,
                       c("Manuscript validation cohorts (327 samples, 18 studies)" = "validation",
                         "Example: one cohort (SRP113329, proliferating vs RAS senescence)" = "example",
                         "Upload my own" = "upload")),
          actionButton("load_example", "Load example (SRP113329)", class = "btn-sm btn-default"),
          br(), br(),
          conditionalPanel(
            "input.src == 'upload'",
            fileInput("expr_file", "Expression CSV (genes x samples)", accept = c(".csv", ".gz", ".tsv", ".txt")),
            info_box("HGNC symbols in the first column, one sample per remaining column. ",
                     "Raw counts, log2-CPM or VST are best — each ratio compares two genes ",
                     "within one sample. ", strong("FPKM/TPM work, approximately"),
                     ": gene-length normalisation shifts every ratio by a fixed amount, which can ",
                     "flip ratios close to zero. On GSE113957 (143 samples deposited as FPKM, ",
                     "also available as counts) 85% of ratio signs and 136 of 143 calls agreed."),
            fileInput("meta_file", "Sample metadata CSV (optional)", accept = ".csv"),
            helpText("A SampleID column plus anything you want to group, filter or score by ",
                     "(e.g. Condition, Dataset, Truth)."),
            div(style = "font-size:12px;", "Examples: ",
                downloadLink("dl_ex_expr", "expression"), " · ",
                downloadLink("dl_ex_meta", "metadata"), " · ",
                downloadLink("dl_ex_genes", "gene list"))
          ),
          h4("2 · Samples"),
          uiOutput("cls_filter"),
          uiOutput("n_samples_msg"),
          h4("3 · Panels"),
          uiOutput("panel_pick"),
          checkboxInput("adv_panels", "Advanced: compare control panels (POOL-180, ALLGENE-180)", FALSE),
          uiOutput("primary_pick"),
          h4("4 · Labels"),
          uiOutput("label_picks"),
          hr(),
          h4("Export"),
          fluidRow(
            column(6, numericInput("exp_w", "Width (in)",  9, 3, 30, 0.5)),
            column(6, numericInput("exp_h", "Height (in)", 7, 3, 40, 0.5))),
          fluidRow(
            column(6, numericInput("exp_dpi", "DPI", 300, 72, 600, 1)),
            column(6, radioButtons("exp_fmt", "Format", c(PNG = "png", PDF = "pdf"), inline = TRUE))),
          downloadButton("dl_plot", "Download current figure", class = "btn-sm"),
          br(), br(),
          downloadButton("dl_results", "Download results, all panels (CSV)", class = "btn-sm")
        ),
        mainPanel(
          width = 9,
          div(class = "ecmx-phone-only ecmx-note",
              "Data, panels and labels are below the figures. Swipe the tabs sideways for more views."),
          uiOutput("coverage_msg"),
          uiOutput("lineage_msg"),
          tabsetPanel(
            id = "viz",
            tabPanel("Phenotype space",    value = "lda",  br(), plotOutput("p_lda", height = "auto")),
            tabPanel("Posteriors",         value = "post", br(), plotOutput("p_post", height = "auto")),
            tabPanel("Centroid distances", value = "dist", br(), plotOutput("p_dist", height = "auto")),
            tabPanel("Compare panels",     value = "cmp",  br(), uiOutput("cmp_msg"),
                     plotOutput("p_cmp", height = "auto"), div(class = "ecmx-scroll", tableOutput("t_agree"))),
            tabPanel("Accuracy",           value = "acc",  br(), uiOutput("acc_body")),
            tabPanel("Panel heatmap",      value = "heat", br(), uiOutput("heat_pheno_pick"),
                     plotOutput("p_heat", height = "auto")),
            tabPanel("Results table",      value = "tbl",  br(), div(class = "ecmx-scroll", tableOutput("t_res")))
          )
        )
      )
    ),

    # ---- Trajectories ------------------------------------------------------------
    tabPanel(
      "Trajectories",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          info_box("Each condition is represented by its ", strong("medoid"),
                   " (a real sample), placed by its distance to two training centroids. ",
                   "Protractor angles carry a bootstrap 95% CI. Uses the primary panel and ",
                   "the samples selected on the Classify tab."),
          uiOutput("traj_ds_pick"),
          uiOutput("traj_group"),
          uiOutput("traj_anchors"),
          uiOutput("traj_order"),
          hr(),
          h5("Protractor"),
          uiOutput("traj_from_to"),
          numericInput("n_boot", "Bootstrap resamples", 2000, 200, 5000, 100)
        ),
        mainPanel(
          width = 9,
          plotOutput("p_traj", height = "auto"),
          hr(),
          fluidRow(column(7, plotOutput("p_prot", height = "auto")),
                   column(5, br(), tableOutput("t_angles")))
        )
      )
    ),

    # ---- Train a new model ----------------------------------------------------------
    tabPanel(
      "Train a new model",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          info_box("Builds a ", strong("separate"), " model from your samples, labels and genes, ",
                   "for this session. The published ADSIG-180 model is never changed; your model ",
                   "is added next to it on the Classify tab."),
          h4("1 · Training data"),
          radioButtons("tr_src", NULL,
                       c("Manuscript training cohort (96 samples, 10 studies)" = "builtin",
                         "Upload my own" = "upload")),
          conditionalPanel("input.tr_src == 'upload'",
            fileInput("tr_counts", "Counts CSV (genes x samples)", accept = c(".csv", ".gz", ".tsv", ".txt")),
            fileInput("tr_meta", "Metadata CSV", accept = ".csv"),
            info_box("Labels can be ", strong("anything you like"),
                     " — the six reference states are not required. Every class needs at ",
                     "least 2 samples. Raw counts are best; log2-CPM also works.")),
          uiOutput("tr_cols"),
          uiOutput("tr_filter"),
          uiOutput("tr_n_msg"),
          h4("2 · Gene universe"),
          radioButtons("tr_genes", NULL,
                       c("Curated matrisome list (SI Table 1, 1,139 genes)" = "curated",
                         "Upload my own gene list" = "upload_list",
                         "Paste gene symbols" = "paste")),
          conditionalPanel("input.tr_genes == 'upload_list'",
            fileInput("tr_genelist", "Gene list CSV", accept = c(".csv", ".txt")),
            helpText("A Gene column. Add a Function column with Assembly / Destruction / ",
                     "Signalling to use the ADSIG design; any list works with POOL.")),
          conditionalPanel("input.tr_genes == 'paste'",
            textAreaInput("tr_genetext", NULL, rows = 5,
                          placeholder = "COL1A1 COL1A2 MMP2 TIMP1 ...")),
          h4("3 · Panel"),
          selectInput("tr_design", "Design",
                      c("ADSIG — Assembly/Destruction + Signalling pairs (ECMExplainer)" = "ADSIG",
                        "POOL — any pair within the gene list" = "POOL",
                        "ALLGENE — any pair within the most variable genes genome-wide" = "ALLGENE")),
          conditionalPanel("input.tr_design == 'ALLGENE'",
            numericInput("tr_topn", "Most variable genes to use", 2000, 200, 5000, 100),
            helpText("ALLGENE ignores the gene list.")),
          numericInput("tr_npp", "Ratios per phenotype", 30, 2, LIM$max_per_pheno, 1),
          selectInput("tr_rank", "Selection ranking",
                      c("TSP delta (manuscript)" = "tsp_delta", "eta-squared" = "eta2")),
          numericInput("tr_minmod", "Minimum ratios per module (ADSIG)", 0, 0, 60, 1),
          numericInput("tr_eta", "eta-squared floor", 0.10, 0, 0.9, 0.01),
          checkboxInput("tr_loso", "Also estimate nested leave-one-study-out accuracy", TRUE),
          helpText("Re-selects and re-fits the panel without each study in turn; needs a Study column."),
          actionButton("do_train", "Train model", class = "btn-primary"),
          br(), br(),
          downloadButton("dl_model", "Download model (.rds)", class = "btn-sm"),
          if (!WEB) tagList(br(), br(),
            fileInput("model_file", "Load a saved model (.rds)", accept = ".rds"))
        ),
        mainPanel(
          width = 9,
          uiOutput("train_status"),
          fluidRow(column(6, h5("Leave-one-sample-out (optimistic: selection saw every sample)"),
                          tableOutput("t_perclass")),
                   column(6, h5("Confusion matrix (LOOCV, rows = predicted)"),
                          tableOutput("t_confusion"))),
          uiOutput("loso_box"),
          hr(),
          plotOutput("p_train_lda", height = "auto"),
          hr(),
          h5("Selected panel"),
          div(class = "ecmx-scroll", style = "max-height:560px; overflow-y:auto;", tableOutput("t_panel"))
        )
      )
    ),

    # ---- ECM Aging Score -----------------------------------------------------------
    tabPanel(
      "ECM Aging Score",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          info_box("Scores any sample — no age information needed. The published score uses ",
                   "96 age-associated ratios with frozen weights (SI Table 4), fitted on healthy ",
                   "dermal fibroblasts aged 1–94 years (GSE113957). Higher score = older-like matrisome. ",
                   "Uses the samples selected on the Classify tab."),
          radioButtons("age_src", "Weights",
                       c("Published (SI Table 4, 96 ratios)" = "si4",
                         "Fit my own from an age column" = "fit")),
          conditionalPanel("input.age_src == 'fit'",
            helpText("Per-ratio age slopes are fitted by OLS on the continuous log-ratio of the ",
                     "primary panel's ratios, filtered at BH FDR, then applied to the sign features."),
            uiOutput("age_col_pick"),
            numericInput("age_fdr", "BH FDR threshold", 0.05, 0.001, 0.5, 0.005),
            actionButton("do_age", "Fit age weights", class = "btn-primary")),
          hr(),
          uiOutput("age_group_pick"),
          conditionalPanel("input.age_src == 'si4'",
            checkboxInput("age_ref", "Show the GSE113957 age reference", TRUE)),
          downloadButton("dl_age_scores", "Download scores (CSV)", class = "btn-sm"), br(), br(),
          downloadButton("dl_betas", "Download weights (CSV)", class = "btn-sm")
        ),
        mainPanel(
          width = 9,
          uiOutput("age_status"),
          plotOutput("p_age", height = "auto"),
          hr(), tableOutput("t_betas")
        )
      )
    ),

    # ---- Manuscript ------------------------------------------------------------------
    tabPanel(
      "Manuscript",
      br(),
      fluidRow(
        column(8,
          info_box(
            "The manuscript's independent validation, recomputed live by this app from the ",
            "bundled cohorts: 327 samples from 18 studies, of which 195 carry a ground-truth ",
            "label (Supplementary Table 6). Every prediction below is identical to the ",
            "published table; the repository's ", code("tests/test_validation.R"), " checks this."),
          h4("Panel summary (SI Table 6, sheet B)"),
          div(class = "ecmx-scroll", tableOutput("ms_summary")),
          fluidRow(column(6, plotOutput("ms_radar", height = "auto")),
                   column(6, plotOutput("ms_acc", height = "auto"))),
          h4("Training-cohort cross-validation (SI Tables 3 and 5)"),
          div(class = "ecmx-scroll", tableOutput("ms_training")),
          h4("Between-panel test (SI Table 6, sheet E)"),
          div(class = "ecmx-scroll", tableOutput("ms_mcnemar")),
          plotOutput("ms_conf", height = "auto"),
          plotOutput("ms_ds", height = "auto"),
          actionButton("ms_open", "Explore these cohorts on the Classify tab", class = "btn-default")
        ),
        column(4,
          div(class = "well",
            h4("Supplementary tables"),
            tagList(lapply(seq_along(SI_FILES), function(i) div(
              style = "margin-bottom:8px;",
              downloadLink(paste0("dl_si", i), names(SI_FILES)[i]), br(),
              span(style = "font-size:12px; color:#5a6b7d;", SI_FILES[[i]])))),
            hr(),
            h4("Figures"),
            tags$ul(class = "ecmx-steps",
              tags$li("Figure 7B — recall per phenotype: the radar on this tab."),
              tags$li("Phenotype space, posteriors, centroid distances and panel heatmaps: ",
                      "Classify tab, any panel."),
              tags$li("Trajectories and protractor angles: Trajectories tab."),
              tags$li("ECM Aging Score: ECM Aging Score tab.")),
            helpText("Every figure downloads as PNG or PDF at the size and resolution you set.")))
      )
    ),

    # ---- About ---------------------------------------------------------------------------
    tabPanel("About", br(), column(9, uiOutput("about_md")))
  )
)

# ==============================================================================
# server
# ==============================================================================
server <- function(input, output, session) {

  session_models <- reactiveVal(list())
  trained <- reactiveVal(NULL)
  betas   <- reactiveVal(NULL)

  all_panels <- reactive(c(PANELS, session_models()))

  # ---- data -----------------------------------------------------------------------------
  example_ids <- function() VAL_META$SampleID[VAL_META$Dataset == EXAMPLE_DS]
  observeEvent(input$load_example, {
    updateRadioButtons(session, "src", selected = "example")
    updateTabsetPanel(session, "viz", selected = "lda")
  })

  expr_raw <- reactive({
    if (identical(input$src, "validation")) return(VAL_EXPR)
    if (identical(input$src, "example")) return(VAL_EXPR[, example_ids(), drop = FALSE])
    req(input$expr_file)
    m <- tryCatch(ecmx_read_expression(input$expr_file$datapath),
                  error = function(e) { notify_err(e); NULL })
    req(m)
    validate(need(ncol(m) <= LIM$max_samples, sprintf(
      "This instance accepts at most %d samples; your file has %d. Run the app locally for unlimited use.",
      LIM$max_samples, ncol(m))))
    if (!ecmx_looks_logged(m)) m <- ecmx_logcpm(m)
    m
  })

  meta_raw <- reactive({
    if (identical(input$src, "validation")) return(VAL_META)
    if (identical(input$src, "example")) return(VAL_META[VAL_META$Dataset == EXAMPLE_DS, , drop = FALSE])
    if (is.null(input$meta_file)) return(NULL)
    tryCatch(read_meta(input$meta_file$datapath), error = function(e) { notify_err(e); NULL })
  })

  output$cls_filter <- renderUI(filter_slots_ui("cf", meta_raw()))
  output$cf_val1_ui <- renderUI(filter_value_ui(input, "cf", 1, meta_raw()))
  output$cf_val2_ui <- renderUI(filter_value_ui(input, "cf", 2, meta_raw()))

  kept_ids <- reactive({
    ids <- colnames(expr_raw()); md <- meta_raw()
    if (is.null(md)) return(ids)
    keep <- apply_filters(input, "cf", md)
    if (length(keep) == nrow(md)) ids else intersect(ids, keep)
  })

  expr_mat <- reactive({
    ids <- kept_ids()
    validate(need(length(ids) > 0, "No samples left after filtering."))
    expr_raw()[, ids, drop = FALSE]
  })

  meta_df <- reactive({
    md <- meta_raw(); ids <- colnames(expr_mat())
    if (is.null(md)) return(data.frame(SampleID = ids, stringsAsFactors = FALSE))
    out <- md[match(ids, md$SampleID), , drop = FALSE]
    out$SampleID <- ids
    out
  })

  output$n_samples_msg <- renderUI({
    n_all <- ncol(expr_raw()); n <- length(kept_ids())
    helpText(sprintf("%d of %d samples selected.", n, n_all))
  })

  # ---- panels -----------------------------------------------------------------------------
  output$panel_pick <- renderUI({
    ap <- all_panels(); adv <- isTRUE(input$adv_panels)
    # The two control panels are the manuscript's comparison, not alternatives to
    # use: they appear only under "Advanced", and arrive ticked when it is opened.
    show <- if (adv) names(ap) else setdiff(names(ap), CONTROL_KEYS)
    labs <- vapply(ap[show], function(m) m$label, character(1))
    sel <- isolate(panel_sel()); if (adv) sel <- union(sel, CONTROL_KEYS)
    tagList(
      checkboxGroupInput("panels", NULL, setNames(show, labs), selected = intersect(sel, show)),
      helpText(HTML(paste0("<b>ADSIG-180</b> is the ECMExplainer model",
                           if (adv) "; <b>POOL-180</b> and <b>ALLGENE-180</b> are the manuscript's controls",
                           ". Models you train or load are added here."))))
  })
  # The ticked panels, kept outside the widget so a newly trained or loaded
  # model can arrive already ticked when the picker is redrawn.
  panel_sel <- reactiveVal("ADSIG180")
  observeEvent(input$panels, panel_sel(input$panels), ignoreNULL = FALSE, ignoreInit = TRUE)
  add_model <- function(m, key) {
    session_models(c(session_models(), setNames(list(m), key)))
    panel_sel(c(panel_sel(), key))
    updateSelectInput(session, "primary", selected = key)
  }

  output$primary_pick <- renderUI({
    ap <- all_panels(); sel <- input$panels
    req(length(sel))
    ch <- setNames(sel, vapply(ap[sel], function(m) m$label, character(1)))
    cur <- isolate(input$primary)
    selectInput("primary", "Primary panel for the figures", ch,
                selected = if (!is.null(cur) && cur %in% sel) cur else sel[1])
  })

  active_panels <- reactive({
    sel <- input$panels
    validate(need(length(sel) > 0, "Tick at least one panel."))
    all_panels()[intersect(sel, names(all_panels()))]
  })
  primary <- reactive({
    ap <- all_panels(); p <- input$primary
    if (is.null(p) || !p %in% names(ap)) p <- names(active_panels())[1]
    ap[[p]]
  })

  scored <- reactive({
    X <- expr_mat(); ap <- active_panels()       # validation messages pass through
    tryCatch(score_panels(ap, X), error = function(e) { notify_err(e); NULL })
  })
  primary_scored <- reactive({
    s <- scored(); req(s)
    k <- primary()$key
    if (!is.null(s[[k]])) s[[k]] else s[[1]]
  })

  # ---- labels -----------------------------------------------------------------------------
  output$label_picks <- renderUI({
    md <- meta_raw()
    if (is.null(md)) return(helpText("Add metadata to group samples or score them against known labels."))
    cols <- setdiff(names(md), "SampleID")
    pick <- function(cands) { h <- cands[cands %in% cols]; if (length(h)) h[1] else "" }
    # 47 conditions across 18 cohorts make an unreadable legend; the validation
    # cohorts open grouped by their known phenotype instead.
    grp <- switch(input$src %||% "", validation = c("Truth_Phenotype", "Condition"),
                  example = "Condition", c("Condition", "Group", "Phenotype", "Dataset", "Study"))
    tagList(
      selectInput("group_col", "Group samples by", c("(none)" = "", col_choices(cols)), selected = pick(grp)),
      selectInput("truth_col", "Known phenotype (to score accuracy)", c("(none)" = "", col_choices(cols)),
                  selected = pick(c("Truth_Phenotype", "Truth", "Phenotype", "Expected_Phenotype"))),
      selectInput("ds_col", "Dataset / study", c("(none)" = "", col_choices(cols)),
                  selected = pick(c("Dataset", "Study", "Cohort"))),
      checkboxInput("show_hulls", "Shade each group's area in the phenotype space", FALSE))
  })

  col_vec <- function(col) {
    md <- meta_df()
    if (is.null(col) || !nzchar(col) || !col %in% names(md)) return(NULL)
    v <- as.character(md[[col]]); v
  }
  group_vec <- reactive({
    g <- col_vec(input$group_col)
    if (is.null(g)) g <- rep("All samples", ncol(expr_mat()))
    g[is.na(g) | !nzchar(g)] <- "Unlabelled"
    g
  })
  truth_vec <- reactive({
    t <- col_vec(input$truth_col); if (is.null(t)) return(NULL)
    t[!nzchar(t)] <- NA
    t
  })
  ds_vec <- reactive(col_vec(input$ds_col))
  group_label <- reactive({
    if (is.null(input$group_col) || !nzchar(input$group_col)) "Your samples" else pretty_col(input$group_col)
  })

  output$coverage_msg <- renderUI({
    s <- scored(); req(s)
    cov <- vapply(s, function(x) min(x$proj$coverage_per_sample), numeric(1))
    labs <- vapply(s, function(x) x$panel$label, character(1))
    low <- cov < 0.85
    div(class = if (any(low)) "ecmx-warn" else "ecmx-note",
        "Ratio coverage, lowest sample (both genes present): ",
        paste(sprintf("%s %.1f%%", labs, 100 * cov), collapse = " · "), ". ",
        if (any(low))
          strong("Below 85% for at least one panel — every manuscript validation cohort had ≥ 90%. Interpret with caution.")
        else "Ratios with a missing gene vote 0, as in the manuscript.")
  })

  # Fibroblast-identity check (Methods 6.1). The bundled cohorts carry panel genes
  # only, so it runs on uploads.
  output$lineage_msg <- renderUI({
    if (!identical(input$src, "upload")) return(NULL)
    ls <- ecmx_lineage_scores(expr_mat())
    if (is.null(ls)) return(div(class = "ecmx-note",
      "Fibroblast identity check skipped: too few lineage markers in your file ",
      "(it needs at least half of the 12 fibroblast markers and three other lineages)."))
    low <- ls$Fibroblast_margin < 0.54
    if (!any(low)) return(div(class = "ecmx-note", sprintf(
      "Fibroblast identity check: fibroblast markers dominate in all %d samples (smallest margin %.2f).",
      nrow(ls), min(ls$Fibroblast_margin))))
    worst <- ls[low, ][order(ls$Fibroblast_margin[low]), ]
    div(class = "ecmx-warn",
        strong(sprintf("Fibroblast identity check: %d of %d samples look less fibroblast-like than any training sample. ",
                       sum(low), nrow(ls))),
        "ECMExplainer was trained on fibroblasts; in all 96 training samples the fibroblast marker score ",
        "exceeded every other lineage by at least 0.54 log2-CPM. ",
        sprintf("Lowest: %s.", paste(sprintf("%s (margin %.2f, top other lineage %s)",
                                             head(worst$SampleID, 5), head(worst$Fibroblast_margin, 5),
                                             gsub("_", " ", head(worst$Top_other, 5))), collapse = "; ")),
        " Interpret their calls with caution.",
        if (any(worst$Top_other == "Smooth_muscle"))
          " A high smooth-muscle score can also mean myofibroblast activation, as in the training Fibrosis samples.")
  })

  # ---- figures ------------------------------------------------------------------------
  # On a phone-width screen (plot narrower than 600 px) every figure gets the
  # compact theme and a height that follows its width; downloads are unchanged.
  plot_w   <- function(id) session$clientData[[paste0("output_", id, "_width")]] %||% 1000
  on_phone <- function(id) plot_w(id) < 600
  fig_output <- function(id, make, h, phone_h = h) {
    output[[id]] <- renderPlot({
      p <- make()
      if (on_phone(id) && inherits(p, "ggplot")) p + ecmx_compact() else p
    }, height = function() {
      if (!on_phone(id)) return(h)
      if (is.function(phone_h)) phone_h(plot_w(id)) else phone_h
    })
  }

  make_lda <- function(point_scale = 1) {
    s <- primary_scored(); m <- s$panel
    validate(need(ncol(m$train_coords) >= 2, "This panel has fewer than 3 classes, so no 2-D space."))
    ecmx_fig_lda(m$train_coords, m$train_phenotype, m$centroids,
                 new_coords = s$proj$coords, new_label = group_vec(),
                 group_title = group_label(),
                 show_groups = isTRUE(input$show_hulls),
                 title = sprintf("%s phenotype space", m$label), point_scale = point_scale)
  }
  fig_lda <- reactive(make_lda())
  fig_post <- reactive({
    s <- primary_scored()
    ecmx_fig_posteriors(s$proj$posterior, colnames(expr_mat()), group_vec()) +
      labs(subtitle = s$panel$label)
  })
  fig_dist <- reactive({
    s <- primary_scored()
    ecmx_fig_distance_heatmap(ecmx_distances(s$proj$coords, s$panel$centroids), group_vec()) +
      labs(subtitle = s$panel$label)
  })
  fig_heat <- reactive({
    m <- primary()
    sel <- if (is.null(input$heat_pheno) || "__all__" %in% input$heat_pheno) NULL else input$heat_pheno
    ecmx_fig_panel_heatmap(m$sign_matrix_train, m$train_phenotype, m$features, selecting = sel)
  })
  preds_long <- reactive({
    s <- scored(); req(s)
    do.call(rbind, lapply(s, function(x) data.frame(
      SampleID = colnames(expr_mat()), Group = group_vec(), Panel = x$panel$label,
      Predicted = as.character(x$proj$class), stringsAsFactors = FALSE)))
  })
  fig_cmp <- reactive(ecmx_fig_panel_agreement(preds_long()))

  acc_tables <- reactive({
    t <- truth_vec(); req(t)
    s <- scored(); req(s)
    validate(need(any(!is.na(t)), "The known-phenotype column is empty for these samples."))
    validation_tables(s, t, ds_vec())
  })
  fig_acc_radar <- reactive(ecmx_fig_recall_radar(acc_tables()$recall))
  fig_acc_conf  <- reactive(ecmx_fig_confusion(acc_tables()$confusion))
  fig_acc_ds    <- reactive({ d <- acc_tables()$per_dataset; req(d); ecmx_fig_dataset_accuracy(d) })

  output$heat_pheno_pick <- renderUI({
    f <- primary()$features
    ph <- sort(unique(stats::na.omit(f$Selecting_Phenotype)))
    n  <- table(f$Selecting_Phenotype)
    lbl <- setNames(c("__all__", ph),
                    c(sprintf("All (%d ratios, thinned)", nrow(f)),
                      sprintf("%s (%d ratios)", ph, as.integer(n[ph]))))
    selectizeInput("heat_pheno", "Ratios selecting for", choices = lbl,
                   selected = "__all__", multiple = TRUE)
  })

  fig_output("p_lda", function() make_lda(if (on_phone("p_lda")) 0.6 else 1), 580, function(w) round(w * 0.66 + 200 + 17 * max(0, length(unique(group_vec())) - 6)))
  fig_output("p_post", function() fig_post(), 560)
  fig_output("p_dist", function() fig_dist(), 440)
  fig_output("p_heat", function() fig_heat(), 820)
  fig_output("p_cmp", function() fig_cmp(), 460)

  output$cmp_msg <- renderUI({
    s <- scored(); req(s)
    if (length(s) < 2) return(info_box("Tick two or more panels to compare them."))
    NULL
  })
  output$t_agree <- renderTable({
    p <- preds_long(); req(length(unique(p$Panel)) > 1)
    w <- reshape(p[, c("SampleID", "Panel", "Predicted")], idvar = "SampleID",
                 timevar = "Panel", direction = "wide")
    pc <- grep("^Predicted\\.", names(w))
    agree <- apply(w[, pc, drop = FALSE], 1, function(x) length(unique(x)) == 1)
    g <- group_vec()[match(w$SampleID, colnames(expr_mat()))]
    tab <- aggregate(agree, by = list(Group = g), FUN = function(x) c(n = length(x), all = sum(x)))
    data.frame(Group = tab$Group, Samples = tab$x[, "n"],
               `All panels agree` = tab$x[, "all"],
               `Agreement (%)` = round(100 * tab$x[, "all"] / tab$x[, "n"], 1), check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  output$acc_body <- renderUI({
    if (is.null(truth_vec())) return(info_box(
      "Choose a ", strong("known phenotype"), " column under 4 · Labels to score the panels. ",
      "Samples with an empty label are skipped. The manuscript validation cohorts come with one."))
    tagList(
      helpText("Accuracy counts a sample as correct when the predicted phenotype equals its label. ",
               "Balanced accuracy averages recall over the labelled phenotypes."),
      div(class = "ecmx-scroll", tableOutput("t_acc")),
      fluidRow(column(6, plotOutput("p_acc_radar", height = "auto")),
               column(6, plotOutput("p_acc_ci", height = "auto"))),
      plotOutput("p_acc_conf", height = "auto"),
      div(class = "ecmx-scroll", tableOutput("t_mcnemar")),
      if (!is.null(ds_vec())) plotOutput("p_acc_ds", height = "auto"))
  })
  output$t_acc <- renderTable({
    s <- acc_tables()$summary
    fmt_tbl(s[, c("Panel", "N_scored", "Micro_accuracy_pct", "Micro_CI95_low", "Micro_CI95_high",
                  "Balanced_accuracy_pct", "Control_specificity_pct", "Stress_sensitivity_pct",
                  "N_datasets", "Mean_perDataset_accuracy_pct")])
  }, striped = TRUE, spacing = "xs", na = "–")
  output$t_mcnemar <- renderTable({ m <- acc_tables()$mcnemar; req(m); fmt_tbl(m) },
                                  striped = TRUE, spacing = "xs")
  fig_output("p_acc_radar", function() fig_acc_radar(), 430)
  fig_output("p_acc_ci", function() ecmx_fig_accuracy(acc_tables()$summary), 430)
  fig_output("p_acc_conf", function() fig_acc_conf(), 380)
  fig_output("p_acc_ds", function() fig_acc_ds(), 540)

  current_fig <- reactive({
    switch(input$viz %||% "lda",
           lda = fig_lda(), post = fig_post(), dist = fig_dist(), heat = fig_heat(),
           cmp = fig_cmp(), acc = fig_acc_radar(), NULL)
  })
  output$dl_plot <- downloadHandler(
    filename = function() sprintf("ECMExplainer_%s_%s.%s", input$viz %||% "figure",
                                  Sys.Date(), input$exp_fmt),
    content = function(file) {
      p <- current_fig()
      validate(need(!is.null(p), "This tab has no figure to export."))
      w <- max(3, min(30, input$exp_w %||% 9))
      h <- max(3, min(40, input$exp_h %||% 7))
      d <- max(72, min(600, input$exp_dpi %||% 300))
      ggsave(file, p, width = w, height = h, dpi = d,
             device = if (identical(input$exp_fmt, "pdf")) grDevices::cairo_pdf else "png",
             limitsize = FALSE, bg = "white")
    })

  results_tbl <- reactive({
    s <- scored(); req(s)
    ids <- colnames(expr_mat())
    base <- data.frame(SampleID = ids, Group = group_vec(), stringsAsFactors = FALSE)
    if (!is.null(truth_vec())) base$Known_phenotype <- truth_vec()
    for (x in s) {
      k <- x$panel$key; cf <- x$conf
      base[[paste0("Pred_", k)]]    <- as.character(x$proj$class)
      base[[paste0("MaxPost_", k)]] <- round(cf$MaxPost, 4)
      base[[paste0("Margin_", k)]]  <- round(cf$Margin, 4)
      base[[paste0("Entropy_", k)]] <- round(cf$Entropy, 4)
      base[[paste0("Closest_", k)]] <- cf$Closest
      base[[paste0("MinDist_", k)]] <- round(cf$MinDist, 4)
      base[[paste0("Mahal_", k)]]   <- round(cf$Mahal, 3)
      base[[paste0("OOD95_", k)]]   <- cf$OOD95
      post <- round(as.data.frame(x$proj$posterior), 4)
      names(post) <- paste0("P_", k, "_", names(post))
      base <- cbind(base, post)
    }
    base
  })
  output$t_res <- renderTable({
    r <- results_tbl()
    keep <- c("SampleID", "Group", intersect("Known_phenotype", names(r)),
              grep("^(Pred|MaxPost|Margin|Closest|MinDist)_", names(r), value = TRUE))
    if ("Known_phenotype" %in% names(r))
      r$Known_phenotype[is.na(r$Known_phenotype)] <- "not scored"
    r[, keep, drop = FALSE]
  }, striped = TRUE, spacing = "xs", digits = 3, na = "—")
  output$dl_results <- downloadHandler(
    filename = function() paste0("ECMExplainer_results_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(results_tbl(), f, row.names = FALSE))

  # ---- trajectories -------------------------------------------------------------------------
  output$traj_group <- renderUI({
    md <- meta_raw()
    validate(need(!is.null(md), "Add sample metadata to build trajectories."))
    cols <- setdiff(names(md), "SampleID")
    selectInput("tr_group_col", "Condition column", col_choices(cols),
                selected = intersect(c("Condition", "Phenotype", cols), cols)[1])
  })
  output$traj_anchors <- renderUI({
    ph <- rownames(primary()$centroids)
    tagList(
      selectInput("anchor_a", "Anchor A", ph,
                  selected = if ("Proliferating" %in% ph) "Proliferating" else ph[1]),
      selectInput("anchor_b", "Anchor B", ph,
                  selected = if ("Senescence" %in% ph) "Senescence" else ph[min(2, length(ph))]))
  })
  # One dataset at a time by default: the 18 validation cohorts carry about 40
  # conditions, which make an unreadable trajectory.
  output$traj_ds_pick <- renderUI({
    d <- ds_vec(); if (is.null(d)) return(NULL)
    ds <- sort(unique(stats::na.omit(d)))
    if (length(ds) < 2) return(NULL)
    selectInput("traj_ds", "Dataset", c("All datasets" = "__all__", ds),
                selected = if ("SRP117883" %in% ds) "SRP117883" else ds[1])
  })
  traj_idx <- reactive({
    d <- ds_vec(); n <- ncol(expr_mat())
    if (is.null(d) || is.null(input$traj_ds) || input$traj_ds == "__all__") return(seq_len(n))
    which(d == input$traj_ds)
  })
  traj_groups <- reactive({
    req(input$tr_group_col)
    col_vec(input$tr_group_col)[traj_idx()]
  })
  output$traj_order <- renderUI({
    g <- unique(stats::na.omit(traj_groups())); req(length(g))
    selectizeInput("traj_seq", "Trajectory order (pick in sequence)",
                   choices = g, selected = g, multiple = TRUE)
  })
  output$traj_from_to <- renderUI({
    g <- unique(stats::na.omit(traj_groups())); req(length(g) >= 2)
    tagList(selectInput("prot_from", "From (origin)", g, selected = g[1]),
            selectInput("prot_to",   "To (endpoint)", g, selected = g[2]))
  })
  fig_output("p_traj", function() {
    s <- primary_scored(); req(input$anchor_a, input$anchor_b)
    validate(need(input$anchor_a != input$anchor_b, "Choose two different anchors."))
    cen <- s$panel$centroids
    req(all(c(input$anchor_a, input$anchor_b) %in% rownames(cen)))
    sep <- sqrt(sum((cen[input$anchor_a, ] - cen[input$anchor_b, ])^2))
    ecmx_fig_trajectory(ecmx_distances(s$proj$coords[traj_idx(), , drop = FALSE], cen), traj_groups(),
                        input$anchor_a, input$anchor_b,
                        order_levels = input$traj_seq, centroid_sep = sep)
  }, 480, function(w) round(w * 0.85 + 200 + 17 * max(0, length(unique(stats::na.omit(traj_groups()))) - 6)))
  angles_tbl <- reactive({
    s <- primary_scored(); req(input$prot_from, input$prot_to)
    g <- traj_groups(); xy <- s$proj$coords[traj_idx(), , drop = FALSE]
    from <- xy[which(g == input$prot_from), , drop = FALSE]
    to   <- xy[which(g == input$prot_to),   , drop = FALSE]
    validate(need(nrow(from) >= 2 && nrow(to) >= 2,
                  "The protractor needs at least 2 samples in each group."))
    ecmx_trajectory_angles(from, to, s$panel$centroids, n_boot = input$n_boot)
  })
  fig_output("p_prot", function() ecmx_fig_protractor(angles_tbl(), sprintf("%s → %s", input$prot_from, input$prot_to)), 460)
  output$t_angles <- renderTable({
    a <- angles_tbl()
    data.frame(Phenotype = a$Phenotype, `Angle (deg)` = round(a$Angle, 1),
               `95% CI` = ifelse(is.na(a$Lower), "-", sprintf("%.0f – %.0f", a$Lower, a$Upper)),
               check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  # ---- retrain ----------------------------------------------------------------------------------
  tr_data <- reactive({
    if (identical(input$tr_src, "builtin")) {
      return(list(counts = readRDS(TRAIN_CNT),
                  meta = read.csv(TRAIN_MET, stringsAsFactors = FALSE), builtin = TRUE))
    }
    req(input$tr_counts, input$tr_meta)
    cnt <- tryCatch(ecmx_read_expression(input$tr_counts$datapath),
                    error = function(e) { notify_err(e); NULL })
    md <- tryCatch(read_meta(input$tr_meta$datapath), error = function(e) { notify_err(e); NULL })
    req(cnt, md)
    list(counts = cnt, meta = md, builtin = FALSE)
  })
  output$tr_cols <- renderUI({
    md <- tr_data()$meta; cols <- setdiff(names(md), "SampleID")
    fluidRow(
      column(6, selectInput("tr_pheno_col", "Label column", cols,
                            selected = intersect(c("Phenotype", cols), cols)[1])),
      column(6, selectInput("tr_study_col", "Study column", c("(none)" = "", cols),
                            selected = intersect(c("Study", "Dataset", ""), c(cols, ""))[1])))
  })
  output$tr_filter <- renderUI(filter_slots_ui("tf", tr_data()$meta))
  output$tf_val1_ui <- renderUI(filter_value_ui(input, "tf", 1, tr_data()$meta))
  output$tf_val2_ui <- renderUI(filter_value_ui(input, "tf", 2, tr_data()$meta))
  tr_selected <- reactive({
    d <- tr_data(); keep <- apply_filters(input, "tf", d$meta)
    ids <- intersect(colnames(d$counts), keep)
    req(input$tr_pheno_col)
    md <- d$meta[match(ids, d$meta$SampleID), , drop = FALSE]
    lab <- as.character(md[[input$tr_pheno_col]])
    ok <- !is.na(lab) & nzchar(lab)
    md <- md[ok, , drop = FALSE]
    data.frame(SampleID = md$SampleID, Phenotype = lab[ok],
               Study = if (!is.null(input$tr_study_col) && nzchar(input$tr_study_col))
                 as.character(md[[input$tr_study_col]]) else NA_character_,
               stringsAsFactors = FALSE)
  })
  output$tr_n_msg <- renderUI({
    s <- tr_selected(); n <- table(s$Phenotype)
    helpText(sprintf("%d samples, %d classes: %s", nrow(s), length(n),
                     paste0(names(n), " (", as.integer(n), ")", collapse = ", ")))
  })

  tr_genelist <- reactive({
    switch(input$tr_genes,
           curated = GENELIST,
           upload_list = { req(input$tr_genelist); read_genelist(input$tr_genelist$datapath) },
           paste = { req(nzchar(input$tr_genetext %||% "")); read_genelist(text = input$tr_genetext) })
  })

  train_task <- ExtendedTask$new(function(args) {
    if (ASYNC) {
      promises::future_promise({
        source(file.path(app_dir, "R/ecm_core.R"))
        run_training(args)
      }, globals = list(args = args, app_dir = APP_DIR, run_training = run_training),
         packages = "MASS", seed = TRUE)
    } else {
      promises::promise_resolve(run_training(args))
    }
  })

  observeEvent(input$do_train, {
    tryCatch({
      d <- tr_data(); s <- tr_selected()
      if (!d$builtin && nrow(s) > LIM$max_train)
        stop(sprintf("This instance trains on at most %d samples; you selected %d. Run locally for unlimited training.",
                     LIM$max_train, nrow(s)))
      if (nrow(s) < 4) stop("Select at least 4 samples to train.")
      meta <- s; if (all(is.na(meta$Study))) meta$Study <- NULL
      args <- list(
        prep = list(counts = d$counts[, s$SampleID, drop = FALSE], meta = meta,
                    genelist = if (input$tr_design == "ALLGENE") NULL else tr_genelist(),
                    design = input$tr_design,
                    top_n = as.integer(input$tr_topn %||% 2000L)),
        train = list(design = input$tr_design, n_per_pheno = as.integer(input$tr_npp),
                     eta_min = input$tr_eta, rank_by = input$tr_rank,
                     min_per_module = as.integer(input$tr_minmod)),
        loso = isTRUE(input$tr_loso) && !is.null(meta$Study))
      train_task$invoke(args)
    }, error = notify_err)
  })

  observeEvent(train_task$result(), {
    fit <- tryCatch(train_task$result(), error = function(e) { notify_err(e); NULL })
    req(fit)
    key <- sprintf("custom_%d", length(session_models()) + 1L)
    fit$label <- sprintf("Your %s-%d (%s)", fit$design, fit$panel_size, format(Sys.time(), "%H:%M"))
    fit <- prepare_panel(fit, fit$label, key)
    trained(fit)
    add_model(fit, key)
    showNotification(sprintf("Trained %s: LOOCV %.1f%%%s. Added to the Classify tab.",
                             fit$label, fit$loocv_acc,
                             if (!is.null(fit$loso)) sprintf(", nested LOSO %.1f%%", fit$loso$accuracy) else ""),
                     type = "message", duration = 10)
  })

  observeEvent(input$model_file, {
    tryCatch({
      m <- readRDS(input$model_file$datapath)
      need_fields <- c("pca", "lda", "features", "n_pcs", "train_phenotype")
      if (!is.list(m) || !all(need_fields %in% names(m)))
        stop("This file is not an ECMExplainer model (it lacks ", paste(need_fields, collapse = ", "), ").")
      key <- sprintf("custom_%d", length(session_models()) + 1L)
      lab <- sprintf("Loaded: %s", tools::file_path_sans_ext(input$model_file$name))
      add_model(prepare_panel(m, lab, key), key)
      showNotification(paste(lab, "added to the Classify tab."), type = "message")
    }, error = notify_err)
  })

  output$train_status <- renderUI({
    if (identical(train_task$status(), "running"))
      return(info_box(strong("Training… "), "ALLGENE and nested LOSO take longest ",
                      "(about 5 s and 30 s on a laptop). The rest of the app stays usable."))
    m <- trained()
    if (is.null(m)) return(info_box(
      "Choose data, a gene universe and a design, then train. The model is added to the ",
      strong("Classify"), " tab, where you can test it on the validation cohorts or your own data. ",
      "With the manuscript cohort and default settings, ADSIG, POOL and ALLGENE reproduce the ",
      "published ADSIG-180, POOL-180 and ALLGENE-180 panels exactly."))
    same <- names(Filter(function(p) setequal(p$features$Ratio_ID, m$features$Ratio_ID), PANELS))
    tagList(
      div(class = "ecmx-note",
          strong(sprintf("%s — %d ratios, %d PCs, %d classes, n = %d",
                         m$label, m$panel_size, m$n_pcs, length(m$phenotypes), m$n_training)),
          br(), sprintf("LOOCV %.2f%% · macro-F1 %.3f · ranking %s · min per module %d",
                        m$loocv_acc, m$macro_f1, m$rank_by, m$min_per_module)),
      if (length(same)) div(class = "ecmx-ok", sprintf(
        "Identical ratio set to the published %s.", PANELS[[same[1]]]$label)))
  })
  output$loso_box <- renderUI({
    m <- trained(); req(m)
    if (is.null(m$loso)) return(helpText("Nested leave-one-study-out not run (needs a Study column with 2+ studies)."))
    l <- m$loso
    tagList(
      h5("Nested leave-one-study-out (honest: each study held out of selection and fitting)"),
      div(class = "ecmx-note", sprintf(
        "Accuracy %.2f%% · restricted to phenotypes seen in ≥ 3 studies %.2f%% · balanced %.2f%% · macro-F1 %.3f",
        l$accuracy, l$restricted, l$balanced, l$macro_f1)),
      renderTable(l$per_class, striped = TRUE, spacing = "xs"))
  })
  output$t_perclass  <- renderTable({ req(trained()); trained()$per_class },
                                    striped = TRUE, spacing = "xs")
  output$t_confusion <- renderTable({ req(trained()); as.data.frame.matrix(trained()$confusion) },
                                    rownames = TRUE, striped = TRUE, spacing = "xs")
  output$t_panel <- renderTable({
    req(trained()); f <- trained()$features
    f[order(f$Selecting_Phenotype, -f$TSP_Delta),
      c("Ratio_ID", "Module", "Selecting_Phenotype", "Direction", "TSP_Delta", "eta_sq")]
  }, striped = TRUE, spacing = "xs", digits = 3)
  fig_output("p_train_lda", function() {
    m <- trained(); req(m)
    validate(need(ncol(m$train_coords) >= 2, "Need at least 3 classes to draw a 2-D LDA space."))
    ecmx_fig_lda(m$train_coords, m$train_phenotype, m$centroids,
                 title = sprintf("%s training space", m$label))
  }, 480, function(w) round(w * 0.66 + 170))
  output$dl_model <- downloadHandler(
    filename = function() sprintf("ECMExplainer_%s_%d_%s.rds", trained()$design,
                                  trained()$panel_size, Sys.Date()),
    content  = function(f) { req(trained()); saveRDS(trained(), f) })

  # ---- aging score ------------------------------------------------------------------------------
  output$age_col_pick <- renderUI({
    md <- meta_df()
    num <- names(md)[vapply(md, is.numeric, logical(1))]
    validate(need(length(num) > 0, "To fit your own weights, add sample metadata with a numeric age column."))
    selectInput("age_col", "Age column", num,
                selected = grep("age", num, ignore.case = TRUE, value = TRUE)[1])
  })
  observeEvent(input$do_age, {
    tryCatch({
      md <- meta_df(); req(input$age_col)
      m <- expr_mat()
      ages <- suppressWarnings(as.numeric(md[[input$age_col]]))
      keep <- !is.na(ages)
      if (sum(keep) < 5) stop("Need at least 5 samples with an age value.")
      b <- ecmx_fit_age_betas(m[, keep, drop = FALSE], ages[keep],
                              ecmx_model_ratios(primary()), fdr = input$age_fdr)
      if (is.null(b) || !nrow(b)) stop("No ratio passed the FDR threshold.")
      betas(b)
      showNotification(sprintf("%d ratios passed FDR < %.3f", nrow(b), input$age_fdr), type = "message")
    }, error = notify_err)
  })
  age_betas <- reactive(if (identical(input$age_src, "fit")) betas() else AGE_WEIGHTS)
  age_scores <- reactive({
    b <- age_betas(); req(b); m <- expr_mat()
    data.frame(SampleID = colnames(m), Group = age_group_vec(),
               ECM_Aging_Score = ecmx_aging_score(m, b), stringsAsFactors = FALSE)
  })
  output$age_group_pick <- renderUI({
    md <- meta_raw(); if (is.null(md)) return(NULL)
    cols <- setdiff(names(md), "SampleID")
    cur <- isolate(input$group_col)
    selectInput("age_group", "Group scores by", c("(none)" = "", col_choices(cols)),
                selected = if (!is.null(cur) && cur %in% cols) cur else "")
  })
  age_group_vec <- reactive({
    g <- col_vec(input$age_group) %||% rep("All samples", ncol(expr_mat()))
    g[is.na(g) | !nzchar(g)] <- "Unlabelled"
    g
  })
  output$age_status <- renderUI({
    b <- age_betas()
    if (is.null(b)) return(info_box("Choose an age column and fit weights, or switch back to the published weights."))
    genes <- unique(trimws(unlist(strsplit(b$Ratio_ID, "/", fixed = TRUE))))
    have <- mean(genes %in% rownames(expr_raw()))
    div(class = if (have < 0.85) "ecmx-warn" else "ecmx-note",
        if (identical(input$age_src, "fit"))
          sprintf("Your weights: %d age-associated ratios at BH FDR < %.3f. ", nrow(b), input$age_fdr)
        else sprintf("Published weights: %d ratios (SI Table 4). ", nrow(b)),
        sprintf("%.0f%% of their genes are in your data; ratios with a missing gene vote 0.", 100 * have))
  })
  fig_output("p_age", function() {
    a <- age_scores()
    ref <- if (!identical(input$age_src, "fit") && isTRUE(input$age_ref)) AGE_REF
    ecmx_fig_aging_score(a$ECM_Aging_Score, a$Group, reference = ref)
  }, 460)
  output$t_betas <- renderTable({
    b <- age_betas(); req(b); b <- b[order(-abs(b$beta)), ]
    fdr <- if (!is.null(b$fdr)) b$fdr else b$FDR
    head(data.frame(Ratio = b$Ratio_ID, Weight = formatC(b$beta, format = "g", digits = 3),
                    FDR = formatC(fdr, format = "g", digits = 2), check.names = FALSE), 25)
  }, striped = TRUE, spacing = "xs")
  output$dl_age_scores <- downloadHandler(
    filename = function() paste0("ECM_aging_scores_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(age_scores(), f, row.names = FALSE))
  output$dl_betas <- downloadHandler(
    filename = function() paste0("ECM_age_weights_", Sys.Date(), ".csv"),
    content  = function(f) { req(age_betas()); write.csv(age_betas(), f, row.names = FALSE) })

  # ---- example files ---------------------------------------------------------------------------
  output$dl_ex_expr <- downloadHandler("ECMExplainer_example_expression.csv",
    function(f) file.copy("data/examples/example_expression.csv", f))
  output$dl_ex_meta <- downloadHandler("ECMExplainer_example_metadata.csv",
    function(f) file.copy("data/examples/example_metadata.csv", f))
  output$dl_ex_genes <- downloadHandler("ECMExplainer_curated_genelist.csv",
    function(f) write.csv(GENELIST[, c("Gene", "Function")], f, row.names = FALSE))

  # ---- manuscript ---------------------------------------------------------------------------
  output$ms_summary <- renderTable({
    s <- MS_TABLES$summary
    fmt_tbl(data.frame(Panel = s$Panel, `Scored` = s$N_scored,
                       `Balanced acc. (%)` = s$Balanced_accuracy_pct,
                       `Micro acc. (%)` = s$Micro_accuracy_pct,
                       `95% CI` = sprintf("%.1f–%.1f", s$Micro_CI95_low, s$Micro_CI95_high),
                       `Mean per-dataset (%)` = s$Mean_perDataset_accuracy_pct,
                       `Control specificity (%)` = s$Control_specificity_pct,
                       `Stress sensitivity (%)` = s$Stress_sensitivity_pct,
                       check.names = FALSE))
  }, striped = TRUE, spacing = "xs")
  fig_output("ms_radar", function() ecmx_fig_recall_radar(MS_TABLES$recall) +
                                  labs(title = "Figure 7B · Recall per phenotype"), 440)
  fig_output("ms_acc", function() ecmx_fig_accuracy(MS_TABLES$summary), 440)
  fig_output("ms_conf", function() ecmx_fig_confusion(MS_TABLES$confusion), 380)
  fig_output("ms_ds", function() ecmx_fig_dataset_accuracy(MS_TABLES$per_dataset), 560)
  output$ms_mcnemar <- renderTable(fmt_tbl(MS_TABLES$mcnemar), striped = TRUE, spacing = "xs")
  output$ms_training <- renderTable(fmt_tbl(PANEL_CV, 2), striped = TRUE, spacing = "xs")
  lapply(seq_along(SI_FILES), function(i) {
    output[[paste0("dl_si", i)]] <- downloadHandler(
      filename = function() names(SI_FILES)[i],
      content = function(f) file.copy(file.path("supplementary", names(SI_FILES)[i]), f))
  })
  observeEvent(input$ms_open, {
    updateRadioButtons(session, "src", selected = "validation")
    updateCheckboxInput(session, "adv_panels", value = TRUE)
    updateTabsetPanel(session, "tabs", selected = "Classify")
    updateTabsetPanel(session, "viz", selected = "acc")
  })

  output$about_md <- renderUI({
    md <- unlist(lapply(c("ABOUT.md", "METHODS.md"), function(f)
      if (file.exists(f)) readLines(f, warn = FALSE)))
    # METHODS.md opens with its own H1; under ABOUT.md's "Methods" heading, drop it.
    md <- md[!(md == "# Methods and provenance")]
    HTML(commonmark::markdown_html(paste(md, collapse = "\n")))
  })
}

shinyApp(ui, server)
