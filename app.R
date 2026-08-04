# ==============================================================================
# ECMExplainer — Shiny application
#
# Positions fibroblast transcriptomes in the ECMExplainer functional matrisome
# phenotype space, using the frozen ADSIG-180 model from the manuscript.
#
# Run locally (no limits):        shiny::runApp()
# Run as a public instance:       ECMX_MODE=web shiny::runApp()
#
# Web mode caps upload size, sample count and retraining size so a shared
# instance stays responsive. Nothing scientific differs between the modes — the
# model, the projection and the figures are identical.
# ==============================================================================
suppressPackageStartupMessages({
  library(shiny); library(ggplot2); library(MASS)
})
source("R/ecm_core.R")
source("R/ecm_figures.R")

ECMX_MODE <- tolower(Sys.getenv("ECMX_MODE", "local"))
WEB <- identical(ECMX_MODE, "web")
LIM <- list(
  upload_mb     = as.numeric(Sys.getenv("ECMX_MAX_UPLOAD_MB", if (WEB) "20"  else "2000")),
  max_samples   = as.integer(Sys.getenv("ECMX_MAX_SAMPLES",   if (WEB) "60"  else "100000")),
  max_train     = as.integer(Sys.getenv("ECMX_MAX_TRAIN",     if (WEB) "60"  else "100000")),
  max_per_pheno = as.integer(Sys.getenv("ECMX_MAX_PER_PHENO", if (WEB) "30"  else "200"))
)
options(shiny.maxRequestSize = LIM$upload_mb * 1024^2)

MODEL     <- readRDS("data/TrainedModel_ADSIG_180.rds")
GENELIST  <- read.csv("data/ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
TRAIN_CNT <- "data/salmon_counts_Combined_FINAL.rds"
TRAIN_MET <- "data/salmon_metadata_Combined_FINAL.csv"
HAS_TRAIN <- file.exists(TRAIN_CNT)

TRAIN_COORDS <- local({
  pcs <- as.data.frame(MODEL$pca$x[, seq_len(MODEL$n_pcs), drop = FALSE])
  colnames(pcs) <- paste0("PC", seq_len(MODEL$n_pcs))
  predict(MODEL$lda, pcs)$x
})
CENTROIDS <- ecmx_centroids(MODEL)

notify_err <- function(e) showNotification(conditionMessage(e), type = "error", duration = 12)
info_box   <- function(...) div(class = "ecmx-note", ...)

# base R gained %||% in 4.4.0; define it if we are on something older.
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

model_features <- function(m) {
  f <- m$features
  if (is.null(f$Ratio_ID) && !is.null(f$Ratio)) f$Ratio_ID <- f$Ratio
  if (is.null(f$Module)   && !is.null(f$Category)) f$Module <- f$Category
  if (is.null(f$eta_sq)) f$eta_sq <- NA_real_
  f
}

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  title = "ECMExplainer",
  # Styling is pure CSS: no JavaScript, no web fonts, no extra R packages, so it
  # adds nothing to load time and nothing to the dependency surface. Palette
  # carried over from the original ECMExplainer app — navy/blue gradient wordmark,
  # #2d6a9f accent, #3c8dbc card rule, #f7f9fc page.
  tags$head(tags$style(HTML("
    :root {
      --ecmx-navy:#1a3a5c; --ecmx-blue:#2d6a9f; --ecmx-light:#3c8dbc;
      --ecmx-page:#f7f9fc; --ecmx-rule:#dbe4ee;
    }
    body { background:var(--ecmx-page); }
    .container-fluid { max-width:1680px; }

    /* ---- Header ---- */
    .ecmx-head {
      background:linear-gradient(135deg, var(--ecmx-navy) 0%, var(--ecmx-blue) 100%);
      color:#fff; padding:16px 22px; margin:-15px -15px 18px -15px;
      border-bottom:3px solid var(--ecmx-light);
    }
    .ecmx-head .wordmark { font-size:1.6em; font-weight:700; letter-spacing:.05em; }
    .ecmx-head .authors {
      font-size:.72em; font-weight:400; font-style:italic;
      color:#e4ebf2; margin-left:8px; white-space:nowrap;
    }
    .ecmx-head .sub { color:#e4ebf2; font-size:12.5px; margin-top:5px; }

    /* ---- Cards ---- */
    .well {
      background:#fff; border:1px solid var(--ecmx-rule);
      border-top:3px solid var(--ecmx-light); border-radius:4px;
      box-shadow:0 1px 2px rgba(26,58,92,.06);
    }
    .well h4 { color:var(--ecmx-navy); font-weight:700; font-size:14px;
               letter-spacing:.02em; margin-top:14px; }
    .well h4:first-child { margin-top:0; }
    .well h5 { color:var(--ecmx-blue); font-weight:700; font-size:13px; }
    .tab-content { background:#fff; border:1px solid var(--ecmx-rule);
                   border-top:none; border-radius:0 0 4px 4px; padding:6px 14px 14px; }
    .nav-tabs { border-bottom:1px solid var(--ecmx-rule); }
    .nav-tabs > li > a { color:#5a6b7d; font-weight:600; font-size:13px; }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      color:var(--ecmx-blue); border-top:2px solid var(--ecmx-blue);
    }

    /* ---- Callouts ---- */
    .ecmx-note { background:#f0f7ff; border:1px solid #b8d4f0;
                 border-left:4px solid var(--ecmx-blue); padding:9px 12px;
                 font-size:12.5px; margin:10px 0; border-radius:3px; color:#33475b; }
    .ecmx-warn { background:#fff8e1; border:1px solid #f5dfa8;
                 border-left:4px solid #f39c12; padding:9px 12px;
                 font-size:12.5px; margin:10px 0; border-radius:3px; color:#6b5320; }

    /* ---- Controls ---- */
    .btn-primary { background:var(--ecmx-blue); border-color:var(--ecmx-navy); }
    .btn-primary:hover, .btn-primary:focus { background:var(--ecmx-navy); border-color:var(--ecmx-navy); }
    .btn-default { border-color:#c3d1e0; color:var(--ecmx-navy); }
    .btn-default:hover { border-color:var(--ecmx-blue); color:var(--ecmx-blue); }
    label { color:#33475b; font-weight:600; font-size:12.5px; }
    .help-block, .shiny-input-container > .help-block { font-size:11.5px; color:#616e78; }
    .shiny-output-error-validation { color:#8a6d1f; font-size:13px; }
    table.data, .table { font-size:12px; }
    .table > thead > tr > th { color:var(--ecmx-navy); border-bottom:2px solid var(--ecmx-rule); }
  "))),

  div(class = "ecmx-head",
      div(span(class = "wordmark", "ECMExplainer"),
          span(class = "authors", "— Eleanor Sheekey, Roberto A. Avelar & Hironobu Fujiwara")),
      div(class = "sub",
          "A functional matrisome framework for decoding cellular stress states and aging trajectories")),

  if (WEB) div(class = "ecmx-warn",
    strong("Public demo instance. "),
    sprintf("Uploads are capped at %g MB / %d samples, and retraining at %d samples. ",
            LIM$upload_mb, LIM$max_samples, LIM$max_train),
    "Uploaded files are held in memory for your session only and are never written to disk. ",
    "For unrestricted use, run the app locally — see the Methods tab."),

  tabsetPanel(
    id = "tabs",

    tabPanel(
      "Classify",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("1 · Expression data"),
          radioButtons("src", NULL,
                       c("Bundled demo cohort (24 samples)" = "demo",
                         "Upload my own" = "upload")),
          conditionalPanel(
            "input.src == 'upload'",
            fileInput("expr_file", "Expression CSV", accept = c(".csv", ".gz")),
            info_box("Genes in the first column, one sample per remaining column. ",
                     "Raw counts or log2-CPM both work — the app detects which. ",
                     strong("Length-normalised units (TPM/FPKM/RPKM) are not valid"),
                     ", because per-gene length factors do not cancel inside the ratio."),
            fileInput("meta_file", "Sample metadata CSV (optional)", accept = ".csv"),
            helpText("SampleID plus any grouping columns.")
          ),
          h4("2 · Model"),
          uiOutput("model_pick"),
          hr(),
          uiOutput("group_pick"),
          conditionalPanel("input.viz == 'heat'", uiOutput("heat_pheno_pick")),
          hr(),
          h4("3 · Export"),
          fluidRow(
            column(6, numericInput("exp_w", "Width (in)",  9, 3, 30, 0.5)),
            column(6, numericInput("exp_h", "Height (in)", 7, 3, 40, 0.5))),
          numericInput("exp_dpi", "DPI", 300, 72, 600, 1),
          radioButtons("exp_fmt", NULL, c(PNG = "png", PDF = "pdf"), inline = TRUE),
          downloadButton("dl_plot", "Download current figure", class = "btn-sm"),
          br(), br(),
          downloadButton("dl_results", "Download results (CSV)", class = "btn-sm")
        ),
        mainPanel(
          width = 9,
          uiOutput("coverage_msg"),
          tabsetPanel(
            id = "viz",
            tabPanel("Phenotype space",    value = "lda",  br(),
                     plotOutput("p_lda",  height = "560px")),
            tabPanel("Posteriors",         value = "post", br(),
                     plotOutput("p_post", height = "440px")),
            tabPanel("Centroid distances", value = "dist", br(),
                     plotOutput("p_dist", height = "420px")),
            tabPanel("Panel heatmap",      value = "heat", br(),
                     plotOutput("p_heat", height = "760px")),
            tabPanel("Results table",      value = "tbl",  br(), tableOutput("t_res"))
          )
        )
      )
    ),

    tabPanel(
      "Trajectories",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          info_box("Methods 5.6–5.7: each condition is represented by its ",
                   strong("medoid"), " (a real sample, so the plotted point respects ",
                   "the LDA triangle inequality), and protractor angles carry a ",
                   "2,000-resample bootstrap CI."),
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
          plotOutput("p_traj", height = "480px"),
          hr(),
          fluidRow(column(7, plotOutput("p_prot", height = "460px")),
                   column(5, br(), tableOutput("t_angles")))
        )
      )
    ),

    tabPanel(
      "Retrain",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Training data"),
          radioButtons("tr_src", NULL,
                       c("Manuscript cohort (96 samples)" = "builtin",
                         "Upload my own" = "upload")),
          conditionalPanel("input.tr_src == 'upload'",
            fileInput("tr_counts", "Counts CSV", accept = c(".csv", ".gz")),
            fileInput("tr_meta", "Metadata CSV (SampleID, Phenotype, Study)", accept = ".csv"),
            info_box("Phenotype may contain ", strong("any labels you like"),
                     " — the six reference states are not required. Every class ",
                     "needs at least 2 samples.")),
          hr(),
          h4("Panel"),
          selectInput("tr_design", "Design",
                      c("ADSIG (A/D + SIG/SIG — the manuscript model)" = "ADSIG",
                        "POOL (unstructured control)" = "POOL")),
          numericInput("tr_npp", "Ratios per phenotype", 30, 3, LIM$max_per_pheno, 1),
          selectInput("tr_rank", "Selection ranking",
                      c("TSP delta (default)" = "tsp_delta",
                        "eta-squared" = "eta2")),
          helpText(paste("TSP delta ranks by one-versus-rest separation and reproduces",
                         "the published panel; eta-squared ranks by overall",
                         "between-phenotype variance.")),
          numericInput("tr_minmod", "Minimum ratios per module", 0, 0, 30, 1),
          helpText("Guarantees both A/D and SIG are represented. 0 = no floor."),
          numericInput("tr_eta", "eta-squared floor", 0.10, 0, 0.9, 0.01),
          hr(),
          actionButton("do_train", "Train model", class = "btn-primary"),
          br(), br(),
          downloadButton("dl_model", "Download model bundle (.rds)", class = "btn-sm")
        ),
        mainPanel(
          width = 9,
          uiOutput("train_status"),
          fluidRow(column(6, tableOutput("t_perclass")),
                   column(6, tableOutput("t_confusion"))),
          hr(),
          plotOutput("p_train_lda", height = "480px"),
          hr(),
          h5("Selected panel"), tableOutput("t_panel")
        )
      )
    ),

    tabPanel(
      "ECM Aging Score",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          info_box("Methods 5.9: per-ratio age slopes are fitted by OLS on the ",
                   "continuous log-ratio, filtered at BH FDR < 0.05, then applied to ",
                   "the sign features. Higher score = older, by construction."),
          uiOutput("age_col_pick"),
          numericInput("age_fdr", "BH FDR threshold", 0.05, 0.001, 0.5, 0.005),
          actionButton("do_age", "Fit age weights", class = "btn-primary"),
          hr(),
          uiOutput("age_group_pick"),
          downloadButton("dl_betas", "Download weights (CSV)", class = "btn-sm")
        ),
        mainPanel(
          width = 9,
          uiOutput("age_status"),
          plotOutput("p_age", height = "460px"),
          hr(), tableOutput("t_betas")
        )
      )
    ),

    tabPanel("Methods & provenance", br(), column(9, uiOutput("methods_md")))
  )
)

# ==============================================================================
# server
# ==============================================================================
server <- function(input, output, session) {

  user_model <- reactiveVal(NULL)
  trained    <- reactiveVal(NULL)
  betas      <- reactiveVal(NULL)

  active_model <- reactive({
    if (identical(input$model_src, "session") && !is.null(user_model())) user_model() else MODEL
  })

  output$model_pick <- renderUI({
    ch <- c("ADSIG-180 (manuscript)" = "default")
    if (!is.null(user_model())) ch <- c(ch, "Model trained this session" = "session")
    radioButtons("model_src", NULL, ch)
  })

  expr_mat <- reactive({
    if (identical(input$src, "demo")) {
      m <- ecmx_read_expression("data/demo_counts.csv.gz")
    } else {
      req(input$expr_file)
      m <- tryCatch(ecmx_read_expression(input$expr_file$datapath),
                    error = function(e) { notify_err(e); NULL })
      req(m)
      validate(need(ncol(m) <= LIM$max_samples, sprintf(
        "This instance accepts at most %d samples; your file has %d. Run the app locally for unlimited use.",
        LIM$max_samples, ncol(m))))
    }
    if (!ecmx_looks_logged(m)) m <- ecmx_logcpm(m)
    m
  })

  meta_df <- reactive({
    if (identical(input$src, "demo"))
      return(read.csv("data/demo_metadata.csv", stringsAsFactors = FALSE))
    if (is.null(input$meta_file)) return(NULL)
    read.csv(input$meta_file$datapath, stringsAsFactors = FALSE)
  })

  projection <- reactive({
    m <- expr_mat(); req(m)
    tryCatch(ecmx_project(active_model(), m), error = function(e) { notify_err(e); NULL })
  })

  group_vec <- reactive({
    req(projection())
    ids <- colnames(expr_mat()); md <- meta_df(); col <- input$group_col
    if (is.null(md) || is.null(col) || !nzchar(col) || !(col %in% names(md)))
      return(setNames(rep("All samples", length(ids)), ids))
    g <- as.character(md[[col]][match(ids, md$SampleID)])
    g[is.na(g)] <- "unassigned"
    setNames(g, ids)
  })

  output$group_pick <- renderUI({
    md <- meta_df()
    if (is.null(md)) return(helpText("Upload metadata to group samples by condition."))
    cols <- setdiff(names(md), "SampleID")
    tagList(
      selectInput("group_col", "Group samples by", c("(none)" = "", cols),
                  selected = if ("Phenotype" %in% cols) "Phenotype" else ""),
      checkboxInput("show_hulls", "Outline groups on the phenotype space", TRUE),
      info_box("Grouping drives the ", strong("facets"), " on Posteriors, the ",
               strong("rows"), " on Centroid distances, the ", strong("colours and "),
               strong("outlines"), " here, and the Group column in the results table. ",
               "It does not change the projection itself — every sample sits where the ",
               "frozen model puts it."),
      uiOutput("group_size_msg")
    )
  })

  # State plainly how the current column partitions the samples. A group of one or
  # two samples gets no hull (a convex hull needs three points), which otherwise
  # looks like the outline setting failing to apply.
  output$group_size_msg <- renderUI({
    g <- group_vec(); req(length(g))
    n <- table(g)
    div(class = "ecmx-note",
        sprintf("%d group%s: %d outlined, %d too small to outline (needs ≥3 samples). %s",
                length(n), if (length(n) == 1) "" else "s",
                sum(n >= 3), sum(n < 3),
                paste0(names(n), " (", as.integer(n), ")", collapse = ", ")))
  })

  output$coverage_msg <- renderUI({
    p <- projection(); req(p)
    n <- length(ecmx_model_ratios(active_model()))
    div(class = if (p$coverage >= 0.85) "ecmx-note" else "ecmx-warn",
        sprintf("Ratio coverage %.1f%% (%d of %d panel ratios had both genes present). ",
                100 * p$coverage, round(p$coverage * n), n),
        if (p$coverage < 0.85)
          strong("Below the 85% floor met by all 21 validation cohorts — interpret with caution.")
        else "Ratios with a missing member are neutralised to 0, per the deployment rule in Methods 5.5.")
  })

  # Figures are built in reactives, not inside renderPlot, so the download
  # handler exports exactly the object on screen rather than rebuilding it.
  group_label <- reactive({
    if (is.null(input$group_col) || !nzchar(input$group_col)) "Your samples" else input$group_col
  })

  fig_lda <- reactive({
    p <- projection(); req(p)
    ecmx_fig_lda(TRAIN_COORDS, MODEL$train_phenotype, CENTROIDS,
                 new_coords = p$coords, new_label = unname(group_vec()),
                 group_title = group_label(),
                 show_groups = isTRUE(input$show_hulls %||% TRUE),
                 title = "Projection into the frozen ECM phenotype space")
  })
  fig_post <- reactive({
    p <- projection(); req(p)
    ecmx_fig_posteriors(p$posterior, colnames(expr_mat()), unname(group_vec()))
  })
  fig_dist <- reactive({
    p <- projection(); req(p)
    ecmx_fig_distance_heatmap(ecmx_distances(p$coords, CENTROIDS), unname(group_vec()))
  })
  fig_heat <- reactive({
    m <- active_model()
    sel <- if (is.null(input$heat_pheno) || "__all__" %in% input$heat_pheno) NULL else input$heat_pheno
    ecmx_fig_panel_heatmap(m$sign_matrix_train, m$train_phenotype,
                           model_features(m), selecting = sel)
  })

  output$heat_pheno_pick <- renderUI({
    m <- active_model(); f <- model_features(m)
    ph <- sort(unique(stats::na.omit(f$Selecting_Phenotype)))
    n  <- table(f$Selecting_Phenotype)
    lbl <- setNames(c("__all__", ph),
                    c(sprintf("All (%d ratios, thinned)", nrow(f)),
                      sprintf("%s (%d ratios)", ph, as.integer(n[ph]))))
    selectizeInput("heat_pheno", "Heatmap: ratios selecting for", choices = lbl,
                   selected = "__all__", multiple = TRUE)
  })

  output$p_lda  <- renderPlot(fig_lda())
  output$p_post <- renderPlot(fig_post())
  output$p_dist <- renderPlot(fig_dist())
  output$p_heat <- renderPlot(fig_heat())

  current_fig <- reactive({
    switch(input$viz %||% "lda",
           lda = fig_lda(), post = fig_post(), dist = fig_dist(), heat = fig_heat(),
           NULL)
  })

  output$dl_plot <- downloadHandler(
    filename = function() sprintf("ECMExplainer_%s_%s.%s", input$viz %||% "figure",
                                  Sys.Date(), input$exp_fmt),
    content = function(file) {
      p <- current_fig()
      validate(need(!is.null(p), "This tab has no figure to export."))
      # Clamp so a stray value cannot ask for a 30,000-pixel canvas.
      w <- max(3, min(30, input$exp_w %||% 9))
      h <- max(3, min(40, input$exp_h %||% 7))
      d <- max(72, min(600, input$exp_dpi %||% 300))
      ggsave(file, p, width = w, height = h, dpi = d,
             device = if (identical(input$exp_fmt, "pdf")) grDevices::cairo_pdf else "png",
             limitsize = FALSE)
    })

  results_tbl <- reactive({
    p <- projection(); req(p)
    d <- ecmx_distances(p$coords, CENTROIDS)
    data.frame(SampleID = colnames(expr_mat()), Group = unname(group_vec()),
               Predicted = as.character(p$class),
               MaxPosterior = round(apply(p$posterior, 1, max), 3),
               round(as.data.frame(p$posterior), 3),
               setNames(round(as.data.frame(d), 3), paste0("dist_", colnames(d))),
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  output$t_res <- renderTable(results_tbl(), striped = TRUE, spacing = "xs", digits = 3)
  output$dl_results <- downloadHandler(
    filename = function() paste0("ECMExplainer_results_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(results_tbl(), f, row.names = FALSE))

  # -- trajectories ------------------------------------------------------------
  output$traj_group <- renderUI({
    md <- meta_df()
    validate(need(!is.null(md), "Upload sample metadata to build trajectories."))
    selectInput("tr_group_col", "Condition column", setdiff(names(md), "SampleID"),
                selected = if ("Phenotype" %in% names(md)) "Phenotype" else NULL)
  })

  output$traj_anchors <- renderUI({
    ph <- rownames(CENTROIDS)
    tagList(
      selectInput("anchor_a", "Anchor A", ph,
                  selected = if ("Proliferating" %in% ph) "Proliferating" else ph[1]),
      selectInput("anchor_b", "Anchor B", ph,
                  selected = if ("Senescence" %in% ph) "Senescence" else ph[min(2, length(ph))]))
  })

  traj_groups <- reactive({
    md <- meta_df(); req(md, input$tr_group_col)
    as.character(md[[input$tr_group_col]][match(colnames(expr_mat()), md$SampleID)])
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

  output$p_traj <- renderPlot({
    p <- projection(); req(p, input$anchor_a, input$anchor_b)
    validate(need(input$anchor_a != input$anchor_b, "Choose two different anchors."))
    sep <- sqrt(sum((CENTROIDS[input$anchor_a, ] - CENTROIDS[input$anchor_b, ])^2))
    ecmx_fig_trajectory(ecmx_distances(p$coords, CENTROIDS), traj_groups(),
                        input$anchor_a, input$anchor_b,
                        order_levels = input$traj_seq, centroid_sep = sep)
  })

  angles_tbl <- reactive({
    p <- projection(); req(p, input$prot_from, input$prot_to)
    g <- traj_groups()
    from <- p$coords[which(g == input$prot_from), , drop = FALSE]
    to   <- p$coords[which(g == input$prot_to),   , drop = FALSE]
    validate(need(nrow(from) >= 2 && nrow(to) >= 2,
                  "The protractor needs at least 2 samples in each group (Methods 5.6)."))
    ecmx_trajectory_angles(from, to, CENTROIDS, n_boot = input$n_boot)
  })

  output$p_prot <- renderPlot(
    ecmx_fig_protractor(angles_tbl(), sprintf("%s → %s", input$prot_from, input$prot_to)))

  output$t_angles <- renderTable({
    a <- angles_tbl()
    data.frame(Phenotype = a$Phenotype, `Angle (deg)` = round(a$Angle, 1),
               `95% CI` = ifelse(is.na(a$Lower), "-", sprintf("%.0f – %.0f", a$Lower, a$Upper)),
               check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  # -- retrain -----------------------------------------------------------------
  observeEvent(input$do_train, {
    tryCatch({
      builtin <- identical(input$tr_src, "builtin")
      if (builtin) {
        if (!HAS_TRAIN) stop("The 96-sample cohort is not bundled in this build.")
        cnt <- readRDS(TRAIN_CNT); mt <- read.csv(TRAIN_MET, stringsAsFactors = FALSE)
      } else {
        req(input$tr_counts, input$tr_meta)
        cnt <- ecmx_read_expression(input$tr_counts$datapath)
        mt  <- read.csv(input$tr_meta$datapath, stringsAsFactors = FALSE)
      }
      # The bundled cohort is a known, vetted 96 samples — exempt from the upload
      # cap, so the public instance can still demonstrate the published training.
      if (!builtin && ncol(cnt) > LIM$max_train)
        stop(sprintf("This instance trains on at most %d samples; you supplied %d. Run locally for unlimited training.",
                     LIM$max_train, ncol(cnt)))
      if (!all(c("SampleID", "Phenotype") %in% names(mt)))
        stop("Metadata needs at least SampleID and Phenotype columns.")

      withProgress(message = "Training", value = 0, {
        fit <- ecmx_train(cnt, mt, GENELIST, design = input$tr_design,
                          n_per_pheno = as.integer(input$tr_npp),
                          eta_min = input$tr_eta, rank_by = input$tr_rank,
                          min_per_module = as.integer(input$tr_minmod),
                          progress = function(f, m) setProgress(value = f, detail = m))
        trained(fit); user_model(fit)
        showNotification(sprintf("Trained: %d ratios, LOOCV %.1f%%",
                                 fit$panel_size, fit$loocv_acc),
                         type = "message", duration = 8)
      })
    }, error = notify_err)
  })

  output$train_status <- renderUI({
    m <- trained()
    if (is.null(m)) return(info_box(
      "Train a model to see its cross-validated performance. The built-in cohort with ",
      strong("ADSIG / 30 per phenotype / TSP-delta"), " reproduces the published ",
      "published ADSIG-180 panel exactly: LOOCV 92.71%, macro-F1 0.921, ",
      "76 A/D + 104 SIG."))
    div(class = "ecmx-note",
        strong(sprintf("%s-%d — LOOCV %.2f%%, macro-F1 %.3f, %d PCs, %d classes, n=%d",
                       m$design, m$panel_size, m$loocv_acc, m$macro_f1, m$n_pcs,
                       length(m$phenotypes), m$n_training)),
        br(), sprintf("Ranking: %s · min per module: %d · gene-reuse cap: %d",
                      m$rank_by, m$min_per_module, m$max_appear),
        br(), "This model is now selectable on the Classify tab.")
  })

  output$t_perclass  <- renderTable({ req(trained()); trained()$per_class },
                                    striped = TRUE, spacing = "xs")
  output$t_confusion <- renderTable({ req(trained()); as.data.frame.matrix(trained()$confusion) },
                                    rownames = TRUE, striped = TRUE, spacing = "xs")
  output$t_panel <- renderTable({
    req(trained()); f <- model_features(trained())
    f[order(-f$TSP_Delta), c("Ratio_ID", "Module", "Selecting_Phenotype",
                             "Direction", "TSP_Delta", "eta_sq")]
  }, striped = TRUE, spacing = "xs", digits = 3)

  output$p_train_lda <- renderPlot({
    m <- trained(); req(m)
    pcs <- as.data.frame(m$pca$x[, seq_len(m$n_pcs), drop = FALSE])
    colnames(pcs) <- paste0("PC", seq_len(m$n_pcs))
    co <- predict(m$lda, pcs)$x
    validate(need(ncol(co) >= 2, "Need at least 3 classes to draw a 2-D LDA space."))
    ecmx_fig_lda(co, m$train_phenotype, ecmx_centroids(m),
                 title = sprintf("%s-%d training space", m$design, m$panel_size))
  })

  output$dl_model <- downloadHandler(
    filename = function() sprintf("ECMExplainer_%s_%d.rds", trained()$design, trained()$panel_size),
    content  = function(f) { req(trained()); saveRDS(trained(), f) })

  # -- aging score -------------------------------------------------------------
  output$age_col_pick <- renderUI({
    md <- meta_df()
    validate(need(!is.null(md), "Upload metadata containing an age column."))
    num <- names(md)[vapply(md, is.numeric, logical(1))]
    validate(need(length(num) > 0, "No numeric column found to use as age."))
    selectInput("age_col", "Age column", num,
                selected = grep("age", num, ignore.case = TRUE, value = TRUE)[1])
  })

  observeEvent(input$do_age, {
    tryCatch({
      md <- meta_df(); req(md, input$age_col)
      m <- expr_mat()
      ages <- suppressWarnings(as.numeric(md[[input$age_col]][match(colnames(m), md$SampleID)]))
      keep <- !is.na(ages)
      if (sum(keep) < 5) stop("Need at least 5 samples with an age value.")
      b <- ecmx_fit_age_betas(m[, keep, drop = FALSE], ages[keep],
                              ecmx_model_ratios(active_model()), fdr = input$age_fdr)
      if (is.null(b) || !nrow(b)) stop("No ratio passed the FDR threshold.")
      betas(b)
      showNotification(sprintf("%d ratios passed FDR < %.3f", nrow(b), input$age_fdr),
                       type = "message")
    }, error = notify_err)
  })

  output$age_group_pick <- renderUI({
    md <- meta_df(); if (is.null(md)) return(NULL)
    selectInput("age_group", "Group scores by", setdiff(names(md), "SampleID"))
  })

  output$age_status <- renderUI({
    b <- betas()
    if (is.null(b)) return(info_box(
      "Fit age weights to compute the score. In the manuscript these are fitted on ",
      "GSE113957 healthy dermal fibroblasts (1–94 y) and then frozen."))
    div(class = "ecmx-note",
        sprintf("%d age-associated ratios at BH FDR < %.3f.", nrow(b), input$age_fdr))
  })

  output$p_age <- renderPlot({
    b <- betas(); req(b); m <- expr_mat()
    md <- meta_df()
    g <- if (!is.null(md) && !is.null(input$age_group) && input$age_group %in% names(md))
      as.character(md[[input$age_group]][match(colnames(m), md$SampleID)]) else rep("All", ncol(m))
    ecmx_fig_aging_score(ecmx_aging_score(m, b), g)
  })

  output$t_betas <- renderTable({
    b <- betas(); req(b); b <- b[order(-abs(b$beta)), ]
    head(data.frame(Ratio = b$Ratio_ID, beta = round(b$beta, 4),
                    p = signif(b$p, 3), FDR = signif(b$fdr, 3)), 25)
  }, striped = TRUE, spacing = "xs")

  output$dl_betas <- downloadHandler(
    filename = function() paste0("ECM_age_weights_", Sys.Date(), ".csv"),
    content  = function(f) { req(betas()); write.csv(betas(), f, row.names = FALSE) })

  output$methods_md <- renderUI({
    if (!file.exists("METHODS.md")) return(p("METHODS.md not found."))
    HTML(commonmark::markdown_html(paste(readLines("METHODS.md", warn = FALSE), collapse = "\n")))
  })
}

shinyApp(ui, server)
