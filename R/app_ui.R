# ==============================================================================
# app_ui.R  –  ECMExplainER Shiny App UI (golem-style)
# ==============================================================================

app_ui <- function(request) {
  shinydashboard::dashboardPage(
    skin = "blue",

    # ---- Header ----
    shinydashboard::dashboardHeader(
      title = shiny::tags$span(
        shiny::tags$img(
          src    = "ecm_logo_placeholder.png",
          height = "30px",
          style  = "margin-right:8px; vertical-align:middle; display:none;"  # hidden until logo added
        ),
        shiny::tags$span("ECMExplainER", style = "font-size:1.25em; font-weight:700; letter-spacing:0.04em;"),
        shiny::tags$span(
          " \u2014 By Eleanor Sheekey & Roberto A. Avelar",
          style = paste0(
            "font-size:0.62em; font-weight:400; font-style:italic;",
            "opacity:0.85; vertical-align:middle; margin-left:6px;"
          )
        )
      ),
      titleWidth = 520
    ),

    # ---- Sidebar ----
    shinydashboard::dashboardSidebar(
      width = 220,
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Training Wizard", tabName = "train_wizard", icon = shiny::icon("flask")),
        shinydashboard::menuItem("Test Wizard",     tabName = "upload",       icon = shiny::icon("upload")),
        shinydashboard::menuItem("Training Plots",  tabName = "train_plots",  icon = shiny::icon("chart-bar")),
        shinydashboard::menuItem("Predictions",     tabName = "preds",        icon = shiny::icon("dna")),
        shinydashboard::menuItem("Downloads",       tabName = "dl",           icon = shiny::icon("download")),
        shinydashboard::menuItem("About",           tabName = "about",        icon = shiny::icon("info-circle"))
      )
    ),

    # ---- Body ----
    shinydashboard::dashboardBody(

      # Waiter (full-screen progress overlay)
      waiter::useWaiter(),

      # Custom CSS
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          /* ---- Layout ---- */
          .content-wrapper { background-color: #f7f9fc; }
          .box { border-top: 3px solid #3c8dbc; }

          /* ---- Header branding ---- */
          .main-header .logo {
            width: 520px !important;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            background: linear-gradient(135deg, #1a3a5c 0%, #2d6a9f 100%) !important;
          }
          .main-header .navbar {
            margin-left: 520px !important;
          }
          .main-header .logo .ecm-title {
            font-size: 1.22em;
            font-weight: 700;
            letter-spacing: 0.05em;
            color: #ffffff;
          }
          .main-header .logo .ecm-authors {
            font-size: 0.63em;
            font-weight: 400;
            font-style: italic;
            color: rgba(255,255,255,0.82);
            margin-left: 5px;
          }

          /* ---- Sidebar to match header width ---- */
          .main-sidebar { padding-top: 50px; }

          /* ---- Misc ---- */
          .ratio-badge { font-size: 1.1em; font-weight: bold; color: #2d6a9f; }
          .info-text   { font-size: 0.92em; color: #555; }

          /* ---- Wizard stepper ---- */
          .wizard-steps {
            display: flex; align-items: center; margin-bottom: 22px;
            padding: 14px 20px; background: #fff;
            border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.10);
          }
          .wizard-step {
            display: flex; flex-direction: column; align-items: center;
            flex: 1; position: relative;
          }
          .wizard-step:not(:last-child)::after {
            content: ''; position: absolute;
            top: 18px; left: 60%; width: 80%; height: 3px;
            background: #d0dce8;
          }
          .wizard-step.done::after   { background: #2d6a9f; }
          .wizard-step.active::after { background: #d0dce8; }
          .step-circle {
            width: 36px; height: 36px; border-radius: 50%;
            display: flex; align-items: center; justify-content: center;
            font-size: 1em; font-weight: 700; border: 3px solid #d0dce8;
            background: #fff; color: #aaa; z-index: 1;
          }
          .wizard-step.done .step-circle {
            background: #2d6a9f; border-color: #2d6a9f; color: #fff;
          }
          .wizard-step.active .step-circle {
            background: #fff; border-color: #2d6a9f; color: #2d6a9f;
          }
          .step-label {
            font-size: 0.72em; margin-top: 5px; color: #888; text-align: center;
          }
          .wizard-step.active .step-label { color: #2d6a9f; font-weight: 600; }
          .wizard-step.done  .step-label  { color: #2d6a9f; }

          /* ---- Filter chip ---- */
          .filter-panel { background:#eef4fb; border-radius:6px; padding:12px 16px; margin-bottom:10px; }

          /* ---- Clickable wizard steps ---- */
          .wizard-step.done .step-circle { cursor: pointer; }
          .wizard-step.done .step-circle:hover { opacity: 0.75; }
          .wizard-step.done .step-label  { cursor: pointer; }
        "))
      ),

      shinydashboard::tabItems(

        # ================================================================
        # TAB 0: Training Wizard
        # ================================================================
        shinydashboard::tabItem(
          tabName = "train_wizard",

          shiny::uiOutput("train_wizard_steps_ui"),

          # --------------------------------------------------
          # TRAIN STEP 1 — Upload training data
          # --------------------------------------------------
          shiny::conditionalPanel(
            condition = "output.train_wizard_step == 1",

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("upload"), " Step 1 – Upload Training Data"),
                width       = 6, status = "primary", solidHeader = TRUE,

                shiny::p("Upload raw count matrix (genes × samples) for training.",
                         " First column = gene symbols; column headers = sample IDs.",
                         class = "info-text"),
                shiny::fileInput("train_counts_file", "Raw Count Matrix (.csv)",
                                 accept = ".csv", buttonLabel = "Browse\u2026"),
                shiny::radioButtons(
                  "train_assay_type", "Assay Type",
                  choices  = c("RNA-seq" = "RNA-seq"),
                  selected = "RNA-seq", inline = TRUE
                ),
                shiny::hr(),
                shiny::uiOutput("train_norm_options_ui")
              ),

              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("table"), " Step 1 – Training Metadata"),
                width       = 6, status = "primary", solidHeader = TRUE,

                shiny::p(shiny::tags$b("Required columns:"),
                         " Sample_ID, Phenotype (the class to predict).",
                         " Optional: Study_ID, AssayType (for mixed data).",
                         class = "info-text"),
                shiny::fileInput("train_meta_file", "Sample Metadata (.csv)",
                                 accept = ".csv", buttonLabel = "Browse\u2026"),
                shiny::uiOutput("train_batch_col_ui"),
                shiny::uiOutput("train_phenotype_col_ui")
              )
            ),

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("dna"), " Step 1 – ECM Dictionary"),
                width       = 6, status = "info", solidHeader = TRUE,
                collapsible = TRUE, collapsed = FALSE,

                shiny::p("CSV with at least columns",
                         shiny::code("Gene"), "and", shiny::code("ECMExplained_Function"),
                         "(values: ASSEMBLY, DESTRUCTION, SIGNALLING).",
                         class = "info-text"),
                shiny::radioButtons(
                  "train_dict_source", "ECM Dictionary",
                  choices  = c("Use built-in dictionary (default)" = "builtin",
                               "Upload custom dictionary"           = "upload"),
                  selected = "builtin"
                ),
                shiny::conditionalPanel(
                  condition = "input.train_dict_source == 'upload'",
                  shiny::fileInput("train_dict_file", "ECM Dictionary (.csv)",
                                   accept = ".csv", buttonLabel = "Browse\u2026")
                )
              ),

              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("sliders-h"), " Step 1 – Feature Selection Options"),
                width       = 6, status = "info", solidHeader = TRUE,
                collapsible = TRUE, collapsed = TRUE,

                shiny::selectInput(
                  "train_species",
                  shiny::tags$span("Organism / species",
                    shiny::tags$small(" — recorded in model bundle for reference", style = "color:#888;")),
                  choices  = c("Human (Homo sapiens)"        = "Human",
                               "Mouse (Mus musculus)"        = "Mouse",
                               "Rat (Rattus norvegicus)"     = "Rat",
                               "Other / custom"              = "Other"),
                  selected = "Human"
                ),
                shiny::conditionalPanel(
                  condition = "input.train_species == 'Other'",
                  shiny::textInput("train_species_custom", "Specify organism",
                                   placeholder = "e.g. Zebrafish")
                ),
                shiny::hr(style = "margin:8px 0;"),
                shiny::numericInput("train_max_ratios", "Max ratios per panel (AD / SIG)",
                                    value = 50L, min = 5L, max = 200L, step = 5L),
                shiny::numericInput("train_candidate_pool",
                                    shiny::tags$span("Candidate pool per panel",
                                      shiny::tags$small(" — screened before 50/50 split",
                                                        style = "color:#888;")),
                                    value = 100L, min = 50L, max = 500L, step = 25L),
                shiny::numericInput("train_eta_sq_min", "Minimum η² for ratio",
                                    value = 0.30, min = 0.05, max = 0.95, step = 0.05),
                shiny::numericInput("train_synergy_delta", "Minimum synergy Δη²",
                                    value = 0.05, min = 0.00, max = 0.50, step = 0.01),

                shiny::hr(style = "margin:10px 0;"),
                shiny::tags$label("Ratio selection method", class = "control-label"),
                shinyWidgets::radioGroupButtons(
                  inputId  = "train_selection_method",
                  label    = NULL,
                  choices  = c(
                    "V19 \u03b7\u00b2 greedy"       = "eta_sq",
                    "V3 TSP-\u0394 stratified (k-TSP)" = "tsp_delta"
                  ),
                  selected  = "eta_sq",
                  justified = TRUE,
                  size      = "sm"
                ),
                shiny::conditionalPanel(
                  condition = "input.train_selection_method == 'tsp_delta'",
                  shiny::br(),
                  shiny::numericInput(
                    "train_n_per_pheno",
                    shiny::tags$span("Ratios per phenotype per module",
                      shiny::tags$small(" — V3 default: 12 \u00d7 6 = 72/module",
                                        style = "color:#888;")),
                    value = 12L, min = 2L, max = 50L, step = 1L
                  )
                ),
                shiny::hr(style = "margin:10px 0;"),
                shiny::tags$label("Add a custom phenotype label", class = "control-label"),
                shiny::tags$small(
                  shiny::tags$em(" — type a new label name; samples with this exact",
                                 " value in the Phenotype column will be included"),
                  style = "display:block; color:#777; margin-bottom:6px;"
                ),
                shiny::fluidRow(
                  shiny::column(8,
                    shiny::textInput("train_new_phenotype_name",
                                     label       = NULL,
                                     placeholder = "e.g. Apoptosis")
                  ),
                  shiny::column(4,
                    shiny::br(),
                    shinyWidgets::actionBttn(
                      "train_add_phenotype_btn", "Add",
                      style = "bordered", color = "primary", size = "sm"
                    )
                  )
                ),
                shiny::uiOutput("train_custom_phenotypes_ui"),
                shiny::tags$label("Max gene appearances per panel",
                                  class = "control-label"),
                shiny::tags$small(
                  shiny::tags$em(" — 1 = maximum feature diversity (V19 default); higher = more lenient"),
                  style = "display:block; color:#777; margin-bottom:4px;"
                ),
                shiny::fluidRow(
                  shiny::column(8,
                    shiny::sliderInput("train_max_appearances",
                                       label  = NULL,
                                       min    = 1L, max = 5L, value = 1L, step = 1L,
                                       ticks  = TRUE)
                  ),
                  shiny::column(4,
                    shiny::checkboxInput("train_max_appearances_unlimited",
                                         "No limit", value = FALSE)
                  )
                ),

                shiny::hr(style = "margin:10px 0;"),
                shiny::tags$label("Restore settings from a previous run", class = "control-label"),
                shiny::tags$small(
                  shiny::tags$em(" — upload a params .json or model .rds to restore all settings"),
                  style = "display:block; color:#777; margin-bottom:6px;"
                ),
                shiny::fileInput(
                  "train_params_upload",
                  label    = NULL,
                  accept   = c(".json", ".rds"),
                  buttonLabel = "Upload log\u2026",
                  placeholder = "params .json or model .rds"
                ),
                shiny::uiOutput("train_params_restore_ui")
              )
            ),

            shiny::fluidRow(
              shiny::column(12,
                shinyWidgets::actionBttn(
                  "train_step1_next", "Next: Review & Run  \u2192",
                  style = "gradient", color = "primary", size = "md"
                )
              )
            )
          ),

          # --------------------------------------------------
          # TRAIN STEP 2 — Review & Run
          # --------------------------------------------------
          shiny::conditionalPanel(
            condition = "output.train_wizard_step == 2",

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("play"), " Step 2 – Run Training"),
                width       = 12, status = "success", solidHeader = TRUE,

                shiny::uiOutput("train_data_summary_ui"),
                shiny::br(),

                shiny::fluidRow(
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      "train_step2_back", "\u2190 Back",
                      style = "bordered", color = "default", size = "md"
                    )
                  ),
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      "run_train_btn", " Run Training",
                      icon  = shiny::icon("flask"),
                      style = "gradient", color = "warning", size = "lg"
                    )
                  ),
                  shiny::column(6,
                    shiny::uiOutput("train_run_status_ui")
                  )
                ),
                shiny::br(),
                shinyWidgets::progressBar(
                  id = "train_pb", value = 0, total = 100,
                  title = "", display_pct = TRUE, status = "warning"
                )
              )
            )
          ),

          # --------------------------------------------------
          # TRAIN STEP 3 — Results & Download
          # --------------------------------------------------
          shiny::conditionalPanel(
            condition = "output.train_wizard_step == 3",

            shiny::fluidRow(
              shinydashboard::box(
                title = shiny::tags$span(shiny::icon("check-circle"), " Training Complete"),
                width = 12, status = "success", solidHeader = TRUE,

                shiny::uiOutput("train_results_summary_ui"),

                # LOOCV accuracy callout
                shiny::uiOutput("train_loocv_ui"),

                shiny::hr(),
                shiny::h5("Selected Feature Ratios"),
                DT::DTOutput("train_features_table"),

                shiny::hr(),

                # ---- Save model panel ----
                shiny::div(
                  style = "background:#f0f7ff; border:1px solid #b8d4f0; border-radius:6px; padding:16px 18px; margin-bottom:14px;",
                  shiny::h5(shiny::icon("save"), " Save Trained Model", style = "margin-top:0;"),
                  shiny::fluidRow(
                    shiny::column(5,
                      shiny::textInput(
                        "train_model_name",
                        shiny::tags$span("Model name", shiny::tags$small(" — embedded in the bundle and used as the filename", style = "color:#888; font-weight:400;")),
                        value       = "",
                        placeholder = "e.g. Human_Fibroblast_v1"
                      )
                    ),
                    shiny::column(3,
                      shiny::br(),
                      shiny::downloadButton("dl_train_model_bundle",
                                            "Download Model (.rds)",
                                            class = "btn btn-warning",
                                            style = "margin-top:6px;")
                    ),
                    shiny::column(2,
                      shiny::br(),
                      shiny::downloadButton("dl_train_params_json",
                                            "Params (.json)",
                                            class = "btn btn-default",
                                            style = "margin-top:6px;")
                    ),
                    shiny::column(2,
                      shiny::br(),
                      shiny::downloadButton("dl_train_features_csv",
                                            "Ratios (.csv)",
                                            class = "btn btn-default",
                                            style = "margin-top:6px;")
                    )
                  ),
                  shiny::hr(style = "margin:10px 0;"),
                  shiny::fluidRow(
                    shiny::column(12,
                      shiny::actionButton(
                        "train_make_default_btn",
                        shiny::tagList(shiny::icon("star"), " Make this a default model in the app"),
                        class = "btn btn-success",
                        style = "margin-bottom:6px;"
                      ),
                      shiny::uiOutput("train_make_default_ui")
                    )
                  ),
                  shiny::p(
                    shiny::icon("info-circle"),
                    " The model .rds and the params .json both record every training setting.",
                    " Upload either file in the Test Wizard to automatically restore all parameters.",
                    " Click",
                    shiny::strong(" Make this a default model"),
                    " to save it to the app's training_data/ folder so it appears in the",
                    " Test Wizard default model dropdown for all future sessions.",
                    class = "info-text", style = "margin-bottom:0;"
                  )
                ),

                shiny::fluidRow(
                  shiny::column(4,
                    shinyWidgets::actionBttn(
                      "train_step3_back", "\u2190 Back to Upload Data",
                      style = "bordered", color = "default", size = "md"
                    )
                  )
                ),

                shiny::hr(),
                shiny::fluidRow(
                  shinydashboard::box(
                    title = "Training LDA Space", width = 6,
                    status = "info", solidHeader = FALSE,
                    shiny::plotOutput("plot_lda_train_wizard", height = "400px")
                  ),
                  shinydashboard::box(
                    title = "Training Global PCA", width = 6,
                    status = "info", solidHeader = FALSE,
                    shiny::plotOutput("plot_pca_train_wizard2", height = "400px")
                  )
                )
              )
            )
          )
        ),

        # ================================================================
        # TAB 1: Test Wizard — Steps 1-4
        # ================================================================
        shinydashboard::tabItem(
          tabName = "upload",

          # ---- Step indicator ----
          shiny::uiOutput("wizard_steps_ui"),

          # ==================================================
          # STEP 1 — Configure training source, conditions, upload data
          # ==================================================
          shiny::conditionalPanel(
            condition = "output.wizard_step == 1",

            # ── Training source + projection mode ──────────────────────────
            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("database"), " Training Model Source"),
                width       = 6, status = "primary", solidHeader = TRUE,

                shiny::radioButtons(
                  "test_model_source",
                  label   = NULL,
                  choices = c(
                    "Use a default model"                        = "default",
                    "Upload my own model bundle (.rds)"          = "upload",
                    "Use model trained in this session"          = "session"
                  ),
                  selected = "default"
                ),

                # ── Default model selector ──────────────────────────────
                shiny::conditionalPanel(
                  condition = "input.test_model_source == 'default'",
                  shiny::selectInput(
                    "test_default_model",
                    label   = shiny::tags$span(shiny::icon("star"), " Select default model"),
                    choices = NULL   # populated server-side from ECM_DEFAULT_MODELS
                  )
                ),

                # ── Upload model ────────────────────────────────────────
                shiny::conditionalPanel(
                  condition = "input.test_model_source == 'upload'",
                  shiny::fileInput("test_model_file", "Model Bundle (.rds)",
                                   accept = ".rds", buttonLabel = "Browse\u2026")
                ),

                # ── Session-trained notice ──────────────────────────────
                shiny::conditionalPanel(
                  condition = "input.test_model_source == 'session'",
                  shiny::tags$p(
                    shiny::icon("info-circle"),
                    " Uses the model you just trained in the Training Wizard.",
                    " Run the Training Wizard first if no model is shown below.",
                    class = "info-text"
                  )
                ),

                shiny::hr(style = "margin:8px 0;"),
                # Unified info panel — shows for all 3 sources
                shiny::uiOutput("active_model_info_ui")
              ),

              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("project-diagram"), " Projection Mode"),
                width       = 6, status = "primary", solidHeader = TRUE,

                shiny::radioButtons(
                  "test_projection_mode",
                  label   = NULL,
                  choices = c(
                    "Mode A \u2014 Unanchored (no reference group, default)" = "unanchored",
                    "Mode B \u2014 Anchored (align using reference group)"   = "anchored"
                  ),
                  selected = "unanchored"
                ),

                shiny::conditionalPanel(
                  condition = "input.test_projection_mode == 'anchored'",
                  shiny::hr(),
                  shiny::p(shiny::icon("info-circle"),
                           " Mode B aligns your samples to the training atlas using one or more",
                           " known reference groups (e.g. Proliferating controls → atlas Proliferating,",
                           " Senescence controls → atlas Senescence).",
                           " For each anchor pair, the centroid of your reference samples is shifted",
                           " onto the corresponding atlas centroid.",
                           " The final shift is the sample-count–weighted average across all pairs,",
                           " and is then applied to all your samples before classification.",
                           " Conditions not listed as anchor pairs are projected without explicit alignment.",
                           class = "info-text"),
                  shiny::br(),
                  shiny::uiOutput("test_anchor_col_ui"),
                  shiny::uiOutput("anchor_pairs_ui")
                )
              )
            ),

            # ── Dynamic condition selector (populated from active training) ─
            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("filter"), " Select Training Conditions"),
                width       = 8,
                status      = "primary",
                solidHeader = TRUE,
                shiny::p(
                  "Choose which phenotype classes to include in the LDA model.",
                  " The list updates automatically when you switch training source.",
                  shiny::tags$b(" Minimum two conditions required."),
                  class = "info-text"
                ),
                shiny::br(),
                shiny::uiOutput("condition_selector_ui")
              ),
              shinydashboard::box(
                title       = "Condition Summary",
                width       = 4,
                status      = "info",
                solidHeader = FALSE,
                shiny::uiOutput("condition_summary_ui")
              )
            ),

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("upload"), " Upload Count Matrix"),
                width       = 6,
                status      = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,

                shiny::p(shiny::tags$b("Required:"), " CSV with genes as rows, samples as columns.",
                         " First column = gene symbols. Column headers = sample IDs.",
                         class = "info-text"),
                shiny::fileInput("counts_file", "Raw Count Matrix (.csv)",
                                 accept = ".csv", buttonLabel = "Browse\u2026"),
                shiny::radioButtons(
                  "assay_type_choice", "Assay Type",
                  choices  = c("RNA-seq" = "RNA-seq"),
                  selected = "RNA-seq", inline = TRUE
                ),
                shiny::hr(),
                shiny::selectInput(
                  inputId  = "dup_gene_strategy",
                  label    = shiny::tags$span(
                    "Duplicate gene handling",
                    shiny::tags$small(" (ratio genes only)", style = "color:#888; font-weight:400;")
                  ),
                  choices  = c(
                    "Remove all duplicate genes"       = "remove",
                    "Average counts across duplicates" = "mean",
                    "Keep highest-expressed duplicate" = "max"
                  ),
                  selected = "remove"
                ),
                shiny::hr(),
                shiny::div(
                  style = "background:#fff8e1; border-left:4px solid #f39c12; padding:8px 10px; border-radius:3px; margin-top:6px;",
                  shiny::checkboxInput(
                    inputId = "counts_already_cpm",
                    label   = shiny::tags$span(
                      shiny::icon("exclamation-circle", style = "color:#e67e22;"),
                      " My data is already CPM-normalised (skip library-size normalisation)",
                      style = "font-size:0.9em;"
                    ),
                    value = FALSE
                  ),
                  shiny::uiOutput("cpm_skip_warning_ui"),
                  shiny::hr(style = "margin:6px 0;"),
                  shiny::uiOutput("scale_counts_ui")
                )
              ),

              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("table"), " Upload Sample Metadata"),
                width       = 6,
                status      = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,

                shiny::p(shiny::tags$b("Required column:"), " Sample_ID (must match count matrix columns).",
                         " Optional: Study_ID, AssayType, and any other columns to regress out.",
                         class = "info-text"),
                shiny::fileInput("meta_file", "Sample Metadata (.csv)",
                                 accept = ".csv", buttonLabel = "Browse\u2026"),
                shiny::uiOutput("test_group_col_ui"),
                shiny::uiOutput("batch_col_ui"),
                shiny::uiOutput("protect_cols_ui")
              )
            ),

            # ---- Sample selector (appears after metadata is uploaded) ----
            shiny::uiOutput("sample_selector_box_ui"),

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("play"), " Run ECMExplainER"),
                width       = 12,
                status      = "success",
                solidHeader = TRUE,

                shiny::p(
                  "Training on: ", shiny::uiOutput("selected_conditions_chips", inline = TRUE),
                  class = "info-text"
                ),

                # ---- Prediction feature mode ----
                shiny::fluidRow(
                  shiny::column(12,
                    shiny::div(
                      style = "background:#f8f9fa; border-left:4px solid #3c8dbc; padding:10px 14px; margin-bottom:12px; border-radius:3px;",
                      shiny::tags$b("Prediction features:"),
                      shiny::radioButtons(
                        inputId  = "prediction_mode",
                        label    = NULL,
                        choices  = c(
                          "ECM ratios (log\u2082A \u2212 log\u2082B) \u2014 recommended" = "ratios",
                          "Binary TSP (sign A > B, 0/1) \u2014 batch-robust"           = "binary_tsp",
                          "All ECM genes (individual log\u2082CPM)"                    = "ecm_genes",
                          "All genes (full transcriptome)"                             = "all_genes"
                        ),
                        selected = "ratios",
                        inline   = FALSE
                      )
                    )
                  )
                ),

                shiny::fluidRow(
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      inputId = "step2_back",
                      label   = "\u2190 Back",
                      style   = "bordered",
                      color   = "default",
                      size    = "md"
                    )
                  ),
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      inputId = "run_btn",
                      label   = " Run ECMExplainER",
                      icon    = shiny::icon("play"),
                      style   = "gradient",
                      color   = "primary",
                      size    = "lg"
                    )
                  ),
                  shiny::column(6,
                    shiny::uiOutput("run_status_ui")
                  )
                ),
                shiny::br(),
                shinyWidgets::progressBar(
                  id = "pb", value = 0, total = 100, title = "",
                  display_pct = TRUE, status = "info"
                ),
                shiny::br(),
                shiny::uiOutput("ratio_info_ui")
              )
            )
          ),

          # ==================================================
          # STEP 2 — Training PCA
          # ==================================================
          shiny::conditionalPanel(
            condition = "output.wizard_step == 2",

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("chart-bar"), " Step 2 – Training PCA Space"),
                width       = 12,
                status      = "info",
                solidHeader = TRUE,

                shiny::p(
                  "Your training set (filtered to selected conditions) projected into PCA space.",
                  " Use the controls below to colour by any metadata column.",
                  class = "info-text"
                ),

                shiny::fluidRow(
                  shiny::column(4,
                    shiny::selectInput(
                      inputId  = "pca_train_colour_col",
                      label    = "Colour by",
                      choices  = c("Phenotype" = "Phenotype"),
                      selected = "Phenotype"
                    )
                  ),
                  shiny::column(4,
                    shiny::selectInput("step3_pc_x", "X axis", choices = "PC1", selected = "PC1")
                  ),
                  shiny::column(4,
                    shiny::selectInput("step3_pc_y", "Y axis", choices = "PC2", selected = "PC2")
                  )
                ),

                shiny::plotOutput("plot_pca_train_wizard", height = "500px"),

                shiny::br(),
                shiny::fluidRow(
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      "step3_back", "\u2190 Back to Configure",
                      style = "bordered", color = "default", size = "md"
                    )
                  ),
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      "step3_next", "Next: Project Test Data  \u2192",
                      style = "gradient", color = "primary", size = "md"
                    )
                  )
                )
              )
            )
          ),

          # ==================================================
          # STEP 3 — Test projection
          # ==================================================
          shiny::conditionalPanel(
            condition = "output.wizard_step == 3",

            shiny::fluidRow(
              shinydashboard::box(
                title       = shiny::tags$span(shiny::icon("dna"), " Step 3 – Test Sample Projection"),
                width       = 12,
                status      = "warning",
                solidHeader = TRUE,

                shiny::p(
                  "Your samples (triangles) projected into the trained LDA space.",
                  " Filter which samples are plotted using the controls below.",
                  class = "info-text"
                ),

                shiny::div(class = "filter-panel",
                  shiny::fluidRow(
                    shiny::column(4,
                      shiny::uiOutput("step4_filter_col_ui")
                    ),
                    shiny::column(4,
                      shiny::uiOutput("step4_filter_val_ui")
                    ),
                    shiny::column(4,
                      shiny::br(),
                      shiny::checkboxInput(
                        "step4_show_all_train",
                        "Show all training samples",
                        value = TRUE
                      )
                    )
                  )
                ),

                shiny::plotOutput("plot_lda_combined_wizard", height = "520px"),
                shiny::div(
                  style = "text-align:right; margin-top:4px;",
                  shiny::downloadButton("dl_lda_combined_wizard", "Download Plot",
                                        class = "btn-sm btn-default")
                ),

                shiny::br(),
                shiny::fluidRow(
                  shiny::column(3,
                    shinyWidgets::actionBttn(
                      "step4_back", "\u2190 Back to Training PCA",
                      style = "bordered", color = "default", size = "md"
                    )
                  ),
                  shiny::column(6,
                    shiny::uiOutput("run_status_ui2")
                  )
                )
              )
            ),

            shiny::fluidRow(
              shinydashboard::box(
                title = "Diagnostic Report",
                width = 12, status = "primary", solidHeader = TRUE,
                shiny::div(
                  style = "text-align:right; margin-bottom:6px;",
                  shiny::downloadButton("dl_report_wizard", "Download CSV",
                                        class = "btn-sm btn-default")
                ),
                shiny::p("Select rows to highlight those samples on the LDA plot above.",
                         class = "info-text"),
                DT::DTOutput("report_table_wizard")
              )
            )
          )
        ),

        # ================================================================
        # TAB 2: Training Plots
        # ================================================================
        shinydashboard::tabItem(
          tabName = "train_plots",

          shiny::fluidRow(
            shinydashboard::box(
              title = "Training Global PCA",
              width = 6, status = "info", solidHeader = TRUE,
              shiny::plotOutput("plot_pca_train", height = "480px")
            ),
            shinydashboard::box(
              title = "Training LDA",
              width = 6, status = "info", solidHeader = TRUE,
              shiny::plotOutput("plot_lda_train", height = "480px")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Training PCA – colour by metadata column",
              width = 12, status = "info", solidHeader = TRUE,
              collapsible = TRUE,

              shiny::fluidRow(
                shiny::column(4,
                  shiny::selectInput(
                    inputId  = "pca_train_colour_col",
                    label    = "Colour by",
                    choices  = c("(none – hide plot)" = ""),
                    selected = ""
                  )
                )
              ),
              shiny::uiOutput("plot_pca_train_meta_ui")
            )
          )
        ),

        # ================================================================
        # TAB 3: Predictions
        # ================================================================
        shinydashboard::tabItem(
          tabName = "preds",

          shiny::fluidRow(
            shinydashboard::box(
              title = "Phenotype Probability Barplot",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsible = TRUE,

              # ── Controls row ──────────────────────────────────────────
              shiny::fluidRow(
                shiny::column(4,
                  shiny::uiOutput("probs_filter_col_ui")
                ),
                shiny::column(4,
                  shiny::uiOutput("probs_filter_val_ui")
                ),
                shiny::column(4,
                  shiny::uiOutput("probs_facet_cols_ui")
                )
              ),
              shiny::hr(style = "margin:6px 0 10px 0;"),

              # ── Plot (height grows with facets) ───────────────────────
              shiny::uiOutput("plot_probs_ui"),
              shiny::div(
                style = "text-align:right; margin-top:4px;",
                shiny::downloadButton("dl_probs_plot_preds", "Download Plot (.png)",
                                      class = "btn-sm btn-default")
              )
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "LDA: Training + Test Overlay",
              width = 12, status = "warning", solidHeader = TRUE,
              shiny::plotOutput("plot_lda_combined", height = "480px")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = shiny::tags$span(shiny::icon("th"), " Sign Matrix Heatmap"),
              width = 12, status = "info", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              shiny::tags$small(
                "Sign values {−1, 0, +1} = blue / cream / red. ",
                "Rows = ECM ratio pairs, columns = samples ordered by phenotype.",
                style = "color:#555;"
              ),
              shiny::br(),
              shiny::fluidRow(
                shiny::column(4,
                  shiny::selectInput(
                    "sign_heatmap_source",
                    "Samples to show",
                    choices  = c("Test samples only"         = "test",
                                 "Training samples only"     = "train",
                                 "Both (test + training)"    = "both"),
                    selected = "test"
                  )
                ),
                shiny::column(4,
                  shiny::checkboxInput("sign_heatmap_cluster_rows",
                                       "Cluster ratio rows", value = TRUE)
                ),
                shiny::column(4,
                  shiny::br(),
                  shiny::downloadButton("dl_sign_heatmap", "Download (.png)",
                                        class = "btn-sm btn-default")
                )
              ),
              shiny::plotOutput("plot_sign_heatmap", height = "500px")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = shiny::tags$span(shiny::icon("map"), " Distance to Training Centroids"),
              width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              shiny::tags$small(
                "Euclidean distance of each test sample to each phenotype centroid in LDA space. ",
                "Red = close, blue = far. Predictions correspond to the closest centroid.",
                style = "color:#555;"
              ),
              shiny::br(),
              shiny::fluidRow(
                shiny::column(4,
                  shiny::checkboxInput("dist_heatmap_aggregate",
                                       "Aggregate by predicted phenotype (median)", value = FALSE)
                ),
                shiny::column(4,
                  shiny::br(),
                  shiny::downloadButton("dl_dist_heatmap", "Download (.png)",
                                        class = "btn-sm btn-default")
                )
              ),
              shiny::plotOutput("plot_dist_heatmap", height = "420px")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Ratio PCA – user samples",
              width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE,

              shiny::fluidRow(
                shiny::column(4,
                  shiny::uiOutput("ratio_pca_colour_ui")
                ),
                shiny::column(4,
                  shiny::checkboxInput(
                    inputId = "ratio_pca_show_training",
                    label   = "Overlay training samples",
                    value   = FALSE
                  )
                )
              ),
              shiny::uiOutput("plot_ratio_pca_ui")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Diagnostic Report Table",
              width = 12, status = "primary", solidHeader = TRUE,
              DT::DTOutput("report_table")
            )
          )
        ),

        # ================================================================
        # TAB 4: Downloads
        # ================================================================
        shinydashboard::tabItem(
          tabName = "dl",

          shiny::fluidRow(
            shinydashboard::box(
              title = "Download Results",
              width = 6, status = "success", solidHeader = TRUE,

              shiny::p("All downloads become available after a successful run.", class = "info-text"),
              shiny::br(),

              shiny::downloadButton("dl_report",         "Diagnostic Report (.csv)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_pca_train",      "Training PCA Plot (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_lda_train",      "Training LDA Plot (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_pca_train_meta", "Training PCA (metadata colour) (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_probs_plot",     "Probability Barplot (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_lda_combined",   "Combined LDA Plot (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_sign_heatmap_dl","Sign Matrix Heatmap (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_dist_heatmap_dl","Distance to Centroids Heatmap (.png)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_ratio_pca",      "Ratio PCA – User Samples (.png)")
            ),

            shinydashboard::box(
              title = "Download Intermediate Data",
              width = 6, status = "success", solidHeader = TRUE,

              shiny::downloadButton("dl_norm_counts",    "Normalised Counts (.csv)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_ratio_mat",      "ECM Ratio Matrix (.csv)"),
              shiny::br(), shiny::br(),
              shiny::downloadButton("dl_avail_features", "Available Features (.csv)")
            )
          )
        ),

        # ================================================================
        # TAB 5: About
        # ================================================================
        shinydashboard::tabItem(
          tabName = "about",

          shiny::fluidRow(
            shinydashboard::box(
              title = "About ECMExplainER", width = 12,
              status = "primary", solidHeader = FALSE,

              shiny::h4("Overview"),
              shiny::p("ECMExplainER is a machine-learning diagnostic tool that predicts",
                       "fibroblast phenotypic states (Proliferating, Senescence, Quiescence,",
                       "Fibrosis, Inflammation, OxidativeStress) from bulk",
                       "RNA-seq transcriptomics using ECM gene expression ratios.",
                       "Phenotypic classes are derived from the loaded training data."),
              shiny::h4("Pipeline Summary"),
              shiny::tags$ol(
                shiny::tags$li("Upload raw count matrix (genes × samples) and sample metadata."),
                shiny::tags$li("Counts are normalised to log\u2082(CPM + 1)."),
                shiny::tags$li("Batch effects (e.g., study) are regressed out via ComBat."),
                shiny::tags$li("Available ECM ratios are detected from the uploaded gene set."),
                shiny::tags$li("The classifier is retrained on available ratios (dual PCA + LDA)."),
                shiny::tags$li("New samples are projected and phenotype probabilities are reported.")
              ),
              shiny::h4("Training Data"),
              shiny::p("The default training set comprises 82 human fibroblast samples from 6",
                       "published studies, spanning 7 phenotypic states. Training data can be",
                       "replaced by updating files in the app's", shiny::code("training_data/"), "directory."),
              shiny::h4("Notes"),
              shiny::tags$ul(
                shiny::tags$li("Gene symbols must be HGNC human gene symbols."),
                shiny::tags$li("The metadata file must include a Sample_ID column matching the count matrix column headers."),
                shiny::tags$li("Only ratios whose constituent genes are present in your data will be used.")
              )
            )
          )
        )
      )
    )
  )
}
