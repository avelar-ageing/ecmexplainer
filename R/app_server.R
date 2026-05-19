# ==============================================================================
# app_server.R  –  ECMExplainER Shiny App Server (golem-style)
# ==============================================================================

app_server <- function(input, output, session) {

  # ── Full reactive traceback to console whenever any reactive errors ──────────
  options(shiny.error = function() {
    cat("\n========== REACTIVE ERROR ==========\n")
    traceback(4)
    cat("=====================================\n\n")
  })
  options(shiny.fullstacktrace = TRUE)

  # ---- Private palette helper -------------------------------------------------
  # Delegates to .pal_for_levels() defined in mod_plots.R (sourced first).
  # Kept as a thin wrapper here so call-sites in this file stay readable.
  .cond_pal <- function(levels_vec) .pal_for_levels(levels_vec)

  # ---- Resolve training data directory ----------------------------------------
  # Works both when run via app.R at shiny_app/ root and from golem::run_app()
  training_data_dir <- {
    candidates <- c(
      file.path(getwd(), "training_data"),
      system.file("training_data", package = "ECMExplainERApp"),
      file.path(dirname(getwd()), "training_data")
    )
    found <- candidates[dir.exists(candidates)]
    if (length(found) == 0)
      stop("Cannot locate training_data/ directory. Ensure it is present next to app.R.")
    found[1]
  }

  # ---- Load training data once at startup ----
  training <- tryCatch(
    load_training_data(training_data_dir),
    error = function(e) {
      shiny::showNotification(paste("Failed to load training data:", e$message),
                              type = "error", duration = NULL)
      NULL
    }
  )

  # ---- Reactive values to hold pipeline results ----
  rv <- shiny::reactiveValues(
    norm_counts       = NULL,
    user_meta         = NULL,
    ratio_res         = NULL,
    train_obj         = NULL,
    train_meta_filt   = NULL,
    pred_obj          = NULL,
    run_success       = FALSE,
    wizard_step       = 1L,        # 1=conditions, 2=upload, 3=train PCA, 4=projection
    # Anchored-mode extra results
    pred_obj_aligned  = NULL       # Mode B aligned prediction
  )

  # ---- Reactive values for the training wizard ----
  rv_train <- shiny::reactiveValues(
    train_obj   = NULL,   # trained model object
    features    = NULL,   # selected ratio features data.frame
    train_meta  = NULL,   # metadata used for training
    train_expr  = NULL,   # log2CPM matrix used for training
    ecm_dict    = NULL,   # ecm dictionary used
    step        = 1L      # 1=upload, 2=review+run, 3=results
  )

  # ================================================================
  # Wizard step outputs (read by conditionalPanel in UI)
  # ================================================================
  output$wizard_step <- shiny::reactive({ rv$wizard_step })
  shiny::outputOptions(output, "wizard_step", suspendWhenHidden = FALSE)

  # Step indicator HTML — completed steps are clickable links
  output$wizard_steps_ui <- shiny::renderUI({
    step   <- rv$wizard_step
    labels <- c("Configure & Upload", "Training PCA", "Test Projection")
    circles <- lapply(seq_along(labels), function(i) {
      cls       <- if (i < step) "wizard-step done" else if (i == step) "wizard-step active" else "wizard-step"
      icon_html <- if (i < step) shiny::icon("check") else as.character(i)
      circle_content <- shiny::div(class = "step-circle", icon_html)
      label_content  <- shiny::div(class = "step-label", labels[[i]])
      # Done steps get a clickable action link to jump back
      if (i < step) {
        shiny::div(class = cls,
          shiny::tags$a(
            href    = "#",
            onclick = sprintf("Shiny.setInputValue('wizard_jump', %d, {priority:'event'}); return false;", i),
            style   = "text-decoration:none; color:inherit;",
            circle_content
          ),
          shiny::tags$a(
            href    = "#",
            onclick = sprintf("Shiny.setInputValue('wizard_jump', %d, {priority:'event'}); return false;", i),
            style   = "text-decoration:none; color:inherit;",
            label_content
          )
        )
      } else {
        shiny::div(class = cls, circle_content, label_content)
      }
    })
    shiny::div(class = "wizard-steps", circles)
  })

  # Jump to any previous step when a completed circle is clicked
  shiny::observeEvent(input$wizard_jump, {
    target <- as.integer(input$wizard_jump)
    if (!is.na(target) && target >= 1L && target < rv$wizard_step) {
      rv$wizard_step <- target
    }
  })

  # ================================================================
  # TRAINING WIZARD — step indicator + navigation
  # ================================================================

  output$train_wizard_step <- shiny::reactive({ rv_train$step })
  shiny::outputOptions(output, "train_wizard_step", suspendWhenHidden = FALSE)

  output$train_wizard_steps_ui <- shiny::renderUI({
    step   <- rv_train$step
    labels <- c("Upload Data", "Review & Run", "Results")
    circles <- lapply(seq_along(labels), function(i) {
      cls  <- if (i < step) "wizard-step done" else if (i == step) "wizard-step active" else "wizard-step"
      icon_html <- if (i < step) shiny::icon("check") else as.character(i)
      circle_content <- shiny::div(class = "step-circle", icon_html)
      label_content  <- shiny::div(class = "step-label", labels[[i]])
      # Done steps are clickable to jump back
      if (i < step) {
        shiny::div(class = cls,
          shiny::tags$a(
            href    = "#",
            onclick = sprintf("Shiny.setInputValue('train_wizard_jump', %d, {priority:'event'}); return false;", i),
            style   = "text-decoration:none; color:inherit;",
            circle_content
          ),
          shiny::tags$a(
            href    = "#",
            onclick = sprintf("Shiny.setInputValue('train_wizard_jump', %d, {priority:'event'}); return false;", i),
            style   = "text-decoration:none; color:inherit;",
            label_content
          )
        )
      } else {
        shiny::div(class = cls, circle_content, label_content)
      }
    })
    shiny::div(class = "wizard-steps", circles)
  })

  shiny::observeEvent(input$train_wizard_jump, {
    target <- as.integer(input$train_wizard_jump)
    if (!is.na(target) && target >= 1L && target < rv_train$step) {
      rv_train$step <- target
    }
  })

  shiny::observeEvent(input$train_step1_next, {
    if (is.null(input$train_counts_file)) {
      shiny::showNotification("Please upload a training count matrix.", type = "warning")
      return()
    }
    if (is.null(input$train_meta_file)) {
      shiny::showNotification("Please upload training metadata.", type = "warning")
      return()
    }
    rv_train$step <- 2L
  })

  shiny::observeEvent(input$train_step2_back, { rv_train$step <- 1L })
  shiny::observeEvent(input$train_step3_back, { rv_train$step <- 1L })

  # Custom phenotype labels — accumulate user additions
  rv_custom_phenotypes <- shiny::reactiveVal(character(0))

  shiny::observeEvent(input$train_add_phenotype_btn, {
    new_label <- trimws(input$train_new_phenotype_name)
    if (nchar(new_label) == 0) return()
    current <- rv_custom_phenotypes()
    if (new_label %in% current) {
      shiny::showNotification(paste0("'", new_label, "' is already in the list."), type = "warning")
      return()
    }
    rv_custom_phenotypes(c(current, new_label))
    shiny::updateTextInput(session, "train_new_phenotype_name", value = "")
    shiny::showNotification(
      paste0("Phenotype '", new_label, "' added. Samples labelled with this value will be included."),
      type = "message", duration = 5
    )
  })

  output$train_custom_phenotypes_ui <- shiny::renderUI({
    lbls <- rv_custom_phenotypes()
    if (length(lbls) == 0) return(NULL)
    shiny::tags$div(
      shiny::tags$small("Custom phenotypes:", style = "color:#555; font-weight:600;"),
      shiny::tags$ul(
        lapply(lbls, function(l) shiny::tags$li(l, style = "font-size:0.85em;")),
        style = "margin:2px 0 0 12px; padding:0;"
      )
    )
  })

  # Training metadata — reactive reader
  rv_train_meta_raw <- shiny::reactive({
    shiny::req(input$train_meta_file)
    tryCatch(
      read.csv(input$train_meta_file$datapath, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  })

  # Training normalisation options — scale checkbox is greyed out when CPM box is ticked
  output$train_norm_options_ui <- shiny::renderUI({
    already_cpm <- isTRUE(input$train_counts_already_cpm)
    scale_val   <- if (already_cpm) FALSE else isTRUE(input$train_scale_counts)
    shiny::div(
      style = "background:#fff8e1; border-left:4px solid #f39c12; padding:8px 10px; border-radius:3px;",
      shiny::checkboxInput(
        "train_counts_already_cpm",
        shiny::tags$span(shiny::icon("exclamation-circle", style = "color:#e67e22;"),
                         " Data already CPM-normalised (skip library-size normalisation)",
                         style = "font-size:0.9em;"),
        value = already_cpm
      ),
      shiny::div(
        style = if (already_cpm) "opacity:0.45; pointer-events:none;" else "",
        shiny::checkboxInput(
          "train_scale_counts",
          shiny::tags$span(shiny::icon("compress-alt", style = "color:#e67e22;"),
                           " Scale counts before normalisation (recount3 / coverage counts)",
                           style = "font-size:0.9em;"),
          value = scale_val
        )
      ),
      if (already_cpm)
        shiny::tags$p("\u2139 Scaling unavailable \u2014 data is already CPM-normalised.",
                      style = "font-size:0.82em; color:#b0a070; font-style:italic; margin:2px 0 0 22px;")
      else
        shiny::tags$p("Tick scale if counts come from recount3 / coverage-sum pipelines.",
                      style = "font-size:0.82em; color:#7f6000; margin:2px 0 0 22px;")
    )
  })

  # Batch column picker for training metadata
  output$train_batch_col_ui <- shiny::renderUI({
    meta_raw <- rv_train_meta_raw()
    if (is.null(meta_raw)) return(NULL)
    cols <- setdiff(colnames(meta_raw), "Sample_ID")
    shinyWidgets::pickerInput(
      inputId  = "train_batch_cols",
      label    = "Regress out (batch columns)",
      choices  = cols,
      selected = if ("Study_ID" %in% cols) "Study_ID" else NULL,
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(
        actionsBox        = TRUE,
        liveSearch        = TRUE,
        noneSelectedText  = "None (no batch correction)",
        selectedTextFormat = "count > 1"
      )
    )
  })

  # Phenotype column picker
  output$train_phenotype_col_ui <- shiny::renderUI({
    meta_raw <- rv_train_meta_raw()
    if (is.null(meta_raw)) return(NULL)
    cols <- setdiff(colnames(meta_raw), "Sample_ID")
    shiny::selectInput(
      "train_phenotype_col",
      label   = shiny::tags$span("Phenotype column",
                                  shiny::tags$small(" — class to predict", style = "color:#888;")),
      choices  = cols,
      selected = if ("Phenotype" %in% cols) "Phenotype" else cols[1]
    )
  })

  # Summary shown on step 2
  output$train_data_summary_ui <- shiny::renderUI({
    meta_raw <- rv_train_meta_raw()
    if (is.null(meta_raw)) return(shiny::p("No metadata loaded.", class = "info-text"))

    n_samples   <- nrow(meta_raw)
    ph_col      <- input$train_phenotype_col
    ph_summary  <- if (!is.null(ph_col) && ph_col %in% colnames(meta_raw)) {
      tbl <- table(meta_raw[[ph_col]])
      paste(paste0(names(tbl), " (n=", tbl, ")"), collapse = ", ")
    } else "—"

    max_app_raw <- input$train_max_appearances
    unlimited   <- isTRUE(input$train_max_appearances_unlimited)
    max_app_str <- if (unlimited) "no limit"
                   else paste0("max ", max_app_raw, " per gene")

    shiny::tagList(
      shiny::tags$ul(
        shiny::tags$li(sprintf("Samples in metadata: %d", n_samples)),
        shiny::tags$li(sprintf("Phenotype classes: %s", ph_summary)),
        shiny::tags$li(sprintf("Max ratios per panel: %d",
                               as.integer(input$train_max_ratios %||% 50))),
        shiny::tags$li(sprintf("Candidate pool per panel: %d",
                               as.integer(input$train_candidate_pool %||% 100))),
        shiny::tags$li(sprintf("Diversity constraint: %s", max_app_str)),
        shiny::tags$li(sprintf("Minimum η²: %.2f  |  Minimum Δη²: %.2f",
                               input$train_eta_sq_min %||% 0.30,
                               input$train_synergy_delta %||% 0.05))
      )
    )
  })

  # ── Training run pipeline ────────────────────────────────────────────────
  shiny::observeEvent(input$run_train_btn, {

    if (is.null(input$train_counts_file)) {
      shiny::showNotification("Please upload a count matrix.", type = "error"); return()
    }
    if (is.null(input$train_meta_file)) {
      shiny::showNotification("Please upload metadata.", type = "error"); return()
    }

    w <- waiter::Waiter$new(html = waiter::spin_flower(), color = waiter::transparent(0.7))
    w$show()
    on.exit(w$hide(), add = TRUE)

    withCallingHandlers({

      # 1. Read counts
      shinyWidgets::updateProgressBar(session, "train_pb", value = 5,
                                      title = "Reading count matrix…")
      train_counts_raw <- tryCatch({
        df <- read.csv(input$train_counts_file$datapath,
                       header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
        gene_col <- df[[1]]; df <- df[, -1, drop = FALSE]
        mat <- as.matrix(df); storage.mode(mat) <- "numeric"
        rownames(mat) <- gene_col; mat
      }, error = function(e) {
        shiny::showNotification(paste("Count read error:", e$message), type = "error"); NULL
      })
      if (is.null(train_counts_raw)) return()

      # 2. Read metadata
      shinyWidgets::updateProgressBar(session, "train_pb", value = 10,
                                      title = "Reading metadata…")
      train_meta <- tryCatch(
        read.csv(input$train_meta_file$datapath, stringsAsFactors = FALSE),
        error = function(e) { shiny::showNotification(paste("Metadata error:", e$message), type = "error"); NULL }
      )
      if (is.null(train_meta)) return()
      if (!"Sample_ID" %in% colnames(train_meta)) {
        shiny::showNotification("Metadata must have 'Sample_ID' column.", type = "error"); return()
      }

      ph_col <- input$train_phenotype_col
      if (is.null(ph_col) || !ph_col %in% colnames(train_meta)) {
        shiny::showNotification(
          paste("Phenotype column", shQuote(ph_col), "not found in metadata."), type = "error"); return()
      }

      # Align counts and metadata
      shared <- intersect(colnames(train_counts_raw), train_meta$Sample_ID)
      if (length(shared) == 0) {
        shiny::showNotification("No matching Sample_IDs between counts and metadata.", type = "error"); return()
      }
      train_counts_raw <- train_counts_raw[, shared, drop = FALSE]
      train_meta       <- train_meta[match(shared, train_meta$Sample_ID), , drop = FALSE]

      # 3. Load ECM dictionary
      shinyWidgets::updateProgressBar(session, "train_pb", value = 15,
                                      title = "Loading ECM dictionary…")
      ecm_dict <- if (!is.null(input$train_dict_file) && input$train_dict_source == "upload") {
        tryCatch(
          read.csv(input$train_dict_file$datapath, stringsAsFactors = FALSE),
          error = function(e) { shiny::showNotification(paste("Dict error:", e$message), type = "error"); NULL }
        )
      } else {
        tryCatch(
          read.csv(file.path(training_data_dir, "ECMExplained_Master_List_Expanded_Fixed.csv"),
                   stringsAsFactors = FALSE),
          error = function(e) { shiny::showNotification(paste("Built-in dict error:", e$message), type = "error"); NULL }
        )
      }
      if (is.null(ecm_dict)) return()
      if (!all(c("Gene", "ECMExplained_Function") %in% colnames(ecm_dict))) {
        shiny::showNotification("ECM dictionary must have 'Gene' and 'ECMExplained_Function' columns.", type = "error"); return()
      }

      # 4. Scale (recount3) if requested
      if (isTRUE(input$train_scale_counts)) {
        train_counts_raw <- transform_counts(train_counts_raw)
        shiny::showNotification("Training counts scaled (÷ colSums/1e6).", type = "message", duration = 5)
      }

      # 5. Normalise
      shinyWidgets::updateProgressBar(session, "train_pb", value = 25,
                                      title = "Normalising training counts…")

      norm_res_train <- tryCatch(
        normalize_counts(
          counts_matrix   = train_counts_raw,
          skip_rnaseq_cpm = isTRUE(input$train_counts_already_cpm)
        ),
        error = function(e) { shiny::showNotification(paste("Normalisation error:", e$message), type = "error"); NULL }
      )
      if (is.null(norm_res_train)) return()
      train_log_cpm <- norm_res_train$log_cpm

      # 6. Batch correction
      bc_col_train     <- NULL   # initialise before conditional — avoids undefined-variable error
      batch_cols_train <- input$train_batch_cols
      if (!is.null(batch_cols_train) && length(batch_cols_train) > 0) {
        shinyWidgets::updateProgressBar(session, "train_pb", value = 35,
                                        title = "Applying batch correction…")
        bc_col_train <- if (length(batch_cols_train) == 1) {
          batch_cols_train
        } else {
          cpd <- apply(train_meta[, batch_cols_train, drop = FALSE], 1, paste, collapse = "_")
          train_meta$CompoundBatch__ <- cpd
          "CompoundBatch__"
        }

        train_log_cpm <- tryCatch({
          # Protect phenotype from being regressed out
          bc_res <- apply_batch_correction(
            log_cpm_mat  = train_log_cpm,
            metadata_df  = train_meta,
            batch_col    = bc_col_train,
            protect_cols = ph_col
          )
          bc_res$corrected_mat
        }, error = function(e) {
          shiny::showNotification(paste("Batch correction failed (skipped):", e$message), type = "warning")
          train_log_cpm
        })
      }

      # 7. Feature selection — V19 η² greedy OR V3 TSP-Δ stratified
      shinyWidgets::updateProgressBar(session, "train_pb", value = 55,
                                      title = "Running synergy filter (this may take a moment)…")
      max_app <- if (isTRUE(input$train_max_appearances_unlimited)) NULL else
                   as.integer(input$train_max_appearances %||% 1L)

      sel_method <- input$train_selection_method %||% "eta_sq"

      selected_features <- if (sel_method == "tsp_delta") {
        tryCatch(
          run_feature_selection_tsp(
            logcpm_mat    = train_log_cpm,
            ecm_dict      = ecm_dict,
            phenotype_vec = factor(train_meta[[ph_col]]),
            n_per_pheno   = as.integer(input$train_n_per_pheno %||% 12L),
            eta_sq_min    = as.numeric(input$train_eta_sq_min  %||% 0.10),
            var_min       = 0.005
          ),
          error = function(e) {
            shiny::showNotification(paste("TSP-Delta selection error:", e$message), type = "error")
            NULL
          }
        )
      } else {
        tryCatch(
          run_feature_selection(
            train_expr          = train_log_cpm,
            ecm_dict            = ecm_dict,
            phenotype_vec       = factor(train_meta[[ph_col]]),
            max_ratios          = as.integer(input$train_max_ratios %||% 50),
            candidate_pool_size = as.integer(input$train_candidate_pool %||% 100L),
            max_appearances     = max_app,
            eta_sq_min          = as.numeric(input$train_eta_sq_min %||% 0.30),
            synergy_delta       = as.numeric(input$train_synergy_delta %||% 0.05)
          ),
          error = function(e) {
            shiny::showNotification(paste("Feature selection error:", e$message), type = "error")
            NULL
          }
        )
      }
      if (is.null(selected_features)) return()

      shiny::showNotification(
        sprintf("%d features selected (%d AD, %d SIG).",
                nrow(selected_features),
                sum(selected_features$Panel == "AD"),
                sum(selected_features$Panel == "SIG")),
        type = "message", duration = 8
      )

      # 8. Compute training ratios + train model
      shinyWidgets::updateProgressBar(session, "train_pb", value = 75,
                                      title = "Training PCA + LDA model…")
      train_ratio_res <- tryCatch(
        compute_ratios(train_log_cpm, selected_features),
        error = function(e) { shiny::showNotification(paste("Ratio computation error:", e$message), type = "error"); NULL }
      )
      if (is.null(train_ratio_res)) return()

      new_train_obj <- tryCatch(
        train_ecmexplainer(
          ratio_mat      = train_ratio_res$ratio_mat,
          avail_features = train_ratio_res$avail_features,
          phenotype_vec  = factor(train_meta[[ph_col]][
            match(colnames(train_ratio_res$ratio_mat), train_meta$Sample_ID)
          ]),
          metadata_df    = train_meta,
          batch_col      = bc_col_train,
          protect_cols   = ph_col
        ),
        error = function(e) { shiny::showNotification(paste("Training error:", e$message), type = "error"); NULL }
      )
      if (is.null(new_train_obj)) return()

      # Store results
      # Ensure 'Phenotype' column always exists in train_meta so that
      # plot_lda_training() / plot_pca_training() (which hard-code that column
      # name) work regardless of what the user named their phenotype column.
      if (!is.null(ph_col) && ph_col != "Phenotype" && ph_col %in% colnames(train_meta)) {
        train_meta$Phenotype <- train_meta[[ph_col]]
      }

      rv_train$train_obj  <- new_train_obj
      rv_train$features   <- selected_features
      rv_train$train_meta <- train_meta
      rv_train$train_expr <- train_log_cpm
      rv_train$ecm_dict   <- ecm_dict

      shinyWidgets::updateProgressBar(session, "train_pb", value = 100, title = "Done!")
      rv_train$step <- 3L
      shiny::showNotification("Training complete! Model ready.", type = "message", duration = 6)

    }, warning = function(w) {
      shiny::showNotification(paste("Warning:", conditionMessage(w)), type = "warning", duration = 8)
      invokeRestart("muffleWarning")
    })
  })

  # Training run status
  output$train_run_status_ui <- shiny::renderUI({
    if (rv_train$step < 3L) return(NULL)
    n_feat <- if (!is.null(rv_train$features)) nrow(rv_train$features) else 0L
    shiny::tags$p(shiny::icon("check-circle", style = "color:#27ae60;"),
                  sprintf(" Training complete. %d features selected.", n_feat),
                  style = "font-weight:700; color:#27ae60;")
  })

  # Training results summary
  output$train_results_summary_ui <- shiny::renderUI({
    shiny::req(rv_train$train_obj, rv_train$features, rv_train$train_meta)
    n_samp <- ncol(rv_train$train_expr)
    n_feat <- nrow(rv_train$features)
    n_ad   <- sum(rv_train$features$Panel == "AD")
    n_sig  <- sum(rv_train$features$Panel == "SIG")
    classes <- levels(rv_train$train_obj$phenotype_vec)
    shiny::tagList(
      shiny::fluidRow(
        shiny::column(3, shiny::tags$div(class = "ratio-badge", n_samp),
                      shiny::p("Training samples", class = "info-text")),
        shiny::column(3, shiny::tags$div(class = "ratio-badge", n_feat),
                      shiny::p("Selected ratios", class = "info-text")),
        shiny::column(3, shiny::tags$div(class = "ratio-badge",
                                         sprintf("%d AD / %d SIG", n_ad, n_sig)),
                      shiny::p("Panel breakdown", class = "info-text")),
        shiny::column(3, shiny::tags$div(class = "ratio-badge", length(classes)),
                      shiny::p(paste("Classes:", paste(classes, collapse = ", ")),
                               class = "info-text"))
      )
    )
  })

  # LOOCV accuracy + batch correction level badge
  output$train_loocv_ui <- shiny::renderUI({
    shiny::req(rv_train$train_obj)
    loocv    <- rv_train$train_obj$loocv_acc
    bc_level <- rv_train$train_obj$batch_correction_level %||% "NONE"
    if (is.null(loocv) || is.na(loocv)) return(NULL)
    shiny::div(
      style = paste0(
        "background:#eafaf1; border-left:4px solid #27ae60;",
        " padding:10px 14px; border-radius:4px; margin:10px 0;"
      ),
      shiny::tags$span(
        shiny::icon("check-double", style = "color:#27ae60;"),
        sprintf(" LOOCV Accuracy: %.1f%%", loocv),
        style = "font-size:1.15em; font-weight:bold; color:#1e8449;"
      ),
      shiny::tags$span(
        sprintf("  |  Batch correction: %s", bc_level),
        style = "font-size:0.9em; color:#555; margin-left:10px;"
      )
    )
  })

  # Training features table
  output$train_features_table <- DT::renderDT({
    shiny::req(rv_train$features)
    DT::datatable(
      rv_train$features[, c("Panel", "Gene_A", "Gene_B", "eta_sq",
                             "max_single", "synergy_delta"), drop = FALSE],
      rownames = FALSE, class = "stripe hover compact",
      options  = list(pageLength = 15, scrollX = TRUE, dom = "frtip")
    ) |> DT::formatRound(c("eta_sq", "max_single", "synergy_delta"), digits = 3)
  })

  # Training wizard LDA + PCA plots
  output$plot_lda_train_wizard <- shiny::renderPlot({
    shiny::req(rv_train$train_obj, rv_train$train_meta)
    plot_lda_training(rv_train$train_obj, rv_train$train_meta)
  })
  output$plot_pca_train_wizard2 <- shiny::renderPlot({
    shiny::req(rv_train$train_obj, rv_train$train_meta)
    ph_col <- input$train_phenotype_col %||% "Phenotype"
    plot_pca_training_meta(rv_train$train_obj, rv_train$train_meta, ph_col)
  })

  # ----------------------------------------------------------------
  # Helper: build the training_params record from current UI inputs.
  # Called at save time so the record is always up to date.
  # ----------------------------------------------------------------
  .build_training_params <- function() {
    species_val <- input$train_species %||% "Human"
    if (identical(species_val, "Other"))
      species_val <- paste0("Other: ", trimws(input$train_species_custom %||% ""))

    list(
      model_name              = trimws(input$train_model_name %||% ""),
      species                 = species_val,
      max_ratios_per_panel    = as.integer(input$train_max_ratios        %||% 50L),
      candidate_pool          = as.integer(input$train_candidate_pool    %||% 100L),
      eta_sq_min              = as.numeric(input$train_eta_sq_min        %||% 0.30),
      synergy_delta           = as.numeric(input$train_synergy_delta     %||% 0.05),
      max_appearances         = if (isTRUE(input$train_max_appearances_unlimited)) NA_integer_
                                else as.integer(input$train_max_appearances %||% 1L),
      max_appearances_unlimited = isTRUE(input$train_max_appearances_unlimited),
      counts_already_cpm      = isTRUE(input$train_counts_already_cpm),
      scale_counts            = isTRUE(input$train_scale_counts),
      dict_source             = input$train_dict_source %||% "builtin",
      batch_cols              = input$train_batch_cols  %||% character(0),
      phenotype_col           = input$train_phenotype_col %||% "Phenotype",
      loocv_acc               = rv_train$train_obj$loocv_acc,
      batch_correction_level  = rv_train$train_obj$batch_correction_level %||% "NONE",
      n_training_samples      = if (!is.null(rv_train$train_expr)) ncol(rv_train$train_expr) else NA_integer_,
      n_features_selected     = if (!is.null(rv_train$features)) nrow(rv_train$features) else NA_integer_,
      phenotype_classes       = if (!is.null(rv_train$train_obj)) levels(rv_train$train_obj$phenotype_vec) else character(0),
      created                 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ecmexplainer_version    = "V19"
    )
  }

  # Download model bundle
  output$dl_train_model_bundle <- shiny::downloadHandler(
    filename = function() {
      nm <- trimws(input$train_model_name %||% "")
      if (nchar(nm) == 0) nm <- paste0("ECMExplainER_ModelBundle_", Sys.Date())
      paste0(nm, ".rds")
    },
    content  = function(file) {
      shiny::req(rv_train$train_obj, rv_train$features)
      params <- .build_training_params()
      bundle <- list(
        model_name      = params$model_name,
        training_params = params,
        train_obj       = rv_train$train_obj,
        features        = rv_train$features,
        train_meta      = rv_train$train_meta,
        train_expr      = rv_train$train_expr,
        ecm_dict        = rv_train$ecm_dict,
        version         = paste0("ECMExplainER-custom-", Sys.Date())
      )
      saveRDS(bundle, file)
    }
  )

  # Download training params as JSON log file
  output$dl_train_params_json <- shiny::downloadHandler(
    filename = function() {
      nm <- trimws(input$train_model_name %||% "")
      if (nchar(nm) == 0) nm <- paste0("ECMExplainER_TrainingParams_", Sys.Date())
      paste0(nm, "_params.json")
    },
    content  = function(file) {
      shiny::req(rv_train$train_obj, rv_train$features)
      params <- .build_training_params()
      if (!requireNamespace("jsonlite", quietly = TRUE))
        stop("Package 'jsonlite' is required to export params as JSON.")
      jsonlite::write_json(params, path = file, pretty = TRUE, auto_unbox = TRUE)
    }
  )

  # Download selected features CSV
  output$dl_train_features_csv <- shiny::downloadHandler(
    filename = function() {
      nm <- trimws(input$train_model_name %||% "")
      if (nchar(nm) == 0) nm <- paste0("ECMExplainER_SelectedFeatures_", Sys.Date())
      paste0(nm, "_features.csv")
    },
    content  = function(file) {
      shiny::req(rv_train$features)
      write.csv(rv_train$features, file, row.names = FALSE)
    }
  )

  # ----------------------------------------------------------------
  # Restore training parameters from an uploaded .json or .rds file.
  # Triggered whenever the user selects a file in train_params_upload.
  # ----------------------------------------------------------------
  .apply_training_params <- function(params) {
    # Numeric inputs
    if (!is.null(params$max_ratios_per_panel))
      shiny::updateNumericInput(session, "train_max_ratios",     value = as.integer(params$max_ratios_per_panel))
    if (!is.null(params$candidate_pool))
      shiny::updateNumericInput(session, "train_candidate_pool", value = as.integer(params$candidate_pool))
    if (!is.null(params$eta_sq_min))
      shiny::updateNumericInput(session, "train_eta_sq_min",     value = as.numeric(params$eta_sq_min))
    if (!is.null(params$synergy_delta))
      shiny::updateNumericInput(session, "train_synergy_delta",  value = as.numeric(params$synergy_delta))

    # Max appearances slider + unlimited checkbox
    unlimited <- isTRUE(params$max_appearances_unlimited) || identical(params$max_appearances, NA) ||
                 identical(params$max_appearances, NA_integer_)
    shiny::updateCheckboxInput(session, "train_max_appearances_unlimited", value = unlimited)
    if (!unlimited && !is.null(params$max_appearances) && !is.na(params$max_appearances))
      shiny::updateSliderInput(session, "train_max_appearances",
                               value = min(5L, max(1L, as.integer(params$max_appearances))))

    # Species
    if (!is.null(params$species)) {
      known <- c("Human", "Mouse", "Rat")
      if (params$species %in% known) {
        shiny::updateSelectInput(session, "train_species", selected = params$species)
      } else {
        shiny::updateSelectInput(session, "train_species", selected = "Other")
        shiny::updateTextInput(session, "train_species_custom",
                               value = sub("^Other:\\s*", "", params$species))
      }
    }

    # Model name
    if (!is.null(params$model_name) && nchar(trimws(params$model_name)) > 0)
      shiny::updateTextInput(session, "train_model_name", value = params$model_name)

    # Checkboxes
    if (!is.null(params$counts_already_cpm))
      shiny::updateCheckboxInput(session, "train_counts_already_cpm", value = isTRUE(params$counts_already_cpm))
    if (!is.null(params$scale_counts))
      shiny::updateCheckboxInput(session, "train_scale_counts", value = isTRUE(params$scale_counts))
  }

  # Reactive: extract training_params from whatever file was uploaded (JSON or RDS)
  rv_restored_params <- shiny::reactive({
    shiny::req(input$train_params_upload)
    path <- input$train_params_upload$datapath
    ext  <- tolower(tools::file_ext(input$train_params_upload$name))
    tryCatch({
      if (ext == "json") {
        if (!requireNamespace("jsonlite", quietly = TRUE))
          stop("jsonlite required to read .json files.")
        jsonlite::read_json(path, simplifyVector = TRUE)
      } else if (ext == "rds") {
        bundle <- readRDS(path)
        bundle$training_params
      } else {
        NULL
      }
    }, error = function(e) {
      shiny::showNotification(paste("Could not read params file:", e$message), type = "error")
      NULL
    })
  })

  shiny::observeEvent(rv_restored_params(), {
    params <- rv_restored_params()
    if (is.null(params)) return()
    .apply_training_params(params)
    shiny::showNotification(
      shiny::tags$span(shiny::icon("check-circle"), " Training parameters restored from uploaded file."),
      type = "message", duration = 6
    )
  })

  # Small status UI shown below the upload input
  output$train_params_restore_ui <- shiny::renderUI({
    params <- rv_restored_params()
    if (is.null(params)) return(NULL)
    nm  <- params$model_name %||% "(unnamed)"
    sp  <- params$species    %||% "?"
    mr  <- params$max_ratios_per_panel %||% "?"
    shiny::div(
      style = "background:#eaf7ea; border-left:4px solid #27ae60; padding:8px 12px; border-radius:3px; margin-top:4px;",
      shiny::icon("check", style = "color:#27ae60;"),
      shiny::tags$b(" Parameters restored"),
      shiny::tags$br(),
      sprintf("Model: %s  |  Species: %s  |  Ratios/panel: %s", nm, sp, mr)
    )
  })

  # Step 1 → 2
  # Step 2 → 1 (back)
  shiny::observeEvent(input$step2_back, { rv$wizard_step <- 1L })

  # Step 2 → 3 (after successful run)
  shiny::observeEvent(input$step3_back, { rv$wizard_step <- 1L })
  shiny::observeEvent(input$step3_next, {
    shiny::req(rv$run_success)
    rv$wizard_step <- 3L
  })

  # Step 3 → 2 (back)
  shiny::observeEvent(input$step4_back, { rv$wizard_step <- 2L })

  # ----------------------------------------------------------------
  # Dynamic condition selector — populated from whichever training
  # data is currently active (built-in OR custom).
  # Re-renders whenever test_model_source or the custom model file changes.
  # ----------------------------------------------------------------
  output$condition_selector_ui <- shiny::renderUI({
    at <- tryCatch(resolve_active_training(), error = function(e) training)
    if (is.null(at) || is.null(at$train_meta)) {
      return(shiny::p("Training data not yet loaded.", class = "info-text"))
    }
    pheno_col <- if ("Phenotype" %in% colnames(at$train_meta)) "Phenotype" else colnames(at$train_meta)[2]
    avail_conds <- sort(unique(at$train_meta[[pheno_col]]))

    cond_col <- .cond_pal(avail_conds)

    n_samp_per <- table(at$train_meta[[pheno_col]])
    choices_named <- setNames(avail_conds,
      paste0(avail_conds, " (n=", n_samp_per[avail_conds], ")"))

    shiny::tagList(
      shiny::p(
        shiny::icon("info-circle"),
        sprintf(" %d phenotype classes available in the selected training data.",
                length(avail_conds)),
        class = "info-text"
      ),
      shinyWidgets::checkboxGroupButtons(
        inputId   = "selected_conditions",
        label     = NULL,
        choices   = choices_named,
        selected  = avail_conds,
        checkIcon = list(yes = shiny::icon("check"), no = shiny::icon("times")),
        status    = "primary",
        size      = "normal"
      ),
      shiny::br(),
      shiny::uiOutput("condition_validation_ui")
    )
  })

  # Condition validation message
  output$condition_validation_ui <- shiny::renderUI({
    conds <- input$selected_conditions
    n <- length(conds)
    if (n < 2) {
      shiny::tags$p(shiny::icon("exclamation-triangle"),
                    " Please select at least 2 conditions.",
                    style = "color:#c0392b; font-weight:600;")
    } else {
      shiny::tags$p(shiny::icon("check-circle"),
                    sprintf(" %d conditions selected.", n),
                    style = "color:#27ae60; font-weight:600;")
    }
  })

  # Condition summary card — also dynamic
  output$condition_summary_ui <- shiny::renderUI({
    at  <- tryCatch(resolve_active_training(), error = function(e) training)
    sel <- input$selected_conditions
    if (is.null(at) || is.null(at$train_meta) || is.null(sel)) return(NULL)

    pheno_col <- if ("Phenotype" %in% colnames(at$train_meta)) "Phenotype" else colnames(at$train_meta)[2]
    avail_conds <- sort(unique(at$train_meta[[pheno_col]]))

    cond_col <- .cond_pal(avail_conds)

    chips <- lapply(avail_conds, function(c) {
      col <- if (!is.null(cond_col[[c]])) cond_col[[c]] else "#999"
      if (c %in% sel) {
        shiny::span(c, style = sprintf(
          "display:inline-block; margin:3px 4px; padding:3px 10px; border-radius:12px;
           background:%s; color:#fff; font-size:0.82em; font-weight:600;", col))
      } else {
        shiny::span(c, style =
          "display:inline-block; margin:3px 4px; padding:3px 10px; border-radius:12px;
           background:#ddd; color:#888; font-size:0.82em; text-decoration:line-through;")
      }
    })
    n_train <- sum(at$train_meta[[pheno_col]] %in% sel)
    shiny::tagList(
      shiny::div(style = "line-height:2;", chips),
      shiny::p(sprintf("(%d training samples selected)", n_train), class = "info-text")
    )
  })

  # Chips shown in run box
  output$selected_conditions_chips <- shiny::renderUI({
    sel <- input$selected_conditions
    if (is.null(sel) || length(sel) == 0) return(NULL)
    pal  <- .cond_pal(sel)
    chips <- lapply(sel, function(c) {
      col <- pal[[c]] %||% "#2d6a9f"
      shiny::span(c, style = sprintf(
        "display:inline-block; margin:2px 3px; padding:2px 8px; border-radius:10px;
         background:%s; color:#fff; font-size:0.8em;", col))
    })
    shiny::tagList(chips)
  })

  # ================================================================
  # UI helpers
  # ================================================================

  # Info box for loaded custom model bundle
  # NOTE: superseded by active_model_info_ui (unified for all 3 model sources).
  # Kept as a no-op so any residual uiOutput("test_model_info_ui") calls still work.
  output$test_model_info_ui <- shiny::renderUI({ NULL })

  # When a bundle .rds is uploaded in the Test Wizard, restore the training
  # parameters into the Training Wizard inputs so re-training is easy.
  shiny::observeEvent(input$test_model_file, {
    shiny::req(input$test_model_file)
    bundle <- tryCatch(readRDS(input$test_model_file$datapath), error = function(e) NULL)
    params <- bundle$training_params
    if (is.null(params)) return()
    .apply_training_params(params)
    shiny::showNotification(
      shiny::tags$span(
        shiny::icon("magic"),
        " Training parameters restored from model bundle — check Training Wizard to re-run."
      ),
      type = "message", duration = 7
    )
  })

  # Group/condition column picker for test metadata (which column labels sample groups)
  output$test_group_col_ui <- shiny::renderUI({
    meta_raw <- rv_meta_raw()
    if (is.null(meta_raw)) return(NULL)
    # Offer all non-Sample_ID columns; try to guess a sensible default
    cols <- setdiff(colnames(meta_raw), "Sample_ID")
    if (length(cols) == 0) return(NULL)
    default_col <- if ("Condition" %in% cols) "Condition" else
                   if ("Phenotype" %in% cols) "Phenotype" else
                   if ("Group"     %in% cols) "Group"     else cols[1]
    shiny::selectInput(
      inputId  = "test_group_col",
      label    = shiny::tags$span(
        shiny::icon("tag"),
        " Sample group / condition column",
        shiny::tags$small(" (used for plot colouring and condition chips)",
                          style = "color:#888; font-weight:400;")
      ),
      choices  = cols,
      selected = default_col
    )
  })

  # Mode B anchor column picker (populated from test metadata)
  output$test_anchor_col_ui <- shiny::renderUI({
    meta_raw <- rv_meta_raw()
    if (is.null(meta_raw)) return(NULL)
    cols <- setdiff(colnames(meta_raw), "Sample_ID")
    shiny::selectInput(
      "test_anchor_col", "Reference group column",
      choices  = c("(choose)" = "", cols),
      selected = ""
    )
  })

  # ── Mode B: multi-pair anchor system ──────────────────────────────────────
  # rv_anchor_pairs: list of pairs, each pair = list(vals = char vec, target = string)
  # Pairs are added/removed dynamically; the first pair is always present.

  rv_anchor_pairs <- shiny::reactiveVal(list(
    list(vals = character(0), target = "")
  ))

  # Reset pairs when the anchor column changes (values in old column are meaningless)
  shiny::observeEvent(input$test_anchor_col, {
    rv_anchor_pairs(list(
      list(vals = character(0), target = "")
    ))
  }, ignoreInit = TRUE)

  # Add a new (empty) anchor pair
  shiny::observeEvent(input$anchor_pair_add, {
    current <- rv_anchor_pairs()
    current[[length(current) + 1L]] <- list(vals = character(0), target = "")
    rv_anchor_pairs(current)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # Remove pair i: we register up to 20 possible remove buttons once at startup.
  # The button only exists in the DOM when pairs > 1, so clicks only fire then.
  lapply(seq_len(20L), function(idx) {
    btn_id <- paste0("anchor_pair_remove_", idx)
    shiny::observeEvent(input[[btn_id]], {
      cur <- rv_anchor_pairs()
      if (idx <= length(cur) && length(cur) > 1L) {
        cur[[idx]] <- NULL
        rv_anchor_pairs(cur)
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
  })

  # Render the full pairs UI
  output$anchor_pairs_ui <- shiny::renderUI({
    shiny::req(input$test_anchor_col, input$test_anchor_col != "")
    meta_raw <- rv_meta_raw()
    if (is.null(meta_raw) || !input$test_anchor_col %in% colnames(meta_raw)) return(NULL)

    col_vals     <- sort(unique(as.character(meta_raw[[input$test_anchor_col]])))
    pairs        <- rv_anchor_pairs()

    active_training <- isolate(resolve_active_training())
    atlas_classes <- if (!is.null(active_training$train_obj)) {
      levels(active_training$train_obj$phenotype_vec)
    } else if (!is.null(active_training$train_meta) &&
               "Phenotype" %in% colnames(active_training$train_meta)) {
      sort(unique(as.character(active_training$train_meta$Phenotype)))
    } else {
      character(0)
    }

    pair_uis <- lapply(seq_along(pairs), function(i) {
      pair <- pairs[[i]]
      cur_vals   <- if (length(pair$vals) > 0) pair$vals else NULL
      cur_target <- if (!is.null(pair$target) && pair$target %in% atlas_classes)
                      pair$target
                    else if (length(atlas_classes) > 0)
                      atlas_classes[1]
                    else
                      NULL

      shiny::wellPanel(
        style = "padding:8px 10px; margin-bottom:6px; background:#f8f9fa; border:1px solid #dee2e6;",
        shiny::fluidRow(
          shiny::column(
            5,
            shinyWidgets::pickerInput(
              inputId  = paste0("anchor_pair_vals_", i),
              label    = shiny::tags$span(
                shiny::tags$b(paste0("Pair ", i, ":")),
                " my reference value(s)",
                style = "font-size:0.88em;"
              ),
              choices  = col_vals,
              selected = cur_vals,
              multiple = TRUE,
              options  = shinyWidgets::pickerOptions(
                actionsBox       = TRUE,
                liveSearch       = TRUE,
                noneSelectedText = "Select group(s)",
                container        = "body"
              )
            )
          ),
          shiny::column(
            5,
            shiny::selectInput(
              inputId  = paste0("anchor_pair_target_", i),
              label    = shiny::tags$span(
                "\u2192 atlas phenotype",
                style = "font-size:0.88em;"
              ),
              choices  = atlas_classes,
              selected = cur_target
            )
          ),
          shiny::column(
            2,
            shiny::br(),
            if (length(pairs) > 1L)
              shiny::actionButton(
                inputId = paste0("anchor_pair_remove_", i),
                label   = shiny::icon("times"),
                class   = "btn btn-danger btn-sm",
                style   = "margin-top:4px;"
              )
            else
              shiny::tags$span()   # nothing to show for last remaining pair
          )
        )
      )
    })

    shiny::tagList(
      shiny::tags$label("Anchor pair(s): reference values \u2192 atlas phenotype",
                        style = "font-weight:600; font-size:0.9em; margin-bottom:4px; display:block;"),
      shiny::tags$small(
        shiny::icon("info-circle"),
        " Each pair maps one or more of your sample groups to a known atlas phenotype.",
        " Groups not listed here are projected without explicit alignment.",
        style = "color:#666; display:block; margin-bottom:8px;"
      ),
      shiny::tagList(pair_uis),
      shiny::actionButton(
        "anchor_pair_add",
        label = shiny::tagList(shiny::icon("plus"), " Add anchor pair"),
        class = "btn btn-default btn-sm",
        style = "margin-top:4px;"
      )
    )
  })

  # Dynamically populate batch column picker from uploaded metadata
  output$batch_col_ui <- shiny::renderUI({
    req(input$meta_file)
    meta_raw <- tryCatch(
      read.csv(input$meta_file$datapath, nrows = 1),
      error = function(e) NULL
    )
    if (is.null(meta_raw)) return(NULL)
    cols <- setdiff(colnames(meta_raw), "Sample_ID")
    shinyWidgets::pickerInput(
      inputId  = "batch_cols",
      label    = "Regress out (batch effect columns)",
      choices  = cols,
      selected = if ("Study_ID" %in% cols) "Study_ID" else NULL,
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(
        actionsBox        = TRUE,
        liveSearch        = TRUE,
        selectedTextFormat = "count > 1"
      )
    )
  })

  # Protect columns — biological covariates to preserve during ComBat
  output$protect_cols_ui <- shiny::renderUI({
    req(input$meta_file)
    # Only show when at least one batch column is selected
    if (is.null(input$batch_cols) || length(input$batch_cols) == 0) return(NULL)
    meta_raw <- tryCatch(
      read.csv(input$meta_file$datapath, nrows = 1),
      error = function(e) NULL
    )
    if (is.null(meta_raw)) return(NULL)
    # Exclude Sample_ID and any already selected batch columns
    avail_cols <- setdiff(colnames(meta_raw), c("Sample_ID", input$batch_cols))
    if (length(avail_cols) == 0) return(NULL)
    shinyWidgets::pickerInput(
      inputId  = "protect_cols",
      label    = shiny::tags$span(
        "Protect (biological covariates)",
        shiny::tags$small(" — preserved during batch correction", style = "color:#888;")
      ),
      choices  = avail_cols,
      selected = NULL,
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(
        actionsBox        = TRUE,
        liveSearch        = TRUE,
        noneSelectedText  = "None (intercept only)",
        selectedTextFormat = "count > 1"
      )
    )
  })

  # ================================================================
  # Training source resolver
  # Returns the active training bundle (built-in OR custom session OR loaded .rds)
  # ================================================================

  resolve_active_training <- shiny::reactive({
    mode <- input$test_model_source %||% "default"

    if (mode == "default") {
      # Load from the registered default models list (built-in + session extras)
      chosen_key <- input$test_default_model %||% names(ECM_DEFAULT_MODELS)[1]
      entry      <- ECM_DEFAULT_MODELS[[chosen_key]] %||%
                    tryCatch(rv_extra_defaults()[[chosen_key]], error = function(e) NULL)
      if (is.null(entry)) return(training)          # fallback: V19 foundation
      return(.load_default_bundle(entry, training_data_dir))
    }

    if (mode == "upload") {
      shiny::req(input$test_model_file)
      bundle <- tryCatch(readRDS(input$test_model_file$datapath), error = function(e) NULL)
      if (!is.null(bundle)) return(bundle)
    }

    if (mode == "session") {
      # Priority: model trained in this session
      if (!is.null(rv_train$train_obj) && !is.null(rv_train$features)) {
        return(list(
          train_obj       = rv_train$train_obj,
          features        = rv_train$features,
          train_meta      = rv_train$train_meta,
          train_expr      = rv_train$train_expr,
          ecm_dict        = rv_train$ecm_dict,
          model_name      = trimws(input$train_model_name %||% "Session model"),
          training_params = list(
            model_name = trimws(input$train_model_name %||% "Session model"),
            species    = input$train_species %||% "Human",
            loocv_acc  = rv_train$train_obj$loocv_acc
          ),
          version = paste0("session-", Sys.Date())
        ))
      }
      shiny::showNotification(
        "No session-trained model found. Please run the Training Wizard first.",
        type = "warning", duration = 5
      )
    }

    # Final fallback: built-in training data
    training
  })

  # ── Auto-detect feature_method from active model and update prediction_mode ─
  shiny::observe({
    at <- tryCatch(resolve_active_training(), error = function(e) NULL)
    if (is.null(at)) return()
    detected_method <- at$feature_method %||%
                       at$training_params$feature_method %||%
                       NULL
    if (!is.null(detected_method)) {
      shiny::updateRadioButtons(session, "prediction_mode", selected = detected_method)
    }
  })

  # ── Populate the default model dropdown on startup ─────────────────────────
  shiny::observe({
    shiny::updateSelectInput(session, "test_default_model",
                             choices  = .default_model_choices(),
                             selected = names(ECM_DEFAULT_MODELS)[1])
  })

  # ── Reactive: session-level extra defaults (added via "Make Default" button) ─
  rv_extra_defaults <- shiny::reactiveVal(list())

  # Placeholder for the "Make Default" status UI — filled in by the observer
  output$train_make_default_ui <- shiny::renderUI({ NULL })

  # ── Unified active model info panel ────────────────────────────────────────
  # Replaces the old test_model_info_ui (upload-only).  Fires for all 3 modes.
  output$active_model_info_ui <- shiny::renderUI({
    mode <- input$test_model_source %||% "default"

    at <- tryCatch(resolve_active_training(), error = function(e) NULL)
    if (is.null(at)) return(shiny::tags$p("No model loaded.", style = "color:#c0392b;"))

    # Extract info from bundle or training list
    p        <- at$training_params %||% list()
    nm       <- at$model_name %||% p$model_name %||% at$version %||% "Unknown model"
    if (nchar(trimws(nm)) == 0) nm <- "Unnamed model"
    sp       <- p$species   %||% "Human"
    loocv    <- p$loocv_acc %||% at$train_obj$loocv_acc %||% NULL
    n_samp   <- if (!is.null(at$train_expr)) ncol(at$train_expr) else
                if (!is.null(at$train_meta)) nrow(at$train_meta) else "?"
    n_feat   <- if (!is.null(at$features)) nrow(at$features) else "?"
    classes  <- if (!is.null(at$train_obj)) {
      tryCatch(levels(at$train_obj$phenotype_vec), error = function(e) NULL)
    } else if (!is.null(at$train_meta) && "Phenotype" %in% colnames(at$train_meta)) {
      sort(unique(at$train_meta$Phenotype))
    } else NULL

    # Study breakdown from train_meta
    studies <- if (!is.null(at$train_meta) && "Study_ID" %in% colnames(at$train_meta)) {
      studies_tbl <- table(at$train_meta$Study_ID)
      paste(paste0(names(studies_tbl), " (n=", studies_tbl, ")"), collapse = ", ")
    } else NULL

    mode_label <- switch(mode,
      default = shiny::tags$span(shiny::icon("star"),    " Default model", style = "color:#8e44ad; font-size:0.8em;"),
      upload  = shiny::tags$span(shiny::icon("upload"),  " Uploaded bundle", style = "color:#2d6a9f; font-size:0.8em;"),
      session = shiny::tags$span(shiny::icon("flask"),   " Session-trained model", style = "color:#27ae60; font-size:0.8em;")
    )

    shiny::div(
      style = "background:#eef4fb; border-left:4px solid #2d6a9f; border-radius:6px; padding:10px 14px; margin-top:4px;",
      mode_label, shiny::br(),
      shiny::tags$b(shiny::icon("database"), " ", nm), shiny::br(),
      shiny::tags$small(
        sprintf("%s samples  \u2022  %s features", n_samp, n_feat),
        style = "color:#555;"
      ),
      if (!is.null(classes)) shiny::tagList(
        shiny::br(),
        shiny::tags$small(
          paste("Classes:", paste(classes, collapse = ", ")),
          style = "color:#555;"
        )
      ) else NULL,
      if (!is.null(loocv)) shiny::tagList(
        shiny::br(),
        shiny::tags$small(sprintf("Species: %s  \u2022  LOOCV: %.1f%%", sp, loocv),
                          style = "color:#555;")
      ) else shiny::tagList(
        shiny::br(),
        shiny::tags$small(sprintf("Species: %s", sp), style = "color:#555;")
      ),
      if (!is.null(studies)) shiny::tagList(
        shiny::br(),
        shiny::tags$small(paste("Studies:", studies), style = "color:#777; font-size:0.85em;")
      ) else NULL
    )
  })

  # ── "Make this a default model" — saves bundle to training_data/ + registers ─
  shiny::observeEvent(input$train_make_default_btn, {
    shiny::req(rv_train$train_obj, rv_train$features)

    nm <- trimws(input$train_model_name %||% "")
    if (nchar(nm) == 0) {
      shiny::showNotification("Please enter a model name before making it a default.",
                              type = "warning"); return()
    }

    params <- .build_training_params()
    bundle <- list(
      model_name      = nm,
      training_params = params,
      train_obj       = rv_train$train_obj,
      features        = rv_train$features,
      train_meta      = rv_train$train_meta,
      train_expr      = rv_train$train_expr,
      ecm_dict        = rv_train$ecm_dict,
      version         = paste0("custom-", Sys.Date())
    )

    # Save to training_data/
    safe_nm   <- gsub("[^A-Za-z0-9_-]", "_", nm)
    dest_file <- paste0("ECMExplainer_Default_", safe_nm, ".rds")
    dest_path <- file.path(training_data_dir, dest_file)
    tryCatch({
      saveRDS(bundle, dest_path)
    }, error = function(e) {
      shiny::showNotification(paste("Could not write to training_data/:", e$message),
                              type = "error"); return()
    })

    # Register in session-level extras (does NOT modify the source file on disk —
    # persists for the lifetime of this R session only; to make permanent, add to
    # default_models.R manually with bundle_file = dest_file).
    extras <- rv_extra_defaults()
    label  <- paste0(nm, " (custom, ", Sys.Date(), ")")
    extras[[label]] <- list(
      label       = label,
      bundle_file = dest_file,
      description = paste("Custom model trained", Sys.Date())
    )
    rv_extra_defaults(extras)

    # Update the dropdown
    all_choices <- c(.default_model_choices(),
                     setNames(names(extras), names(extras)))
    shiny::updateSelectInput(session, "test_default_model",
                             choices  = all_choices,
                             selected = label)

    output$train_make_default_ui <- shiny::renderUI({
      shiny::div(
        style = "background:#eaf7ea; border-left:4px solid #27ae60; padding:8px 12px; border-radius:3px; margin-top:6px;",
        shiny::icon("check-circle", style = "color:#27ae60;"),
        shiny::tags$b(sprintf(" '%s' saved as a default model.", nm)), shiny::br(),
        shiny::tags$small(sprintf("File: training_data/%s", dest_file), style = "color:#555;"), shiny::br(),
        shiny::tags$small(
          "It is now available in the Test Wizard default model dropdown for this session.",
          " To make it permanent across restarts, add an entry to R/default_models.R.",
          style = "color:#777;"
        )
      )
    })

    shiny::showNotification(
      shiny::tags$span(shiny::icon("star"), sprintf(" '%s' registered as a default model.", nm)),
      type = "message", duration = 7
    )
  })

  # ── resolve_active_training must also check session extras ─────────────────
  # (Patch: after "Make Default" adds to rv_extra_defaults, lookup must cover both)
  .get_registry_entry <- function(key) {
    entry <- ECM_DEFAULT_MODELS[[key]]
    if (!is.null(entry)) return(entry)
    rv_extra_defaults()[[key]]
  }

  # ================================================================
  # Sample selector — reactive metadata store + table UI
  # ================================================================

  # Store the full parsed metadata whenever a new file is uploaded
  rv_meta_raw <- shiny::reactive({
    shiny::req(input$meta_file)
    tryCatch(
      read.csv(input$meta_file$datapath, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  })

  # Render the whole sample-selector box (only when meta is loaded)
  output$sample_selector_box_ui <- shiny::renderUI({
    meta <- rv_meta_raw()
    if (is.null(meta) || !"Sample_ID" %in% colnames(meta)) return(NULL)

    all_cols    <- setdiff(colnames(meta), "Sample_ID")
    filter_cols <- all_cols[vapply(all_cols, function(col)
      length(unique(meta[[col]])) <= 50, logical(1))]

    shinydashboard::box(
      title       = shiny::tags$span(shiny::icon("users"), " Select Samples to Analyse"),
      width       = 12,
      status      = "warning",
      solidHeader = TRUE,
      collapsible = TRUE,

      shiny::p(
        "Use the dropdowns to filter the browser, then",
        shiny::strong("tick rows"), "and click",
        shiny::strong("Add Selected."),
        " In the right table, tick rows and click",
        shiny::strong("Remove Selected,"), "or use",
        shiny::strong("Clear All."),
        " Only samples in the Selected Samples table will be analysed.",
        class = "info-text"
      ),

      # ── LEFT PANEL: browser + pickers ─────────────────────────────────────
      shiny::fluidRow(
        shiny::column(
          width = 7,
          shiny::h5(shiny::icon("filter"), " Sample Browser",
                    style = "font-weight:700; margin-bottom:6px;"),

          # Per-column pickers
          if (length(filter_cols) > 0) {
            shiny::div(
              class = "filter-panel",
              shiny::tags$script(shiny::HTML(
                "$(document).on('changed.bs.select', '[id^=\"col_filter_\"]', function() {",
                "  var n = $(this).val() ? $(this).val().length : 0;",
                "  Shiny.setInputValue(this.id + '_n', n, {priority: 'event'});",
                "});"
              )),
              shiny::fluidRow(
                lapply(filter_cols, function(col) {
                  vals <- ifelse(is.na(unique(meta[[col]])), "NA",
                                 as.character(unique(meta[[col]])))
                  shiny::column(
                    width = max(2, min(4, floor(12 / length(filter_cols)))),
                    shinyWidgets::pickerInput(
                      inputId  = paste0("col_filter_", col),
                      label    = col,
                      choices  = vals,
                      selected = vals,
                      multiple = TRUE,
                      options  = shinyWidgets::pickerOptions(
                        actionsBox         = TRUE,
                        liveSearch         = TRUE,
                        selectedTextFormat = "count > 2",
                        size               = 8
                      )
                    )
                  )
                })
              )
            )
          },

          # Reset-filters + Add Selected buttons
          shiny::div(
            style = "margin-bottom:6px; display:flex; gap:6px; align-items:center;",
            shinyWidgets::actionBttn(
              "reset_filters", "Reset Filters",
              style = "bordered", color = "default", size = "sm",
              icon  = shiny::icon("undo")
            ),
            shinyWidgets::actionBttn(
              "add_selected_samples", "Add Selected \u2192",
              style = "bordered", color = "primary", size = "sm",
              icon  = shiny::icon("plus")
            ),
            shinyWidgets::actionBttn(
              "add_all_filtered", "Add All Filtered \u2192",
              style = "bordered", color = "warning", size = "sm",
              icon  = shiny::icon("angle-double-right")
            )
          ),

          # Browser table — multi-select rows
          DT::DTOutput("sample_browser_table", height = "300px")
        ),

        # ── RIGHT PANEL: selected samples ──────────────────────────────────
        shiny::column(
          width = 5,
          shiny::h5(shiny::icon("check-circle"), " Selected Samples",
                    style = "font-weight:700; margin-bottom:6px;"),

          shiny::fluidRow(
            shiny::column(6, shiny::uiOutput("sample_count_ui")),
            shiny::column(6,
              shiny::div(
                style = "text-align:right; padding-top:4px; display:flex; gap:6px; justify-content:flex-end;",
                shinyWidgets::actionBttn(
                  "remove_selected_samples", "Remove Selected",
                  style = "bordered", color = "warning", size = "sm",
                  icon  = shiny::icon("minus")
                ),
                shinyWidgets::actionBttn(
                  "clear_selected_samples", "Clear All",
                  style = "bordered", color = "danger", size = "sm",
                  icon  = shiny::icon("trash")
                )
              )
            )
          ),
          shiny::br(),

          # Selected samples table — multi-select rows
          DT::DTOutput("sample_selected_table", height = "300px"),

          shiny::br(),
          shiny::hr(style = "margin:8px 0;"),
          shiny::h6(shiny::icon("filter"), " Low-Expression Filter",
                    style = "font-weight:700; margin-bottom:4px; color:#555;"),
          shiny::p(
            "Applied to raw counts ",
            shiny::strong("after"), " sample selection, ",
            shiny::strong("before"), " log-normalisation.",
            " Genes not meeting both thresholds are dropped.",
            class = "info-text", style = "font-size:0.85em;"
          ),
          shiny::fluidRow(
            shiny::column(6,
              shiny::numericInput(
                inputId = "min_cpm_filter",
                label   = shiny::tags$span("Min CPM", shiny::tags$small(" (per sample)", style="color:#888;")),
                value   = 1,
                min     = 0,
                step    = 0.5
              )
            ),
            shiny::column(6,
              shiny::numericInput(
                inputId = "min_pct_samples",
                label   = shiny::tags$span("Min % samples passing", shiny::tags$small(" (0–100)", style="color:#888;")),
                value   = 30,
                min     = 0,
                max     = 100,
                step    = 5
              )
            )
          )
        )
      )
    )
  })

  # ==========================================================================
  # SAMPLE SELECTOR — two-table approach
  #
  #   picker_filter_ids  — rows currently visible in the browser (from pickers)
  #   selected_ids_rv    — reactiveVal: the user's confirmed selection
  #
  #   Browser table  : multi-row selectable
  #     "Add Selected"     → add DT-highlighted rows to selected_ids_rv
  #     "Add All Filtered" → add everything in picker_filter_ids
  #   Selected table : multi-row selectable
  #     "Remove Selected"  → drop DT-highlighted rows from selected_ids_rv
  #     "Clear All"        → selected_ids_rv = character(0)
  #   "Reset Filters"      → reset all pickers to full selection
  # ==========================================================================

  to_chr <- function(x) ifelse(is.na(x), "NA", as.character(x))

  selected_ids_rv <- shiny::reactiveVal(character(0))

  # On metadata load: auto-select all samples and reset pickers
  shiny::observeEvent(rv_meta_raw(), {
    meta <- rv_meta_raw()
    if (is.null(meta) || !"Sample_ID" %in% colnames(meta)) return()
    selected_ids_rv(meta$Sample_ID)
    all_cols    <- setdiff(colnames(meta), "Sample_ID")
    filter_cols <- all_cols[vapply(all_cols, function(col)
      length(unique(meta[[col]])) <= 50, logical(1))]
    for (col in filter_cols) {
      vals <- ifelse(is.na(unique(meta[[col]])), "NA", as.character(unique(meta[[col]])))
      shinyWidgets::updatePickerInput(session,
        inputId = paste0("col_filter_", col), selected = vals)
    }
  })

  # What the pickers currently show (browser filter)
  picker_filter_ids <- shiny::reactive({
    meta <- rv_meta_raw()
    shiny::req(meta, "Sample_ID" %in% colnames(meta))
    all_cols    <- setdiff(colnames(meta), "Sample_ID")
    filter_cols <- all_cols[vapply(all_cols, function(col)
      length(unique(meta[[col]])) <= 50, logical(1))]
    keep <- rep(TRUE, nrow(meta))
    for (col in filter_cols) {
      sv       <- input[[paste0("col_filter_", col)]]
      n_sel    <- input[[paste0("col_filter_", col, "_n")]]
      all_vals <- to_chr(unique(meta[[col]]))
      if (!is.null(n_sel) && n_sel == 0) { keep[] <- FALSE; break }
      if (is.null(sv))            next
      if (setequal(sv, all_vals)) next
      if (length(sv) == 0)       { keep[] <- FALSE; break }
      keep <- keep & (to_chr(meta[[col]]) %in% sv)
    }
    meta$Sample_ID[keep]
  })

  # Helper: get the Sample_IDs corresponding to DT row indices
  browser_ids_for_rows <- function(rows) {
    meta <- rv_meta_raw()
    if (is.null(meta) || length(rows) == 0) return(character(0))
    ids  <- picker_filter_ids()
    df   <- meta[meta$Sample_ID %in% ids, , drop = FALSE]
    df$Sample_ID[rows]
  }

  selected_ids_for_rows <- function(rows) {
    meta <- rv_meta_raw()
    if (is.null(meta) || length(rows) == 0) return(character(0))
    ids <- selected_ids_rv()
    df  <- meta[meta$Sample_ID %in% ids, , drop = FALSE]
    df$Sample_ID[rows]
  }

  # Add rows highlighted in browser table
  shiny::observeEvent(input$add_selected_samples, {
    rows     <- input$sample_browser_table_rows_selected
    new_ids  <- browser_ids_for_rows(rows)
    if (length(new_ids) == 0) {
      shiny::showNotification("No rows selected in the browser table.", type = "warning", duration = 3)
      return()
    }
    existing <- selected_ids_rv()
    added    <- new_ids[!new_ids %in% existing]
    combined <- unique(c(existing, new_ids))
    selected_ids_rv(combined)
    shiny::showNotification(
      sprintf("%d new sample(s) added. Total: %d.", length(added), length(combined)),
      type = "message", duration = 3
    )
  })

  # Add everything currently visible in the browser
  shiny::observeEvent(input$add_all_filtered, {
    new_ids  <- picker_filter_ids()
    existing <- selected_ids_rv()
    added    <- new_ids[!new_ids %in% existing]
    combined <- unique(c(existing, new_ids))
    selected_ids_rv(combined)
    shiny::showNotification(
      sprintf("%d new sample(s) added. Total: %d.", length(added), length(combined)),
      type = "message", duration = 3
    )
  })

  # Remove rows highlighted in selected table
  shiny::observeEvent(input$remove_selected_samples, {
    rows    <- input$sample_selected_table_rows_selected
    rm_ids  <- selected_ids_for_rows(rows)
    if (length(rm_ids) == 0) {
      shiny::showNotification("No rows selected in the selected-samples table.", type = "warning", duration = 3)
      return()
    }
    remaining <- selected_ids_rv()[!selected_ids_rv() %in% rm_ids]
    selected_ids_rv(remaining)
    shiny::showNotification(
      sprintf("%d sample(s) removed. Total: %d.", length(rm_ids), length(remaining)),
      type = "warning", duration = 3
    )
  })

  # Clear All
  shiny::observeEvent(input$clear_selected_samples, {
    selected_ids_rv(character(0))
  })

  # Reset Filters → all pickers back to full selection
  shiny::observeEvent(input$reset_filters, {
    meta <- rv_meta_raw()
    shiny::req(meta)
    all_cols    <- setdiff(colnames(meta), "Sample_ID")
    filter_cols <- all_cols[vapply(all_cols, function(col)
      length(unique(meta[[col]])) <= 50, logical(1))]
    for (col in filter_cols) {
      vals <- ifelse(is.na(unique(meta[[col]])), "NA", as.character(unique(meta[[col]])))
      shinyWidgets::updatePickerInput(session,
        inputId = paste0("col_filter_", col), selected = vals)
    }
  })

  # sel_ids_rv = confirmed selection (pipeline-facing)
  sel_ids_rv <- shiny::reactive({ selected_ids_rv() })

  # ── Browser table (left) — multi-select rows ─────────────────────────────
  output$sample_browser_table <- DT::renderDT({
    meta <- rv_meta_raw()
    shiny::req(meta, "Sample_ID" %in% colnames(meta))
    ids <- picker_filter_ids()
    df  <- meta[meta$Sample_ID %in% ids, , drop = FALSE]
    DT::datatable(
      df,
      selection = list(mode = "multiple", target = "row"),
      rownames  = FALSE,
      class     = "stripe hover compact",
      options   = list(
        pageLength = 10,
        scrollX    = TRUE,
        dom        = "frtip",
        scrollY    = "240px",
        scroller   = TRUE
      )
    )
  })

  # ── Selected samples table (right) — multi-select rows ───────────────────
  output$sample_selected_table <- DT::renderDT({
    meta <- rv_meta_raw()
    shiny::req(meta, "Sample_ID" %in% colnames(meta))
    ids <- selected_ids_rv()
    df  <- meta[meta$Sample_ID %in% ids, , drop = FALSE]
    DT::datatable(
      df,
      selection = list(mode = "multiple", target = "row"),
      rownames  = FALSE,
      class     = "stripe hover compact",
      options   = list(
        pageLength = 10,
        scrollX    = TRUE,
        dom        = "frtip",
        scrollY    = "240px",
        scroller   = TRUE
      )
    )
  })

  # ── Sample count badge ───────────────────────────────────────────────────
  output$sample_count_ui <- shiny::renderUI({
    meta  <- rv_meta_raw()
    total <- if (!is.null(meta)) nrow(meta) else 0
    n_sel <- length(selected_ids_rv())
    colour <- if (n_sel == 0) "#c0392b" else if (n_sel < total) "#e67e22" else "#27ae60"
    shiny::tags$p(shiny::tags$span(
      sprintf("%d / %d selected", n_sel, total),
      style = sprintf("font-size:1.05em; font-weight:700; color:%s;", colour)
    ))
  })

  # ── Warning shown inline when 'already CPM' is ticked ───────────────────
  output$cpm_skip_warning_ui <- shiny::renderUI({
    if (!isTRUE(input$counts_already_cpm)) return(NULL)
    shiny::p(
      shiny::icon("info-circle"),
      " Library-size normalisation will be skipped. Your values will be treated as CPM ",
      "and passed directly to log\u2082(CPM + 1). The low-expression filter will still ",
      "be applied using the CPM values as-is.",
      style = "font-size:0.82em; color:#7d6608; margin:4px 0 0 24px;"
    )
  })

  # ── Scale data checkbox — disabled + unchecked when 'already CPM' is ticked
  output$scale_counts_ui <- shiny::renderUI({
    already_cpm <- isTRUE(input$counts_already_cpm)
    # When CPM box is checked, force scale_counts_flag to FALSE so the pipeline
    # never sees a stale TRUE from a previously checked state.
    cb_value <- if (already_cpm) FALSE else isTRUE(input$scale_counts_flag)
    shiny::tagList(
      shiny::div(
        style = if (already_cpm) "opacity:0.45; pointer-events:none;" else "",
        shiny::checkboxInput(
          inputId = "scale_counts_flag",
          label   = shiny::tags$span(
            shiny::icon("compress-alt", style = "color:#e67e22;"),
            " Scale data before normalisation",
            style = "font-size:0.9em;"
          ),
          value = cb_value
        )
      ),
      shiny::tags$p(
        style = paste0(
          "font-size:0.82em; margin:2px 0 0 22px; ",
          if (already_cpm) "color:#b0a070; font-style:italic;" else "color:#7f6000;"
        ),
        if (already_cpm)
          "Not available — data is already CPM-normalised."
        else
          shiny::tagList(
            "Tick this if your counts come from ",
            shiny::tags$strong("recount3"),
            " or similar pipelines that store raw base-level coverage sums.",
            " Counts will be divided by (column sums \u00f7 10\u2076) to convert",
            " coverage integers to read-count equivalents before the CPM step."
          )
      )
    )
  })

  # ── Pipeline-facing accessor ─────────────────────────────────────────────
  selected_sample_ids <- shiny::reactive({ sel_ids_rv() })

  # Run status display
  run_status_content <- shiny::reactive({
    tryCatch({
    if (!rv$run_success) return(NULL)
    n_avail <- rv$ratio_res$n_avail
    n_total <- rv$ratio_res$n_total
    shiny::tagList(
      shiny::tags$span(class = "ratio-badge",
        sprintf("\u2713 %d / %d ECM ratios available in your data", n_avail, n_total)
      )
    )
    }, error = function(e) { cat("\n[CRASH] run_status_content: ", conditionMessage(e), "\n"); traceback(3); NULL })
  })

  output$run_status_ui  <- shiny::renderUI({ run_status_content() })
  output$run_status_ui2 <- shiny::renderUI({ run_status_content() })

  output$ratio_info_ui <- shiny::renderUI({
    if (!rv$run_success) return(NULL)
    avail <- rv$ratio_res$avail_features
    ad_n  <- sum(avail$Panel %in% c("A/D", "AD"))
    sig_n <- sum(avail$Panel %in% c("SIG", "Signalling", "SIGNALLING"))
    shiny::p(sprintf("A/D ratios: %d   |   SIG ratios: %d", ad_n, sig_n),
             class = "info-text")
  })

  # ================================================================
  # Main run pipeline
  # ================================================================
  shiny::observeEvent(input$run_btn, {

    if (is.null(input$counts_file)) {
      shiny::showNotification(
        shiny::tags$span(shiny::icon("exclamation-triangle"), " Please upload a count matrix before running."),
        type = "error", duration = 5
      )
      return()
    }
    if (is.null(input$meta_file)) {
      shiny::showNotification(
        shiny::tags$span(shiny::icon("exclamation-triangle"), " Please upload sample metadata before running."),
        type = "error", duration = 5
      )
      return()
    }

    # Resolve which training bundle to use
    active_training <- isolate(resolve_active_training())
    if (is.null(active_training)) {
      shiny::showNotification("Training data not available.", type = "error")
      return()
    }

    # ---- Waiter: full-screen progress ----
    w <- waiter::Waiter$new(
      html  = waiter::spin_flower(),
      color = waiter::transparent(0.7)
    )
    w$show()
    on.exit(w$hide(), add = TRUE)

    # ── DEBUG: outer tryCatch to catch any error that slips past inner handlers ──
    tryCatch(

    isolate(withCallingHandlers({

      message("[DEBUG] run_btn fired — entering pipeline")

      # Step 1 – Read user counts
      shinyWidgets::updateProgressBar(session, "pb", value = 5,
                                      title = "Reading uploaded count matrix…")
      counts_raw <- tryCatch({
        # Read WITHOUT row.names to avoid "duplicate row.names" error on parse
        df <- read.csv(input$counts_file$datapath,
                       header           = TRUE,
                       check.names      = FALSE,
                       stringsAsFactors = FALSE)
        # First column is gene symbols
        gene_col  <- df[[1]]
        df        <- df[, -1, drop = FALSE]
        mat       <- as.matrix(df)
        storage.mode(mat) <- "numeric"
        rownames(mat) <- gene_col
        mat
      }, error = function(e) {
        shiny::showNotification(paste("Count matrix read error:", e$message), type = "error")
        NULL
      })
      if (is.null(counts_raw)) return()
      message("[DEBUG] Step 1 OK: counts read, ", nrow(counts_raw), " genes")

      # Step 1b – Filter to only genes relevant for ECM ratios BEFORE touching duplicates.
      #           This keeps the matrix small and means duplicates only matter for
      #           the handful of genes actually used in the 100 ratios.
      ratio_genes <- unique(c(active_training$features$Gene_A, active_training$features$Gene_B))
      keep_rows   <- rownames(counts_raw) %in% ratio_genes
      counts_filt <- counts_raw[keep_rows, , drop = FALSE]

      # Step 1c – Resolve duplicate gene symbols (among ratio genes only)
      dup_genes <- unique(rownames(counts_filt)[duplicated(rownames(counts_filt))])
      if (length(dup_genes) > 0) {
        strategy <- input$dup_gene_strategy
        msg <- sprintf(
          "%d duplicate gene symbol(s) found in ratio genes (%s). Strategy: %s.",
          length(dup_genes),
          paste(head(dup_genes, 5), collapse = ", "),
          switch(strategy,
            mean   = "averaging counts",
            max    = "keeping highest-expressed",
            remove = "removing all duplicates"
          )
        )
        shiny::showNotification(msg, type = "warning", duration = 10)

        if (strategy == "remove") {
          counts_filt <- counts_filt[!rownames(counts_filt) %in% dup_genes, , drop = FALSE]
        } else {
          all_genes  <- rownames(counts_filt)
          uniq_genes <- unique(all_genes)
          resolved   <- lapply(uniq_genes, function(g) {
            rows <- counts_filt[all_genes == g, , drop = FALSE]
            if (nrow(rows) == 1) return(rows)
            if (strategy == "mean") {
              matrix(colMeans(rows, na.rm = TRUE), nrow = 1, dimnames = list(g, colnames(rows)))
            } else {  # max: row with highest mean expression
              matrix(rows[which.max(rowMeans(rows, na.rm = TRUE)), ], nrow = 1,
                     dimnames = list(g, colnames(rows)))
            }
          })
          counts_filt <- do.call(rbind, resolved)
        }
      }

      # Step 2 – Read user metadata
      shinyWidgets::updateProgressBar(session, "pb", value = 15,
                                      title = "Reading metadata…")
      user_meta <- tryCatch(
        read.csv(input$meta_file$datapath, stringsAsFactors = FALSE),
        error = function(e) {
          shiny::showNotification(paste("Metadata read error:", e$message), type = "error")
          NULL
        }
      )
      if (is.null(user_meta)) return()

      # Validate metadata
      if (!"Sample_ID" %in% colnames(user_meta)) {
        shiny::showNotification("Metadata must have a 'Sample_ID' column.", type = "error")
        return()
      }

      # ---- Apply sample selection from the interactive table ----
      sel_ids <- isolate(selected_sample_ids())
      if (length(sel_ids) == 0) {
        shiny::showNotification(
          "No samples selected. Please select at least one sample in the sample table.",
          type = "error"
        )
        return()
      }
      # Restrict metadata and counts to the user-selected samples
      user_meta   <- user_meta[user_meta$Sample_ID %in% sel_ids, , drop = FALSE]
      counts_filt <- counts_filt[, colnames(counts_filt) %in% sel_ids, drop = FALSE]

      # Align counts columns to metadata order
      shared_samples <- intersect(colnames(counts_filt), user_meta$Sample_ID)
      if (length(shared_samples) == 0) {
        shiny::showNotification(
          "No matching Sample_IDs between count matrix columns and selected metadata samples.",
          type = "error"
        )
        return()
      }
      counts_filt <- counts_filt[, shared_samples, drop = FALSE]
      user_meta   <- user_meta[match(shared_samples, user_meta$Sample_ID), , drop = FALSE]

      # ── Scale raw counts (recount3 / coverage-scaled data) ───────────────
      # Divides each sample's counts by (colSums / 1e6) to convert
      # base-level coverage integers to read-count equivalents.
      # Must run BEFORE the low-expression CPM filter and before normalisation.
      if (isTRUE(input$scale_counts_flag)) {
        counts_filt <- transform_counts(counts_filt)
        shiny::showNotification(
          shiny::tags$span(
            shiny::icon("compress-alt"),
            " Counts scaled (\u00f7 colSums/10\u2076) for recount3 / coverage-count compatibility."
          ),
          type = "message", duration = 7
        )
      }

      # ── Low-expression filter ────────────────────────────────────────────
      # Applied to raw counts AFTER sample selection, BEFORE log-normalisation.
      # CPM is computed from the selected samples' raw counts.
      # The filter is based only on RNA-seq columns (microarray values are not
      # raw counts and should not be tested against a CPM threshold).
      shinyWidgets::updateProgressBar(session, "pb", value = 20,
                                      title = "Applying low-expression filter…")
      {
        min_cpm         <- max(0, as.numeric(input$min_cpm_filter),  na.rm = TRUE)
        min_pct         <- max(0, min(100, as.numeric(input$min_pct_samples), na.rm = TRUE))
        already_cpm     <- isTRUE(input$counts_already_cpm)

        # All columns are RNA-seq (microarray removed in V19)
        filt_cols <- seq_len(ncol(counts_filt))

        n_genes_before <- nrow(counts_filt)

        if (length(filt_cols) > 0 && min_cpm > 0) {
          rna_block <- counts_filt[, filt_cols, drop = FALSE]

          if (already_cpm) {
            # Data already in CPM — use values directly for the filter
            cpm_block <- rna_block
          } else {
            lib_sizes <- colSums(rna_block)
            lib_sizes[lib_sizes == 0] <- 1
            cpm_block <- sweep(rna_block, 2, lib_sizes / 1e6, "/")
          }

          min_samples_needed <- ceiling(min_pct / 100 * length(filt_cols))
          keep_genes <- rowSums(cpm_block >= min_cpm) >= min_samples_needed
          counts_filt <- counts_filt[keep_genes, , drop = FALSE]

          n_removed <- n_genes_before - nrow(counts_filt)
          if (n_removed > 0) {
            shiny::showNotification(
              shiny::tags$span(
                shiny::icon("filter"),
                sprintf(
                  " Low-expression filter: removed %d gene(s), kept %d (≥ %.1f CPM in ≥ %.0f%% of RNA-seq samples).",
                  n_removed, nrow(counts_filt), min_cpm, min_pct
                )
              ),
              type = "message", duration = 8
            )
          }
        }
      }

      # Step 4 – Normalise (use the filtered, deduplicated matrix)
      shinyWidgets::updateProgressBar(session, "pb", value = 25,
                                      title = "Normalising counts…")

      already_cpm_flag <- isTRUE(input$counts_already_cpm)
      if (already_cpm_flag) {
        shiny::showNotification(
          shiny::tags$span(
            shiny::icon("exclamation-circle", style = "color:#e67e22;"),
            " Library-size normalisation skipped (data treated as pre-normalised CPM). ",
            "Applying log\u2082(CPM + 1) only."
          ),
          type = "warning", duration = 7
        )
      }

      norm_res <- tryCatch(
        normalize_counts(
          counts_matrix   = counts_filt,
          skip_rnaseq_cpm = already_cpm_flag
        ),
        error = function(e) {
          shiny::showNotification(paste("Normalisation error:", e$message), type = "error")
          NULL
        }
      )
      if (is.null(norm_res)) return()
      norm_mat <- norm_res$log_cpm
      message("[DEBUG] Step 4 OK: normalised, ", nrow(norm_mat), " genes x ", ncol(norm_mat), " samples")

      # Step 5 – Batch correction on user data (if columns selected)
      batch_cols_sel <- input$batch_cols
      if (!is.null(batch_cols_sel) && length(batch_cols_sel) > 0) {

        shinyWidgets::updateProgressBar(session, "pb", value = 35,
                                        title = "Applying batch correction (ComBat)…")
        # For multi-column batch: combine into single compound batch variable
        if (length(batch_cols_sel) == 1) {
          bc_col <- batch_cols_sel
        } else {
          compound_batch <- apply(user_meta[, batch_cols_sel, drop = FALSE], 1, paste, collapse = "_")
          user_meta$CompoundBatch__ <- compound_batch
          bc_col <- "CompoundBatch__"
        }

        norm_mat <- tryCatch({
          bc_res <- apply_batch_correction(
            log_cpm_mat  = norm_mat,
            metadata_df  = user_meta,
            batch_col    = bc_col,
            protect_cols = if (length(input$protect_cols) > 0) input$protect_cols else NULL
          )
          bc_res$corrected_mat
        }, error = function(e) {
          shiny::showNotification(
            paste("Batch correction failed (skipped):", e$message),
            type = "warning"
          )
          norm_mat
        })
      }

      rv$norm_counts <- norm_mat
      rv$user_meta   <- user_meta

      # ── Resolve prediction mode ───────────────────────────────────────────
      pred_mode <- input$prediction_mode %||% "ratios"

      # Step 6 – Detect available features (ratios, ECM genes, or all genes)
      if (pred_mode == "ratios") {
        # ── Default: ratio mode ───────────────────────────────────────────
        shinyWidgets::updateProgressBar(session, "pb", value = 50,
                                        title = "Detecting available ECM ratios…")
        ratio_res <- tryCatch(
          compute_ratios(norm_mat, active_training$features),
          error = function(e) {
            shiny::showNotification(paste("Ratio detection error:", e$message), type = "error")
            NULL
          }
        )
        if (is.null(ratio_res)) return()
        rv$ratio_res <- ratio_res
        message("[DEBUG] Step 6 OK: ratios detected, ", ratio_res$n_avail, "/", ratio_res$n_total, " available")

        shinyWidgets::updateProgressBar(session, "pb", value = 65,
                                        title = "Retraining ECMExplainER on available ratios…")

        # Filter training set to selected conditions
        selected_conds   <- input$selected_conditions
        train_meta_filt  <- active_training$train_meta[
          active_training$train_meta$Phenotype %in% selected_conds, , drop = FALSE
        ]
        train_expr_filt  <- active_training$train_expr[, train_meta_filt$Sample_ID, drop = FALSE]

        train_ratio_res <- tryCatch(
          compute_ratios(train_expr_filt, ratio_res$avail_features),
          error = function(e) {
            shiny::showNotification(paste("Training ratio computation error:", e$message),
                                    type = "error")
            NULL
          }
        )
        if (is.null(train_ratio_res)) return()

        train_obj <- tryCatch(
          train_ecmexplainer(
            ratio_mat      = train_ratio_res$ratio_mat,
            avail_features = train_ratio_res$avail_features,
            phenotype_vec  = factor(
              train_meta_filt$Phenotype[
                match(colnames(train_ratio_res$ratio_mat), train_meta_filt$Sample_ID)
              ]
            ),
            metadata_df    = train_meta_filt,
            batch_col      = if ("Study_ID" %in% colnames(train_meta_filt)) "Study_ID" else NULL,
            protect_cols   = "Phenotype"
          ),
          error = function(e) {
            shiny::showNotification(paste("Training error:", e$message), type = "error")
            NULL
          }
        )
        if (is.null(train_obj)) return()
        rv$train_obj       <- train_obj
        rv$train_meta_filt <- train_meta_filt
        message("[DEBUG] Step 7 OK: trained (ratio mode), conditions=",
                paste(levels(train_obj$phenotype_vec), collapse = ","),
                " n_pcs_ad=", train_obj$n_pcs_ad, " n_pcs_sig=", train_obj$n_pcs_sig)

        shinyWidgets::updateProgressBar(session, "pb", value = 82,
                                        title = "Projecting and predicting new samples…")
        pred_obj <- tryCatch(
          predict_ecmexplainer(
            test_expr_mat  = norm_mat,
            train_obj      = train_obj,
            avail_features = ratio_res$avail_features
          ),
          error = function(e) {
            shiny::showNotification(paste("Prediction error:", e$message), type = "error")
            NULL
          }
        )
        if (is.null(pred_obj)) return()
        rv$pred_obj         <- pred_obj
        rv$pred_obj_aligned <- NULL
        message("[DEBUG] Step 8 OK: predicted (ratio mode)")

      } else if (pred_mode == "binary_tsp") {
        # ── Binary k-TSP mode ─────────────────────────────────────────────
        shinyWidgets::updateProgressBar(session, "pb", value = 50,
                                        title = "Encoding binary TSP features…")

        # features table must have Gene_A / Gene_B / Ratio_ID columns
        feat_tbl <- active_training$features
        if (is.null(feat_tbl) || !all(c("Gene_A","Gene_B") %in% colnames(feat_tbl))) {
          shiny::showNotification(
            "Binary TSP requires a ratio-pair feature table with Gene_A / Gene_B columns.",
            type = "error"
          )
          return()
        }

        # Filter training meta to selected conditions
        selected_conds   <- input$selected_conditions
        train_meta_filt  <- active_training$train_meta[
          active_training$train_meta$Phenotype %in% selected_conds, , drop = FALSE
        ]
        train_expr_filt  <- active_training$train_expr[, train_meta_filt$Sample_ID, drop = FALSE]

        # Convert to log2(CPM+1) if the matrix contains raw counts.
        # Heuristic: if the 99th-percentile value is > 50, treat as raw counts.
        train_logcpm_filt <- if (quantile(train_expr_filt, 0.99) > 50) {
          message("[BinaryTSP] Detected raw counts — converting to log2(CPM+1).")
          to_logcpm_ktsp(train_expr_filt)
        } else {
          message("[BinaryTSP] Training data appears already normalised — using as-is.")
          train_expr_filt
        }

        shinyWidgets::updateProgressBar(session, "pb", value = 65,
                                        title = "Training k-TSP classifier…")
        train_obj <- tryCatch(
          train_ecmexplainer(
            ratio_mat      = NULL,
            avail_features = feat_tbl,
            phenotype_vec  = factor(
              train_meta_filt$Phenotype[
                match(colnames(train_expr_filt), train_meta_filt$Sample_ID)
              ]
            ),
            logcpm_mat     = train_logcpm_filt,
            feature_method = "binary_tsp",
            metadata_df    = train_meta_filt,
            batch_col      = NULL,   # binary encoding is batch-robust
            protect_cols   = "Phenotype"
          ),
          error = function(e) {
            shiny::showNotification(paste("k-TSP training error:", e$message), type = "error")
            NULL
          }
        )
        if (is.null(train_obj)) return()

        rv$train_obj       <- train_obj
        rv$train_meta_filt <- train_meta_filt

        # Build a synthetic ratio_res for downstream UI compatibility
        rv$ratio_res <- list(
          ratio_mat      = NULL,
          avail_features = feat_tbl,
          n_avail        = nrow(feat_tbl),
          n_total        = nrow(feat_tbl)
        )

        shinyWidgets::updateProgressBar(session, "pb", value = 82,
                                        title = "Predicting with k-TSP classifier…")
        pred_obj <- tryCatch(
          predict_ecmexplainer(
            test_expr_mat  = norm_mat,
            train_obj      = train_obj,
            avail_features = feat_tbl
          ),
          error = function(e) {
            shiny::showNotification(paste("k-TSP prediction error:", e$message), type = "error")
            NULL
          }
        )
        if (is.null(pred_obj)) return()
        rv$pred_obj         <- pred_obj
        rv$pred_obj_aligned <- NULL
        message("[DEBUG] Step 8 OK: predicted (binary_tsp mode)")

      } else {
        # ── Gene mode (ecm_genes or all_genes) ───────────────────────────
        shinyWidgets::updateProgressBar(session, "pb", value = 50,
                                        title = sprintf("Preparing gene-level features (%s)…", pred_mode))

        # Filter genes to ECM subset or use all genes
        # ecm_master may only be present in the built-in training bundle
        ecm_dict_genes <- if (!is.null(active_training$ecm_master) &&
                               "Gene" %in% colnames(active_training$ecm_master)) {
          unique(active_training$ecm_master$Gene)
        } else if (!is.null(training$ecm_master) && "Gene" %in% colnames(training$ecm_master)) {
          unique(training$ecm_master$Gene)
        } else {
          character(0)
        }
        if (pred_mode == "ecm_genes") {
          if (length(ecm_dict_genes) == 0) {
            shiny::showNotification(
              "ECM gene dictionary not available — falling back to all-genes mode.",
              type = "warning"
            )
            user_genes  <- rownames(norm_mat)
            train_genes <- rownames(active_training$train_expr)
          } else {
            user_genes  <- intersect(rownames(norm_mat),                    ecm_dict_genes)
            train_genes <- intersect(rownames(active_training$train_expr),  ecm_dict_genes)
          }
        } else {   # all_genes
          user_genes  <- rownames(norm_mat)
          train_genes <- rownames(active_training$train_expr)
        }

        # Keep only genes present in BOTH user and training data
        shared_genes <- intersect(user_genes, train_genes)
        if (length(shared_genes) == 0) {
          shiny::showNotification(
            sprintf("No shared genes between user data and training data in '%s' mode.", pred_mode),
            type = "error"
          )
          return()
        }
        shiny::showNotification(
          sprintf("[%s] Using %d shared genes for PCA + LDA.", pred_mode, length(shared_genes)),
          type = "message", duration = 6
        )

        # Build fake avail_features where Gene_B = Gene_A (ratio = 0, but we will
        # bypass ratio computation via precomputed_ratio_mat in predict_ecmexplainer)
        gene_features <- data.frame(
          Panel  = "AD",
          Gene_A = shared_genes,
          Gene_B = shared_genes,
          stringsAsFactors = FALSE
        )
        rownames(gene_features) <- NULL

        # Gene expression sub-matrices (genes × samples)
        user_gene_mat  <- norm_mat[shared_genes,                    , drop = FALSE]
        train_gene_mat <- active_training$train_expr[shared_genes, , drop = FALSE]

        # Build ratio_mat for train_ecmexplainer(): here "ratio_mat" IS the gene matrix
        # (features × samples). Row names = gene names (no " / " separator).
        train_ratio_mat_gene <- train_gene_mat   # genes x training_samples

        # Filter training metadata to selected conditions
        selected_conds   <- input$selected_conditions
        train_meta_filt  <- active_training$train_meta[
          active_training$train_meta$Phenotype %in% selected_conds, , drop = FALSE
        ]
        train_ratio_mat_filt <- train_ratio_mat_gene[, train_meta_filt$Sample_ID, drop = FALSE]

        shinyWidgets::updateProgressBar(session, "pb", value = 65,
                                        title = sprintf("Training PCA + LDA on %d genes…", length(shared_genes)))

        train_obj <- tryCatch(
          train_ecmexplainer(
            ratio_mat      = train_ratio_mat_filt,
            avail_features = gene_features,
            phenotype_vec  = factor(
              train_meta_filt$Phenotype[
                match(colnames(train_ratio_mat_filt), train_meta_filt$Sample_ID)
              ]
            ),
            metadata_df    = train_meta_filt,
            batch_col      = if ("Study_ID" %in% colnames(train_meta_filt)) "Study_ID" else NULL,
            protect_cols   = "Phenotype"
          ),
          error = function(e) {
            shiny::showNotification(paste("Training error (gene mode):", e$message), type = "error")
            NULL
          }
        )
        if (is.null(train_obj)) return()

        # Build a synthetic ratio_res so the rest of the app has a valid rv$ratio_res
        gene_ratio_res <- list(
          ratio_mat     = user_gene_mat,
          avail_features = gene_features,
          n_avail       = nrow(gene_features),
          n_total       = nrow(gene_features)
        )
        rv$ratio_res       <- gene_ratio_res
        rv$train_obj       <- train_obj
        rv$train_meta_filt <- train_meta_filt
        message("[DEBUG] Step 7 OK: trained (", pred_mode, " mode), ",
                length(shared_genes), " genes, conditions=",
                paste(levels(train_obj$phenotype_vec), collapse = ","))

        shinyWidgets::updateProgressBar(session, "pb", value = 82,
                                        title = "Projecting and predicting new samples…")

        pred_obj <- tryCatch(
          predict_ecmexplainer(
            test_expr_mat         = norm_mat,
            train_obj             = train_obj,
            avail_features        = gene_features,
            precomputed_ratio_mat = user_gene_mat
          ),
          error = function(e) {
            shiny::showNotification(paste("Prediction error (gene mode):", e$message), type = "error")
            NULL
          }
        )
        if (is.null(pred_obj)) return()
        rv$pred_obj         <- pred_obj
        rv$pred_obj_aligned <- NULL
        message("[DEBUG] Step 8 OK: predicted (", pred_mode, " mode)")
      }  # end gene mode branch

      # ── Mode B: Anchored projection (multi-pair centroid shift in PC space) ─
      if (!is.null(input$test_projection_mode) &&
          input$test_projection_mode == "anchored") {

        message("[DEBUG] Mode B entered, anchor_col='", input$test_anchor_col, "'")

        anchor_col <- input$test_anchor_col
        ok_col     <- !is.null(anchor_col) && anchor_col != "" &&
                       anchor_col %in% colnames(user_meta)

        if (!ok_col) {
          shiny::showNotification(
            "Mode B: anchor column not specified or not found in metadata — skipping alignment.",
            type = "warning"
          )
        } else {
          # Collect valid pairs directly from input (bypass reactive rv_anchor_pairs
          # to avoid any reactive-state issues inside a non-reactive context).
          # Scan up to 20 possible pair slots (same upper bound as the static remove observers).
          valid_pairs <- list()
          for (i in seq_len(20L)) {
            p_vals   <- input[[paste0("anchor_pair_vals_",   i)]]
            p_target <- input[[paste0("anchor_pair_target_", i)]]
            # Stop scanning once both inputs are NULL (no more pair slots rendered)
            if (is.null(p_vals) && is.null(p_target)) break
            if (!is.null(p_vals) && length(p_vals) > 0 &&
                !is.null(p_target) && nzchar(as.character(p_target)[1L])) {
              valid_pairs[[length(valid_pairs) + 1L]] <- list(
                vals   = as.character(p_vals),
                target = as.character(p_target)[1L]
              )
            }
          }

          if (length(valid_pairs) == 0L) {
            shiny::showNotification(
              "Mode B: no anchor pairs configured with values and a target — skipping alignment.",
              type = "warning"
            )
          } else {
            message("[DEBUG] Mode B: ", length(valid_pairs), " valid pair(s) collected")
            for (dbg_i in seq_along(valid_pairs))
              message("[DEBUG]   pair ", dbg_i, ": vals=",
                      paste(valid_pairs[[dbg_i]]$vals, collapse=","),
                      " -> target='", valid_pairs[[dbg_i]]$target, "'")
            tryCatch({
              val_pcs     <- pred_obj$val_pcs       # rows = test samples
              train_pc_df <- train_obj$train_pcs    # rows = training samples
              train_pheno <- as.character(train_obj$phenotype_vec)
              message("[DEBUG] Mode B tryCatch: val_pcs class=", class(val_pcs)[1],
                      " dim=", nrow(val_pcs), "x", ncol(val_pcs),
                      " train_pc_df class=", class(train_pc_df)[1],
                      " dim=", nrow(train_pc_df), "x", ncol(train_pc_df),
                      " n_pheno=", length(train_pheno))

              # Compute per-pair shift weighted by number of reference samples;
              # final shift = weighted average across all valid pairs.
              n_dims        <- ncol(val_pcs)
              total_weight  <- 0L
              shift_accum   <- numeric(n_dims)
              pair_summaries <- character(0)

              for (pair in valid_pairs) {
                message("[DEBUG] pair loop: vals=", paste(pair$vals, collapse=","),
                        " target='", pair$target, "'")
                ref_ids <- user_meta$Sample_ID[
                  as.character(user_meta[[anchor_col]]) %in% pair$vals
                ]
                message("[DEBUG]   ref_ids n=", length(ref_ids))
                if (length(ref_ids) == 0L) {
                  shiny::showNotification(
                    sprintf("Mode B pair: no test samples matched '%s' — this pair skipped.",
                            paste(pair$vals, collapse = ", ")),
                    type = "warning"
                  )
                  next
                }

                atlas_rows <- train_pheno == pair$target
                message("[DEBUG]   atlas_rows sum=", sum(atlas_rows))
                if (sum(atlas_rows) == 0L) {
                  shiny::showNotification(
                    sprintf("Mode B pair: no training samples with phenotype '%s' — this pair skipped.",
                            pair$target),
                    type = "warning"
                  )
                  next
                }

                message("[DEBUG]   computing atlas_centroid")
                atlas_centroid <- colMeans(as.matrix(train_pc_df[atlas_rows,  , drop = FALSE]))
                message("[DEBUG]   computing user_centroid, user_ref_rows=",
                        sum(rownames(val_pcs) %in% ref_ids))
                user_ref_rows  <- rownames(val_pcs) %in% ref_ids
                user_centroid  <- colMeans(as.matrix(val_pcs[user_ref_rows, , drop = FALSE]))
                message("[DEBUG]   pair_shift computed OK")
                pair_shift     <- atlas_centroid - user_centroid

                pair_w        <- length(ref_ids)
                shift_accum   <- shift_accum + pair_shift * pair_w
                total_weight  <- total_weight + pair_w

                pair_summaries <- c(pair_summaries,
                  sprintf("%s \u2192 atlas '%s' (n=%d)",
                          paste(pair$vals, collapse = "/"), pair$target, pair_w))
              }
              message("[DEBUG] pair loop complete, total_weight=", total_weight)

              if (total_weight == 0L) {
                shiny::showNotification(
                  "Mode B: all anchor pairs had zero matched samples — skipping alignment.",
                  type = "warning"
                )
              } else {
                # Weighted-average shift
                shift <- shift_accum / total_weight
                message("[DEBUG] post-loop: shift length=", length(shift),
                        " class=", class(shift)[1])

                # Apply to ALL test samples
                val_pcs_aligned <- sweep(val_pcs, 2, shift, "+")
                message("[DEBUG] sweep OK: val_pcs_aligned class=", class(val_pcs_aligned)[1],
                        " dim=", nrow(val_pcs_aligned), "x", ncol(val_pcs_aligned))

                # Re-run LDA on aligned PC scores
                message("[DEBUG] calling predict(lda_model) ...")
                preds_aligned <- predict(train_obj$lda_model, val_pcs_aligned)
                message("[DEBUG] predict OK: class=", class(preds_aligned)[1],
                        " names=", paste(names(preds_aligned), collapse=","))

                message("[DEBUG] extracting posterior ...")
                post_aligned <- as.data.frame(round(preds_aligned$posterior * 100, 2))
                message("[DEBUG] posterior OK: dim=", nrow(post_aligned), "x", ncol(post_aligned))

                message("[DEBUG] building final_rep_aligned ...")
                final_rep_aligned <- data.frame(
                  Sample_ID           = rownames(val_pcs_aligned),
                  Predicted_Phenotype = as.character(preds_aligned$class),
                  stringsAsFactors    = FALSE
                )
                final_rep_aligned <- cbind(final_rep_aligned, post_aligned)
                message("[DEBUG] final_rep_aligned OK: dim=", nrow(final_rep_aligned), "x", ncol(final_rep_aligned))

                message("[DEBUG] extracting preds_aligned$x: class=", class(preds_aligned$x)[1],
                        " dim=", paste(dim(preds_aligned$x), collapse="x"))
                aligned_pred_obj <- list(
                  final_report   = final_rep_aligned,
                  val_pcs        = as.data.frame(val_pcs_aligned),
                  val_lda_coords = preds_aligned$x,
                  anchor_col     = anchor_col,
                  anchor_pairs   = valid_pairs,
                  shift_vector   = shift
                )
                message("[DEBUG] aligned_pred_obj built OK")

                # Overwrite rv$pred_obj so ALL downstream plots / tables use the
                # aligned result automatically.  Also store the aligned version
                # separately in rv$pred_obj_aligned for any future Mode-B-aware code.
                message("[DEBUG] assigning rv$pred_obj ...")
                rv$pred_obj         <- aligned_pred_obj
                message("[DEBUG] rv$pred_obj assigned OK")
                rv$pred_obj_aligned <- aligned_pred_obj
                message("[DEBUG] rv$pred_obj_aligned assigned OK")

                message("[DEBUG] building notification string ...")
                notif_msg <- paste0(
                  "Mode B alignment applied (", length(pair_summaries), " pair(s)): ",
                  paste(pair_summaries, collapse = "; "), "."
                )
                message("[DEBUG] notification string OK: ", substr(notif_msg, 1, 60))
                shiny::showNotification(notif_msg, type = "message", duration = 10)
                message("[DEBUG] showNotification OK")
              }
            }, error = function(e) {
              emsg <- if (inherits(e, "condition")) conditionMessage(e) else paste(as.character(e), collapse=" ")
              cat("\n========== MODE B ERROR ==========\n")
              cat("Message: ", emsg, "\n")
              cat("Call:    ", deparse(conditionCall(e)), "\n")
              cat("==================================\n\n")
              shiny::showNotification(
                paste("Mode B alignment failed (unanchored result used):", emsg),
                type = "warning"
              )
            })
          }
        }
      }

      shinyWidgets::updateProgressBar(session, "pb", value = 100,
                                      title = "Done!")
      message("[DEBUG] updateProgressBar OK")
      rv$run_success <- TRUE
      message("[DEBUG] rv$run_success set TRUE")
      rv$wizard_step <- 2L   # → Training PCA step
      message("[DEBUG] rv$wizard_step set 2L")
      shiny::showNotification("ECMExplainER run complete!", type = "message", duration = 5)
      message("[DEBUG] final showNotification OK — observer DONE")

    }, warning = function(w) {
      shiny::showNotification(paste("Warning:", conditionMessage(w)),
                              type = "warning", duration = 8)
      invokeRestart("muffleWarning")
    })) # end withCallingHandlers + isolate

    , error = function(e) {
      # ── DEBUG: print full call stack to R console so the crash site is visible ──
      cat("\n========== RUN OBSERVER ERROR ==========\n")
      cat("Message: ", conditionMessage(e), "\n")
      cat("Call:    ", deparse(conditionCall(e)), "\n")
      cat("Traceback:\n")
      tb <- sys.calls()
      for (i in seq_along(tb)) cat(sprintf("  [%d] %s\n", i, deparse(tb[[i]])[1]))
      cat("========================================\n\n")
      shiny::showNotification(
        paste0("ERROR (see R console for traceback): ", conditionMessage(e)),
        type = "error", duration = 15
      )
    }) # end tryCatch

  }) # end observeEvent(input$run_btn)

  # ================================================================
  # Plots – Training
  # ================================================================

  output$plot_pca_train <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$train_meta_filt)
    plot_pca_training(rv$train_obj, rv$train_meta_filt)
  })

  output$plot_lda_train <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$train_meta_filt)
    plot_lda_training(rv$train_obj, rv$train_meta_filt)
  })

  # Wizard Step 3 PCA: update PC axis choices dynamically
  shiny::observe({
    shiny::req(rv$train_obj)
    n_pcs <- ncol(rv$train_obj$global_pca$x)
    pc_choices <- paste0("PC", seq_len(n_pcs))
    shiny::updateSelectInput(session, "step3_pc_x", choices = pc_choices, selected = "PC1")
    shiny::updateSelectInput(session, "step3_pc_y", choices = pc_choices, selected = "PC2")
  })

  output$plot_pca_train_wizard <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$train_meta_filt)
    col_by <- if (!is.null(input$pca_train_colour_col) && input$pca_train_colour_col != "")
                input$pca_train_colour_col else "Phenotype"
    pc_x <- if (!is.null(input$step3_pc_x)) input$step3_pc_x else "PC1"
    pc_y <- if (!is.null(input$step3_pc_y)) input$step3_pc_y else "PC2"
    plot_pca_training_meta(rv$train_obj, rv$train_meta_filt, col_by,
                           pc_x = pc_x, pc_y = pc_y)
  })

  # Populate "colour by" dropdown from training metadata columns
  shiny::observe({
    shiny::req(rv$train_meta_filt)
    meta_cols <- setdiff(colnames(rv$train_meta_filt), "Sample_ID")
    # Remove columns that are byte-for-byte identical to the internal 'Phenotype' alias,
    # so that users never see a duplicate "Condition / Phenotype" pair in the dropdown.
    if ("Phenotype" %in% meta_cols) {
      ph_vals   <- as.character(rv$train_meta_filt$Phenotype)
      meta_cols <- meta_cols[!vapply(meta_cols, function(col) {
        col != "Phenotype" &&
          identical(as.character(rv$train_meta_filt[[col]]), ph_vals)
      }, logical(1))]
    }
    shiny::updateSelectInput(session, "pca_train_colour_col",
                             choices  = meta_cols,
                             selected = if ("Phenotype" %in% meta_cols) "Phenotype" else meta_cols[1])
  })

  # Conditionally render the training PCA coloured by metadata (Training Plots tab)
  output$plot_pca_train_meta_ui <- shiny::renderUI({
    shiny::req(rv$train_obj)
    if (is.null(input$pca_train_colour_col) || input$pca_train_colour_col == "")
      return(shiny::p("Select a column above to display this plot.", class = "info-text"))
    shiny::plotOutput("plot_pca_train_meta", height = "480px")
  })

  output$plot_pca_train_meta <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$train_meta_filt, input$pca_train_colour_col,
               input$pca_train_colour_col != "")
    plot_pca_training_meta(rv$train_obj, rv$train_meta_filt, input$pca_train_colour_col)
  })

  # ================================================================
  # Plots – Predictions
  # ================================================================

  # ==========================================================================
  # PREDICTIONS TAB — Phenotype Probability Barplot controls
  # ==========================================================================

  # ── Filter: column selector ──────────────────────────────────────────────
  output$probs_filter_col_ui <- shiny::renderUI({
    shiny::req(rv$user_meta)
    cols <- setdiff(colnames(rv$user_meta), c("Sample_ID", "CompoundBatch__"))
    grp  <- isolate(input$test_group_col)
    default_col <- if (!is.null(grp) && grp %in% cols) grp else ""
    shiny::selectInput(
      "probs_filter_col", "Filter samples by",
      choices  = c("(show all)" = "", cols),
      selected = default_col
    )
  })

  # ── Filter: value picker (multi) ─────────────────────────────────────────
  output$probs_filter_val_ui <- shiny::renderUI({
    shiny::req(rv$user_meta, input$probs_filter_col, input$probs_filter_col != "")
    vals <- sort(unique(as.character(rv$user_meta[[input$probs_filter_col]])))
    shinyWidgets::pickerInput(
      "probs_filter_val", "Keep values",
      choices  = vals,
      selected = vals,
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(
        actionsBox         = TRUE,
        liveSearch         = TRUE,
        selectedTextFormat = "count > 2",
        size               = 8
      )
    )
  })

  # ── Facet: column picker (multi) ─────────────────────────────────────────
  output$probs_facet_cols_ui <- shiny::renderUI({
    shiny::req(rv$user_meta)
    cols <- setdiff(colnames(rv$user_meta), c("Sample_ID", "CompoundBatch__"))
    shinyWidgets::pickerInput(
      "probs_facet_cols", "Facet by",
      choices  = cols,
      selected = NULL,
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(
        actionsBox         = TRUE,
        liveSearch         = TRUE,
        noneSelectedText   = "None (no facet)",
        selectedTextFormat = "count > 1",
        size               = 8
      )
    )
  })

  # ── Reactive: filtered + metadata-annotated report df ───────────────────
  probs_plot_data <- shiny::reactive({
    shiny::req(rv$pred_obj, rv$user_meta)
    tryCatch({
    report   <- rv$pred_obj$final_report
    meta     <- rv$user_meta

    # Apply sample filter
    col_sel <- input$probs_filter_col
    val_sel <- input$probs_filter_val
    if (!is.null(col_sel) && col_sel != "" && !is.null(val_sel) &&
        col_sel %in% colnames(meta)) {
      keep_ids <- meta$Sample_ID[as.character(meta[[col_sel]]) %in% val_sel]
      report   <- report[report$Sample_ID %in% keep_ids, , drop = FALSE]
    }

    # Attach metadata columns for faceting
    facet_cols <- input$probs_facet_cols
    if (!is.null(facet_cols) && length(facet_cols) > 0) {
      keep_meta_cols <- intersect(facet_cols, colnames(meta))
      if (length(keep_meta_cols) > 0) {
        report <- merge(
          report,
          meta[, c("Sample_ID", keep_meta_cols), drop = FALSE],
          by = "Sample_ID", all.x = TRUE
        )
      }
    }

    report
    }, error = function(e) {
      cat("\n[CRASH] probs_plot_data: ", conditionMessage(e), "\n"); traceback(3); NULL
    })
  })

  # ── Dynamic plot height: grows with number of facet panels ───────────────
  probs_plot_height <- shiny::reactive({
    shiny::req(rv$pred_obj)
    facet_cols <- input$probs_facet_cols
    if (is.null(facet_cols) || length(facet_cols) == 0 || is.null(rv$user_meta))
      return(420L)

    # Estimate number of facet panels = product of unique levels across facet cols
    n_panels <- prod(vapply(facet_cols, function(col) {
      if (!col %in% colnames(rv$user_meta)) return(1L)
      length(unique(rv$user_meta[[col]]))
    }, integer(1)))

    # Stack vertically: ~300px per row of panels (2 cols per row)
    n_rows <- ceiling(n_panels / 2)
    as.integer(min(max(420L, n_rows * 300L), 2400L))
  })

  output$plot_probs_ui <- shiny::renderUI({
    h <- probs_plot_height()
    shiny::plotOutput("plot_probs", height = paste0(h, "px"))
  })

  output$plot_probs <- shiny::renderPlot({
    shiny::req(rv$pred_obj)
    tryCatch({
      plot_data  <- probs_plot_data()
      facet_cols <- input$probs_facet_cols
      facet_cols <- if (!is.null(facet_cols) && length(facet_cols) > 0) facet_cols else NULL
      plot_prediction_probs(rv$pred_obj, plot_df = plot_data, facet_cols = facet_cols)
    }, error = function(e) { cat("\n[CRASH] plot_probs: ", conditionMessage(e), "\n"); traceback(3); NULL })
  })

  output$plot_lda_combined <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$pred_obj, rv$train_meta_filt)
    tryCatch({
      plot_lda_combined(rv$train_obj, rv$pred_obj, rv$train_meta_filt)
    }, error = function(e) { cat("\n[CRASH] plot_lda_combined: ", conditionMessage(e), "\n"); traceback(3); NULL })
  })

  # ---- Sign matrix heatmap ----
  output$plot_sign_heatmap <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$pred_obj)
    tryCatch({
      # Determine which sign matrix to show
      source_sel <- input$sign_heatmap_source %||% "test"
      train_obj  <- rv$train_obj
      pred_obj   <- rv$pred_obj

      # Build sign matrix depending on source selector
      sign_mat <- if (source_sel == "train" && !is.null(train_obj$sign_mat)) {
        train_obj$sign_mat
      } else if (source_sel == "both" && !is.null(pred_obj$sign_mat_val) &&
                 !is.null(train_obj$sign_mat)) {
        cbind(train_obj$sign_mat, pred_obj$sign_mat_val)
      } else {
        pred_obj$sign_mat_val
      }
      if (is.null(sign_mat)) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = "Sign matrix not available\n(requires binary_tsp / TSP-Delta model)")
        return()
      }

      # Phenotype annotation
      if (source_sel == "train") {
        pheno_vec <- setNames(as.character(train_obj$train_phenotype),
                              colnames(sign_mat))
      } else if (source_sel == "both" &&
                 !is.null(train_obj$sign_mat) && !is.null(pred_obj$sign_mat_val)) {
        pheno_vec <- c(
          setNames(as.character(train_obj$train_phenotype), colnames(train_obj$sign_mat)),
          setNames(
            if (!is.null(rv$user_meta)) {
              rv$user_meta$Phenotype[match(colnames(pred_obj$sign_mat_val),
                                           rv$user_meta$Sample_ID)]
            } else rep("Test", ncol(pred_obj$sign_mat_val)),
            colnames(pred_obj$sign_mat_val)
          )
        )
      } else {
        pheno_vec <- if (!is.null(rv$user_meta)) {
          setNames(rv$user_meta$Phenotype[
            match(colnames(sign_mat), rv$user_meta$Sample_ID)
          ], colnames(sign_mat))
        } else NULL
      }
      pheno_vec[is.na(pheno_vec)] <- "Unknown"

      plot_sign_heatmap(sign_mat,
                        pheno_vec    = pheno_vec,
                        cluster_rows = isTRUE(input$sign_heatmap_cluster_rows))
    }, error = function(e) {
      cat("\n[CRASH] plot_sign_heatmap: ", conditionMessage(e), "\n"); traceback(3); NULL
    })
  })

  # ---- Distance to centroids heatmap ----
  output$plot_dist_heatmap <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$pred_obj)
    tryCatch({
      train_obj <- rv$train_obj
      pred_obj  <- rv$pred_obj
      if (is.null(pred_obj$val_lda_coords)) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = "LDA coordinates not available for this model type.")
        return()
      }
      aggregate_by_cond <- isTRUE(input$dist_heatmap_aggregate)
      cond_vec <- if (aggregate_by_cond && !is.null(pred_obj$prediction))
        pred_obj$prediction else NULL

      plot_distance_heatmap(pred_obj$val_lda_coords, train_obj,
                            condition_vec = cond_vec)
    }, error = function(e) {
      cat("\n[CRASH] plot_dist_heatmap: ", conditionMessage(e), "\n"); traceback(3); NULL
    })
  })

  # ---- Download handlers for new heatmaps ----
  output$dl_sign_heatmap <- shiny::downloadHandler(
    filename = function() paste0("sign_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      p <- isolate(output$plot_sign_heatmap())
      if (!is.null(p))
        ggplot2::ggsave(file, p, width = 14, height = 8, dpi = 180)
    }
  )
  output$dl_dist_heatmap <- shiny::downloadHandler(
    filename = function() paste0("dist_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      p <- isolate(output$plot_dist_heatmap())
      if (!is.null(p))
        ggplot2::ggsave(file, p, width = 10, height = 7, dpi = 180)
    }
  )

  # ---- Step 4: filter UI ----
  output$step4_filter_col_ui <- shiny::renderUI({
    shiny::req(rv$user_meta)
    cols <- setdiff(colnames(rv$user_meta), c("Sample_ID", "CompoundBatch__"))
    grp  <- isolate(input$test_group_col)
    default_col <- if (!is.null(grp) && grp %in% cols) grp else ""
    shiny::selectInput(
      "step4_filter_col", "Filter test samples by column",
      choices  = c("(show all)" = "", cols),
      selected = default_col
    )
  })

  output$step4_filter_val_ui <- shiny::renderUI({
    shiny::req(rv$user_meta, input$step4_filter_col, input$step4_filter_col != "")
    vals <- unique(rv$user_meta[[input$step4_filter_col]])
    shinyWidgets::pickerInput(
      "step4_filter_val", "Show only",
      choices  = as.character(vals),
      selected = as.character(vals),
      multiple = TRUE,
      options  = shinyWidgets::pickerOptions(actionsBox = TRUE, liveSearch = TRUE,
                                             selectedTextFormat = "count > 2")
    )
  })

  # Helper: return filtered pred_obj based on step4 selections
  filtered_pred_obj <- shiny::reactive({
    shiny::req(rv$pred_obj)
    tryCatch({
    pred <- rv$pred_obj

    col_sel <- input$step4_filter_col
    val_sel <- input$step4_filter_val

    # ── DEBUG ──
    message("[DEBUG] filtered_pred_obj fired: col_sel='", col_sel %||% "NULL",
            "' val_lda class=", class(pred$val_lda_coords)[1],
            " val_lda dim=", paste(dim(pred$val_lda_coords), collapse="x"))

    if (!is.null(col_sel) && col_sel != "" && !is.null(val_sel) &&
        !is.null(rv$user_meta) && col_sel %in% colnames(rv$user_meta)) {

      keep_ids <- rv$user_meta$Sample_ID[
        as.character(rv$user_meta[[col_sel]]) %in% val_sel
      ]
      keep_idx <- which(pred$final_report$Sample_ID %in% keep_ids)

      message("[DEBUG] filtered_pred_obj: keep_idx length=", length(keep_idx),
              " val_lda_coords class=", class(pred$val_lda_coords)[1])

      # Always apply the filter — even if keep_idx is empty (0 samples selected).
      # An empty result is correct: the LDA plot will show only training data.
      pred$final_report   <- pred$final_report[keep_idx, , drop = FALSE]
      pred$val_lda_coords <- pred$val_lda_coords[keep_idx, , drop = FALSE]
      if (!is.null(pred$val_pcs))
        pred$val_pcs <- pred$val_pcs[keep_idx, , drop = FALSE]
    }
    pred
    }, error = function(e) { cat("\n[CRASH] filtered_pred_obj: ", conditionMessage(e), "\n"); traceback(3); rv$pred_obj })
  })

  # ---- Step 4 wizard combined LDA plot ----
  output$plot_lda_combined_wizard <- shiny::renderPlot({
    shiny::req(rv$train_obj, rv$train_meta_filt)
    tryCatch({
    pred_to_use <- filtered_pred_obj()
    show_train  <- isTRUE(input$step4_show_all_train)

    # Highlighted samples = rows selected in the report table
    sel_rows  <- input$report_table_wizard_rows_selected
    hl_ids <- if (length(sel_rows) > 0) {
      pred_to_use$final_report$Sample_ID[sel_rows]
    } else character(0)

    train_meta_arg <- if (show_train) rv$train_meta_filt else NULL
    plot_lda_combined(rv$train_obj, pred_to_use, train_meta_arg,
                      highlight_ids = hl_ids)
    }, error = function(e) { cat("\n[CRASH] plot_lda_combined_wizard: ", conditionMessage(e), "\n"); traceback(3); NULL })
  })

  # ---- Step 4 wizard report table ----
  output$report_table_wizard <- DT::renderDT({
    shiny::req(rv$pred_obj)
    pred <- filtered_pred_obj()
    DT::datatable(
      pred$final_report,
      selection = list(mode = "multiple", target = "row"),
      options  = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
                      buttons = c("copy", "csv", "excel")),
      rownames = FALSE, class = "stripe hover compact"
    ) |>
      DT::formatRound(
        columns = setdiff(colnames(pred$final_report),
                          c("Sample_ID", "Predicted_Phenotype")),
        digits  = 2
      )
  })

  # ================================================================
  # Ratio PCA – user samples
  # ================================================================

  # Populate colour dropdown from user metadata columns (after run)
  output$ratio_pca_colour_ui <- shiny::renderUI({
    shiny::req(rv$user_meta)
    meta_cols <- setdiff(colnames(rv$user_meta), c("Sample_ID", "CompoundBatch__"))
    # Default to the user-selected group column (from step 1), falling back to first col
    default_col <- {
      grp <- isolate(input$test_group_col)
      if (!is.null(grp) && grp %in% meta_cols) grp else
        if (length(meta_cols) > 0) meta_cols[1] else ""
    }
    shiny::selectInput(
      inputId  = "ratio_pca_colour_col",
      label    = "Colour by",
      choices  = c("(none – hide plot)" = "", meta_cols),
      selected = default_col
    )
  })

  # Conditionally render ratio PCA plot
  output$plot_ratio_pca_ui <- shiny::renderUI({
    shiny::req(rv$pred_obj)
    col_sel <- input$ratio_pca_colour_col
    if (is.null(col_sel) || col_sel == "")
      return(shiny::p("Select a metadata column above to display this plot.", class = "info-text"))
    shiny::plotOutput("plot_ratio_pca", height = "500px")
  })

  output$plot_ratio_pca <- shiny::renderPlot({
    shiny::req(rv$ratio_res, rv$user_meta, rv$pred_obj)
    col_sel <- input$ratio_pca_colour_col
    shiny::req(!is.null(col_sel) && col_sel != "")

    show_train <- isTRUE(input$ratio_pca_show_training)

    # Build training ratio matrix for overlay (same available features)
    train_ratio_mat <- if (show_train) {
      tryCatch({
        at <- resolve_active_training()
        tr2 <- compute_ratios(at$train_expr, rv$ratio_res$avail_features)
        tr2$ratio_mat
      }, error = function(e) NULL)
    } else NULL

    plot_ratio_pca_user(
      ratio_res       = rv$ratio_res,
      user_meta       = rv$user_meta,
      colour_by       = col_sel,
      train_ratio_mat = train_ratio_mat,
      train_meta      = rv$train_meta_filt,
      show_training   = show_train,
      train_obj       = rv$train_obj
    )
  })

  # ================================================================
  # Prediction table
  # ================================================================

  output$report_table <- DT::renderDT({
    shiny::req(rv$pred_obj)
    DT::datatable(
      rv$pred_obj$final_report,
      options = list(
        pageLength = 20,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel")
      ),
      rownames = FALSE,
      class    = "stripe hover compact"
    ) |>
      DT::formatRound(
        columns = setdiff(colnames(rv$pred_obj$final_report),
                          c("Sample_ID", "Predicted_Phenotype")),
        digits  = 2
      )
  })

  # ================================================================
  # Downloads – Figures
  # ================================================================

  output$dl_pca_train <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_TrainingPCA_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$train_obj, rv$train_meta_filt)
      ggplot2::ggsave(file,
                      plot   = plot_pca_training(rv$train_obj, rv$train_meta_filt),
                      width  = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  output$dl_lda_train <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_TrainingLDA_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$train_obj, rv$train_meta_filt)
      ggplot2::ggsave(file,
                      plot   = plot_lda_training(rv$train_obj, rv$train_meta_filt),
                      width  = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  output$dl_probs_plot <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_Probabilities_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$pred_obj)
      ggplot2::ggsave(file,
                      plot   = plot_prediction_probs(rv$pred_obj),
                      width  = max(8, 0.5 * nrow(rv$pred_obj$final_report)),
                      height = 6, dpi = 200, bg = "white")
    }
  )

  # ── Predictions-tab barplot download (respects current filter + facet) ──
  output$dl_probs_plot_preds <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_Probabilities_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$pred_obj)
      plot_data  <- probs_plot_data()
      facet_cols <- input$probs_facet_cols
      facet_cols <- if (!is.null(facet_cols) && length(facet_cols) > 0) facet_cols else NULL
      n_panels   <- if (!is.null(facet_cols) && !is.null(rv$user_meta)) {
        prod(vapply(facet_cols, function(col) {
          if (!col %in% colnames(rv$user_meta)) 1L
          else length(unique(rv$user_meta[[col]]))
        }, integer(1)))
      } else 1L
      n_samples  <- length(unique(plot_data$Sample_ID))
      plot_width  <- max(8, 0.4 * ceiling(n_samples / max(1L, n_panels)) + 4 * min(2L, n_panels))
      plot_height <- max(5, ceiling(n_panels / 2) * 4)
      ggplot2::ggsave(
        file,
        plot   = plot_prediction_probs(rv$pred_obj, plot_df = plot_data,
                                        facet_cols = facet_cols),
        width  = plot_width,
        height = plot_height,
        dpi    = 200,
        bg     = "white"
      )
    }
  )

  output$dl_lda_combined <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_LDA_Combined_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$train_obj, rv$pred_obj, rv$train_meta_filt)
      ggplot2::ggsave(file,
                      plot   = plot_lda_combined(rv$train_obj, rv$pred_obj,
                                                  rv$train_meta_filt),
                      width  = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  output$dl_pca_train_meta <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_TrainingPCA_meta_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$train_obj, rv$train_meta_filt,
                 input$pca_train_colour_col, input$pca_train_colour_col != "")
      ggplot2::ggsave(file,
                      plot  = plot_pca_training_meta(rv$train_obj, rv$train_meta_filt,
                                                      input$pca_train_colour_col),
                      width = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  output$dl_ratio_pca <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_RatioPCA_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$ratio_res, rv$user_meta, input$ratio_pca_colour_col,
                 input$ratio_pca_colour_col != "")
      show_train <- isTRUE(input$ratio_pca_show_training)
      train_ratio_mat <- if (show_train) {
        tryCatch({
          at <- isolate(resolve_active_training())
          tr <- compute_ratios(at$train_expr, rv$ratio_res$avail_features)
          tr$ratio_mat
        }, error = function(e) NULL)
      } else NULL
      ggplot2::ggsave(file,
                      plot  = plot_ratio_pca_user(rv$ratio_res, rv$user_meta,
                                                   input$ratio_pca_colour_col,
                                                   train_ratio_mat, rv$train_meta_filt,
                                                   show_train,
                                                   train_obj = isolate(rv$train_obj)),
                      width = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  # ================================================================
  # Downloads – CSVs
  # ================================================================

  output$dl_report <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_DiagnosticReport_", Sys.Date(), ".csv"),
    content  = function(file) {
      shiny::req(rv$pred_obj)
      write.csv(rv$pred_obj$final_report, file, row.names = FALSE)
    }
  )

  output$dl_norm_counts <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_NormalisedCounts_", Sys.Date(), ".csv"),
    content  = function(file) {
      shiny::req(rv$norm_counts)
      write.csv(as.data.frame(rv$norm_counts), file)
    }
  )

  output$dl_ratio_mat <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_RatioMatrix_", Sys.Date(), ".csv"),
    content  = function(file) {
      shiny::req(rv$ratio_res)
      write.csv(as.data.frame(rv$ratio_res$ratio_mat), file)
    }
  )

  output$dl_avail_features <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_AvailableFeatures_", Sys.Date(), ".csv"),
    content  = function(file) {
      shiny::req(rv$ratio_res)
      write.csv(rv$ratio_res$avail_features, file, row.names = FALSE)
    }
  )

  # ---- Wizard step 4: download LDA plot (respects current filter + highlights) ----
  output$dl_lda_combined_wizard <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_LDA_Wizard_", Sys.Date(), ".png"),
    content  = function(file) {
      shiny::req(rv$train_obj, rv$train_meta_filt)
      pred_to_use    <- filtered_pred_obj()
      show_train     <- isTRUE(input$step4_show_all_train)
      train_meta_arg <- if (show_train) rv$train_meta_filt else NULL
      sel_rows       <- input$report_table_wizard_rows_selected
      hl_ids <- if (length(sel_rows) > 0) pred_to_use$final_report$Sample_ID[sel_rows] else character(0)
      ggplot2::ggsave(file,
        plot   = plot_lda_combined(rv$train_obj, pred_to_use, train_meta_arg,
                                   highlight_ids = hl_ids),
        width  = 10, height = 7, dpi = 200, bg = "white")
    }
  )

  # ---- Wizard step 4: download diagnostic report CSV (respects current filter) ----
  output$dl_report_wizard <- shiny::downloadHandler(
    filename = function() paste0("ECMExplainER_DiagnosticReport_", Sys.Date(), ".csv"),
    content  = function(file) {
      shiny::req(rv$pred_obj)
      write.csv(filtered_pred_obj()$final_report, file, row.names = FALSE)
    }
  )
}
