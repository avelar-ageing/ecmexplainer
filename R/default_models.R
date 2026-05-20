# ==============================================================================
# default_models.R
#
# Registry of built-in (default) ECMExplainer model bundles.
#
# Each entry is a named list with:
#   label       — display name shown in the UI dropdown
#   bundle_file — filename inside training_data/ that is the .rds model bundle
#                 (NULL = no pre-built bundle; load from foundation files instead)
#   description — one-line description shown in the info panel
#
# HOW TO ADD A NEW DEFAULT MODEL
# ─────────────────────────────────────────────────────────────────────────────
# 1. Train your model in the Training Wizard and save the bundle (.rds).
# 2. Copy the .rds into shiny_app/training_data/
# 3. Add an entry to ECM_DEFAULT_MODELS below.
# 4. Restart the app — your model will appear in the "Default Model" dropdown.
#
# Note: the special value bundle_file = NULL tells the app to reconstruct the
# active-training list from the V19 foundation files (the original behaviour).
# ==============================================================================

ECM_DEFAULT_MODELS <- list(

  "ECMExplainer V19 — Human fibroblasts (83 samples, 6 phenotypes)" = list(
    label       = "ECMExplainer V19 — Human fibroblasts (83 samples, 6 phenotypes)",
    bundle_file = "ECMExplorer_v19_TrainedModel.rds",
    description = paste(
      "Reference model trained on 83 human fibroblast bulk RNA-seq samples",
      "from 8 GEO studies (GSE63577, GSE74324, GSE75643, GSE102674,",
      "GSE117337, GSE117444, GSE221693, GSE284342).",
      "Phenotypes: Fibrosis, Inflammation, OxidativeStress,",
      "Proliferating, Quiescence, Senescence.",
      "LOOCV: 90.4%  |  Species: Human  |  Ratios: 10 AD + 9 SIG"
    )
  ),

  "ECMExplainer V3 k-TSP — Human fibroblasts (96 samples, binary TSP)" = list(
    label       = "ECMExplainer V3 k-TSP — Human fibroblasts (96 samples, binary TSP)",
    bundle_file = NULL,
    # data now lives alongside V19 files in training_data/ — portable across machines
    ktsp_data_dir = NULL,           # NULL = use training_data_dir
    counts_file   = "salmon_counts_Combined_FINAL.rds",
    meta_file     = "salmon_metadata_Combined_FINAL.csv",
    gene_list_file = NULL,          # gene list not needed; ratios come from V19 CSV
    feature_method = "binary_tsp",
    description = paste(
      "Binary k-TSP model trained on 96 human fibroblast samples (10 GEO studies).",
      "Features: sign(log\u2082CPM_A \u2212 log\u2082CPM_B) for each ECM ratio pair.",
      "Phenotypes: Fibrosis, Hypoxia, Inflammation, Proliferating, Quiescence, Senescence.",
      "Batch-robust encoding — no ComBat required.",
      "Species: Human  |  Method: binary_tsp  |  Ratios: V19 feature set (10 AD + 9 SIG)"
    )
  )

)

# ---------------------------------------------------------------------------
# .load_default_bundle()
#
# Given a registry entry, return the full active-training list:
#   $train_obj, $features, $train_meta, $train_expr, $ecm_dict
#
# If bundle_file is NULL (legacy mode), load from V19 foundation files.
# If the bundle contains all keys, return it directly.
# If the bundle only contains train_obj (old slim format), supplement from
# the foundation files.
# ---------------------------------------------------------------------------
.load_default_bundle <- function(entry, training_data_dir) {

  # ---- Case 0: k-TSP entry — raw counts + metadata ----
  if (!is.null(entry$counts_file) && !is.null(entry$feature_method) &&
      entry$feature_method == "binary_tsp") {
    # ktsp_data_dir: explicit override path; if NULL fall back to training_data_dir
    kdir <- entry$ktsp_data_dir %||% training_data_dir
    counts_path <- file.path(kdir, entry$counts_file)
    meta_path   <- file.path(kdir, entry$meta_file)
    # gene_list_file is optional for binary_tsp (ratios come from V19 CSV)
    gene_path   <- if (!is.null(entry$gene_list_file))
                     file.path(kdir, entry$gene_list_file)
                   else NULL

    if (!file.exists(counts_path))
      stop("[DefaultModels] k-TSP counts file not found: ", counts_path)
    if (!file.exists(meta_path))
      stop("[DefaultModels] k-TSP metadata file not found: ", meta_path)

    counts_mat <- readRDS(counts_path)           # genes × samples raw counts
    meta_df    <- read.csv(meta_path, row.names = NULL, stringsAsFactors = FALSE)
    gene_list  <- if (!is.null(gene_path) && file.exists(gene_path))
                    read.csv(gene_path, stringsAsFactors = FALSE)
                  else NULL

    # ---- Harmonise metadata column names to app standard ----
    # k-TSP CSV uses SampleID/Study; the app expects Sample_ID/Study_ID
    if ("SampleID" %in% colnames(meta_df) && !"Sample_ID" %in% colnames(meta_df))
      names(meta_df)[names(meta_df) == "SampleID"] <- "Sample_ID"
    if ("Study"    %in% colnames(meta_df) && !"Study_ID"  %in% colnames(meta_df))
      names(meta_df)[names(meta_df) == "Study"]    <- "Study_ID"

    # ---- Load V19 ratio pairs as the feature table ----
    # The V19 ratios (Gene_A, Gene_B, Category) are used for sign encoding.
    # This gives the same high-performing ratio pairs selected for the V19 model,
    # applied with binary sign encoding to the larger 96-sample training set.
    v19_ratios_path <- file.path(
      dirname(dirname(kdir)),                       # back to workspace root
      "shiny_app", "training_data", "ECMExplorer_v19_Final_Ratios.csv"
    )
    # Fallback search: look relative to training_data_dir if kdir path fails
    if (!file.exists(v19_ratios_path)) {
      v19_ratios_path <- file.path(training_data_dir, "ECMExplorer_v19_Final_Ratios.csv")
    }

    features_df <- NULL
    if (file.exists(v19_ratios_path)) {
      feat_raw <- read.csv(v19_ratios_path, stringsAsFactors = FALSE)
      # Normalise column names: Category → Module
      if ("Category" %in% colnames(feat_raw) && !"Module" %in% colnames(feat_raw))
        names(feat_raw)[names(feat_raw) == "Category"] <- "Module"
      # Add Ratio_ID if missing
      if (!"Ratio_ID" %in% colnames(feat_raw))
        feat_raw$Ratio_ID <- paste0(feat_raw$Gene_A, " / ", feat_raw$Gene_B)
      features_df <- feat_raw[, c("Ratio_ID", "Gene_A", "Gene_B", "Module"), drop = FALSE]
    } else {
      message("[DefaultModels] V19 ratios CSV not found; features will be NULL.",
              " Binary TSP requires a feature table — provide via Upload model.")
    }

    return(list(
      train_obj        = NULL,            # not pre-trained; train on demand
      features         = features_df,
      train_meta       = meta_df,
      train_expr       = counts_mat,      # raw counts; normalise in server
      ecm_dict         = gene_list,
      model_name       = entry$label,
      feature_method   = entry$feature_method %||% "binary_tsp",
      training_params  = list(feature_method = entry$feature_method %||% "binary_tsp"),
      version          = "V3-kTSP"
    ))
  }

  # ---- Case 1: no pre-built bundle — reconstruct from foundation files ----
  if (is.null(entry$bundle_file)) {
    return(load_training_data(training_data_dir))
  }

  bundle_path <- file.path(training_data_dir, entry$bundle_file)
  if (!file.exists(bundle_path)) {
    warning("[DefaultModels] Bundle file not found: ", bundle_path,
            " — falling back to foundation files.")
    return(load_training_data(training_data_dir))
  }

  bundle <- tryCatch(
    readRDS(bundle_path),
    error = function(e) {
      warning("[DefaultModels] Could not read bundle: ", e$message,
              " — falling back to foundation files.")
      NULL
    }
  )
  if (is.null(bundle)) return(load_training_data(training_data_dir))

  # ---- Case 2: full bundle (all keys present) — use directly ----
  has_all <- all(c("train_obj", "features", "train_meta", "train_expr") %in% names(bundle))
  if (has_all) return(bundle)

  # ---- Case 3: slim bundle (only train_obj, e.g. old V19 TrainedModel.rds) ----
  # Supplement with foundation data for meta / expr.
  message("[DefaultModels] Bundle lacks train_meta/train_expr; supplementing from foundation files.")
  foundation <- load_training_data(training_data_dir)

  list(
    train_obj  = bundle$lda_model %||% bundle$train_obj %||% bundle,
    features   = bundle$features %||% foundation$features,
    train_meta = bundle$train_meta %||% foundation$train_meta,
    train_expr = bundle$train_expr %||% foundation$train_expr,
    ecm_dict   = bundle$ecm_dict  %||% foundation$ecm_master,
    model_name = bundle$model_name %||% entry$label,
    training_params = bundle$training_params %||% NULL,
    version    = bundle$version %||% "V19"
  )
}

# ---------------------------------------------------------------------------
# .default_model_choices()
# Returns a named character vector suitable for selectInput(choices = ...)
# ---------------------------------------------------------------------------
.default_model_choices <- function() {
  setNames(names(ECM_DEFAULT_MODELS), names(ECM_DEFAULT_MODELS))
}
