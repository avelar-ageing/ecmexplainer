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
