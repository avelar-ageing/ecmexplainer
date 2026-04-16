# ==============================================================================
# MODULE: Data Loading
# Loads and validates the training data from the training_data/ directory.
#
# V19 file names (preferred):
#   - ECMExplorer_v19_Foundation.rds       (genes x samples RAW logCPM matrix)
#   - ECMExplorer_v19_Metadata.csv         (Sample_ID, Study_ID, Phenotype, CellType, ...)
#   - ECMExplorer_v19_Final_Ratios.csv     (Category, Gene_A, Gene_B)
#   - ECMExplained_Master_List_Expanded_Fixed.csv (optional; for feature re-selection)
#
# V12 fallback (backward compatibility):
#   - ECMExplorer_v12_logCPM_corrected_82s.rds
#   - ECMExplorer_v12_FullMetadata_82samples.csv
#   - ECMExplainer_V2_Final_100_Features.csv
# ==============================================================================

#' Resolve a training file path, preferring V19 names over V12 fallback names.
#'
#' @param data_dir character; directory to search
#' @param v19_name character; preferred V19 filename
#' @param v12_name character; fallback V12 filename
#' @param required logical; if TRUE, stop if neither file is found
#' @return character; resolved absolute path, or NA if not required and absent
.resolve_training_file <- function(data_dir, v19_name, v12_name, required = TRUE) {
  v19_path <- file.path(data_dir, v19_name)
  v12_path <- file.path(data_dir, v12_name)

  if (file.exists(v19_path)) {
    message("[DataLoader] Using V19 file: ", v19_name)
    return(v19_path)
  }
  if (!is.null(v12_name) && file.exists(v12_path)) {
    message("[DataLoader] V19 file '", v19_name, "' not found — falling back to V12: ", v12_name)
    return(v12_path)
  }
  if (required)
    stop("[DataLoader] Required training file not found.\n",
         "  Looked for: ", v19_path, "\n",
         "  Fallback:   ", v12_path)
  NA_character_
}

#' Load and validate all training data from a directory.
#'
#' Prefers V19 file names; silently falls back to V12 files when V19 are absent.
#'
#' @param data_dir character; path to the directory containing the training files
#' @return list with $train_expr, $train_meta, $features, $ecm_master
load_training_data <- function(data_dir) {

  # ---- Resolve file paths (V19 preferred, V12 fallback) ----
  expr_file     <- .resolve_training_file(data_dir,
                     "ECMExplorer_v19_Foundation.rds",
                     "ECMExplorer_v12_logCPM_corrected_82s.rds")
  meta_file     <- .resolve_training_file(data_dir,
                     "ECMExplorer_v19_Metadata.csv",
                     "ECMExplorer_v12_FullMetadata_82samples.csv")
  features_file <- .resolve_training_file(data_dir,
                     "ECMExplorer_v19_Final_Ratios.csv",
                     "ECMExplainer_V2_Final_100_Features.csv")
  ecm_file      <- file.path(data_dir, "ECMExplained_Master_List_Expanded_Fixed.csv")

  # ---- Load expression matrix ----
  mat_obj    <- readRDS(expr_file)
  train_expr <- if (inherits(mat_obj, "DGEList")) mat_obj$counts else as.matrix(mat_obj)

  # ---- Load metadata ----
  train_meta <- read.csv(meta_file, stringsAsFactors = FALSE)

  # Align metadata rows to expression matrix columns
  idx <- match(colnames(train_expr), train_meta$Sample_ID)
  if (any(is.na(idx)))
    warning("[DataLoader] Some expression matrix samples not found in metadata. Dropping.")
  train_meta <- train_meta[na.omit(idx), , drop = FALSE]
  train_expr <- train_expr[, !is.na(idx), drop = FALSE]

  # ---- Load features ----
  features <- read.csv(features_file, stringsAsFactors = FALSE)
  if (!all(c("Gene_A", "Gene_B") %in% colnames(features)))
    stop("[DataLoader] Features CSV must have columns Gene_A and Gene_B.")

  # Normalise column naming: V19 uses "Category"; V12/old files used "Panel".
  # Downstream code expects "Panel", so add/alias as needed.
  if ("Category" %in% colnames(features) && !"Panel" %in% colnames(features)) {
    features$Panel <- features$Category
  } else if (!"Panel" %in% colnames(features)) {
    message("[DataLoader] No Panel/Category column found in features file – defaulting all to AD.")
    features$Panel <- "AD"
  }

  # Normalise Panel values: "A/D" -> "AD" for consistency
  features$Panel <- gsub("^A/D$", "AD", features$Panel, ignore.case = FALSE)

  # ---- Load ECM master list (optional) ----
  ecm_master <- if (file.exists(ecm_file)) {
    read.csv(ecm_file, stringsAsFactors = FALSE)
  } else {
    message("[DataLoader] ECM Master List not found – feature re-selection will be unavailable.")
    NULL
  }

  list(
    train_expr = train_expr,
    train_meta = train_meta,
    features   = features,
    ecm_master = ecm_master
  )
}
