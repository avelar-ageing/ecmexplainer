# ==============================================================================
# mod_binary_tsp.R  —  Binary k-TSP (sign-encoded ratio) training and prediction
# ==============================================================================
#
# Implements the ECMExplainer V3 binary method:
#   raw counts → log2(CPM+1) → sign(logCPM_A − logCPM_B) ∈ {−1, 0, +1}
#   → split-module PCA (AD + SIG) → concatenate PCs → LDA
#
# KEY INVARIANCE PROPERTY
# The sign of a log-ratio is invariant to library size, normalisation method,
# and study-level batch offsets — so NO batch correction is required or applied.
#
# BUNDLE STRUCTURE  (encoding = "sign")
#   $encoding           "sign"
#   $max_appear         integer (gene uniqueness constraint used at selection)
#   $n_training         integer
#   $phenotypes         character vector of phenotype levels
#   $train_phenotype    factor[n_training] — true labels
#   $features           data.frame: Ratio_ID, Gene_A, Gene_B, Module ("AD"/"SIG")
#   $sign_matrix_train  [n_ratios × n_training]  rownames = Ratio_ID
#   $mat_ad_train_corr  [n_training × n_AD_ratios]  colnames = Ratio_ID
#   $mat_sig_train_corr [n_training × n_SIG_ratios] colnames = Ratio_ID
#   $pca_ad / $pca_sig  prcomp objects (input: sign-encoded [samples × ratios])
#   $n_pcs_ad / $n_pcs_sig integers
#   $lda_model          lda object  (predictor columns: AD_PC1…, SIG_PC1…)
#   $loocv_acc          numeric scalar (%)
#   $loocv_class        factor[n_training] — LOOCV predicted classes
#   $norm_method        "logcpm"  (always)
#   $batch_correction_level  "NONE"
#   $species            character (e.g. "Human")
#   $study_col          character column name used for study in metadata
# ==============================================================================

# ---- Dependency guard --------------------------------------------------------
if (!requireNamespace("MASS", quietly = TRUE))
  stop("Package 'MASS' is required for mod_binary_tsp.R.")


# ==============================================================================
# to_logcpm_ktsp()
#   Convert raw integer count matrix [genes × samples] to log2(CPM + 1).
#   Identical to the k-TSP manuscript helper; kept local to avoid namespace clash.
# ==============================================================================
to_logcpm_ktsp <- function(mat) {
  log2(sweep(mat, 2, colSums(mat) / 1e6, "/") + 1)
}


# ==============================================================================
# compute_sign_ratios()
#   Build sign-encoded {-1, 0, +1} matrix from a logCPM matrix.
#
#   logcpm_mat  [genes × samples]
#   features_df data.frame with columns Gene_A, Gene_B, Ratio_ID (or Ratio)
#
#   Returns [samples × ratios] integer matrix, colnames = Ratio_ID.
#   Missing genes → sign forced to 0 with a warning.
# ==============================================================================
compute_sign_ratios <- function(logcpm_mat, features_df) {
  # Normalise column name (Ratio vs Ratio_ID)
  if (!"Ratio_ID" %in% colnames(features_df) && "Ratio" %in% colnames(features_df))
    features_df$Ratio_ID <- features_df$Ratio

  n_samp  <- ncol(logcpm_mat)
  n_ratio <- nrow(features_df)
  all_genes <- unique(c(features_df$Gene_A, features_df$Gene_B))
  n_miss    <- length(setdiff(all_genes, rownames(logcpm_mat)))
  if (n_miss > 0)
    message(sprintf("[BinaryTSP] %d / %d ratio genes absent — sign forced to 0.",
                    n_miss, length(all_genes)))

  ratio_mat <- vapply(seq_len(n_ratio), function(i) {
    ga <- features_df$Gene_A[i]; gb <- features_df$Gene_B[i]
    if (!ga %in% rownames(logcpm_mat) || !gb %in% rownames(logcpm_mat))
      return(rep(0L, n_samp))
    as.integer(sign(logcpm_mat[ga, ] - logcpm_mat[gb, ]))
  }, integer(n_samp))

  colnames(ratio_mat) <- features_df$Ratio_ID
  rownames(ratio_mat) <- colnames(logcpm_mat)
  ratio_mat
}


# ==============================================================================
# compute_continuous_ratios()
#   Continuous logCPM_A − logCPM_B for biology analytics (NOT the classifier input).
#   Returns [samples × ratios] numeric matrix.
# ==============================================================================
compute_continuous_ratios <- function(logcpm_mat, features_df) {
  if (!"Ratio_ID" %in% colnames(features_df) && "Ratio" %in% colnames(features_df))
    features_df$Ratio_ID <- features_df$Ratio

  n_samp  <- ncol(logcpm_mat)
  n_ratio <- nrow(features_df)

  ratio_mat <- vapply(seq_len(n_ratio), function(i) {
    ga <- features_df$Gene_A[i]; gb <- features_df$Gene_B[i]
    if (!ga %in% rownames(logcpm_mat) || !gb %in% rownames(logcpm_mat))
      return(rep(NA_real_, n_samp))
    as.numeric(logcpm_mat[ga, ] - logcpm_mat[gb, ])
  }, numeric(n_samp))

  colnames(ratio_mat) <- features_df$Ratio_ID
  rownames(ratio_mat) <- colnames(logcpm_mat)
  ratio_mat
}


# ==============================================================================
# .get_npc_ktsp()  [internal]
#   Choose the minimum PCs that capture >= 95% cumulative variance.
# ==============================================================================
.get_npc_ktsp <- function(pca) {
  cv  <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  idx <- which(cv >= 0.95)
  max(1L, if (length(idx) == 0L) ncol(pca$x) else idx[1L])
}


# ==============================================================================
# .pad_to_pca_ktsp()  [internal]
#   Pad a [samples × subset_cols] matrix to the full column layout expected
#   by a stored prcomp object. Fills absent columns with 0.
# ==============================================================================
.pad_to_pca_ktsp <- function(mat, pca_obj) {
  train_feats <- names(pca_obj$center)
  out <- matrix(0, nrow = nrow(mat), ncol = length(train_feats),
                dimnames = list(rownames(mat), train_feats))
  shared <- intersect(colnames(mat), train_feats)
  if (length(shared) > 0) out[, shared] <- mat[, shared, drop = FALSE]
  out
}


# ==============================================================================
# train_ktsp()
#
# Train a binary k-TSP ECMExplainer model on a pre-selected feature set.
#
# logcpm_mat    [genes × samples]  — already log2(CPM+1)
# features_df   data.frame with columns Gene_A, Gene_B, Ratio_ID (or Ratio),
#               Module ("AD" or "SIG")
# phenotype_vec factor[n_samples]
# species       character (e.g. "Human")
# study_col     character; column name for Study in metadata_df (if any)
# metadata_df   optional data.frame; used only to store study_col info in bundle
# max_appear    integer; gene uniqueness used at selection time (informational)
#
# Returns a named list (bundle) with $encoding = "sign".
# ==============================================================================
train_ktsp <- function(logcpm_mat,
                       features_df,
                       phenotype_vec,
                       species      = "Human",
                       study_col    = "Study_ID",
                       metadata_df  = NULL,
                       max_appear   = 1L) {

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for train_ktsp().")

  phenotype_vec <- as.factor(phenotype_vec)
  n_samples     <- ncol(logcpm_mat)
  n_ratios      <- nrow(features_df)

  # Normalise Ratio_ID column name
  if (!"Ratio_ID" %in% colnames(features_df) && "Ratio" %in% colnames(features_df))
    features_df$Ratio_ID <- features_df$Ratio
  if (!"Module" %in% colnames(features_df) && "Category" %in% colnames(features_df))
    features_df$Module <- features_df$Category

  stopifnot(
    all(c("Ratio_ID", "Gene_A", "Gene_B", "Module") %in% colnames(features_df))
  )

  message("[BinaryTSP] Building sign matrix (", n_ratios, " ratios × ",
          n_samples, " samples)…")

  # ---- Sign-encode all ratios ------------------------------------------------
  sign_mat <- t(compute_sign_ratios(logcpm_mat, features_df))
  # sign_mat is now [ratios × samples]; ensure rownames match Ratio_ID
  rownames(sign_mat) <- features_df$Ratio_ID

  # ---- Drop invariant ratios (all same sign — no PCA variance) ---------------
  n_unique <- apply(sign_mat, 1, function(x) length(unique(x)))
  keep     <- n_unique > 1L
  if (!all(keep)) {
    n_drop <- sum(!keep)
    message(sprintf("[BinaryTSP] Dropping %d invariant ratio(s) (constant sign).", n_drop))
    sign_mat    <- sign_mat[keep, , drop = FALSE]
    features_df <- features_df[keep, , drop = FALSE]
  }

  ratio_ids <- features_df$Ratio_ID
  ad_ids    <- features_df$Ratio_ID[features_df$Module == "AD"]
  sig_ids   <- features_df$Ratio_ID[features_df$Module == "SIG"]

  # [samples × ratios] matrices for PCA
  mat_ad_sign  <- t(sign_mat[ad_ids,  , drop = FALSE])   # samples × AD_ratios
  mat_sig_sign <- t(sign_mat[sig_ids, , drop = FALSE])   # samples × SIG_ratios

  # ---- Split-module PCA -------------------------------------------------------
  message("[BinaryTSP] Split-module PCA (AD=", length(ad_ids),
          " ratios, SIG=", length(sig_ids), " ratios)…")
  pca_ad  <- prcomp(mat_ad_sign,  center = TRUE, scale. = TRUE)
  pca_sig <- prcomp(mat_sig_sign, center = TRUE, scale. = TRUE)

  n_pcs_ad  <- .get_npc_ktsp(pca_ad)
  n_pcs_sig <- .get_npc_ktsp(pca_sig)
  message(sprintf("[BinaryTSP] AD PCs: %d  |  SIG PCs: %d", n_pcs_ad, n_pcs_sig))

  # ---- Concatenated PC matrix for LDA ----------------------------------------
  train_pcs <- data.frame(
    pca_ad$x[,  seq_len(n_pcs_ad),  drop = FALSE],
    pca_sig$x[, seq_len(n_pcs_sig), drop = FALSE]
  )
  colnames(train_pcs) <- c(paste0("AD_PC",  seq_len(n_pcs_ad)),
                            paste0("SIG_PC", seq_len(n_pcs_sig)))

  # ---- LDA -------------------------------------------------------------------
  message("[BinaryTSP] Training LDA…")
  lda_model  <- MASS::lda(phenotype_vec ~ ., data = train_pcs)
  lda_cv     <- suppressWarnings(
                  MASS::lda(train_pcs, grouping = phenotype_vec, CV = TRUE))
  loocv_acc  <- mean(lda_cv$class == phenotype_vec) * 100
  message(sprintf("[BinaryTSP] LOOCV accuracy: %.2f%%", loocv_acc))

  # ---- Bundle ----------------------------------------------------------------
  list(
    encoding            = "sign",
    max_appear          = as.integer(max_appear),
    n_training          = n_samples,
    phenotypes          = levels(phenotype_vec),
    train_phenotype     = phenotype_vec,
    species             = species,
    study_col           = study_col,
    norm_method         = "logcpm",
    batch_correction_level = "NONE",

    features            = features_df,
    sign_matrix_train   = sign_mat,       # [ratios × samples]
    mat_ad_train_corr   = mat_ad_sign,    # [samples × AD_ratios]
    mat_sig_train_corr  = mat_sig_sign,   # [samples × SIG_ratios]

    pca_ad              = pca_ad,
    n_pcs_ad            = n_pcs_ad,
    pca_sig             = pca_sig,
    n_pcs_sig           = n_pcs_sig,

    lda_model           = lda_model,
    loocv_acc           = loocv_acc,
    loocv_class         = lda_cv$class,

    # Legacy aliases so unified code paths work
    ad_rows             = ad_ids,
    sig_rows            = sig_ids
  )
}


# ==============================================================================
# predict_ktsp()
#
# Project new samples through a stored sign-encoded model bundle.
#
# logcpm_mat  [genes × test_samples]  — already log2(CPM+1)
# train_obj   list returned by train_ktsp() (or loaded from .rds)
#
# Returns list:
#   $final_report   data.frame (Sample_ID, Predicted_Phenotype, posterior %)
#   $val_pcs        data.frame [samples × AD_PC… + SIG_PC…]  (for LDA plot)
#   $val_lda_coords matrix [samples × LDs]
#   $sign_mat_val   [samples × ratios] integer — concordance use
# ==============================================================================
predict_ktsp <- function(logcpm_mat, train_obj) {

  features_df <- train_obj$features
  pca_ad      <- train_obj$pca_ad
  pca_sig     <- train_obj$pca_sig
  n_pcs_ad    <- train_obj$n_pcs_ad
  n_pcs_sig   <- train_obj$n_pcs_sig
  ad_ids      <- train_obj$ad_rows  %||% features_df$Ratio_ID[features_df$Module == "AD"]
  sig_ids     <- train_obj$sig_rows %||% features_df$Ratio_ID[features_df$Module == "SIG"]

  # ---- Sign-encode test data -------------------------------------------------
  sign_mat_val <- compute_sign_ratios(logcpm_mat, features_df)
  # sign_mat_val [samples × ratios]

  mat_ad_val  <- sign_mat_val[, ad_ids,  drop = FALSE]   # samples × AD_ratios
  mat_sig_val <- sign_mat_val[, sig_ids, drop = FALSE]   # samples × SIG_ratios

  # ---- Project into training PCA space (no data leakage) --------------------
  proj_ad  <- predict(pca_ad,  newdata = .pad_to_pca_ktsp(mat_ad_val,  pca_ad))[,
                seq_len(n_pcs_ad),  drop = FALSE]
  proj_sig <- predict(pca_sig, newdata = .pad_to_pca_ktsp(mat_sig_val, pca_sig))[,
                seq_len(n_pcs_sig), drop = FALSE]

  # ---- Concatenate PCs -------------------------------------------------------
  val_pcs <- data.frame(proj_ad, proj_sig)
  colnames(val_pcs) <- c(paste0("AD_PC",  seq_len(n_pcs_ad)),
                          paste0("SIG_PC", seq_len(n_pcs_sig)))
  rownames(val_pcs) <- colnames(logcpm_mat)

  # Align to exact LDA predictor columns (guards against PC count differences)
  lda_cols <- rownames(train_obj$lda_model$scaling)
  missing  <- setdiff(lda_cols, colnames(val_pcs))
  if (length(missing) > 0) {
    extra <- matrix(0, nrow = nrow(val_pcs), ncol = length(missing),
                    dimnames = list(rownames(val_pcs), missing))
    val_pcs <- cbind(val_pcs, as.data.frame(extra))
  }
  val_pcs <- val_pcs[, lda_cols, drop = FALSE]

  # ---- LDA prediction --------------------------------------------------------
  preds          <- predict(train_obj$lda_model, val_pcs)
  val_lda_coords <- preds$x
  post_probs     <- as.data.frame(round(preds$posterior * 100, 2))

  final_report <- data.frame(
    Sample_ID           = colnames(logcpm_mat),
    Predicted_Phenotype = as.character(preds$class),
    stringsAsFactors    = FALSE
  )
  final_report <- cbind(final_report, post_probs)

  list(
    final_report    = final_report,
    val_pcs         = val_pcs,
    val_lda_coords  = val_lda_coords,
    sign_mat_val    = sign_mat_val
  )
}


# ==============================================================================
# compute_ktsp_concordance()
#
# Per-sample concordance: fraction of ratio signs that match the modal training
# sign for the PREDICTED phenotype.  Used as a quality metric alongside the
# posterior probabilities.
#
# sign_mat_val  [samples × ratios]  from predict_ktsp()$sign_mat_val
# winning_pheno character[n]        predicted class per sample
# train_obj     list from train_ktsp()
#
# Returns data.frame(Sample_ID, Predicted_Phenotype, Concordance [0..1])
# ==============================================================================
compute_ktsp_concordance <- function(sign_mat_val, winning_pheno, train_obj) {
  sm        <- train_obj$sign_matrix_train   # [ratios × n_training]
  pheno_vec <- as.character(train_obj$train_phenotype)
  pheno_lvls <- train_obj$phenotypes

  # Modal sign per ratio per phenotype → [ratios × phenotypes]
  modal_sign <- vapply(pheno_lvls, function(p) {
    idx <- which(pheno_vec == p)
    apply(sm[, idx, drop = FALSE], 1, function(x) {
      tbl <- table(x)
      as.integer(names(tbl)[which.max(tbl)])
    })
  }, integer(nrow(sm)))
  rownames(modal_sign) <- rownames(sm)
  colnames(modal_sign) <- pheno_lvls

  n_s <- nrow(sign_mat_val)
  conc <- vapply(seq_len(n_s), function(i) {
    pred_p <- winning_pheno[i]
    if (is.na(pred_p) || !pred_p %in% colnames(modal_sign)) return(NA_real_)
    shared  <- intersect(colnames(sign_mat_val), rownames(modal_sign))
    if (length(shared) == 0) return(NA_real_)
    v   <- sign_mat_val[i, shared]
    ref <- modal_sign[shared, pred_p]
    nz  <- ref != 0
    if (sum(nz) == 0) return(NA_real_)
    mean(v[nz] == ref[nz])
  }, numeric(1))

  data.frame(
    Sample_ID           = rownames(sign_mat_val),
    Predicted_Phenotype = winning_pheno,
    Concordance         = conc,
    stringsAsFactors    = FALSE
  )
}
