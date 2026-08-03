APP <- "/home/ro/Downloads/Rcode_ECMExplainer-20260309T165711Z-3-001/Rcode_ECMExplainer/shiny_app"

# Source all R modules
for (f in list.files(file.path(APP, "R"), pattern="\\.R$", full.names=TRUE)) {
  cat("Sourcing:", basename(f), "\n")
  source(f)
}

cat("\n--- Testing backend pipeline ---\n")

# 1. Load training data
td <- load_training_data(file.path(APP, "training_data"))
cat("train_expr dims:", dim(td$train_expr), "\n")
cat("train_meta rows:", nrow(td$train_meta), "\n")
cat("features rows:", nrow(td$features), "\n")

# 2. Use 5 random training samples as test data (already logCPM)
set.seed(42)
test_cols <- sample(colnames(td$train_expr), 5)
test_expr_raw <- td$train_expr[, test_cols, drop=FALSE]

# 3. Normalise (simulate raw counts -> logCPM path)
raw_counts_sim <- round(pmax((2^test_expr_raw - 1) * 1e6 / matrix(colSums(2^test_expr_raw - 1), nrow=nrow(test_expr_raw), ncol=5, byrow=TRUE), 0))
norm_res <- normalize_counts(raw_counts_sim)
cat("Norm mat dims:", dim(norm_res$log_cpm), "\n")

# 4. Ratio detection
ratio_res <- compute_ratios(test_expr_raw, td$features)
cat("Available ratios:", ratio_res$n_avail, "/", ratio_res$n_total, "\n")

# 5. Compute ratios on training data using same available features
train_ratio_res <- compute_ratios(td$train_expr, ratio_res$avail_features)

# 6. Train
pheno_vec <- factor(td$train_meta$Phenotype[match(colnames(train_ratio_res$ratio_mat), td$train_meta$Sample_ID)])
train_obj <- train_ecmexplainer(train_ratio_res$ratio_mat, train_ratio_res$avail_features, pheno_vec)
cat("LDA model classes:", paste(train_obj$lda_model$lev, collapse=", "), "\n")
cat("n_pcs_ad:", train_obj$n_pcs_ad, "  n_pcs_sig:", train_obj$n_pcs_sig, "\n")

# 7. Predict
pred_obj <- predict_ecmexplainer(test_expr_raw, train_obj, ratio_res$avail_features)
cat("\n=== Predictions ===\n")
print(pred_obj$final_report[, 1:3])

cat("\n--- All smoke tests passed ---\n")
