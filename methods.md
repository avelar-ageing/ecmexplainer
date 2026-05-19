# ECMExplainer k-TSP — Methods Traceability

Each numbered step maps a section of the manuscript Methods to the implementing script(s).

---

## 1. Training data construction — *n* = 96 samples, 10 GEO studies, 6 phenotypes

**Manuscript:** "Raw counts or author-normalised expression matrices were loaded and converted, where required, to log₂(CPM + 1)…"

| Script | Key lines |
|--------|-----------|
| `001_sampleloading_trainingset_MANUSCRIPT.R` | Combines old Salmon batch + 62 new samples from `txi_ensembl111_ALL_62_FIXED.rds`; applies phenotype mapping (GSE281533 → Hypoxia/Proliferating; GSE295622 → Inflammation/Proliferating); excludes OxidativeStress, HypoxiaFibrosis, and GSE284342; outputs `salmon_counts_Combined_FINAL.rds` + `salmon_metadata_Combined_FINAL.csv` (n = 96) |
| `shiny_app/R/default_models.R` | `ECM_DEFAULT_MODELS` registry entry for k-TSP model; `.load_default_bundle()` loads these two files and harmonises column names (`SampleID` → `Sample_ID`, `Study` → `Study_ID`) |
| `shiny_app/R/app_server.R` | `to_logcpm_ktsp()` call in the binary_tsp run branch converts raw counts → log₂(CPM+1) before training |

---

## 2. Normalisation — log₂(CPM + 1)

**Manuscript:** "raw counts … converted, where required, to log₂(CPM + 1)"

| Script | Key function |
|--------|-------------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `to_logcpm(mat) = log2(sweep(mat, 2, colSums(mat)/1e6, "/") + 1)` |
| `shiny_app/R/mod_binary_tsp.R` | `to_logcpm_ktsp()` — identical; kept local to avoid namespace collision |
| `shiny_app/R/mod_normalize.R` | `normalize_counts()` used for test data in the Test Wizard; outputs `$log_cpm` |

---

## 3. ECM gene filtering — variance > 0.005

**Manuscript:** "the matrix was reduced to the intersection of validation genes and the union of the 144 panel genes … var > 0.005"

| Script | Key lines |
|--------|-----------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `train_expr <- train_expr[apply(train_expr, 1, var) > 0.005, ]` |
| `shiny_app/R/mod_feature_select.R` | `run_feature_selection()` intersects ECM genes from `ecm_dict`, then downstream eta-sq computation implicitly drops zero-variance genes |

---

## 4. Gene pair classification — ASSEMBLY / DESTRUCTION / SIGNALLING

**Manuscript:** "144 sign-encoded matrisome gene ratios partitioned into 72 Assembly/Destruction and 72 Signalling pairs"

| Script | Key lines |
|--------|-----------|
| `ECMExplainer_GeneList.csv` | Master list (1 062 rows); column `ECMExplained_Function` = `ASSEMBLY`, `DESTRUCTION`, or `SIGNALLING` |
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `a_genes`, `d_genes`, `sig_genes` constructed from `ecm_uniq` after `!duplicated(ecm_master$Gene)` |
| `shiny_app/R/mod_feature_select.R` | `run_feature_selection()` splits on `ecm_dict$ECMExplained_Function`; ASSEMBLY × DESTRUCTION → AD panel; SIGNALLING × SIGNALLING → SIG panel |

---

## 5. Candidate pool — synergy η² filter (loose threshold η² > 0.10)

**Manuscript:** "η² candidate pool (loose, > 0.10)"

| Script | Key function |
|--------|-------------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `compute_eta_sq_synergy()` — vectorised outer-product SS computation; `ad_ranked` / `sig_ranked` filtered at `eta_sq > 0.10` and `synergy_delta > 0` |
| `shiny_app/R/mod_feature_select.R` | `compute_eta_sq_synergy()` — same algorithm; `run_feature_selection()` filters at `eta_sq_min` (default 0.30 for V19; set to 0.10 for k-TSP stratified path) |

---

## 6. Feature selection — phenotype-stratified TSP-Δ greedy (12 pairs/phenotype/module, max_appear = 1)

**Manuscript:** "TSP-Δ = |P(sign > 0 | in-phenotype) − P(sign > 0 | out-of-phenotype)| … stratified greedy selection … each gene in at most one ratio"

| Script | Key function |
|--------|-------------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `tsp_delta()` — computes per-phenotype TSP-Δ for every candidate ratio; `greedy_stratified()` — fills `N_PER_PHENO = 12` slots per phenotype per module in TSP-Δ rank order with `max_appear = 1` |
| `shiny_app/R/mod_feature_select.R` | `compute_tsp_delta()` + `greedy_stratified_tsp()` (added in app port); called by `run_feature_selection_tsp()` |

---

## 7. Sign encoding — x_{AB} = sign(log₂X_A − log₂X_B) ∈ {−1, 0, +1}

**Manuscript:** "the sign of a within-sample gene-pair log-ratio is invariant to sequencing depth, normalization, and study-level batch"

| Script | Key function |
|--------|-------------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `build_sign_mat(rlist, expr)` — builds ratio × sample sign matrix |
| `shiny_app/R/mod_binary_tsp.R` | `compute_sign_ratios(logcpm_mat, features_df)` — returns [samples × ratios] integer matrix; missing genes → 0 |
| `shiny_app/R/mod_train.R` | Called inside `train_ktsp()` dispatch when `feature_method = "binary_tsp"` |

---

## 8. Invariant ratio removal

**Manuscript:** "constant-sign ratios dropped before PCA"

| Script | Key lines |
|--------|-----------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `keep <- apply(ratio_mat_sign, 1, function(x) length(unique(x)) > 1)` |
| `shiny_app/R/mod_binary_tsp.R` | `n_unique <- apply(sign_mat, 1, function(x) length(unique(x)))` → `keep <- n_unique > 1L` before PCA |

---

## 9. Split-module PCA — AD (96 × 72) and SIG (96 × 72), ≥95% cumulative variance

**Manuscript:** "PCA was applied separately to each module sign matrix A/D (96 × 72) and SIG (96 × 72), retaining components that cumulatively explained ≥ 95% of variance"

| Script | Key lines |
|--------|-----------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `pca_ad <- prcomp(t(ratio_mat_ad), center = TRUE, scale. = TRUE)` (similar pattern for SIG) |
| `shiny_app/R/mod_binary_tsp.R` | `pca_ad <- prcomp(mat_ad_sign, center = TRUE, scale. = TRUE)`; `.get_npc_ktsp()` returns first PC where cumvar ≥ 0.95 |

---

## 10. LDA — concatenated PCs → 5-dimensional discriminant space, LOOCV

**Manuscript:** "AD and SIG PC scores were concatenated … multi-class linear discriminant analysis (LDA, MASS::lda) … LOOCV (MASS::lda, CV = TRUE)"

| Script | Key lines |
|--------|-----------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `lda_model <- lda(phenotype ~ ., data = train_pcs)`; `lda_cv <- lda(train_pcs, grouping = phenotype, CV = TRUE)` |
| `shiny_app/R/mod_binary_tsp.R` | `train_ktsp()` — same two calls; `loocv_acc <- mean(lda_cv$class == phenotype_vec) * 100` |
| `shiny_app/R/mod_train.R` | `train_ecmexplainer()` — standard ratio path also uses MASS::lda + LOOCV |

---

## 11. Validation projection — reuse training PCA rotations (no leakage)

**Manuscript:** "we reused the rotations learned during training … all validation samples occupy the same subspace as the training samples"

| Script | Key lines |
|--------|-----------|
| `004_kTSP_Validation_SRP113329_MANUSCRIPT.R` | `project_sign()` / `predict(pca_ad, newdata = …)` — projects new sign matrix through stored `pca_ad$rotation` |
| `shiny_app/R/mod_binary_tsp.R` | `predict_ktsp()` — `.pad_to_pca_ktsp()` aligns columns; `predict(pca_ad, newdata = …)` then `predict(train_obj$lda_model, val_pcs)` |
| `shiny_app/R/mod_predict.R` | `predict_ecmexplainer()` dispatches to `predict_ktsp()` when `train_obj$encoding == "sign"` |

---

## 12. LDA posterior probabilities and hard labels

**Manuscript:** "For each sample this produced a posterior probability across the six phenotypes … a hard-assigned label corresponding to the highest posterior"

| Script | Key lines |
|--------|-----------|
| `shiny_app/R/mod_binary_tsp.R` | `preds <- predict(train_obj$lda_model, val_pcs)`; `final_report` includes `Predicted_Phenotype` + posterior columns |
| `shiny_app/R/mod_predict.R` | Standard ratio path: same `predict(lda_model, val_pcs)` call with `preds$posterior` |

---

## 13. Centroid geometry and distance heatmaps

**Manuscript:** "Training phenotype centroids μ_p were computed … per-condition median distances were rendered as a 6-column distance heatmap"

| Script | Key lines |
|--------|-----------|
| `004–008_kTSP_Validation_*.R` | Per-validation script computes `centroid_dist <- sapply(pheno_levels, function(p) sqrt(rowSums((lda_coords - centroid[p,])^2)))` |
| `shiny_app/R/mod_viz_heatmap.R` | `plot_distance_heatmap()` — computes Euclidean distances from `val_lda_coords` to training centroids; renders ggplot2 tile heatmap (added in app port) |

---

## 14. Sign matrix heatmap (ComplexHeatmap / ggplot2 tile)

**Manuscript:** "color scale encodes sign values: blue = −1, cream = 0, red = +1"

| Script | Key lines |
|--------|-----------|
| `002_kTSP_V3_clean_MANUSCRIPT.R` | `pheatmap(ratio_mat_sign, …, color = c("blue","cream","red"))` |
| `shiny_app/R/mod_viz_heatmap.R` | `plot_sign_heatmap()` — ggplot2 `geom_tile()` with same 3-colour scale; avoids `pheatmap` PNG-buffer issues at 300 DPI (added in app port) |

---

## 15. Benchmark models — Model A (all genes + limma) and Model B (ECM genes + limma)

**Manuscript:** "Model A … logCPM matrix of the 18,957 var > 0.005 genes … batch-corrected with limma::removeBatchEffect … Model B … restricted to the 820 ECM-annotated genes"

| Script | Key lines |
|--------|-----------|
| `003_competing_models_MANUSCRIPT.R` | Full competing-model training and evaluation |
| `shiny_app/R/mod_batch.R` | `apply_batch_correction()` wraps `limma::removeBatchEffect` |
| `shiny_app/R/app_ui.R` | `prediction_mode` radio: `"ecm_genes"` (Model B equivalent), `"all_genes"` (Model A equivalent), `"binary_tsp"` (ECMExplainer V3) |

---

## 16. Mahalanobis out-of-distribution flagging

**Manuscript:** "squared Mahalanobis distance to every trained phenotype centroid … flagged at 95% confidence when D²_min > qchisq(0.95, p)"

| Script | Key lines |
|--------|-----------|
| `008_CollaboratorValidation_MANUSCRIPT.R` | Per-sample Mahalanobis distance computation + χ² threshold flagging |
| *(app port — not yet implemented)* | Planned: `plot_mahalanobis_flags()` in `mod_viz_heatmap.R` |

---

## 17. ECM Aging Score

**Manuscript:** "72 [ratios] significantly associated with chronological age … integrated into an ECM Aging Score"

| Script | Key lines |
|--------|-----------|
| `008_GSEaging_kTSP_MANUSCRIPT.R` | Linear model of sign values vs age; Aging Score = weighted sum |
| *(app port — not yet implemented)* | Future extension: dedicated Aging Score tab |

---

*File auto-generated — do not edit manually.*
