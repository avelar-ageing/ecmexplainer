# ==============================================================================
# MODULE: Feature Selection (Synergy Filter)
# Inputs:  train_expr      (genes x samples log2CPM matrix)
#          ecm_dict        (data.frame: Gene, ECMExplained_Function)
#          phenotype_vec   (factor; one label per training sample)
#          max_ratios      (integer; max ratios per panel in final output, default 50)
#          candidate_pool  (integer; candidate pool screened per panel before final
#                           50/50 split — V19 default 100; larger = more thorough search)
#          max_appearances (integer; max times a gene may appear across selected
#                           ratios — V19 default 1 for maximum feature diversity;
#                           NULL = no limit)
#          eta_sq_min      (numeric; minimum eta-squared for the ratio, default 0.30)
#          synergy_delta   (numeric; minimum synergy above best single gene, default 0.05)
# Outputs: data.frame with columns:
#            Gene_A, Gene_B, eta_sq, max_single, synergy_delta, Category, Panel
#          (Panel == Category: "AD" or "SIG")
# ==============================================================================

# Null-coalescing operator (sourced early so all R/ files can use it)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Compute eta-squared for every A×B gene pair plus their individual genes.
#'
#' @param A_mat genes_A × samples numeric matrix (rows = candidate Gene_A)
#' @param D_mat genes_D × samples numeric matrix (rows = candidate Gene_B)
#' @param phenotype factor; length == ncol(A_mat) == ncol(D_mat)
#' @return list: $ratio (n_A × n_D matrix), $A (n_A vector), $D (n_D vector),
#'         $A_names, $D_names
compute_eta_sq_synergy <- function(A_mat, D_mat, phenotype) {

  n          <- ncol(A_mat)
  phen_idx   <- split(seq_len(n), phenotype)
  group_ns   <- sapply(phen_idx, length)

  # Group means
  gm_A <- sapply(phen_idx, function(i) rowMeans(A_mat[, i, drop = FALSE]))
  gm_D <- sapply(phen_idx, function(i) rowMeans(D_mat[, i, drop = FALSE]))

  mu_A <- rowMeans(A_mat)
  mu_D <- rowMeans(D_mat)

  dA <- gm_A - mu_A
  dD <- gm_D - mu_D

  ssb_A <- rowSums(sweep(dA^2, 2, group_ns, "*"))
  ssb_D <- rowSums(sweep(dD^2, 2, group_ns, "*"))

  # Between-group SS for the ratio (vectorised outer product)
  cross_term <- sweep(dA, 2, group_ns, "*") %*% t(dD)
  SS_bet <- outer(ssb_A, ssb_D, "+") - 2 * cross_term

  # Total SS for the ratio
  var_A  <- apply(A_mat, 1, var, na.rm = TRUE)
  var_D  <- apply(D_mat, 1, var, na.rm = TRUE)
  A_ctr  <- A_mat - mu_A
  D_ctr  <- D_mat - mu_D
  cov_AD <- (A_ctr %*% t(D_ctr)) / (n - 1)
  SS_tot <- (outer(var_A, var_D, "+") - 2 * cov_AD) * (n - 1)

  eta_sq_ratio <- pmin(pmax(SS_bet / SS_tot, 0), 1)
  eta_sq_A     <- pmin(pmax(ssb_A / (var_A * (n - 1)), 0), 1)
  eta_sq_D     <- pmin(pmax(ssb_D / (var_D * (n - 1)), 0), 1)

  # Replace NaN/Inf with 0
  eta_sq_ratio[!is.finite(eta_sq_ratio)] <- 0
  eta_sq_A[!is.finite(eta_sq_A)] <- 0
  eta_sq_D[!is.finite(eta_sq_D)] <- 0

  list(
    ratio   = eta_sq_ratio,
    A       = eta_sq_A,
    D       = eta_sq_D,
    A_names = rownames(A_mat),
    D_names = rownames(D_mat)
  )
}


#' Greedy diversity-constrained selection of the top-k ratios.
#'
#' No gene may appear in more than \code{max_appear} selected ratios.
#' Default is 1 (V19 pipeline) — maximum feature diversity.
#'
#' @param ranked_df  data.frame with Gene_A, Gene_B, eta_sq (sorted descending)
#' @param n_ratios   integer; target number of ratios to select
#' @param max_appear integer or NULL (NULL = no limit; 1 = V19 default)
#' @return subset of ranked_df (at most n_ratios rows)
select_top_ratios <- function(ranked_df, n_ratios, max_appear = 1L) {
  no_limit   <- is.null(max_appear)
  gene_count <- new.env(hash = TRUE, parent = emptyenv())
  selected   <- vector("list", n_ratios)
  n_sel      <- 0L

  for (i in seq_len(nrow(ranked_df))) {
    if (n_sel >= n_ratios) break
    a  <- ranked_df$Gene_A[i]
    d  <- ranked_df$Gene_B[i]
    ca <- if (exists(a, envir = gene_count)) get(a, envir = gene_count) else 0L
    cd <- if (exists(d, envir = gene_count)) get(d, envir = gene_count) else 0L

    ok <- no_limit || (ca < max_appear && cd < max_appear)
    if (ok) {
      n_sel            <- n_sel + 1L
      selected[[n_sel]] <- ranked_df[i, , drop = FALSE]
      if (!no_limit) {
        assign(a, ca + 1L, envir = gene_count)
        assign(d, cd + 1L, envir = gene_count)
      }
    }
  }
  do.call(rbind, selected[seq_len(n_sel)])
}


#' Run the full synergy-filter feature selection pipeline.
#'
#' Screens \code{candidate_pool_size} candidates per panel then applies greedy
#' diversity selection to produce \code{max_ratios} final features per panel.
#' This two-stage approach (V19) gives the greedy selector a larger search space,
#' reducing the risk that the diversity constraint cuts off high-quality pairs.
#'
#' @param train_expr        genes × samples log2CPM matrix
#' @param ecm_dict          data.frame: Gene, ECMExplained_Function
#'                          (ASSEMBLY / DESTRUCTION / SIGNALLING)
#' @param phenotype_vec     factor; length == ncol(train_expr)
#' @param max_ratios        integer; final ratios per panel (default 50)
#' @param candidate_pool_size integer; candidate ratios screened per panel before
#'                          the greedy step (default 100; V19 value)
#' @param max_appearances   integer or NULL; max gene appearances across selected
#'                          ratios (default 1; NULL = no limit)
#' @param eta_sq_min        numeric; minimum ratio eta-squared (default 0.30)
#' @param synergy_delta     numeric; required improvement over best single gene (default 0.05)
#' @return data.frame with columns Gene_A, Gene_B, eta_sq, max_single,
#'         synergy_delta, Category ("AD" or "SIG"), Panel ("AD" or "SIG")
run_feature_selection <- function(train_expr,
                                  ecm_dict,
                                  phenotype_vec,
                                  max_ratios          = 50L,
                                  candidate_pool_size = 100L,
                                  max_appearances     = 1L,
                                  eta_sq_min          = 0.30,
                                  synergy_delta       = 0.05) {

  stopifnot(is.matrix(train_expr))
  stopifnot(all(c("Gene", "ECMExplained_Function") %in% colnames(ecm_dict)))
  phenotype_vec <- as.factor(phenotype_vec)
  stopifnot(length(phenotype_vec) == ncol(train_expr))

  # Map gene pools (intersection with available genes in train_expr)
  a_genes   <- intersect(ecm_dict$Gene[ecm_dict$ECMExplained_Function == "ASSEMBLY"],
                         rownames(train_expr))
  d_genes   <- intersect(ecm_dict$Gene[ecm_dict$ECMExplained_Function == "DESTRUCTION"],
                         rownames(train_expr))
  sig_genes <- intersect(ecm_dict$Gene[ecm_dict$ECMExplained_Function == "SIGNALLING"],
                         rownames(train_expr))

  if (length(a_genes) == 0 || length(d_genes) == 0)
    stop("[FeatureSelect] No ASSEMBLY or DESTRUCTION genes found in expression matrix.")
  if (length(sig_genes) == 0)
    message("[FeatureSelect] No SIGNALLING genes found — SIG ratios will be skipped.")

  # ---- AD panel ----
  syn_ad  <- compute_eta_sq_synergy(train_expr[a_genes, , drop = FALSE],
                                    train_expr[d_genes, , drop = FALSE],
                                    phenotype_vec)

  ad_df <- data.frame(
    Gene_A     = rep(syn_ad$A_names, times  = length(syn_ad$D_names)),
    Gene_B     = rep(syn_ad$D_names, each   = length(syn_ad$A_names)),
    eta_sq     = as.vector(syn_ad$ratio),
    eta_sq_A   = rep(syn_ad$A, times = length(syn_ad$D_names)),
    eta_sq_D   = rep(syn_ad$D, each  = length(syn_ad$A_names)),
    stringsAsFactors = FALSE
  )
  ad_df$max_single     <- pmax(ad_df$eta_sq_A, ad_df$eta_sq_D)
  ad_df$synergy_delta  <- ad_df$eta_sq - ad_df$max_single

  ad_ranked <- ad_df[ad_df$synergy_delta > synergy_delta &
                     ad_df$eta_sq > eta_sq_min, , drop = FALSE]
  ad_ranked <- ad_ranked[order(-ad_ranked$eta_sq), , drop = FALSE]

  # Two-stage selection (V19): build a candidate pool first, then apply greedy
  # diversity constraint to get the final max_ratios features.
  ad_pool   <- select_top_ratios(ad_ranked, candidate_pool_size, max_appear = NULL)
  ad_top    <- select_top_ratios(ad_pool,   max_ratios,          max_appear = max_appearances)

  # ---- SIG panel ----
  sig_top <- NULL
  if (length(sig_genes) >= 2) {
    syn_sig <- compute_eta_sq_synergy(train_expr[sig_genes, , drop = FALSE],
                                      train_expr[sig_genes, , drop = FALSE],
                                      phenotype_vec)

    sig_df <- data.frame(
      Gene_A     = rep(syn_sig$A_names, times  = length(syn_sig$D_names)),
      Gene_B     = rep(syn_sig$D_names, each   = length(syn_sig$A_names)),
      eta_sq     = as.vector(syn_sig$ratio),
      eta_sq_A   = rep(syn_sig$A, times = length(syn_sig$D_names)),
      eta_sq_D   = rep(syn_sig$D, each  = length(syn_sig$A_names)),
      stringsAsFactors = FALSE
    )
    sig_df$max_single    <- pmax(sig_df$eta_sq_A, sig_df$eta_sq_D)
    sig_df$synergy_delta <- sig_df$eta_sq - sig_df$max_single

    sig_ranked <- sig_df[sig_df$Gene_A < sig_df$Gene_B &   # no reverse duplicates
                         sig_df$synergy_delta > synergy_delta &
                         sig_df$eta_sq > eta_sq_min, , drop = FALSE]
    sig_ranked <- sig_ranked[order(-sig_ranked$eta_sq), , drop = FALSE]

    # Two-stage selection matching V19
    sig_pool   <- select_top_ratios(sig_ranked, candidate_pool_size, max_appear = NULL)
    sig_top    <- select_top_ratios(sig_pool,   max_ratios,          max_appear = max_appearances)
  }

  # ---- Combine ----
  if (!is.null(ad_top) && nrow(ad_top) > 0)
    ad_top$Category <- "AD"
  if (!is.null(sig_top) && nrow(sig_top) > 0)
    sig_top$Category <- "SIG"

  result <- do.call(rbind, list(ad_top, sig_top))

  if (is.null(result) || nrow(result) == 0)
    stop("[FeatureSelect] No ratios passed the synergy filter. ",
         "Consider relaxing eta_sq_min or synergy_delta.")

  # Panel column expected by downstream code
  result$Panel <- result$Category

  # Clean up helper columns
  result[, c("Gene_A", "Gene_B", "eta_sq", "max_single", "synergy_delta",
             "Category", "Panel"), drop = FALSE]
}


# ==============================================================================
# TSP-DELTA STRATIFIED SELECTION  (V3 / k-TSP manuscript method)
# ==============================================================================
#
# This is the V3 pipeline used in the ECMExplainer k-TSP paper:
#   logCPM -> ECM subset (var>0.005) -> η² candidate pool (> 0.10)
#   -> phenotype-stratified TSP-Δ greedy (n_per_pheno/phenotype/module, max_appear=1)
#
# TSP-Δ for phenotype p = |P(sign(A-B) > 0 | in-p) - P(sign(A-B) > 0 | out-p)|
# Higher TSP-Δ → the sign of this ratio better discriminates p from all others.
# ==============================================================================

#' Build a sign matrix [ratios × samples] from a ratio list and logCPM matrix.
#'
#' @param ratio_vec character vector of "Gene_A / Gene_B" strings
#' @param logcpm_mat genes × samples numeric matrix
#' @return integer matrix [length(ratio_vec) × ncol(logcpm_mat)], values in {-1,0,+1}
build_sign_mat_fs <- function(ratio_vec, logcpm_mat) {
  n_r  <- length(ratio_vec)
  n_s  <- ncol(logcpm_mat)
  mat  <- matrix(0L, nrow = n_r, ncol = n_s,
                 dimnames = list(ratio_vec, colnames(logcpm_mat)))
  for (i in seq_len(n_r)) {
    parts <- strsplit(ratio_vec[i], " / ", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    ga <- trimws(parts[1]); gb <- trimws(parts[2])
    if (ga %in% rownames(logcpm_mat) && gb %in% rownames(logcpm_mat))
      mat[i, ] <- as.integer(sign(logcpm_mat[ga, ] - logcpm_mat[gb, ]))
  }
  mat
}


#' Compute per-phenotype TSP-Δ for a set of candidate ratios.
#'
#' @param ranked_df  data.frame with Gene_A, Gene_B, eta_sq
#' @param logcpm_mat genes × samples numeric matrix
#' @param pheno      factor, length == ncol(logcpm_mat)
#' @param panel_label "AD" or "SIG"
#' @return data.frame with columns: Phenotype, Panel, Gene_A, Gene_B,
#'         Ratio, TSP_Delta, Pr_in_p, Pr_out_p, Direction
compute_tsp_delta <- function(ranked_df, logcpm_mat, pheno, panel_label) {
  ratio_vec <- paste0(ranked_df$Gene_A, " / ", ranked_df$Gene_B)
  sign_mat  <- build_sign_mat_fs(ratio_vec, logcpm_mat)
  up_mat    <- sign_mat > 0L

  pheno <- as.factor(pheno)
  rows  <- vector("list", nlevels(pheno) * nrow(ranked_df))
  idx   <- 0L
  for (p in levels(pheno)) {
    in_p   <- pheno == p
    pr_in  <- rowMeans(up_mat[, in_p,  drop = FALSE])
    pr_out <- rowMeans(up_mat[, !in_p, drop = FALSE])
    delta  <- pr_in - pr_out
    for (r in seq_len(nrow(ranked_df))) {
      idx <- idx + 1L
      rows[[idx]] <- list(
        Phenotype  = p,
        Panel      = panel_label,
        Gene_A     = ranked_df$Gene_A[r],
        Gene_B     = ranked_df$Gene_B[r],
        Ratio      = ratio_vec[r],
        TSP_Delta  = abs(delta[r]),
        Pr_in_p    = pr_in[r],
        Pr_out_p   = pr_out[r],
        Direction  = if (delta[r] > 0) "Up" else if (delta[r] < 0) "Down" else "Flat"
      )
    }
  }
  as.data.frame(do.call(rbind, lapply(rows[seq_len(idx)], as.data.frame)),
                stringsAsFactors = FALSE)
}


#' Phenotype-stratified greedy selection in TSP-Δ rank order.
#'
#' Fills \code{n_per_pheno} ratio slots per phenotype, ensuring no gene
#' appears in more than \code{max_appear} ratios across the module.
#'
#' @param delta_df     data.frame from \code{compute_tsp_delta()}
#' @param n_per_pheno  integer; target ratios per phenotype (default 12)
#' @param max_appear   integer; max appearances per gene across selected ratios (default 1)
#' @return data.frame subset of delta_df (one row per selected ratio × phenotype)
greedy_stratified_tsp <- function(delta_df, n_per_pheno = 12L, max_appear = 1L) {
  delta_df <- delta_df[order(-delta_df$TSP_Delta), , drop = FALSE]
  phenos   <- unique(delta_df$Phenotype)
  counts   <- setNames(rep(0L, length(phenos)), phenos)
  claimed  <- character(0)
  gene_ct  <- new.env(hash = TRUE, parent = emptyenv())
  keep     <- logical(nrow(delta_df))

  for (i in seq_len(nrow(delta_df))) {
    r <- delta_df$Ratio[i]
    p <- delta_df$Phenotype[i]
    a <- delta_df$Gene_A[i]
    b <- delta_df$Gene_B[i]
    if (r %in% claimed)             next
    if (counts[p] >= n_per_pheno)   next
    ca <- if (exists(a, envir = gene_ct)) get(a, envir = gene_ct) else 0L
    cb <- if (exists(b, envir = gene_ct)) get(b, envir = gene_ct) else 0L
    if (ca >= max_appear || cb >= max_appear) next
    keep[i]   <- TRUE
    claimed   <- c(claimed, r)
    counts[p] <- counts[p] + 1L
    assign(a, ca + 1L, envir = gene_ct)
    assign(b, cb + 1L, envir = gene_ct)
  }
  delta_df[keep, , drop = FALSE]
}


#' Full V3 k-TSP feature selection pipeline.
#'
#' Runs: variance filter -> synergy η² pool (threshold > eta_sq_min) ->
#' phenotype-stratified TSP-Δ greedy (n_per_pheno slots/phenotype/module,
#' max_appear = 1).
#'
#' @param logcpm_mat   genes × samples log2(CPM+1) matrix
#' @param ecm_dict     data.frame with columns Gene, ECMExplained_Function
#' @param phenotype_vec factor, length == ncol(logcpm_mat)
#' @param n_per_pheno  integer; target ratios per phenotype per module (default 12)
#' @param eta_sq_min   numeric; minimum η² for candidate pool (default 0.10)
#' @param var_min      numeric; minimum per-gene variance to enter pool (default 0.005)
#' @return data.frame: Ratio_ID, Gene_A, Gene_B, Module, TSP_Delta,
#'         Selecting_Phenotype, Direction (one row per selected ratio)
run_feature_selection_tsp <- function(logcpm_mat,
                                      ecm_dict,
                                      phenotype_vec,
                                      n_per_pheno  = 12L,
                                      eta_sq_min   = 0.10,
                                      var_min      = 0.005) {
  stopifnot(is.matrix(logcpm_mat))
  stopifnot(all(c("Gene", "ECMExplained_Function") %in% colnames(ecm_dict)))
  phenotype_vec <- as.factor(phenotype_vec)
  stopifnot(length(phenotype_vec) == ncol(logcpm_mat))

  ecm_uniq  <- ecm_dict[!duplicated(ecm_dict$Gene), ]
  ecm_in    <- intersect(ecm_uniq$Gene, rownames(logcpm_mat))
  expr_ecm  <- logcpm_mat[ecm_in, , drop = FALSE]
  expr_ecm  <- expr_ecm[apply(expr_ecm, 1, var) > var_min, , drop = FALSE]

  a_genes   <- intersect(ecm_uniq$Gene[ecm_uniq$ECMExplained_Function == "ASSEMBLY"],
                         rownames(expr_ecm))
  d_genes   <- intersect(ecm_uniq$Gene[ecm_uniq$ECMExplained_Function == "DESTRUCTION"],
                         rownames(expr_ecm))
  sig_genes <- intersect(ecm_uniq$Gene[ecm_uniq$ECMExplained_Function == "SIGNALLING"],
                         rownames(expr_ecm))

  if (length(a_genes) == 0 || length(d_genes) == 0)
    stop("[FeatureSelectTSP] No ASSEMBLY or DESTRUCTION genes in matrix.")

  # ---- Candidate pool: η² synergy filter ----
  syn_ad <- compute_eta_sq_synergy(expr_ecm[a_genes, , drop = FALSE],
                                   expr_ecm[d_genes, , drop = FALSE],
                                   phenotype_vec)
  ad_df <- data.frame(
    Gene_A = rep(syn_ad$A_names, times = length(syn_ad$D_names)),
    Gene_B = rep(syn_ad$D_names, each  = length(syn_ad$A_names)),
    eta_sq = as.vector(syn_ad$ratio),
    stringsAsFactors = FALSE
  )
  ad_df$max_single    <- pmax(rep(syn_ad$A, times = length(syn_ad$D_names)),
                               rep(syn_ad$D, each  = length(syn_ad$A_names)))
  ad_ranked <- ad_df[ad_df$eta_sq > eta_sq_min, , drop = FALSE]
  ad_ranked <- ad_ranked[order(-ad_ranked$eta_sq), , drop = FALSE]

  # ---- TSP-Delta + stratified greedy: AD ----
  sel_ad <- NULL
  if (nrow(ad_ranked) > 0) {
    delta_ad <- compute_tsp_delta(ad_ranked, logcpm_mat, phenotype_vec, "AD")
    sel_ad   <- greedy_stratified_tsp(delta_ad, n_per_pheno = n_per_pheno)
  }

  # ---- SIG panel ----
  sel_sig <- NULL
  if (length(sig_genes) >= 2) {
    syn_sig <- compute_eta_sq_synergy(expr_ecm[sig_genes, , drop = FALSE],
                                      expr_ecm[sig_genes, , drop = FALSE],
                                      phenotype_vec)
    sig_df <- data.frame(
      Gene_A = rep(syn_sig$A_names, times = length(syn_sig$D_names)),
      Gene_B = rep(syn_sig$D_names, each  = length(syn_sig$A_names)),
      eta_sq = as.vector(syn_sig$ratio),
      stringsAsFactors = FALSE
    )
    # Exclude self-pairs and reverse duplicates
    sig_df <- sig_df[sig_df$Gene_A < sig_df$Gene_B &
                     sig_df$eta_sq > eta_sq_min, , drop = FALSE]
    sig_df <- sig_df[order(-sig_df$eta_sq), , drop = FALSE]
    if (nrow(sig_df) > 0) {
      delta_sig <- compute_tsp_delta(sig_df, logcpm_mat, phenotype_vec, "SIG")
      sel_sig   <- greedy_stratified_tsp(delta_sig, n_per_pheno = n_per_pheno)
    }
  }

  # ---- Combine ----
  result_rows <- list()
  if (!is.null(sel_ad) && nrow(sel_ad) > 0) {
    ad_uniq <- sel_ad[!duplicated(sel_ad$Ratio), ]
    result_rows[["AD"]] <- data.frame(
      Ratio_ID           = ad_uniq$Ratio,
      Gene_A             = ad_uniq$Gene_A,
      Gene_B             = ad_uniq$Gene_B,
      Module             = "AD",
      TSP_Delta          = ad_uniq$TSP_Delta,
      Selecting_Phenotype = ad_uniq$Phenotype,
      Direction          = ad_uniq$Direction,
      stringsAsFactors   = FALSE
    )
  }
  if (!is.null(sel_sig) && nrow(sel_sig) > 0) {
    sig_uniq <- sel_sig[!duplicated(sel_sig$Ratio), ]
    result_rows[["SIG"]] <- data.frame(
      Ratio_ID           = sig_uniq$Ratio,
      Gene_A             = sig_uniq$Gene_A,
      Gene_B             = sig_uniq$Gene_B,
      Module             = "SIG",
      TSP_Delta          = sig_uniq$TSP_Delta,
      Selecting_Phenotype = sig_uniq$Phenotype,
      Direction          = sig_uniq$Direction,
      stringsAsFactors   = FALSE
    )
  }

  result <- do.call(rbind, result_rows)
  if (is.null(result) || nrow(result) == 0)
    stop("[FeatureSelectTSP] No ratios selected. Try lowering eta_sq_min or n_per_pheno.")

  message(sprintf("[FeatureSelectTSP] Selected: AD=%d | SIG=%d | TOTAL=%d",
                  sum(result$Module == "AD"), sum(result$Module == "SIG"), nrow(result)))
  result
}
