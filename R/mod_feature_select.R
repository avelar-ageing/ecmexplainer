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
