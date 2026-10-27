# ==============================================================================
# ecm_core.R — the ECMExplainer method, as described in the manuscript.
#
# Pure functions only: no Shiny, no global state, no file paths. Everything the
# app does scientifically happens here, so it can be tested headlessly
# (tests/test_manuscript_numbers.R asserts the published numbers).
#
# Section numbers below refer to the Methods of Sheekey, Avelar & Fujiwara,
# "ECMExplainer: a functional matrisome framework for decoding cellular stress
# states and aging trajectories".
# ==============================================================================

# ---- normalisation ----------------------------------------------------------

#' log2(CPM + 1), the training transform (Methods 5.1).
ecmx_logcpm <- function(counts) {
  counts <- as.matrix(counts)
  log2(sweep(counts, 2, colSums(counts) / 1e6, "/") + 1)
}

#' Is a matrix already on a log scale? Used to stop users double-logging.
#' Counts matrices have large maxima; log2-CPM tops out around 20.
ecmx_looks_logged <- function(mat) {
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  is.finite(mx) && mx < 100
}

# ---- ratio encoding (Methods 5.3, deployment rule in 5.5) -------------------

#' Split "GENEA / GENEB" ratio IDs into a two-column data frame.
ecmx_parse_ratios <- function(ratio_ids) {
  parts <- strsplit(ratio_ids, "/", fixed = TRUE)
  data.frame(
    Ratio_ID = ratio_ids,
    Gene_A   = trimws(vapply(parts, `[`, character(1), 1)),
    Gene_B   = trimws(vapply(parts, `[`, character(1), 2)),
    stringsAsFactors = FALSE
  )
}

#' Sign-encode a panel of gene-pair ratios for arbitrary samples.
#'
#' Methods 5.5: "The scorable unit is the ratio pair, not the individual gene."
#' A ratio is computed only when BOTH members are present; ratios with a missing
#' member are neutralised to 0, contributing neither a positive nor a negative
#' vote. Symbol matching is case-insensitive to tolerate GEO casing differences.
#'
#' @param expr genes x samples, log-scale (log2 CPM / RPM — any per-sample
#'   scaling cancels inside the sign).
#' @return list(sign = ratios x samples in {-1,0,1}, coverage = fraction of
#'   ratios with both genes present, missing = ratio IDs that were neutralised)
ecmx_sign_matrix <- function(expr, ratio_ids) {
  rp  <- ecmx_parse_ratios(ratio_ids)
  key <- toupper(rownames(expr))
  # First occurrence wins; callers should collapse duplicate symbols upstream.
  idx <- match(toupper(rp$Gene_A), key)
  jdx <- match(toupper(rp$Gene_B), key)

  out <- matrix(0, nrow = length(ratio_ids), ncol = ncol(expr),
                dimnames = list(ratio_ids, colnames(expr)))
  ok <- !is.na(idx) & !is.na(jdx)
  per_sample <- rep(0, ncol(expr))
  if (any(ok)) {
    r <- expr[idx[ok], , drop = FALSE] - expr[jdx[ok], , drop = FALSE]
    per_sample <- colSums(!is.na(r)) / length(ratio_ids)
    r[is.na(r)] <- 0                       # an NA member cannot vote either way
    out[ok, ] <- sign(r)
  }
  # Coverage counts a ratio as present in a sample only when both genes have a
  # value there, so a gene row that is NA for some samples lowers their coverage.
  list(sign     = out,
       coverage = mean(per_sample),
       coverage_per_sample = stats::setNames(per_sample, colnames(expr)),
       missing  = ratio_ids[!ok])
}

# ---- projection through the frozen model (Methods 5.5) ----------------------

#' Canonical ratio order for a trained model: the columns the PCA was fitted on.
ecmx_model_ratios <- function(model) names(model$pca$center)

#' Project samples into the frozen ECMExplainer coordinate system.
#'
#' Methods 5.5: the stored training centre and scale are applied and the vector
#' rotated by the stored training eigenvectors (predict.prcomp on the saved
#' prcomp object), retaining the same fixed number of components used at
#' training, then classified by the stored LDA. Nothing is re-estimated on the
#' new data — no test-time recentering, because sign features carry no
#' cohort-specific location.
ecmx_project <- function(model, expr) {
  ratios <- ecmx_model_ratios(model)
  sm     <- ecmx_sign_matrix(expr, ratios)

  newdata <- as.data.frame(t(sm$sign))          # samples x ratios
  colnames(newdata) <- ratios
  pcs <- stats::predict(model$pca, newdata = newdata)[, seq_len(model$n_pcs), drop = FALSE]
  pcs <- as.data.frame(pcs)
  colnames(pcs) <- paste0("PC", seq_len(model$n_pcs))

  pred <- stats::predict(model$lda, pcs)
  list(
    coords    = pred$x,                 # LDA discriminant coordinates
    posterior = pred$posterior,
    class     = pred$class,
    pcs       = pcs,
    sign      = sm$sign,
    coverage  = sm$coverage,
    coverage_per_sample = sm$coverage_per_sample,
    missing   = sm$missing
  )
}

#' Training phenotype centroids in the frozen LDA space (Methods 5.6).
#'
#' Per-class means of the training LDA scores. Within-class covariance is
#' approximately isotropic in this sphered space, so plain Euclidean distance is
#' equivalent to the Mahalanobis distance the classifier uses.
ecmx_centroids <- function(model) {
  train_pcs <- as.data.frame(model$pca$x[, seq_len(model$n_pcs), drop = FALSE])
  colnames(train_pcs) <- paste0("PC", seq_len(model$n_pcs))
  tr <- stats::predict(model$lda, train_pcs)$x
  ph <- factor(model$train_phenotype)
  t(vapply(levels(ph), function(p) colMeans(tr[ph == p, , drop = FALSE]),
           numeric(ncol(tr))))
}

#' Euclidean distance from each sample to each training centroid (Methods 5.6).
ecmx_distances <- function(coords, centroids) {
  d <- matrix(NA_real_, nrow(coords), nrow(centroids),
              dimnames = list(rownames(coords), rownames(centroids)))
  for (p in rownames(centroids)) {
    diff <- sweep(coords, 2, centroids[p, ], "-")
    d[, p] <- sqrt(rowSums(diff^2))
  }
  d
}

#' Per-sample confidence and out-of-distribution measures, defined as in the
#' manuscript's validation scripts (SI Table 6 sheet O):
#'   MaxPost  largest posterior; Margin = largest minus second largest;
#'   Entropy  Shannon entropy of the posterior, natural log;
#'   Closest / MinDist  nearest training centroid in LDA space and its distance;
#'   Mahal    Mahalanobis distance of the sample's PC scores from the training
#'            PC scores (training mean and covariance);
#'   OOD95    Mahal above the 95th percentile of chi-squared on n_pcs df.
ecmx_confidence <- function(model, proj) {
  post <- as.matrix(proj$posterior)
  top2 <- t(apply(post, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  ent  <- apply(post, 1, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
  d    <- ecmx_distances(proj$coords, ecmx_centroids(model))
  trp  <- model$pca$x[, seq_len(model$n_pcs), drop = FALSE]
  S    <- stats::cov(trp)
  Si   <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  mah  <- stats::mahalanobis(as.matrix(proj$pcs), colMeans(trp), Si)
  data.frame(MaxPost = top2[, 1], Margin = top2[, 1] - top2[, 2], Entropy = ent,
             Closest = colnames(d)[max.col(-d, ties.method = "first")],
             MinDist = apply(d, 1, min), Mahal = mah,
             OOD95 = mah > stats::qchisq(0.95, model$n_pcs),
             row.names = rownames(post), stringsAsFactors = FALSE)
}

#' Medoid of a set of samples: the real sample minimising total distance to the
#' others (Methods 5.6/5.7). Used instead of a mean so every plotted trajectory
#' point is an actual sample and respects the LDA triangle inequality.
ecmx_medoid <- function(coords) {
  if (nrow(coords) == 1L) return(1L)
  which.min(rowSums(as.matrix(stats::dist(coords))))
}

# ---- training (Methods 5.2–5.4) ---------------------------------------------

#' eta-squared for every A-vs-B ratio, computed on the continuous log-ratio.
#'
#' Vectorised exactly as in the manuscript pipeline: the proportion of variance
#' in the ratio explained by phenotype (one-way).
ecmx_eta_sq <- function(A_mat, D_mat, phenotype) {
  n <- ncol(A_mat)
  idx <- split(seq_len(n), phenotype)
  ns  <- vapply(idx, length, integer(1))
  gmA <- vapply(idx, function(i) rowMeans(A_mat[, i, drop = FALSE]), numeric(nrow(A_mat)))
  gmD <- vapply(idx, function(i) rowMeans(D_mat[, i, drop = FALSE]), numeric(nrow(D_mat)))
  muA <- rowMeans(A_mat); muD <- rowMeans(D_mat)
  dA <- gmA - muA;        dD <- gmD - muD
  ssbA <- rowSums(sweep(dA^2, 2, ns, "*"))
  ssbD <- rowSums(sweep(dD^2, 2, ns, "*"))
  cross <- sweep(dA, 2, ns, "*") %*% t(dD)
  SS_bet <- outer(ssbA, ssbD, "+") - 2 * cross
  vA <- apply(A_mat, 1, stats::var); vD <- apply(D_mat, 1, stats::var)
  cov_AD <- ((A_mat - muA) %*% t(D_mat - muD)) / (n - 1)
  SS_tot <- (outer(vA, vD, "+") - 2 * cov_AD) * (n - 1)
  e <- pmin(pmax(SS_bet / SS_tot, 0), 1)
  e[!is.finite(e)] <- 0
  dimnames(e) <- list(rownames(A_mat), rownames(D_mat))
  e
}

#' Candidate pool: all pairs above the eta-squared floor, ranked descending.
#' @param sym TRUE for within-module (SIG/SIG) pairs, which are unordered and
#'   are therefore deduplicated by alphabetical ordering (Methods 5.3).
#' @param orient how an unordered pair is oriented when `sym` is TRUE:
#'   "alpha" keeps Gene_A < Gene_B alphabetically (ADSIG SIG/SIG and POOL, as
#'   003_FinalTrain_TwoPanels.R); "index" keeps the pair whose Gene_A comes
#'   first in the matrix row order (ALLGENE, as allgenes2.r, whose rows are
#'   sorted by decreasing variance, so Gene_A is the more variable gene).
ecmx_rank_pairs <- function(eta, sym, eta_min = 0.10, category = NA_character_,
                            orient = c("alpha", "index")) {
  orient <- match.arg(orient)
  # Work on column-major cell indices and build gene names only for the pairs
  # that pass the floor; the row order is the same as filtering a full frame.
  idx <- which(as.vector(eta) > eta_min)
  i <- (idx - 1L) %% nrow(eta) + 1L
  j <- (idx - 1L) %/% nrow(eta) + 1L
  if (sym) {
    keep <- if (orient == "alpha") rownames(eta)[i] < colnames(eta)[j] else i < j
    idx <- idx[keep]; i <- i[keep]; j <- j[keep]
  }
  o <- order(-eta[idx])
  df <- data.frame(Gene_A = rownames(eta)[i[o]], Gene_B = colnames(eta)[j[o]],
                   eta_sq = eta[idx[o]], stringsAsFactors = FALSE)
  # A zero-row frame is a legitimate outcome (weak labels, high floor); assigning
  # a length-1 column to it would error, so build the column at the right length.
  df$Category <- rep(category, nrow(df))
  df
}

#' TSP delta for every candidate ratio and phenotype (Methods 5.3):
#' P(A>B | in phenotype) - P(A>B | out of phenotype), signed.
#'
#' Returns a ratios x phenotypes matrix in the row order of `ranked`. Counts of
#' A>B are exact integers, so dividing them by the group sizes gives the same
#' doubles as rowMeans() of the 0/1 indicator; the panel is unchanged. Computed
#' in row chunks, so a genome-wide universe (ALLGENE: ~1.8 million candidate
#' pairs) never needs its full sign matrix in memory.
ecmx_tsp_delta <- function(ranked, expr, phenotype, chunk = 25000L) {
  ia <- match(ranked$Gene_A, rownames(expr)); ib <- match(ranked$Gene_B, rownames(expr))
  ind <- stats::model.matrix(~ phenotype - 1)          # samples x phenotypes, 0/1
  n_in <- colSums(ind); n_out <- nrow(ind) - n_in
  out <- matrix(NA_real_, nrow(ranked), ncol(ind),
                dimnames = list(NULL, levels(phenotype)))
  for (s in seq(1L, nrow(ranked), by = chunk)) {
    k  <- seq.int(s, min(s + chunk - 1L, nrow(ranked)))
    up <- (expr[ia[k], , drop = FALSE] - expr[ib[k], , drop = FALSE]) > 0
    storage.mode(up) <- "double"
    in_cnt  <- up %*% ind
    tot     <- rowSums(up)
    out[k, ] <- sweep(in_cnt, 2, n_in, "/") - sweep(tot - in_cnt, 2, n_out, "/")
  }
  out
}

#' Phenotype-stratified greedy panel selection (Methods 5.3).
#'
#' Candidates are the (ratio, phenotype) combinations, visited in descending
#' order of the ranking statistic. Ties keep phenotype-block order, then the
#' eta-squared order of `ranked`, which is exactly the stable sort of the
#' manuscript pipeline (003_FinalTrain_TwoPanels.R). A candidate is taken when
#' its phenotype's quota is not yet full, the ratio has not been taken for
#' another phenotype, and neither gene already appears in `max_appear` ratios.
#'
#' `rank_by` selects the ordering statistic:
#'   "tsp_delta" — one-versus-rest sign-frequency separation for the tagged
#'                 phenotype. Matches the phenotype-stratified quota architecture,
#'                 and reproduces the published ADSIG-180 panel. The default.
#'   "eta2"      — overall between-phenotype variance in the ratio. Selects a
#'                 different panel; eta-squared is otherwise used only as the
#'                 admission filter.
#' `min_per_module` > 0 first reserves that many ratios for each module
#' (Methods 5.3: "each module must contribute at least five ratios"). Off by
#' default so the published panel reproduces exactly.
#'
#' Works on integer indices throughout and stops as soon as every quota is
#' full, so ALLGENE selection takes well under a second.
#' @param ranked candidate pairs (Gene_A, Gene_B, eta_sq, Category), eta-descending
#' @param delta ratios x phenotypes TSP-delta matrix from ecmx_tsp_delta()
#' @return one row per selected ratio, in selection-rank order
ecmx_greedy_select <- function(ranked, delta, n_per_pheno, max_appear = 1L,
                               rank_by = c("tsp_delta", "eta2"),
                               min_per_module = 0L) {
  rank_by <- match.arg(rank_by)
  N <- nrow(ranked); K <- ncol(delta)
  stat <- if (rank_by == "tsp_delta") abs(as.vector(delta)) else rep(ranked$eta_sq, K)

  # Ordering all N x K candidates is the memory peak for a genome-wide universe.
  # Restricting a stable sort to every candidate at or above a value threshold
  # yields exactly the leading part of the full ordering, so try a generous
  # prefix first and fall back to the full ordering only if quotas are unmet.
  head_n <- max(50000L, 400L * n_per_pheno * K)
  run <- function(cand) {
    ord <- cand[order(-stat[cand])]                      # stable
    row <- (ord - 1L) %% N + 1L                          # candidate -> ratio
    phe <- (ord - 1L) %/% N + 1L                         # candidate -> phenotype
    ecmx_greedy_pass(ranked, delta, row, phe, n_per_pheno, max_appear, min_per_module)
  }
  if (length(stat) > 2L * head_n) {
    thr <- -sort(-stat, partial = head_n)[head_n]
    res <- run(which(stat >= thr))
    if (all(res$counts >= n_per_pheno) && is.null(res$module_short)) return(res$sel)
  }
  run(seq_along(stat))$sel
}

# The greedy pass over candidates already in rank order; see ecmx_greedy_select.
ecmx_greedy_pass <- function(ranked, delta, row, phe, n_per_pheno, max_appear,
                             min_per_module) {
  K <- ncol(delta)

  genes <- unique(c(ranked$Gene_A, ranked$Gene_B))
  ga <- match(ranked$Gene_A, genes); gb <- match(ranked$Gene_B, genes)
  used    <- integer(length(genes))
  claimed <- logical(nrow(ranked))
  counts  <- integer(K)
  keep    <- logical(length(row))

  usable <- function(r) !claimed[r] && used[ga[r]] < max_appear && used[gb[r]] < max_appear
  take <- function(i) {
    r <- row[i]
    keep[i] <<- TRUE; claimed[r] <<- TRUE
    counts[phe[i]] <<- counts[phe[i]] + 1L
    used[ga[r]] <<- used[ga[r]] + 1L; used[gb[r]] <<- used[gb[r]] + 1L
  }

  module_short <- NULL
  if (min_per_module > 0L && !all(is.na(ranked$Category))) {
    cat_c <- ranked$Category[row]
    for (mod in unique(stats::na.omit(cat_c))) {
      got <- 0L
      for (i in which(cat_c == mod)) {
        if (got >= min_per_module) break
        if (counts[phe[i]] >= n_per_pheno || !usable(row[i])) next
        take(i); got <- got + 1L
      }
      if (got < min_per_module) module_short <- c(module_short, mod)
    }
  }

  for (i in seq_along(row)) {
    if (all(counts >= n_per_pheno)) break
    if (keep[i] || counts[phe[i]] >= n_per_pheno || !usable(row[i])) next
    take(i)
  }

  i <- which(keep); r <- row[i]; d <- delta[cbind(r, phe[i])]
  list(counts = counts, module_short = module_short, sel = data.frame(
    Phenotype = colnames(delta)[phe[i]],
    Gene_A = ranked$Gene_A[r], Gene_B = ranked$Gene_B[r],
    Ratio = paste0(ranked$Gene_A[r], " / ", ranked$Gene_B[r]),
    eta_sq = ranked$eta_sq[r], Category = ranked$Category[r],
    TSP_Delta = abs(d),
    Direction = ifelse(d > 0, "Up", ifelse(d < 0, "Down", "Flat")),
    stringsAsFactors = FALSE))
}

#' Number of PCs: >=95% cumulative variance, capped (Methods 5.4).
ecmx_n_pcs <- function(pca, cap = 10L) {
  cv  <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  idx <- which(cv >= 0.95)
  npc <- max(1L, if (length(idx) == 0L) ncol(pca$x) else idx[1])
  min(npc, cap, ncol(pca$x))
}

#' Prepare a training cohort: normalise, filter, and fix the gene universe and
#' the candidate pair structure (Methods 5.2). Done once on all samples; the
#' manuscript's nested leave-one-study-out (019_PanelSweep_Fast.r) keeps this
#' universe fixed and re-runs only selection and fitting inside each fold.
#'
#' @param counts gene x sample matrix of raw counts, or log2-CPM (detected as
#'   for classification; the CPM filter is then applied to 2^x - 1)
#' @param meta   data frame with SampleID, Phenotype and (optionally) Study,
#'   in any row order; matched to the count columns by SampleID
#' @param genelist gene universe: a Gene column, plus (for ADSIG) a function
#'   column whose values are ASSEMBLY / DESTRUCTION / SIGNALLING. Any list of
#'   symbols works for POOL. Ignored for ALLGENE.
#' @param design "ADSIG" (A/D cross-module + SIG/SIG within-module, the
#'   manuscript model), "POOL" (all pairs within the gene list, unstructured
#'   control) or "ALLGENE" (all pairs within the `top_n` most variable expressed
#'   genes genome-wide, the genome-wide control of allgenes2.r)
#' @param top_n ALLGENE universe size
ecmx_prepare <- function(counts, meta, genelist = NULL,
                         design = c("ADSIG", "POOL", "ALLGENE"),
                         var_min = 0.005, top_n = 2000L) {
  design <- match.arg(design)
  counts <- as.matrix(counts)
  meta <- meta[match(colnames(counts), meta$SampleID), , drop = FALSE]
  if (anyNA(meta$SampleID))
    stop("Every column of the count matrix needs a matching SampleID in the metadata.")
  phenotype <- factor(meta$Phenotype)
  if (nlevels(phenotype) < 2L) stop("Need at least two phenotype classes to train.")
  if (min(table(phenotype)) < 2L)
    stop("Every phenotype needs at least 2 samples; smallest class has ",
         min(table(phenotype)), ".")

  if (ecmx_looks_logged(counts)) {
    logcpm <- counts; cpm <- 2^counts - 1
  } else {
    logcpm <- ecmx_logcpm(counts)
    cpm <- sweep(counts, 2, colSums(counts) / 1e6, "/")
  }
  # Expression filter: CPM > 1 in at least as many samples as the smallest
  # phenotype class (Methods 5.2).
  logcpm <- logcpm[rowSums(cpm > 1) >= min(table(phenotype)), , drop = FALSE]

  if (design == "ALLGENE") {
    gv <- apply(logcpm, 1, stats::var)
    pool <- names(sort(gv, decreasing = TRUE))[seq_len(min(top_n, sum(gv > 0)))]
    train_expr <- logcpm[pool, , drop = FALSE]
  } else {
    if (is.null(genelist) || !"Gene" %in% names(genelist))
      stop("The gene list needs a column named Gene.")
    genelist <- genelist[!is.na(genelist$Gene) & nzchar(genelist$Gene), , drop = FALSE]
    genelist <- genelist[!duplicated(genelist$Gene), , drop = FALSE]
    fn_col <- intersect(c("Function", "ECMExplained_Function", "ECMExplainer_Function",
                          "Module"), names(genelist))[1]
    genelist$FN <- if (is.na(fn_col)) NA_character_ else toupper(trimws(genelist[[fn_col]]))
    # Accept the British and American spellings of the signalling module.
    genelist$FN[genelist$FN %in% "SIGNALING"] <- "SIGNALLING"
    if (design == "ADSIG" && all(is.na(genelist$FN)))
      stop("ADSIG needs a Function (or Module) column with Assembly / Destruction / ",
           "Signalling values. Use the POOL design for a plain gene list.")
    train_expr <- logcpm[intersect(genelist$Gene, rownames(logcpm)), , drop = FALSE]
  }
  train_expr <- train_expr[apply(train_expr, 1, stats::var) > var_min, , drop = FALSE]
  if (nrow(train_expr) < 20L)
    stop("Only ", nrow(train_expr), " genes survived filtering — check that ",
         "your gene identifiers are HGNC symbols.")

  blocks <- if (design == "ADSIG") {
    g <- function(fn) intersect(genelist$Gene[genelist$FN %in% fn], rownames(train_expr))
    a <- g("ASSEMBLY"); d <- g("DESTRUCTION"); s <- g("SIGNALLING")
    if (!length(a) || !length(d) || !length(s))
      stop("ADSIG needs Assembly, Destruction and Signalling genes to survive filtering.")
    list(list(A = a, B = d, sym = FALSE, cat = "AD", orient = "alpha"),
         list(A = s, B = s, sym = TRUE,  cat = "SIG", orient = "alpha"))
  } else {
    p <- rownames(train_expr)
    list(list(A = p, B = p, sym = TRUE, cat = design,
              orient = if (design == "ALLGENE") "index" else "alpha"))
  }
  list(design = design, train_expr = train_expr, phenotype = phenotype,
       study = if ("Study" %in% names(meta)) as.character(meta$Study) else NULL,
       blocks = blocks, sample_ids = colnames(counts))
}

#' Select a panel and fit PCA -> LDA on a subset of the prepared samples
#' (Methods 5.3-5.4). With `cols` = all samples this is the deployed model; the
#' nested leave-one-study-out passes the training studies of each fold.
ecmx_fit_panel <- function(prep, cols = seq_along(prep$phenotype), n_per_pheno = 30L,
                           eta_min = 0.10, max_appear = 1L, max_pcs = 10L,
                           rank_by = "tsp_delta", min_per_module = 0L, say = NULL) {
  say <- if (is.function(say)) say else function(f, m) NULL
  ex <- prep$train_expr[, cols, drop = FALSE]
  ph <- droplevels(prep$phenotype[cols])

  say(0.30, "Scoring candidate ratios")
  ranked <- do.call(rbind, lapply(prep$blocks, function(b)
    ecmx_rank_pairs(ecmx_eta_sq(ex[b$A, , drop = FALSE], ex[b$B, , drop = FALSE], ph),
                    sym = b$sym, eta_min, b$cat, orient = b$orient)))
  ranked <- ranked[order(-ranked$eta_sq), , drop = FALSE]
  if (!nrow(ranked)) stop("No gene pair exceeded the eta-squared floor of ", eta_min, ".")

  say(0.55, "Selecting the panel")
  sel <- ecmx_greedy_select(ranked, ecmx_tsp_delta(ranked, ex, ph),
                            n_per_pheno, max_appear, rank_by, min_per_module)
  features <- data.frame(
    Gene_A = sel$Gene_A, Gene_B = sel$Gene_B, Ratio_ID = sel$Ratio,
    TSP_Delta = sel$TSP_Delta, eta_sq = sel$eta_sq,
    Selecting_Phenotype = sel$Phenotype, Direction = sel$Direction,
    Module = sel$Category, stringsAsFactors = FALSE
  )

  say(0.75, "Fitting PCA and LDA")
  sm <- ecmx_sign_matrix(ex, features$Ratio_ID)$sign
  informative <- apply(sm, 1, function(x) length(unique(x)) > 1)
  sm <- sm[informative, , drop = FALSE]
  features <- features[informative, , drop = FALSE]
  if (nrow(sm) < 2L) stop("Fewer than 2 informative ratios were selected.")

  pca <- stats::prcomp(t(sm), center = TRUE, scale. = TRUE)
  n_pcs <- ecmx_n_pcs(pca, max_pcs)
  train_pcs <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  colnames(train_pcs) <- paste0("PC", seq_len(n_pcs))
  list(features = features, sign_matrix_train = sm, pca = pca, n_pcs = n_pcs,
       train_pcs = train_pcs, phenotype = ph,
       # Matrix interface: a formula would capture this function's environment
       # (the whole training matrix) in the saved model.
       lda = MASS::lda(train_pcs, grouping = ph))
}

#' Train an ECMExplainer model end to end (Methods 5.2–5.4).
#'
#' Arguments as ecmx_prepare() and ecmx_fit_panel(); `n_per_pheno` ratios are
#' selected per phenotype, so the panel has (phenotypes x n_per_pheno) ratios.
#' @param progress optional function(fraction, message) for UI feedback
ecmx_train <- function(counts, meta, genelist = NULL,
                       design = c("ADSIG", "POOL", "ALLGENE"),
                       n_per_pheno = 30L,
                       eta_min = 0.10, max_appear = 1L, max_pcs = 10L,
                       rank_by = c("tsp_delta", "eta2"), min_per_module = 0L,
                       var_min = 0.005, top_n = 2000L, progress = NULL,
                       prep = NULL) {
  design  <- match.arg(design)
  rank_by <- match.arg(rank_by)
  say <- function(f, m) if (is.function(progress)) progress(f, m)

  say(0.05, "Normalising and selecting the gene universe")
  if (is.null(prep)) prep <- ecmx_prepare(counts, meta, genelist, design, var_min, top_n)
  phenotype <- prep$phenotype
  fit <- ecmx_fit_panel(prep, n_per_pheno = n_per_pheno, eta_min = eta_min,
                        max_appear = max_appear, max_pcs = max_pcs, rank_by = rank_by,
                        min_per_module = min_per_module, say = say)

  say(0.9, "Cross-validating")
  cv <- suppressWarnings(MASS::lda(fit$train_pcs, grouping = phenotype, CV = TRUE))
  loocv <- mean(cv$class == phenotype) * 100
  cm   <- table(Predicted = cv$class, True = phenotype)
  prec <- diag(cm) / pmax(rowSums(cm), 1)
  rec  <- diag(cm) / pmax(colSums(cm), 1)
  f1   <- ifelse((prec + rec) > 0, 2 * prec * rec / (prec + rec), 0)

  say(1, "Done")
  structure(list(
    encoding = "sign", design = design, max_appear = max_appear,
    rank_by = rank_by, min_per_module = min_per_module, eta_min = eta_min,
    n_per_pheno = n_per_pheno, top_n = if (design == "ALLGENE") top_n else NULL,
    n_training = length(phenotype), panel_size = nrow(fit$features),
    phenotypes = levels(phenotype), train_phenotype = phenotype,
    train_study = prep$study, train_samples = prep$sample_ids,
    universe = rownames(prep$train_expr),
    features = fit$features, sign_matrix_train = fit$sign_matrix_train,
    pca = fit$pca, n_pcs = fit$n_pcs, lda = fit$lda,
    loocv_acc = loocv, macro_f1 = mean(f1), loocv_class = cv$class,
    confusion = cm,
    per_class = data.frame(Phenotype = levels(phenotype),
                           Precision = round(prec, 3), Recall = round(rec, 3),
                           F1 = round(f1, 3), row.names = NULL),
    trained_at = Sys.time()
  ), class = "ecmx_model")
}

#' Cores for parallel cross-validation: ECMX_CORES if set, else all but one.
ecmx_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("ECMX_CORES", "")))
  if (is.na(n)) n <- max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)
  max(1L, n)
}

#' Nested leave-one-study-out accuracy (019_PanelSweep_Fast.r). For each study,
#' the panel is re-selected and re-fitted on the other studies only, then the
#' held-out study is classified. Unlike LOOCV, no held-out sample influences
#' which ratios are chosen, so this estimates how the *procedure* generalises
#' to a new study.
#'
#' @param min_studies phenotypes present in fewer studies than this cannot be
#'   predicted for a held-out study by construction; `restricted` reports
#'   accuracy over the other phenotypes only (SI Table 3 tab C)
#' @return list(per_sample, accuracy, restricted, balanced, macro_f1, per_class)
ecmx_loso <- function(prep, n_per_pheno = 30L, min_studies = 3L, progress = NULL,
                      cores = ecmx_cores(), ...) {
  if (is.null(prep$study)) stop("Leave-one-study-out needs a Study column in the metadata.")
  study <- factor(prep$study); ph <- prep$phenotype
  if (nlevels(study) < 2L) stop("Leave-one-study-out needs at least two studies.")
  ex <- prep$train_expr
  # Folds are independent and deterministic, so they run in parallel where the
  # platform can fork; the result is identical to the sequential loop.
  one_fold <- function(s) {
    tr <- which(study != s); te <- which(study == s)
    if (nlevels(droplevels(ph[tr])) < 2L) return(NULL)
    f <- tryCatch(ecmx_fit_panel(prep, tr, n_per_pheno = n_per_pheno, ...),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    sm <- ecmx_sign_matrix(ex[, te, drop = FALSE], f$features$Ratio_ID)$sign
    pc <- stats::predict(f$pca, t(sm))[, seq_len(f$n_pcs), drop = FALSE]
    colnames(pc) <- paste0("PC", seq_len(f$n_pcs))
    list(te = te, pred = as.character(stats::predict(f$lda, as.data.frame(pc))$class))
  }
  if (is.function(progress)) progress(0.1, sprintf("%d folds", nlevels(study)))
  folds <- if (cores > 1L && .Platform$OS.type == "unix")
    parallel::mclapply(levels(study), one_fold, mc.cores = cores)
  else lapply(levels(study), one_fold)
  pred <- rep(NA_character_, length(ph))
  for (r in folds) if (is.list(r)) pred[r$te] <- r$pred
  truth <- as.character(ph)
  n_st <- tapply(as.character(study), truth, function(x) length(unique(x)))
  well <- names(n_st)[n_st >= min_studies]
  ok <- pred == truth
  rec <- tapply(ok, truth, mean, na.rm = TRUE)
  cm <- table(Predicted = factor(pred, levels = levels(ph)), True = ph)
  pr <- diag(cm) / pmax(rowSums(cm), 1); re <- diag(cm) / pmax(colSums(cm), 1)
  list(per_sample = data.frame(SampleID = prep$sample_ids, Study = as.character(study),
                               True = truth, LOSO_pred = pred, stringsAsFactors = FALSE),
       accuracy = 100 * mean(ok, na.rm = TRUE),
       restricted = 100 * mean(ok[truth %in% well], na.rm = TRUE),
       balanced = 100 * mean(rec),
       macro_f1 = mean(ifelse(pr + re > 0, 2 * pr * re / (pr + re), 0)),
       per_class = data.frame(Phenotype = names(rec), Recall = round(100 * rec, 1),
                              Studies = as.integer(n_st[names(rec)]), row.names = NULL),
       confusion = cm)
}

# ---- accuracy against known labels (SI Table 6) ------------------------------

#' Wilson score interval for k successes in n trials, in percent.
ecmx_wilson <- function(k, n, conf = 0.95) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(1 - (1 - conf) / 2); den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  100 * c(max(0, ctr - hw), min(1, ctr + hw))
}

#' Two-sided exact McNemar test: binomial test on the discordant pairs.
ecmx_mcnemar <- function(correct_a, correct_b) {
  only_a <- sum(correct_a & !correct_b); only_b <- sum(!correct_a & correct_b)
  nd <- only_a + only_b
  list(only_a = only_a, only_b = only_b, discordant = nd,
       p = if (nd == 0) 1 else stats::binom.test(only_b, nd, 0.5)$p.value)
}

#' Validation metrics for one panel's predictions, as in SI Table 6 sheets B-D:
#' micro accuracy with Wilson 95% CI, balanced accuracy (unweighted mean of
#' per-phenotype recall), control specificity (recall on `control`), stress
#' sensitivity (exact-class accuracy on every other truth label), and the
#' unweighted mean and SD of per-dataset accuracy. Rows with no truth are skipped.
ecmx_metrics <- function(pred, truth, dataset = NULL, control = "Proliferating") {
  keep <- !is.na(truth) & nzchar(truth) & !is.na(pred)
  pred <- as.character(pred[keep]); truth <- as.character(truth[keep])
  ok <- pred == truth
  rec <- tapply(ok, truth, mean)
  ci <- ecmx_wilson(sum(ok), length(ok))
  per_ds <- if (!is.null(dataset)) tapply(ok, as.character(dataset[keep]), mean) else NULL
  lev <- sort(unique(c(truth, pred)))
  list(
    summary = data.frame(
      N_scored = length(ok), Micro_accuracy_pct = 100 * mean(ok),
      Micro_CI95_low = ci[1], Micro_CI95_high = ci[2],
      Balanced_accuracy_pct = 100 * mean(rec),
      Control_specificity_pct = if (control %in% truth) 100 * rec[[control]] else NA_real_,
      Stress_sensitivity_pct = if (any(truth != control)) 100 * mean(ok[truth != control]) else NA_real_,
      N_datasets = if (is.null(per_ds)) NA_integer_ else length(per_ds),
      Mean_perDataset_accuracy_pct = if (is.null(per_ds)) NA_real_ else 100 * mean(per_ds),
      SD_perDataset_accuracy = if (length(per_ds) > 1) 100 * stats::sd(per_ds) else NA_real_),
    recall = data.frame(Phenotype = names(rec), N = as.integer(table(truth)[names(rec)]),
                        Recall_pct = 100 * as.numeric(rec), row.names = NULL),
    per_dataset = if (is.null(per_ds)) NULL else
      data.frame(Dataset = names(per_ds), N = as.integer(table(dataset[keep])[names(per_ds)]),
                 Accuracy_pct = 100 * as.numeric(per_ds), row.names = NULL),
    confusion = table(Truth = factor(truth, levels = lev), Predicted = factor(pred, levels = lev)),
    correct = stats::setNames(ok, which(keep)))
}

# ---- ECM Aging Score (Methods 5.9) ------------------------------------------

#' Fit per-ratio age-trend slopes by OLS on the continuous log-ratio, keeping
#' ratios that pass a Benjamini-Hochberg FDR threshold.
ecmx_fit_age_betas <- function(expr, ages, ratio_ids, fdr = 0.05) {
  rp  <- ecmx_parse_ratios(ratio_ids)
  key <- toupper(rownames(expr))
  res <- lapply(seq_len(nrow(rp)), function(i) {
    ia <- match(toupper(rp$Gene_A[i]), key); ib <- match(toupper(rp$Gene_B[i]), key)
    if (is.na(ia) || is.na(ib)) return(NULL)
    r <- expr[ia, ] - expr[ib, ]
    if (stats::var(r, na.rm = TRUE) == 0) return(NULL)
    s <- summary(stats::lm(r ~ ages))$coefficients
    if (nrow(s) < 2) return(NULL)
    data.frame(Ratio_ID = rp$Ratio_ID[i], beta = s[2, 1], p = s[2, 4])
  })
  out <- do.call(rbind, res)
  if (is.null(out) || !nrow(out)) return(out)
  out$fdr <- stats::p.adjust(out$p, method = "BH")
  out[out$fdr < fdr, , drop = FALSE]
}

#' Sign-encoded ECM Aging Score: beta-weighted sum of the sign features
#' (Methods 5.9). Higher = older, by construction of beta as an age slope.
ecmx_aging_score <- function(expr, betas) {
  sm <- ecmx_sign_matrix(expr, betas$Ratio_ID)$sign
  as.vector(t(sm) %*% betas$beta)
}

# ---- trajectory geometry (Methods 5.6) --------------------------------------

#' Angle in degrees between an observed trajectory and each phenotype direction,
#' both drawn from a common origin, with a percentile bootstrap CI.
#'
#' @param from,to coordinate matrices for the two groups defining the trajectory
#' @param centroids training centroids (phenotypes x dimensions)
#' @param n_boot bootstrap resamples; the manuscript uses 2000 with a fixed seed
ecmx_trajectory_angles <- function(from, to, centroids, n_boot = 2000L, seed = 1L) {
  ang <- function(v, d) {
    den <- sqrt(sum(v^2)) * sqrt(sum(d^2))
    if (!is.finite(den) || den == 0) return(NA_real_)
    acos(pmin(pmax(sum(v * d) / den, -1), 1)) * 180 / pi
  }
  one <- function(f, t) {
    o <- colMeans(f); v <- colMeans(t) - o
    vapply(rownames(centroids), function(p) ang(v, centroids[p, ] - o), numeric(1))
  }
  obs <- one(from, to)
  if (nrow(from) < 2L || nrow(to) < 2L || n_boot < 1L) {
    return(data.frame(Phenotype = names(obs), Angle = obs,
                      Lower = NA_real_, Upper = NA_real_,
                      Distance = sqrt(rowSums(sweep(centroids, 2, colMeans(from), "-")^2)),
                      row.names = NULL))
  }
  set.seed(seed)
  boot <- replicate(n_boot, one(from[sample(nrow(from), replace = TRUE), , drop = FALSE],
                                to[sample(nrow(to),   replace = TRUE), , drop = FALSE]))
  data.frame(
    Phenotype = names(obs),
    Angle = apply(boot, 1, stats::median, na.rm = TRUE),
    Lower = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
    Upper = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE),
    Distance = sqrt(rowSums(sweep(centroids, 2, colMeans(from), "-")^2)),
    row.names = NULL
  )
}

# ---- input parsing ----------------------------------------------------------

#' Read a user expression table and collapse duplicate gene symbols.
#'
#' Methods 5.5: "where a symbol mapped to multiple rows, the row with the
#' highest variance in that cohort's log2-CPM matrix was retained".
ecmx_read_expression <- function(path, collapse = TRUE) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                        row.names = NULL)
  if (ncol(df) < 2) stop("Expression file needs a gene column plus >=1 sample column.")
  genes <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat <- mat[!is.na(genes) & nzchar(genes), , drop = FALSE]
  if (collapse && anyDuplicated(rownames(mat))) {
    v <- apply(mat, 1, stats::var, na.rm = TRUE); v[is.na(v)] <- -Inf
    ord <- order(rownames(mat), -v)
    mat <- mat[ord, , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  mat
}
