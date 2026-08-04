# ==============================================================================
# ECMExplainer -- 002c FAST v2: POOL + ADSIG, balanced quota, MAX_PCS, NESTED LOSO
# ------------------------------------------------------------------------------
# ROOT-CAUSE FIX for the slow POOL run:
#   The big sign(logCPM_A - logCPM_B) matrix is INVARIANT to which samples are
#   held out (dropping a study only removes COLUMNS). So we build the full-pool
#   sign matrix ONCE over all 96 samples, then every LOSO fold merely SUBSETS
#   columns and recomputes the cheap eta-sq + TSP_Delta from that cached sign
#   matrix. No build_sign / pair enumeration is ever repeated inside a fold.
#   (Old 002_V5 was fast only because it had NO LOSO -- it built the pool once.)
# Still correct LOSO: ranking/selection uses TRAINING columns only -> no leak.
# Design: balanced n/6 quota, unbiased AD:SIG, MAX_PCS cap, max_appear=1.
# ==============================================================================
suppressPackageStartupMessages({ library(MASS); library(dplyr); library(tidyr); library(ggplot2) })
filter <- dplyr::filter; select <- dplyr::select; set.seed(1)
out_dir <- "Salmon_kTSP_FAST_NESTED"; dir.create(out_dir, showWarnings = FALSE)

# ---- CONFIG -----------------------------------------------------------------
N_PHENO      <- 6L
TARGET_SWEEP <- seq(60L, 180L, by = 20L)
MAX_NPP      <- max(TARGET_SWEEP) %/% N_PHENO
MAX_PCS      <- 10L
MAX_APPEAR   <- 1L
ETA_MIN      <- 0.10
MODULE_FLOOR <- 5L
VERBOSE      <- TRUE
# -----------------------------------------------------------------------------

# ---- [1] LOAD + GLOBAL FILTER -----------------------------------------------
salmon_counts <- readRDS("salmon_counts_Combined_FINAL.rds")
meta <- readr::read_csv("salmon_metadata_Combined_FINAL.csv", show_col_types = FALSE)
meta <- meta[match(colnames(salmon_counts), meta$SampleID), ]
stopifnot(all(meta$SampleID == colnames(salmon_counts)))
to_logcpm <- function(mat) log2(sweep(mat, 2, colSums(mat)/1e6, "/") + 1)
all_logcpm <- to_logcpm(salmon_counts)
phenotype <- factor(meta$Phenotype); study <- factor(meta$Study)
stopifnot(nlevels(phenotype) == N_PHENO)
cpm_mat <- sweep(salmon_counts, 2, colSums(salmon_counts)/1e6, "/")
min_samp <- min(table(phenotype))
expressed <- rownames(salmon_counts)[ rowSums(cpm_mat > 1) >= min_samp ]
all_logcpm <- all_logcpm[expressed, , drop = FALSE]
cat("N =", ncol(all_logcpm), "| studies =", nlevels(study), "| genes =", length(expressed), "\n")
cat("Phenotype x Study (confounding check):\n"); print(table(phenotype, study))

ecm_master <- read.csv("ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors = FALSE)
ecm_master <- ecm_master[!duplicated(ecm_master$Gene), ]
fn_col <- intersect(c("Function","ECMExplained_Function","ECMExplainer_Function"), names(ecm_master))[1]
stopifnot("No function-label column found" = !is.na(fn_col))
ecm_master$FN <- toupper(trimws(ecm_master[[fn_col]]))

# Genes for the pool: ECM intersect + var>0.005 on FULL data (matches 002_V5)
ecm_in <- intersect(ecm_master$Gene, rownames(all_logcpm))
train_expr <- all_logcpm[ecm_in, , drop = FALSE]
train_expr <- train_expr[apply(train_expr,1,var) > 0.005, , drop = FALSE]
cat("ECM genes (var>0.005):", nrow(train_expr), "\n")

# ---- [2] BUILD THE FULL SIGN MATRICES ONCE (the expensive step, done 1x) ----
# POOL: all gene pairs (i<j). ADSIG: AD (A x D) U SIG (S x S, i<j).
make_pairs <- function(genes_A, genes_B, sym){
  if (sym){
    idx <- which(upper.tri(matrix(0, length(genes_A), length(genes_A))), arr.ind=TRUE)
    data.frame(Gene_A=genes_A[idx[,1]], Gene_B=genes_A[idx[,2]], stringsAsFactors=FALSE)
  } else {
    expand.grid(Gene_A=genes_A, Gene_B=genes_B, stringsAsFactors=FALSE)
  }
}
build_full_sign <- function(pairs){
  A <- train_expr[pairs$Gene_A, , drop=FALSE]
  B <- train_expr[pairs$Gene_B, , drop=FALSE]
  sm <- sign(A - B)                              # vectorised; ONE pass
  storage.mode(sm) <- "integer"
  rownames(sm) <- paste0(pairs$Gene_A," / ",pairs$Gene_B)
  sm
}
pool_genes <- rownames(train_expr)
ag <- intersect(ecm_master$Gene[ecm_master$FN=="ASSEMBLY"],    pool_genes)
dg <- intersect(ecm_master$Gene[ecm_master$FN=="DESTRUCTION"], pool_genes)
sg <- intersect(ecm_master$Gene[ecm_master$FN=="SIGNALLING"],  pool_genes)

cat("\n[sign] building full-cohort sign matrices ONCE...\n")
pairs_pool <- make_pairs(pool_genes, pool_genes, sym=TRUE)
SM_POOL    <- build_full_sign(pairs_pool); MOD_POOL <- rep("POOL", nrow(SM_POOL))
cat("  POOL sign matrix:", nrow(SM_POOL), "pairs\n")
pairs_ad  <- make_pairs(ag, dg, sym=FALSE)
pairs_sig <- make_pairs(sg, sg, sym=TRUE)
SM_AD  <- build_full_sign(pairs_ad);  SM_SIG <- build_full_sign(pairs_sig)
SM_ADSIG <- rbind(SM_AD, SM_SIG)
MOD_ADSIG <- c(rep("AD", nrow(SM_AD)), rep("SIG", nrow(SM_SIG)))
cat("  ADSIG sign matrix:", nrow(SM_ADSIG), "pairs (AD", nrow(SM_AD), "SIG", nrow(SM_SIG), ")\n")

# ---- [3] CHEAP fold ops operating on the cached sign matrix -----------------
# eta-sq on signs (group-mean separation) + TSP_Delta, on a COLUMN SUBSET only.
rank_and_delta <- function(SM, MOD, cols){
  sub <- SM[, cols, drop=FALSE]
  ph  <- droplevels(phenotype[cols])
  # drop pairs invariant across these samples (also kills self/degenerate)
  vary <- which(matrixStats_rowAny_diff(sub))
  sub <- sub[vary, , drop=FALSE]; mod <- MOD[vary]
  up  <- sub > 0
  lev <- levels(ph)
  # eta_sq proxy on signs: between-group var of P(up) / total (fast, monotone w/ separation)
  pin_mat <- sapply(lev, function(p) rowMeans(up[, ph==p, drop=FALSE]))
  grand   <- rowMeans(up)
  gn      <- as.integer(table(ph))
  ssb     <- rowSums(sweep((pin_mat - grand)^2, 2, gn, "*"))
  sst     <- apply(up, 1, function(x) sum((x - mean(x))^2))
  eta     <- ifelse(sst>0, ssb/sst, 0)
  keep_e  <- which(eta > ETA_MIN)
  sub <- sub[keep_e,,drop=FALSE]; mod <- mod[keep_e]; up <- up[keep_e,,drop=FALSE]
  pin_mat <- pin_mat[keep_e,,drop=FALSE]
  # TSP_Delta per phenotype (max one-vs-rest sign-freq gap), keep best phenotype tag
  pout_mat <- sapply(lev, function(p) rowMeans(up[, ph!=p, drop=FALSE]))
  dmat <- abs(pin_mat - pout_mat)
  best <- max.col(dmat, ties.method="first")
  data.frame(Ratio=rownames(sub),
             Gene_A=sub_geneA(rownames(sub)), Gene_B=sub_geneB(rownames(sub)),
             Phenotype=lev[best], TSP_Delta=dmat[cbind(seq_len(nrow(dmat)),best)],
             Module=mod, stringsAsFactors=FALSE)
}
matrixStats_rowAny_diff <- function(m) apply(m, 1, function(x) length(unique(x)) > 1)
sub_geneA <- function(r) sub(" / .*","", r); sub_geneB <- function(r) sub(".* / ","", r)

greedy_quota <- function(d, n_per, max_appear=MAX_APPEAR, module_floor=0L){
  d <- arrange(d, desc(TSP_Delta)); claimed <- character(0)
  gc <- new.env(hash=TRUE, parent=emptyenv())
  phl <- levels(factor(d$Phenotype)); cnt <- setNames(integer(length(phl)), phl)
  mods <- unique(d$Module); modc <- setNames(integer(length(mods)), mods)
  keep <- logical(nrow(d))
  free <- function(a,b){ ca<-if(exists(a,envir=gc))get(a,envir=gc)else 0L
  cb<-if(exists(b,envir=gc))get(b,envir=gc)else 0L; ca<max_appear&&cb<max_appear }
  grab <- function(i){ a<-d$Gene_A[i]; b<-d$Gene_B[i]
  keep[i]<<-TRUE; claimed<<-c(claimed,d$Ratio[i]); cnt[d$Phenotype[i]]<<-cnt[d$Phenotype[i]]+1L
  modc[d$Module[i]]<<-modc[d$Module[i]]+1L
  assign(a,(if(exists(a,envir=gc))get(a,envir=gc)else 0L)+1L,envir=gc)
  assign(b,(if(exists(b,envir=gc))get(b,envir=gc)else 0L)+1L,envir=gc) }
  ok <- function(i) !(d$Ratio[i] %in% claimed) && free(d$Gene_A[i], d$Gene_B[i])
  if (length(mods)>1 && module_floor>0)
    for (m in mods) for (i in which(d$Module==m)){
      if (modc[m]>=module_floor) break
      if (cnt[d$Phenotype[i]]>=n_per) next
      if (ok(i)) grab(i) }
  for (p in phl) for (i in which(d$Phenotype==p)){ if (cnt[p]>=n_per) break; if (ok(i)) grab(i) }
  d[keep,,drop=FALSE] %>% distinct(Ratio,.keep_all=TRUE)
}
get_npc <- function(p, cap=MAX_PCS){ cv<-cumsum(p$sdev^2)/sum(p$sdev^2); i<-which(cv>=0.95)
npc<-max(1L,if(length(i)==0)ncol(p$x)else i[1]); min(npc,cap,ncol(p$x)) }
per_class <- function(pred,truth,lev){ cm<-table(factor(pred,levels=lev),factor(truth,levels=lev))
prec<-diag(cm)/pmax(colSums(t(cm)),1); rec<-diag(cm)/pmax(colSums(cm),1)
f1<-ifelse((prec+rec)>0,2*prec*rec/(prec+rec),0); data.frame(Phenotype=lev,F1=round(f1,3)) }

# fit on a SET of selected ratio rownames using cached sign matrix + given cols
fit_panel <- function(SM, sel, cols){
  sm <- SM[sel$Ratio, cols, drop=FALSE]
  vary <- matrixStats_rowAny_diff(sm); sm <- sm[vary,,drop=FALSE]; sel <- sel[vary,,drop=FALSE]
  pca <- prcomp(t(sm), center=TRUE, scale.=TRUE); np <- get_npc(pca)
  trp <- as.data.frame(pca$x[,1:np,drop=FALSE]); colnames(trp)<-paste0("PC",1:np)
  list(lda=lda(trp, grouping=droplevels(phenotype[cols])), pca=pca, np=np,
       sel=sel, trp=trp, n=nrow(sm),
       n_ad=sum(sel$Module=="AD"), n_sig=sum(sel$Module=="SIG"), n_pool=sum(sel$Module=="POOL"))
}
predict_panel <- function(fit, SM, cols){
  sm <- SM[fit$sel$Ratio, cols, drop=FALSE]
  pc <- predict(fit$pca, t(sm))[,1:fit$np,drop=FALSE]
  as.character(predict(fit$lda, as.data.frame(pc) %>% setNames(paste0("PC",1:fit$np)))$class)
}

# ---- [4] CACHE rank+delta ONCE per fold (cheap: column subset of cached SM) -
build_caches <- function(SM, MOD, design){
  if (VERBOSE) cat(sprintf("\n[caches:%s] ranking full + per-fold (column subsets)...\n", design))
  full <- rank_and_delta(SM, MOD, seq_len(ncol(SM)))
  folds <- list(); k<-0L
  for (s in levels(study)){
    cols <- which(study != s)
    if (length(unique(phenotype[cols])) < 2) next
    k <- k+1L
    if (VERBOSE) cat(sprintf("  [%s] fold %d (hold %s)\n", design, k, s))
    folds[[s]] <- rank_and_delta(SM, MOD, cols)
  }
  list(full=full, folds=folds)
}

# ---- [5] SWEEP -------------------------------------------------------------
LEV <- levels(phenotype)
run_design <- function(SM, MOD, design){
  ch <- build_caches(SM, MOD, design)
  mfloor <- if (design=="ADSIG") MODULE_FLOOR else 0L
  sweep_rows <- list(); persample_rows <- list()
  for (n_total in TARGET_SWEEP){
    npp <- n_total %/% N_PHENO
    sel_full <- greedy_quota(ch$full, npp, module_floor=mfloor)
    full <- fit_panel(SM, sel_full, seq_len(ncol(SM)))
    loocv_cv <- suppressWarnings(lda(full$trp, grouping=phenotype, CV=TRUE))
    loocv_acc <- mean(loocv_cv$class==phenotype)*100
    loso_pred <- rep(NA_character_, ncol(SM)); names(loso_pred) <- colnames(SM)
    for (s in names(ch$folds)){
      cols_tr <- which(study != s); cols_te <- which(study == s)
      sel_s <- greedy_quota(ch$folds[[s]], npp, module_floor=mfloor)
      fit_s <- tryCatch(fit_panel(SM, sel_s, cols_tr), error=function(e) NULL)
      if (is.null(fit_s)) next
      loso_pred[cols_te] <- tryCatch(predict_panel(fit_s, SM, cols_te),
                                     error=function(e) rep(NA,length(cols_te)))
    }
    loso_acc <- mean(loso_pred==as.character(phenotype), na.rm=TRUE)*100
    f1_loso  <- mean(per_class(loso_pred, phenotype, LEV)$F1)
    cat(sprintf(" [%s] N=%d (ratios=%d AD=%d SIG=%d POOL=%d) PCs=%d | LOOCV=%.1f%% LOSO=%.1f%% F1=%.3f\n",
                design, n_total, full$n, full$n_ad, full$n_sig, full$n_pool, full$np,
                loocv_acc, loso_acc, f1_loso))
    sweep_rows[[length(sweep_rows)+1]] <- data.frame(Design=design, N_target=n_total,
                                                     N_ratios=full$n, PCs=full$np, AD=full$n_ad, SIG=full$n_sig, POOL=full$n_pool,
                                                     LOOCV=round(loocv_acc,2), LOSO=round(loso_acc,2), LOSO_macroF1=round(f1_loso,3))
    persample_rows[[length(persample_rows)+1]] <- data.frame(Design=design, N_target=n_total,
                                                             SampleID=colnames(SM), Study=as.character(study), True=as.character(phenotype),
                                                             LOOCV_pred=as.character(loocv_cv$class), LOSO_pred=loso_pred)
  }
  big <- greedy_quota(ch$full, MAX_NPP, module_floor=mfloor)
  fname <- if (design=="POOL") "Final_Ratios_POOL_BALANCED.csv" else "Final_Ratios_ADSIG_BALANCED.csv"
  write.csv(big %>% select(Gene_A,Gene_B,Ratio,TSP_Delta,Phenotype,Module),
            file.path(out_dir,fname), row.names=FALSE)
  list(sweep=bind_rows(sweep_rows), persample=bind_rows(persample_rows))
}

res_pool  <- run_design(SM_POOL,  MOD_POOL,  "POOL")
res_adsig <- run_design(SM_ADSIG, MOD_ADSIG, "ADSIG")
perf_all  <- bind_rows(res_pool$sweep, res_adsig$sweep)
write.csv(perf_all, file.path(out_dir,"PanelSweep_LOOCV_LOSO.csv"), row.names=FALSE)
write.csv(bind_rows(res_pool$persample,res_adsig$persample),
          file.path(out_dir,"PanelSweep_perSample.csv"), row.names=FALSE)
cat("\n=== SWEEP SUMMARY (balanced n/6, MAX_PCS=",MAX_PCS,", NESTED LOSO) ===\n",sep=""); print(perf_all)

# ---- [6] FIGURE ------------------------------------------------------------
plt <- perf_all %>% pivot_longer(c(LOOCV,LOSO), names_to="Scheme", values_to="Acc")
p <- ggplot(plt, aes(N_ratios, Acc, colour=Design, linetype=Scheme)) +
  geom_line(linewidth=1)+geom_point(size=2.4)+
  scale_colour_manual(values=c(POOL="#333333", ADSIG="#377EB8"))+
  scale_x_continuous(breaks=TARGET_SWEEP)+theme_classic(base_size=14)+
  theme(plot.title=element_text(face="bold",hjust=.5), legend.position="bottom")+
  labs(title=sprintf("LOOCV (solid) vs NESTED LOSO (dashed) | balanced n/6, PC cap=%d",MAX_PCS),
       x="Number of sign-encoded ratios", y="Accuracy (%)")
ggsave(file.path(out_dir,"PanelSweep_LOOCV_vs_LOSO.png"), p, width=9.5,height=5.5,dpi=300)
cat("\nDONE ->", out_dir, "\n")
