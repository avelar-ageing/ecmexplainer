# ==============================================================================
# ECMExplainer -- 003_FinalTrain_TwoPanels.R
# Trains TWO final models chosen from the NESTED-LOSO sweep:
#   * POOL-120  : unrestricted pool, 20/phenotype  (balanced, simplest)
#   * ADSIG-180 : AD U SIG pooled, 30/phenotype     (best LOSO / macro-F1)
# For EACH model: saturation curve, 2D PCA+LDA(+module views), 3D LDA,
#   mean-sign heatmaps, per-class metrics, saved .rds payload.
# Figure styling, palettes and plot fns carried over verbatim from
#   002b-UNIFIED-144.R so output matches previous versions.
# NOTE: final models are trained on ALL 96 samples (LOOCV reported here);
#   the cross-study LOSO estimates live in the sweep table, not here.
# ==============================================================================
suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(pheatmap); library(scatterplot3d)
})
filter <- dplyr::filter; select <- dplyr::select; pheatmap <- pheatmap::pheatmap
set.seed(1)

out_dir <- "Salmon_kTSP_FINAL_TwoPanels"          # <<< NEW output directory
dir.create(out_dir, showWarnings = FALSE)

PUB_COLORS <- c(Fibrosis="#377EB8", Hypoxia="#E41A1C", Inflammation="#4DAF4A",
                Proliferating="#FF7F00", Quiescence="#A65628", Senescence="#333333")
STUDY_SHAPES <- c(GSE102674=15, GSE117337=16, GSE117444=17, GSE139963=18, GSE221693=19,
                  GSE281533=20, GSE295622=2, GSE63577=5, GSE74324=6, GSE75643=8)

MAX_PCS <- 10L      # PC cap carried from the sweep (prevents the 37-PC overfit)
ETA_MIN <- 0.10

# ---- [1] LOAD + FILTER (identical to 002b-UNIFIED-144 / 002_V5) --------------
cat("\n[1] Loading foundation...\n")
salmon_counts <- readRDS("salmon_counts_Combined_FINAL.rds")
meta <- readr::read_csv("salmon_metadata_Combined_FINAL.csv", show_col_types=FALSE)
meta <- meta[match(colnames(salmon_counts), meta$SampleID), ]
stopifnot(all(meta$SampleID == colnames(salmon_counts)))
to_logcpm <- function(mat) log2(sweep(mat,2,colSums(mat)/1e6,"/")+1)
all_genes_raw <- to_logcpm(salmon_counts)
phenotype <- factor(meta$Phenotype); study_factor <- factor(meta$Study)
cat(" Matrix:",nrow(all_genes_raw),"x",ncol(all_genes_raw),
    "| studies:",nlevels(study_factor),"\n"); print(table(phenotype))

cpm_mat <- sweep(salmon_counts,2,colSums(salmon_counts)/1e6,"/")
min_samp <- min(table(phenotype))
expressed <- rownames(salmon_counts)[ rowSums(cpm_mat>1) >= min_samp ]
all_genes_raw <- all_genes_raw[expressed,,drop=FALSE]
cat(" Low-expr filter (CPM>1 in >=",min_samp,"samples): kept",length(expressed),"genes\n")

ecm_master <- read.csv("ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors=FALSE)
ecm_master <- ecm_master[!duplicated(ecm_master$Gene),]
fn_col <- intersect(c("Function","ECMExplained_Function","ECMExplainer_Function"), names(ecm_master))[1]
stopifnot("No function-label column found"= !is.na(fn_col))
ecm_master$FN <- toupper(trimws(ecm_master[[fn_col]]))

ecm_in <- intersect(ecm_master$Gene, rownames(all_genes_raw))
train_expr <- all_genes_raw[ecm_in,,drop=FALSE]
train_expr <- train_expr[apply(train_expr,1,var)>0.005,,drop=FALSE]
a_genes <- intersect(ecm_master$Gene[ecm_master$FN=="ASSEMBLY"],    rownames(train_expr))
d_genes <- intersect(ecm_master$Gene[ecm_master$FN=="DESTRUCTION"], rownames(train_expr))
s_genes <- intersect(ecm_master$Gene[ecm_master$FN=="SIGNALLING"],  rownames(train_expr))
pool_genes <- rownames(train_expr)
cat(" pool:",length(pool_genes)," A:",length(a_genes)," D:",length(d_genes)," S:",length(s_genes),"\n")

# ---- [2] UTILITIES (verbatim from 002b-UNIFIED-144) -------------------------
build_ratio_mat <- function(ratio_list, expr){
  mat <- matrix(NA,length(ratio_list),ncol(expr),dimnames=list(ratio_list,colnames(expr)))
  for(i in seq_along(ratio_list)){ g<-trimws(unlist(strsplit(ratio_list[[i]],"/")))
  mat[i,]<-expr[g[1],]-expr[g[2],] }; mat }
build_sign_mat <- function(ratio_list, expr) sign(build_ratio_mat(ratio_list,expr))
compute_eta_sq_synergy <- function(A_mat,D_mat,phenotype){
  n<-ncol(A_mat); phen_idx<-split(seq_len(n),phenotype); group_ns<-sapply(phen_idx,length)
  gm_A<-sapply(phen_idx,function(i)rowMeans(A_mat[,i,drop=FALSE]))
  gm_D<-sapply(phen_idx,function(i)rowMeans(D_mat[,i,drop=FALSE]))
  mu_A<-rowMeans(A_mat);mu_D<-rowMeans(D_mat);dA<-gm_A-mu_A;dD<-gm_D-mu_D
  ssb_A<-rowSums(sweep(dA^2,2,group_ns,"*"));ssb_D<-rowSums(sweep(dD^2,2,group_ns,"*"))
  cross<-sweep(dA,2,group_ns,"*")%*%t(dD); SS_bet<-outer(ssb_A,ssb_D,"+")-2*cross
  var_A<-apply(A_mat,1,var);var_D<-apply(D_mat,1,var)
  cov_AD<-((A_mat-mu_A)%*%t(D_mat-mu_D))/(n-1); SS_tot<-(outer(var_A,var_D,"+")-2*cov_AD)*(n-1)
  e_r<-pmin(pmax(SS_bet/SS_tot,0),1); e_A<-pmin(pmax(ssb_A/(var_A*(n-1)),0),1)
  e_D<-pmin(pmax(ssb_D/(var_D*(n-1)),0),1)
  e_r[!is.finite(e_r)]<-0; e_A[!is.finite(e_A)]<-0; e_D[!is.finite(e_D)]<-0
  list(ratio=e_r,A=e_A,D=e_D,A_names=rownames(A_mat),D_names=rownames(D_mat)) }
make_ranked_loose <- function(syn, sym){
  df<-data.frame(Gene_A=rep(syn$A_names,times=length(syn$D_names)),
                 Gene_B=rep(syn$D_names,each=length(syn$A_names)),
                 eta_sq=as.vector(syn$ratio),
                 max_single=pmax(rep(syn$A,times=length(syn$D_names)),
                                 rep(syn$D,each=length(syn$A_names))),
                 stringsAsFactors=FALSE)%>%mutate(synergy_delta=eta_sq-max_single)
  if(sym) df<-filter(df,Gene_A<Gene_B)
  filter(df,eta_sq>ETA_MIN)%>%arrange(desc(eta_sq)) }
compute_tsp_delta <- function(ranked_df,expr,pheno,panel_label){
  rlist<-paste0(ranked_df$Gene_A," / ",ranked_df$Gene_B)
  sm<-build_sign_mat(rlist,expr);up<-sm>0;rows<-list()
  for(p in levels(pheno)){ in_p<-pheno==p
  p_in<-rowMeans(up[,in_p,drop=FALSE]);p_out<-rowMeans(up[,!in_p,drop=FALSE]);delta<-p_in-p_out
  rows[[length(rows)+1]]<-data.frame(Phenotype=p,Panel=panel_label,
                                     Gene_A=ranked_df$Gene_A,Gene_B=ranked_df$Gene_B,Ratio=rlist,TSP_Delta=abs(delta),
                                     Direction=ifelse(delta>0,"Up",ifelse(delta<0,"Down","Flat")),stringsAsFactors=FALSE) }
  bind_rows(rows) }
greedy_stratified <- function(delta_df,n_per_pheno,max_appear=1){
  delta_df<-arrange(delta_df,desc(TSP_Delta));claimed<-character(0)
  gene_ct<-new.env(hash=TRUE,parent=emptyenv())
  counts<-setNames(rep(0L,nlevels(factor(delta_df$Phenotype))),levels(factor(delta_df$Phenotype)))
  keep<-logical(nrow(delta_df))
  for(i in seq_len(nrow(delta_df))){ r<-delta_df$Ratio[i];p<-as.character(delta_df$Phenotype[i])
  a<-delta_df$Gene_A[i];b<-delta_df$Gene_B[i]
  if(r%in%claimed)next; if(counts[p]>=n_per_pheno)next
  ca<-if(exists(a,envir=gene_ct))get(a,envir=gene_ct)else 0L
  cb<-if(exists(b,envir=gene_ct))get(b,envir=gene_ct)else 0L
  if(ca>=max_appear||cb>=max_appear)next
  keep[i]<-TRUE;claimed<-c(claimed,r);counts[p]<-counts[p]+1L
  assign(a,ca+1L,envir=gene_ct);assign(b,cb+1L,envir=gene_ct) }
  cat("   slots filled (target",n_per_pheno,"):\n");print(counts)
  delta_df[keep,,drop=FALSE]%>%distinct(Ratio,.keep_all=TRUE) }
get_npc <- function(pca, cap=MAX_PCS){ cv<-cumsum(pca$sdev^2/sum(pca$sdev^2)); idx<-which(cv>=0.95)
npc<-max(1L,if(length(idx)==0)ncol(pca$x)else idx[1]); min(npc,cap,ncol(pca$x)) }

# ---- [3] BUILD THE TWO CANDIDATE RANKED POOLS (full data) -------------------
cat("\n[2] Building candidate pools...\n")
# POOL: unrestricted symmetric
syn_pool   <- compute_eta_sq_synergy(train_expr[pool_genes,], train_expr[pool_genes,], phenotype)
ranked_POOL <- make_ranked_loose(syn_pool, sym=TRUE) %>% mutate(Category="POOL")
# ADSIG: AD (asym) U SIG (sym), pooled, ranked together
syn_ad  <- compute_eta_sq_synergy(train_expr[a_genes,], train_expr[d_genes,], phenotype)
syn_sig <- compute_eta_sq_synergy(train_expr[s_genes,], train_expr[s_genes,], phenotype)
ranked_ADSIG <- bind_rows(make_ranked_loose(syn_ad,FALSE) %>% mutate(Category="AD"),
                          make_ranked_loose(syn_sig,TRUE) %>% mutate(Category="SIG")) %>%
  arrange(desc(eta_sq))
cat("  POOL pairs:",nrow(ranked_POOL)," | ADSIG pairs:",nrow(ranked_ADSIG),
    "(AD",sum(ranked_ADSIG$Category=="AD"),"SIG",sum(ranked_ADSIG$Category=="SIG"),")\n")

# ---- [4] PER-MODEL TRAIN + ALL FIGURES --------------------------------------
th2d <- function() theme_classic(base_size=16)+
  theme(plot.title=element_text(face="bold",size=17,hjust=0.5),
        axis.title=element_text(face="bold"),legend.title=element_text(face="bold"))
mk2d <- function(coords, meta_df, xl, yl, title, var=NULL){
  df <- data.frame(coords[,1:2,drop=FALSE], meta_df); colnames(df)[1:2] <- c("D1","D2")
  df$Phenotype <- factor(df$Phenotype, levels=names(PUB_COLORS))
  xl2 <- if(!is.null(var)) paste0(xl," (",var[1],"%)") else xl
  yl2 <- if(!is.null(var)) paste0(yl," (",var[2],"%)") else yl
  # keep ellipse only for phenotypes with >=4 points (stat_ellipse needs >=3)
  ok <- names(which(table(df$Phenotype) >= 4))
  ell_df <- df[df$Phenotype %in% ok, , drop=FALSE]
  ggplot(df, aes(D1, D2, colour=Phenotype)) +
    stat_ellipse(data=ell_df, aes(fill=Phenotype), geom="polygon",
                 type="norm", level=0.85, alpha=0.12, colour=NA, show.legend=FALSE) +
    stat_ellipse(data=ell_df, aes(colour=Phenotype),
                 type="norm", level=0.85, linewidth=0.7, show.legend=FALSE) +
    geom_point(aes(shape=Study), size=3.4, alpha=0.9, stroke=0.4) +
    scale_colour_manual(values=PUB_COLORS, name="Phenotype", drop=FALSE) +
    scale_fill_manual(values=PUB_COLORS, guide="none", drop=FALSE) +
    scale_shape_manual(values=STUDY_SHAPES, name="Study") +
    guides(colour=guide_legend(order=1, override.aes=list(shape=16, size=3.6, alpha=1)),
           shape =guide_legend(order=2, override.aes=list(colour="grey30"), ncol=1)) +
    th2d() + labs(title=title, x=xl2, y=yl2)
}
HEAT_COLS <- colorRampPalette(c("#313695","white","#a50026"))(100)
render_heat <- function(rl,sel_ph,fname,lab,npp){
  z<-t(apply(sign(build_ratio_mat(rl,train_expr)),1,function(x)tapply(x,phenotype,mean)))
  ra<-data.frame(Selects=sel_ph);rownames(ra)<-rl
  ord<-order(factor(ra$Selects,levels=levels(phenotype)));z<-z[ord,,drop=FALSE];ra<-ra[ord,,drop=FALSE]
  ca<-data.frame(Phenotype=colnames(z));rownames(ca)<-colnames(z)
  while(!is.null(dev.list())) dev.off()
  png(fname,width=10,height=12,units="in",res=300)
  pheatmap::pheatmap(z,annotation_col=ca,annotation_row=ra,
                     annotation_colors=list(Phenotype=PUB_COLORS,Selects=PUB_COLORS),
                     cluster_cols=TRUE,cluster_rows=FALSE,scale="none",color=HEAT_COLS,
                     breaks=seq(-1,1,length.out=101),border_color="white",fontsize_row=7,fontsize_col=12,
                     gaps_row=cumsum(rle(as.character(ra$Selects))$lengths),
                     main=sprintf("%s (mean sign -1..+1 | %d/phenotype | max_appear=1)",lab,npp))
  dev.off(); cat("   OK",basename(fname),"\n") }

train_one <- function(model_tag, ranked_df, n_per_pheno, has_modules){
  PANEL <- n_per_pheno*6L
  cat("\n",strrep("=",70),"\n  TRAIN ",model_tag," (",n_per_pheno,"/pheno, target ",PANEL,")\n",strrep("=",70),"\n",sep="")
  cat_lookup <- setNames(ranked_df$Category, paste0(ranked_df$Gene_A," / ",ranked_df$Gene_B))
  sel <- greedy_stratified(compute_tsp_delta(ranked_df, train_expr, phenotype, model_tag), n_per_pheno, 1)
  sel$Category <- cat_lookup[sel$Ratio]
  stopifnot("gene in >1 ratio"=max(table(c(sel$Gene_A,sel$Gene_B)))<=1)
  final_features <- sel %>% select(Gene_A,Gene_B,Ratio,TSP_Delta,
                                   Selecting_Phenotype=Phenotype,Direction,Category)
  cat("   selected:",nrow(final_features),
      if(has_modules) sprintf("(AD=%d SIG=%d)",sum(final_features$Category=="AD"),
                              sum(final_features$Category=="SIG")) else "","\n")
  ratio_list_all <- final_features$Ratio
  
  # ---- saturation curve (single PCA, LOOCV, PC cap) ----
  calc_step <- function(rl,n_feat){
    mc<-build_sign_mat(rl,train_expr); keep<-apply(mc,1,function(x)length(unique(x))>1)
    if(!all(keep))mc<-mc[keep,,drop=FALSE]
    if(nrow(mc)<2)return(data.frame(Features=n_feat,Accuracy=NA,CumVar=NA,UniqueGenes=NA))
    pca<-prcomp(t(mc),center=TRUE,scale.=TRUE); npc<-min(get_npc(pca),ncol(mc)-1L)
    cv<-suppressWarnings(lda(as.data.frame(pca$x[,1:npc,drop=FALSE]),grouping=phenotype,CV=TRUE))
    acc<-mean(cv$class==phenotype)*100
    Y<-as.numeric(phenotype);SST<-var(Y)*(length(Y)-1);ve<-(SST-sum(lm(Y~t(mc))$residuals^2))/SST*100
    ug<-length(unique(c(sub(" /.*","",rownames(mc)),sub(".* / ","",rownames(mc)))))
    data.frame(Features=n_feat,Accuracy=acc,CumVar=ve,UniqueGenes=ug) }
  perf <- bind_rows(lapply(seq(6,length(ratio_list_all),by=6), function(n) calc_step(ratio_list_all[1:n],n)))
  write.csv(perf, file.path(out_dir,sprintf("Saturation_%s_%d.csv",model_tag,PANEL)), row.names=FALSE)
  p_sat <- ggplot(perf,aes(Features,Accuracy))+geom_line(linewidth=1,color="#333333")+
    geom_point(size=2.5,color="#377EB8")+theme_classic(base_size=14)+
    theme(plot.title=element_text(face="bold",hjust=0.5))+
    labs(title=sprintf("%s -- LOOCV saturation (%d, PC cap %d)",model_tag,PANEL,MAX_PCS),
         x="Number of sign-encoded ratios",y="LOOCV accuracy (%)")
  ggsave(file.path(out_dir,sprintf("Saturation_%s_%d.png",model_tag,PANEL)),p_sat,width=9,height=5.5,dpi=300)
  
  # ---- final sign -> ONE PCA -> LDA ----
  ratio_mat_sign <- sign(build_ratio_mat(final_features$Ratio, train_expr))
  keep_ratio <- apply(ratio_mat_sign,1,function(x)length(unique(x))>1)
  ratio_mat_sign <- ratio_mat_sign[keep_ratio,,drop=FALSE]; final_features <- final_features[keep_ratio,,drop=FALSE]
  stopifnot(identical(rownames(ratio_mat_sign),final_features$Ratio))
  pca_all <- prcomp(t(ratio_mat_sign),center=TRUE,scale.=TRUE); n_pcs <- get_npc(pca_all)
  train_pcs <- as.data.frame(pca_all$x[,1:n_pcs,drop=FALSE]); colnames(train_pcs)<-paste0("PC",1:n_pcs)
  lda_model <- lda(phenotype~.,data=train_pcs)
  lda_coords <- predict(lda_model,train_pcs)$x
  lda_var <- round(100*lda_model$svd^2/sum(lda_model$svd^2),1)
  lda_cv <- suppressWarnings(lda(train_pcs,grouping=phenotype,CV=TRUE))
  loocv <- mean(lda_cv$class==phenotype)*100
  cm <- table(Predicted=lda_cv$class,True=phenotype)
  prec<-diag(cm)/pmax(colSums(t(cm)),1);rec<-diag(cm)/pmax(colSums(cm),1)
  f1<-ifelse((prec+rec)>0,2*prec*rec/(prec+rec),0);macro_f1<-mean(f1)
  cat(sprintf("   LOOCV=%.2f%% macro-F1=%.3f PCs=%d\n",loocv,macro_f1,n_pcs)); print(cm)
  write.csv(data.frame(Phenotype=levels(phenotype),Precision=round(prec,3),Recall=round(rec,3),F1=round(f1,3)),
            file.path(out_dir,sprintf("PerClass_%s_%d.csv",model_tag,PANEL)),row.names=FALSE)
  write.csv(final_features, file.path(out_dir,sprintf("Final_Ratios_%s_%d.csv",model_tag,PANEL)),row.names=FALSE)
  # ── SAVE PER-SAMPLE SIGN MATRIX (Supp Table 3 second tab) ───────────────────
  sign_export <- as.data.frame(t(ratio_mat_sign))
  sign_export$SampleID <- rownames(sign_export)
  sign_export <- sign_export[, c("SampleID", setdiff(names(sign_export), "SampleID"))]
  write.csv(sign_export,
            file.path(out_dir, sprintf("PerSample_SignMatrix_%s_%d.csv", model_tag, PANEL)),
            row.names = FALSE)
  cat(" Per-sample sign matrix saved:",
      sprintf("PerSample_SignMatrix_%s_%d.csv (%d ratios x %d samples)\n",
              model_tag, PANEL, nrow(ratio_mat_sign), ncol(ratio_mat_sign)))
  
  # ── SAVE CONFUSION MATRIX (Supp Table 5 Tab C support) ───────────────────────
  cm_df <- as.data.frame(cm)
  write.csv(cm_df,
            file.path(out_dir, sprintf("ConfusionMatrix_%s_%d.csv", model_tag, PANEL)),
            row.names = FALSE)
  cat(" Confusion matrix saved:", sprintf("ConfusionMatrix_%s_%d.csv\n", model_tag, PANEL))
  
  # ---- 2D figures: global PCA + LDA (+ AD/SIG views if modular) ----
  pca_var <- round(100*pca_all$sdev^2/sum(pca_all$sdev^2),1)
  ph_df <- data.frame(Phenotype=as.character(phenotype),Study=meta$Study)
  pg <- mk2d(pca_all$x,ph_df,"PC1","PC2",sprintf("%s PCA -- all ratios",model_tag),pca_var)
  pl <- mk2d(lda_coords,ph_df,"LD1","LD2",sprintf("%s LDA (LOOCV %.1f%%)",model_tag,loocv),lda_var)
  if(has_modules){
    mad<-t(ratio_mat_sign[final_features$Category=="AD",,drop=FALSE])
    msig<-t(ratio_mat_sign[final_features$Category=="SIG",,drop=FALSE])
    pca_ad<-prcomp(mad,center=TRUE,scale.=TRUE); pca_sig<-prcomp(msig,center=TRUE,scale.=TRUE)
    pa<-mk2d(pca_ad$x,ph_df,"PC1","PC2",sprintf("AD view (%d)",sum(final_features$Category=="AD")),
             round(100*pca_ad$sdev^2/sum(pca_ad$sdev^2),1))
    ps<-mk2d(pca_sig$x,ph_df,"PC1","PC2",sprintf("SIG view (%d)",sum(final_features$Category=="SIG")),
             round(100*pca_sig$sdev^2/sum(pca_sig$sdev^2),1))
    combo<-(pg|pl)/(pa|ps)+plot_layout(guides="collect")
    H<-14
  } else { combo<-(pg|pl)+plot_layout(guides="collect"); H<-7 }
  ggsave(file.path(out_dir,sprintf("2D_Combined_%s_%d.png",model_tag,PANEL)),
         combo + plot_annotation(
           title=sprintf("ECMExplainer %s (N=96, %d ratios)",model_tag,PANEL),
           theme=theme(plot.title=element_text(size=16,face="bold",hjust=0.5))) &
           theme(legend.position="bottom", legend.box="horizontal"),
         width=if(has_modules) 10 else 9, height=if(has_modules) 14 else 4.6, dpi=600)
  
  # ---- 3D LDA with centroids ----
  ld <- lda_coords; ldv <- lda_var
  if(ncol(ld)>=3){
    ph<-factor(phenotype,levels=names(PUB_COLORS)); cols<-PUB_COLORS[as.character(ph)]
    png(file.path(out_dir,sprintf("3D_LDA_%s_%d.png",model_tag,PANEL)),width=9,height=8,units="in",res=300)
    s3d<-scatterplot3d(ld[,1],ld[,2],ld[,3],color=cols,pch=19,cex.symbols=1.5,angle=55,box=TRUE,
                       xlab=sprintf("LD1 (%.1f%%)",ldv[1]),ylab=sprintf("LD2 (%.1f%%)",ldv[2]),zlab=sprintf("LD3 (%.1f%%)",ldv[3]),
                       main=sprintf("%s %d-ratio -- 3D LDA (LOOCV %.1f%%)",model_tag,PANEL,loocv))
    cen<-t(sapply(levels(ph),function(p)colMeans(ld[ph==p,1:3,drop=FALSE])))
    s3d$points3d(cen[,1],cen[,2],cen[,3],pch=8,cex=2.2,col="black")
    legend("topright",legend=names(PUB_COLORS),col=PUB_COLORS,pch=19,bty="n",cex=.9)
    dev.off(); cat("   OK 3D LDA\n")
  } else cat("   (LDA <3 dims; 3D skipped)\n")
  
  # ---- heatmaps ----
  render_heat(final_features$Ratio, final_features$Selecting_Phenotype,
              file.path(out_dir,sprintf("Heatmap_ALL_%s_%d.png",model_tag,PANEL)),
              sprintf("%s panel",model_tag), n_per_pheno)
  if(has_modules){
    render_heat(final_features$Ratio[final_features$Category=="AD"],
                final_features$Selecting_Phenotype[final_features$Category=="AD"],
                file.path(out_dir,sprintf("Heatmap_AD_%s_%d.png",model_tag,PANEL)),"AD ratios",n_per_pheno)
    render_heat(final_features$Ratio[final_features$Category=="SIG"],
                final_features$Selecting_Phenotype[final_features$Category=="SIG"],
                file.path(out_dir,sprintf("Heatmap_SIG_%s_%d.png",model_tag,PANEL)),"SIG ratios",n_per_pheno)
  }
  
  # ---- save model payload (matches 002b contract) ----
  ff<-final_features; names(ff)[names(ff)=="Ratio"]<-"Ratio_ID"; names(ff)[names(ff)=="Category"]<-"Module"
  sm<-ratio_mat_sign; rownames(sm)<-ff$Ratio_ID
  saveRDS(list(encoding="sign",design=model_tag,max_appear=1L,n_training=96L,panel_size=PANEL,
               phenotypes=levels(phenotype),train_phenotype=phenotype,features=ff,
               sign_matrix_train=sm,pca=pca_all,n_pcs=n_pcs,lda=lda_model,
               loocv_acc=loocv,macro_f1=macro_f1,loocv_class=lda_cv$class),
          file.path(out_dir,sprintf("TrainedModel_%s_%d.rds",model_tag,PANEL)))
  cat("   model saved.\n")
  
  
  # ── SAVE LDA COEFFICIENTS ─────────────────────────────────────────────────────
  # lda_model$scaling: rows = PC axes, cols = LD axes  (discriminant loadings)
  # We want the contribution per RATIO, not per PC.
  # Project via PCA rotation: ratio_loadings = pca_all$rotation[, 1:n_pcs] %*% lda_model$scaling
  
  ratio_to_pca <- pca_all$rotation[, 1:n_pcs, drop = FALSE]   # ratios x PCs
  ld_per_ratio  <- ratio_to_pca %*% lda_model$scaling          # ratios x LDs
  
  coef_df <- as.data.frame(ld_per_ratio)
  coef_df$Ratio <- rownames(coef_df)
  
  # Join phenotype metadata
  coef_df <- coef_df %>%
    left_join(final_features %>% select(Ratio, Selecting_Phenotype, Category, TSP_Delta),
              by = "Ratio") %>%
    mutate(Abs_LD1 = abs(LD1),
           Abs_LD2 = if ("LD2" %in% names(.)) abs(LD2) else NA_real_) %>%
    arrange(desc(Abs_LD1))
  
  write.csv(coef_df,
            file.path(out_dir, sprintf("LDA_Coefficients_%s_%d.csv", model_tag, PANEL)),
            row.names = FALSE)
  cat(" LDA coefficients saved:", sprintf("LDA_Coefficients_%s_%d.csv\n", model_tag, PANEL))
  
  # ── PLOT: distributed weight profile (Figure 1H) ─────────────────────────────
  p_weights <- coef_df %>%
    mutate(Ratio = forcats::fct_reorder(Ratio, Abs_LD1)) %>%
    ggplot(aes(x = Abs_LD1, y = Ratio, fill = Selecting_Phenotype)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = PUB_COLORS) +
    labs(title  = sprintf("%s: LD1 weight distribution across %d ratios", model_tag, PANEL),
         x      = "|LD1 coefficient|",
         y      = NULL,
         fill   = "Selecting phenotype") +
    theme_classic(base_size = 9) +
    theme(axis.text.y    = element_text(size = 3.5),
          legend.position = "right",
          plot.title      = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(out_dir, sprintf("Fig1H_LDA_Weights_%s_%d.png", model_tag, PANEL)),
         p_weights, width = 8, height = 12, dpi = 300)
  cat(" Fig1H weight plot saved.\n")
  
  data.frame(Model=model_tag,Panel=PANEL,N_ratios=nrow(final_features),PCs=n_pcs,
             LOOCV=round(loocv,2),macroF1=round(macro_f1,3),
             AD=sum(final_features$Category=="AD"),SIG=sum(final_features$Category=="SIG"),
             POOL=sum(final_features$Category=="POOL"))
}

# ---- [5] RUN BOTH PANELS ----------------------------------------------------
summ <- bind_rows(
  train_one("POOL",  ranked_POOL,  n_per_pheno=20L, has_modules=FALSE),  # POOL-120
  train_one("ADSIG", ranked_ADSIG, n_per_pheno=30L, has_modules=TRUE)    # ADSIG-180
)
write.csv(summ, file.path(out_dir,"FINAL_TwoPanel_Summary.csv"), row.names=FALSE)
cat("\n=== FINAL TWO-PANEL SUMMARY ===\n"); print(summ)
cat("\nDONE ->", out_dir, "\n")




# ── SAVE PER-SAMPLE SIGN MATRIX (Supp Table 3 second tab) ───────────────────
# ratio_mat_sign already exists in this scope: rows = ratios, cols = 96 samples
sign_export <- as.data.frame(t(ratio_mat_sign))
sign_export$SampleID <- rownames(sign_export)
sign_export <- sign_export[, c("SampleID", setdiff(names(sign_export), "SampleID"))]

write.csv(sign_export,
          file.path(out_dir, sprintf("PerSample_SignMatrix_%s_%d.csv", model_tag, PANEL)),
          row.names = FALSE)
cat(" Per-sample sign matrix saved:",
    sprintf("PerSample_SignMatrix_%s_%d.csv (%d ratios x %d samples)\n",
            model_tag, PANEL, nrow(ratio_mat_sign), ncol(ratio_mat_sign)))

# ── SAVE FULL CONFUSION MATRIX (Supp Table 5 Tab C — precision/recall/F1) ───
# cm already exists in this scope: table(Predicted=lda_cv$class, True=phenotype)
cm_df <- as.data.frame(cm)
write.csv(cm_df,
          file.path(out_dir, sprintf("ConfusionMatrix_%s_%d.csv", model_tag, PANEL)),
          row.names = FALSE)

# Also save the per-class P/R/F1 table with a Model/Panel tag added, so all
# four models (once each is trained through this same function) can be
# rbind-ed into one Supp Table 5 tab without re-deriving anything
perclass_export <- data.frame(
  Model = model_tag,
  Panel_Size = PANEL,
  Phenotype = levels(phenotype),
  Precision = round(prec, 3),
  Recall = round(rec, 3),
  F1 = round(f1, 3)
)
write.csv(perclass_export,
          file.path(out_dir, sprintf("PerClass_PRF1_%s_%d.csv", model_tag, PANEL)),
          row.names = FALSE)
cat(" Confusion matrix + per-class P/R/F1 saved for", model_tag, "\n")

















# ============================================================
#  DUMBBELL PANEL — top-3 ratios x 6 phenotypes
#  A vs B effect size (Cohen's d) + dominance index
#  Self-contained: robust file find, gene-complete filtering,
#  no NaN, no duplicate factor levels.
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(readr)

# ---- 1. Robust inputs (survives -4 / -8 download suffixes) ----------------
long_path  <- list.files(pattern = "perSample_geneExpression_long.*\\.csv$",
                         recursive = TRUE, full.names = TRUE)[1]
ratio_path <- list.files(pattern = "Final_Ratios_ADSIG_BALANCED.*\\.csv$",
                         recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(long_path), !is.na(ratio_path))

long  <- read_csv(long_path,  show_col_types = FALSE)   # SampleID Gene log2CPM Phenotype Study
ratio <- read_csv(ratio_path, show_col_types = FALSE)   # Gene_A Gene_B Ratio TSP_Delta Phenotype Module

# ---- 2. Per-phenotype, per-gene mean expression ---------------------------
gene_means <- long %>%
  group_by(Phenotype, Gene) %>%
  summarise(mean_expr = mean(log2CPM, na.rm = TRUE),
            var_expr  = var(log2CPM,  na.rm = TRUE),
            n         = dplyr::n(), .groups = "drop")
genes_avail <- unique(gene_means$Gene)

# ---- 3. Keep only gene-complete pairs, THEN take top-3 per phenotype ------
top3 <- ratio %>%
  filter(Gene_A %in% genes_avail, Gene_B %in% genes_avail) %>%   # drop missing-gene pairs first
  group_by(Phenotype) %>%
  slice_max(order_by = Ratio, n = 3, with_ties = FALSE) %>%       # swap to TSP_Delta to rank by score
  ungroup() %>%
  mutate(pair = paste(Gene_A, Gene_B, sep = " / "))

# Report anything excluded so it's never silent
dropped <- setdiff(unique(c(ratio$Gene_A, ratio$Gene_B)), genes_avail)
if (length(dropped)) message("Genes absent from long table (pairs skipped): ",
                             paste(dropped, collapse = ", "))

# ---- 4. Effect size + dominance index (join-based, no scoping bug) --------
es <- top3 %>%
  left_join(gene_means, by = c("Phenotype", "Gene_A" = "Gene")) %>%
  rename(mean_A = mean_expr, var_A = var_expr, n_A = n) %>%
  left_join(gene_means, by = c("Phenotype", "Gene_B" = "Gene")) %>%
  rename(mean_B = mean_expr, var_B = var_expr, n_B = n) %>%
  filter(!is.na(mean_A), !is.na(mean_B)) %>%
  mutate(
    sd_pool   = sqrt(((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - 2)),
    cohens_d  = (mean_A - mean_B) / sd_pool,
    dominance = (mean_A - mean_B) / (mean_A + mean_B),
    pair_id   = row_number(),
    label     = paste0(pair, "  (d=", round(cohens_d, 2),
                       ", D=", round(dominance, 2), ")  #", pair_id)
  )

# ---- 5. Reshape to dumbbell endpoints -------------------------------------
plot_df <- es %>%
  mutate(label = fct_reorder(label, dominance)) %>%
  pivot_longer(c(mean_A, mean_B), names_to = "member", values_to = "expr") %>%
  mutate(member = recode(member, mean_A = "Gene A", mean_B = "Gene B"))

# ---- 6. Dumbbell plot -----------------------------------------------------
p_dumbbell <- ggplot(plot_df, aes(x = expr, y = label)) +
  geom_line(aes(group = label), colour = "grey70", linewidth = 1) +
  geom_point(aes(colour = member), size = 3) +
  facet_grid(Phenotype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c("Gene A" = "#1b9e77", "Gene B" = "#d95f02"),
                      name = NULL) +
  labs(x = "Mean log2CPM", y = NULL,
       title = "Top-3 ratio pairs per phenotype",
       subtitle = "A vs B effect size (Cohen's d) and dominance index (D)") +
  theme_minimal(base_size = 11) +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position    = "top")

print(p_dumbbell)

# ---- 7. (optional) assemble beside the heatmap + save ---------------------
# library(patchwork)
# combined <- p_heat + p_dumbbell + plot_layout(widths = c(1, 1.1))
# ggsave("Salmon_kTSP_FINAL_TwoPanels.pdf", combined, width = 14, height = 9)



library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(purrr)

targets <- tibble::tribble(
  ~Target,         ~Gene_A,   ~Gene_B,
  "Fibrosis",      "MFAP2",   "HGSNAT",
  "Hypoxia",       "BMP4",    "PLXNA1",
  "Inflammation",  "IL34",    "MUC1",
  "Proliferating", "C1QTNF7", "LGALS9",
  "Quiescence",    "FSTL3",   "PDGFD",
  "Senescence",    "CXCL12",  "IL6")

# 1. Build per-sample Δlog2 for each target ratio, across ALL phenotypes
#    pmap avoids rowwise scoping problems entirely.
ratio_long <- pmap_dfr(targets, function(Target, Gene_A, Gene_B) {
  long %>%
    filter(Gene %in% c(Gene_A, Gene_B)) %>%
    select(SampleID, samp_pheno = Phenotype, Gene, log2CPM) %>%
    pivot_wider(names_from = Gene, values_from = log2CPM) %>%
    transmute(SampleID, samp_pheno,
              delta = .data[[Gene_A]] - .data[[Gene_B]],
              Target = Target,
              pair   = paste0(Target, ": ", Gene_A, " / ", Gene_B))
}) %>%
  mutate(grp = if_else(samp_pheno == Target, "Target phenotype", "Other phenotypes"))

# 2. Collapse to target-vs-rest medians (one row per ratio)
seg <- ratio_long %>%
  group_by(pair, grp) %>%
  summarise(med = median(delta, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = grp, values_from = med) %>%
  mutate(separation = `Target phenotype` - `Other phenotypes`)

# 3. Diagnostic dumbbell — length of the bar = diagnostic strength
p_diag <- ggplot(seg, aes(y = fct_reorder(pair, separation))) +
  geom_segment(aes(x = `Other phenotypes`, xend = `Target phenotype`,
                   yend = pair), colour = "grey65", linewidth = 1.8) +
  geom_point(aes(x = `Other phenotypes`, colour = "Other phenotypes"), size = 3.5) +
  geom_point(aes(x = `Target phenotype`, colour = "Target phenotype"), size = 3.5) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c("Target phenotype"  = "#1b9e77",
                                 "Other phenotypes"   = "#9970ab"), name = NULL) +
  labs(x = expression("median "*Delta*" log"[2]*"  ( A - B )"), y = NULL,
       title = "Diagnostic value of each ratio for its target phenotype",
       subtitle = "Bar length = separation of target from all other states; dashed = sign-flip") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")

print(p_diag)
ggsave("FIG2_Dumbbell_diagnostic.pdf", p_diag, width = 8, height = 4.5)
ggsave("FIG2_Dumbbell_diagnostic.png", p_diag, width = 8, height = 4.5, dpi = 300)








# ============================================================================
# Supplementary Figure 3C — sample-level clustering, ADSIG-180 model
# Three panels: ALL 180 ratios, AD-only subset, SIG-only subset
# Tests whether samples cluster by phenotype (signal) or by study (batch)
# ============================================================================
library(pheatmap)

# ---- Load the ADSIG-180 final ratio panel (saved by train_one) -----------
final_features_adsig <- read.csv(
  file.path(out_dir, "Final_Ratios_ADSIG_180.csv"),
  stringsAsFactors = FALSE
)
stopifnot(nrow(final_features_adsig) == 180)
cat("Loaded ADSIG-180:", nrow(final_features_adsig), "ratios | AD =",
    sum(final_features_adsig$Category == "AD"), "| SIG =",
    sum(final_features_adsig$Category == "SIG"), "\n")

# ---- Core function: per-sample sign matrix, hierarchically clustered -----
render_sample_heat <- function(rl, sel_ph, expr, phenotype, meta,
                               fname, panel_label, legend_title = "matrix_4",
                               cluster_rows = TRUE, cluster_cols = TRUE) {
  
  sm <- sign(build_ratio_mat(rl, expr))               # rows=ratios, cols=samples, values in {-1,0,+1}
  stopifnot(identical(colnames(sm), colnames(expr)))
  
  col_anno <- data.frame(
    Phenotype = as.character(phenotype),
    Study     = as.character(meta$Study)
  )
  rownames(col_anno) <- colnames(sm)
  
  row_anno <- data.frame(Selects = sel_ph)
  rownames(row_anno) <- rl
  
  study_levels <- sort(unique(as.character(meta$Study)))
  anno_colors <- list(
    Phenotype = PUB_COLORS,
    Selects   = PUB_COLORS,
    Study     = setNames(
      colorRampPalette(c("#4B0082","#1F78B4","#33A02C","#FF7F00","#E31A1C",
                         "#B15928","#6A3D9A","#A6CEE3","#FDBF6F","#CAB2D6"))(length(study_levels)),
      study_levels
    )
  )
  
  HEAT_COLS <- colorRampPalette(c("#313695","white","#a50026"))(3)
  
  while (!is.null(dev.list())) dev.off()
  png(fname, width = 14, height = 12, units = "in", res = 300)
  pheatmap(
    sm,
    annotation_col    = col_anno,
    annotation_row    = row_anno,
    annotation_colors = anno_colors,
    cluster_rows      = cluster_rows,
    cluster_cols      = cluster_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    scale             = "none",
    color             = HEAT_COLS,
    breaks            = c(-1, -0.34, 0.34, 1),
    border_color      = NA,
    fontsize_row      = 5,
    fontsize_col      = 5,
    show_colnames     = FALSE,
    legend_breaks     = c(-1, 0, 1),
    legend_labels     = c("-1", "0", "+1"),
    main = sprintf("%s: sample clustering on %d sign-encoded ratios (%s)",
                   panel_label, nrow(sm), legend_title)
  )
  dev.off()
  cat(" OK ->", basename(fname), "\n")
}

# ---- Panel 1: ALL 180 ratios (global ADSIG panel) --------------------------
render_sample_heat(
  rl        = final_features_adsig$Ratio,
  sel_ph    = final_features_adsig$Selecting_Phenotype,
  expr      = train_expr,
  phenotype = phenotype,
  meta      = meta,
  fname     = file.path(out_dir, "SuppFig3C_SampleClustering_ADSIG180_ALL.png"),
  panel_label  = "ADSIG-180 (ALL)",
  legend_title = "matrix_4"
)

# ---- Panel 2: AD-only subset ------------------------------------------------
ad_ratios <- final_features_adsig[final_features_adsig$Category == "AD", ]
render_sample_heat(
  rl        = ad_ratios$Ratio,
  sel_ph    = ad_ratios$Selecting_Phenotype,
  expr      = train_expr,
  phenotype = phenotype,
  meta      = meta,
  fname     = file.path(out_dir, "SuppFig3C_SampleClustering_ADSIG180_AD.png"),
  panel_label  = sprintf("ADSIG-180 (AD only, n=%d)", nrow(ad_ratios)),
  legend_title = "matrix_4"
)

# ---- Panel 3: SIG-only subset -----------------------------------------------
sig_ratios <- final_features_adsig[final_features_adsig$Category == "SIG", ]
render_sample_heat(
  rl        = sig_ratios$Ratio,
  sel_ph    = sig_ratios$Selecting_Phenotype,
  expr      = train_expr,
  phenotype = phenotype,
  meta      = meta,
  fname     = file.path(out_dir, "SuppFig3C_SampleClustering_ADSIG180_SIG.png"),
  panel_label  = sprintf("ADSIG-180 (SIG only, n=%d)", nrow(sig_ratios)),
  legend_title = "matrix_4"
)

cat("\nDone — 3 sample-clustering heatmaps written to", out_dir, "\n")


