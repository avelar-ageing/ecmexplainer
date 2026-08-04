# ==============================================================================
# ECMExplainer -- Script 002c: A/D + SIG/SIG SATURATION CURVES (LOOCV + LOSO)
# ==============================================================================
# Question: how many sign-encoded ratios are actually needed?
# We sweep n = 1..N_PER_PHENO ratios PER PHENOTYPE PER MODULE and, at each step,
# record BOTH:
#   - LOOCV accuracy (optimistic ceiling; fixed panel, LDA-only CV)
#   - LOSO  accuracy (study-out; the honest generalisation estimate)
# for three panels: A/D module, SIG module, and Combined (A/D + SIG).
#
# IMPORTANT: ratio SELECTION uses all 96 samples (so the saturation ORDER is
# fixed and interpretable), but the LOSO estimate at each step REFITS PCA+LDA
# on the training studies only and predicts the held-out study -- so LOSO here
# reflects model transfer for a panel of that size, not selection leakage.
# (Selection-inside-fold is handled in 002b; here we isolate the size question.)
#
# Sources informing design: Saeb 2017 (subject-wise CV), Varoquaux 2017,
# Geman 2004 / Tan 2005 (rank-based TSP invariance).
# Output: Salmon_kTSP_Saturation_LOSO/
# ==============================================================================

suppressPackageStartupMessages({ library(MASS); library(dplyr); library(tidyr); library(ggplot2) })
filter <- dplyr::filter; set.seed(1)
out_dir <- "Salmon_kTSP_Saturation_LOSO"; dir.create(out_dir, showWarnings = FALSE)

MAX_PER_PHENO <- 12   # sweep 1..12 ratios / phenotype / module
ETA_MIN <- 0.10

# ----- load (same contract as 002b) -----
salmon_counts <- readRDS("salmon_counts_Combined_FINAL.rds")
meta <- readr::read_csv("salmon_metadata_Combined_FINAL.csv", show_col_types = FALSE)
meta <- meta[match(colnames(salmon_counts), meta$SampleID), ]
to_logcpm <- function(m) log2(sweep(m,2,colSums(m)/1e6,"/")+1)
all_logcpm <- to_logcpm(salmon_counts)
phenotype <- factor(meta$Phenotype); study <- factor(meta$Study)

ecm <- read.csv("ECMExplainer_FINAL_curated_genelist.csv", stringsAsFactors=FALSE)
ecm <- ecm[!duplicated(ecm$Gene),]
fn_col <- intersect(c("Function","ECMExplained_Function","ECMExplainer_Function"), names(ecm))[1]
ecm$FN <- toupper(trimws(ecm[[fn_col]]))

expr_full <- all_logcpm[apply(all_logcpm,1,var)>0.005,,drop=FALSE]
A <- intersect(ecm$Gene[ecm$FN=="ASSEMBLY"],    rownames(expr_full))
D <- intersect(ecm$Gene[ecm$FN=="DESTRUCTION"], rownames(expr_full))
S <- intersect(ecm$Gene[ecm$FN=="SIGNALLING"],  rownames(expr_full))
cat("A=",length(A)," D=",length(D)," SIG=",length(S),"\n")

# ---- shared helpers (same maths as 002 / 002b) ----
eta_syn <- function(Am,Dm,ph){ n<-ncol(Am); idx<-split(seq_len(n),ph); gn<-sapply(idx,length)
  gA<-sapply(idx,function(i)rowMeans(Am[,i,drop=FALSE])); gD<-sapply(idx,function(i)rowMeans(Dm[,i,drop=FALSE]))
  mA<-rowMeans(Am);mD<-rowMeans(Dm);dA<-gA-mA;dD<-gD-mD
  sA<-rowSums(sweep(dA^2,2,gn,"*"));sD<-rowSums(sweep(dD^2,2,gn,"*"));cr<-sweep(dA,2,gn,"*")%*%t(dD)
  SSb<-outer(sA,sD,"+")-2*cr; vA<-apply(Am,1,var);vD<-apply(Dm,1,var)
  cv<-((Am-mA)%*%t(Dm-mD))/(n-1); SSt<-(outer(vA,vD,"+")-2*cv)*(n-1)
  e<-pmin(pmax(SSb/SSt,0),1);e[!is.finite(e)]<-0
  list(r=e,An=rownames(Am),Dn=rownames(Dm)) }
rk <- function(s,sym){ df<-data.frame(A=rep(s$An,times=length(s$Dn)),B=rep(s$Dn,each=length(s$An)),
  e=as.vector(s$r),stringsAsFactors=FALSE); if(sym) df<-filter(df,A<B)
  filter(df,e>ETA_MIN)%>%arrange(desc(e)) }
bs <- function(rl,ex){ m<-matrix(0L,length(rl),ncol(ex),dimnames=list(rl,colnames(ex)))
  for(i in seq_along(rl)){g<-trimws(strsplit(rl[i]," / ")[[1]]); m[i,]<-as.integer(sign(ex[g[1],]-ex[g[2],]))}; m }
td <- function(r,ex,ph){ rl<-paste0(r$A," / ",r$B); sm<-bs(rl,ex);up<-sm>0
  do.call(rbind,lapply(levels(ph),function(p){pin<-rowMeans(up[,ph==p,drop=FALSE]);pout<-rowMeans(up[,ph!=p,drop=FALSE])
    data.frame(P=p,Ratio=rl,Delta=abs(pin-pout),stringsAsFactors=FALSE)})) }
grd <- function(d,np){ d<-arrange(d,desc(Delta));cl<-character(0)
  g<-new.env(hash=TRUE,parent=emptyenv());ct<-setNames(integer(nlevels(factor(d$P))),levels(factor(d$P)))
  k<-logical(nrow(d))
  for(i in seq_len(nrow(d))){r<-d$Ratio[i];p<-d$P[i];ab<-trimws(strsplit(r," / ")[[1]])
    if(r%in%cl)next; if(ct[p]>=np)next
    ca<-if(exists(ab[1],envir=g))get(ab[1],envir=g)else 0L; cb<-if(exists(ab[2],envir=g))get(ab[2],envir=g)else 0L
    if(ca>=1||cb>=1)next; k[i]<-TRUE;cl<-c(cl,r);ct[p]<-ct[p]+1L
    assign(ab[1],ca+1L,envir=g);assign(ab[2],cb+1L,envir=g)}
  d[k,,drop=FALSE]%>%distinct(Ratio,.keep_all=TRUE) }
npc <- function(p){cv<-cumsum(p$sdev^2)/sum(p$sdev^2);i<-which(cv>=0.95);max(1L,if(length(i)==0)ncol(p$x)else i[1])}

# fixed selection order on all 96
ph_all <- phenotype
ad_sel  <- grd(td(rk(eta_syn(expr_full[A,,drop=FALSE],expr_full[D,,drop=FALSE],ph_all),FALSE),expr_full,ph_all), MAX_PER_PHENO)
sig_sel <- grd(td(rk(eta_syn(expr_full[S,,drop=FALSE],expr_full[S,,drop=FALSE],ph_all),TRUE ),expr_full,ph_all), MAX_PER_PHENO)
# order ratios by rank within each phenotype so step n = "top n per phenotype"
order_by_pheno <- function(sel){ sel%>%group_by(P)%>%mutate(rk=row_number())%>%ungroup()%>%arrange(rk,desc(Delta)) }
ad_o <- order_by_pheno(ad_sel); sig_o <- order_by_pheno(sig_sel)

# build a PCA->LDA on given ratio set & sample cols
fit_lda <- function(ratios, cols){
  ex<-all_logcpm[,cols,drop=FALSE]; sm<-bs(ratios,ex)
  keep<-apply(sm,1,function(x)length(unique(x))>1); sm<-sm[keep,,drop=FALSE]
  ph<-droplevels(phenotype[cols])
  list(sm=sm, ph=ph, cols=cols, ratios=rownames(sm))
}
loocv_of <- function(ratios){
  ex<-all_logcpm; sm<-bs(ratios,ex); keep<-apply(sm,1,function(x)length(unique(x))>1); sm<-sm[keep,,drop=FALSE]
  if(nrow(sm)<2) return(NA)
  p<-prcomp(t(sm),center=TRUE,scale.=TRUE); k<-npc(p)
  cv<-suppressWarnings(lda(as.data.frame(p$x[,1:k,drop=FALSE]),grouping=phenotype,CV=TRUE))
  mean(cv$class==phenotype)*100
}
loso_of <- function(ratios){
  pred<-rep(NA_character_,ncol(all_logcpm))
  for(s in levels(study)){ te<-which(study==s); tr<-which(study!=s)
    smtr<-bs(ratios, all_logcpm[,tr,drop=FALSE]); keep<-apply(smtr,1,function(x)length(unique(x))>1)
    if(sum(keep)<2) next; rr<-rownames(smtr)[keep]
    ptr<-prcomp(t(smtr[keep,,drop=FALSE]),center=TRUE,scale.=TRUE); k<-npc(ptr)
    md<-lda(as.data.frame(ptr$x[,1:k,drop=FALSE]),grouping=droplevels(phenotype[tr]))
    smte<-bs(rr, all_logcpm[,te,drop=FALSE]); proj<-predict(ptr,t(smte))[,1:k,drop=FALSE]
    pred[te]<-as.character(predict(md,as.data.frame(proj))$class) }
  mean(pred==as.character(phenotype),na.rm=TRUE)*100
}

# ---- sweep ----
rows<-list()
for(n in 1:MAX_PER_PHENO){
  ad_n <- ad_o %>% group_by(P) %>% slice_head(n=n) %>% pull(Ratio)
  sig_n<- sig_o%>% group_by(P) %>% slice_head(n=n) %>% pull(Ratio)
  for(panel in c("A/D","SIG","Combined")){
    rl <- switch(panel, "A/D"=ad_n, "SIG"=sig_n, "Combined"=c(ad_n,sig_n))
    rows[[length(rows)+1]] <- data.frame(
      Panel=panel, Per_Phenotype=n, N_Ratios=length(rl),
      LOOCV=loocv_of(rl), LOSO=loso_of(rl), stringsAsFactors=FALSE)
  }
  cat("  per-phenotype =",n,"done\n")
}
sat<-bind_rows(rows)
write.csv(sat, file.path(out_dir,"Salmon_kTSP_Saturation_LOOCV_LOSO.csv"), row.names=FALSE)

long<-sat%>%pivot_longer(c(LOOCV,LOSO),names_to="Scheme",values_to="Acc")
p<-ggplot(long, aes(N_Ratios, Acc, colour=Scheme, linetype=Scheme))+
  geom_line(linewidth=1)+geom_point(size=2)+
  facet_wrap(~Panel, scales="free_x")+
  scale_colour_manual(values=c(LOOCV="#377EB8", LOSO="#E41A1C"))+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold",hjust=0.5),legend.position="bottom")+
  labs(title="A/D and SIG/SIG saturation: LOOCV vs LOSO",
       subtitle="LOSO = leave-one-study-out (honest); LOOCV = optimistic ceiling",
       x="Number of sign-encoded ratios", y="Accuracy (%)")
ggsave(file.path(out_dir,"Salmon_kTSP_Saturation_LOOCV_LOSO.png"), p, width=11, height=4.5, dpi=300)

cat("\n  Saturation summary (Combined panel):\n")
print(sat%>%filter(Panel=="Combined")%>%select(Per_Phenotype,N_Ratios,LOOCV,LOSO), row.names=FALSE)
cat("\nDONE ->", out_dir, "\n")
