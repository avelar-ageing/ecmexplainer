# Methods

## RNA-seq data acquisition

Publicly available RNA-seq datasets were obtained from nine studies deposited in the NCBI Gene Expression Omnibus (GEO) and Sequence Read Archive (SRA): GSE63577, GSE102674, GSE74324, GSE75643, GSE117337, GSE117444, GSE139963, GSE221693 and GSE284342. Together these comprised 88 biological samples spanning eight fibroblast phenotypic states: Proliferating (*n* = 36), Senescence (*n* = 17), Quiescence (*n* = 15), Fibrosis (*n* = 9), Hypoxia (*n* = 3), Hypoxia/Fibrosis (*n* = 3), Oxidative Stress (*n* = 3) and Inflammation (*n* = 2). Cell types included BJ, IMR90, WI-38, HFF and MRC-5 foreskin fibroblasts, primary lung and cardiac fibroblasts and Tig3 lung fibroblasts. Libraries were sequenced on Illumina HiSeq 2000, HiSeq 2500, NextSeq 500 or NovaSeq X Plus platforms, with both single-end and paired-end layouts represented (Supplementary Table 1).

## Read download and quality control

Raw sequencing reads were downloaded from NCBI SRA using fasterq-dump (SRA Toolkit v3.2.1) with the `--split-3` flag to separate paired-end and single-end reads automatically. For samples sequenced across multiple lanes (GSE102674; four SRR accessions per biological sample), all lane-level FASTQ files were retained as separate files and passed jointly to the quantification step, thereby merging technical replicates at the alignment stage rather than post hoc. Downloaded FASTQ files were compressed with pigz prior to quantification and deleted immediately after each sample was processed to minimise disk usage.

## Transcript quantification

Transcript-level abundance estimation was performed using Salmon v1.11.4 in quasi-mapping mode against a pre-built index of the human reference transcriptome (GRCh38/Ensembl release 111 cDNA; *Homo sapiens* GRCh38.cdna.all.fa.gz). The index was constructed with default *k*-mer size (*k* = 31). Quantification was run with automatic library type detection (`--libType A`), sequence-composition bias correction (`--seqBias`), GC-content bias correction (`--gcBias`; paired-end samples only) and selective alignment (`--validateMappings`). For multi-lane samples, all FASTQ files for a given biological replicate were supplied simultaneously via repeated `-1`/`-2` (paired-end) or `--unmatedReads` (single-end) flags in a single Salmon call. Across all 88 samples, the median quasi-mapping rate was 84.9% (range: 63.0–93.8%).

## Gene-level count summarisation

Transcript-level estimated counts from Salmon were collapsed to gene-level counts using tximport v1.38.2 in R v4.5.1. A transcript-to-gene mapping was constructed from the FASTA headers of the Ensembl 111 transcriptome by extracting the `gene_symbol` field, yielding 205,792 transcript–gene mappings across 36,539 genes. Version suffixes were stripped from transcript identifiers (e.g. ENST00000415118.1 → ENST00000415118) to ensure concordance between the FASTA and Salmon output, and the `ignoreTxVersion = TRUE` and `ignoreAfterBar = TRUE` options were enabled. Raw gene-level counts were obtained with `countsFromAbundance = "no"`.

## Filtering and normalisation

Genes detected in fewer than 20% of samples were removed, reducing the gene set from 35,968 to 21,378. Library-size normalisation was performed on a per-study basis by computing counts per million (CPM) within each study independently, followed by a log₂(CPM + 1) transformation. This per-study normalisation was chosen because median library sizes varied substantially across studies (range: 5.0–39.5 million mapped fragments), and pooled CPM calculation would have allowed between-study depth differences to inflate subsequent batch-correction estimates.

## Batch correction

Inter-study technical variation was removed using the parametric empirical Bayes method ComBat from the sva package (v3.58.0). The batch variable was defined as the GEO study of origin (nine levels). A model matrix encoding the biological variable of interest (phenotypic condition; eight levels) was supplied to preserve genuine biological signal during correction. Genes with zero variance within any single batch (3,815 genes) were excluded from the adjustment, yielding a final batch-corrected matrix of 19,208 genes across 88 samples.

## ECM gene ratio feature engineering

An ECM master gene list was used to define three functional categories: Assembly, Destruction and Signalling. For each pair of genes spanning Assembly × Destruction (A/D pairs) and Signalling × Signalling (SIG pairs), log₂-expression ratios were computed from the batch-corrected expression matrix. The discriminatory power of each ratio was quantified using a synergy-aware η² (eta-squared) statistic that measures the proportion of between-group variance in the ratio that exceeds the maximum η² of either constituent gene alone. Ratios were ranked by η² and filtered to retain those with η² > 0.30 and synergy δ > 0.05. A greedy selection algorithm was applied to limit each gene to a maximum of one appearance across all selected ratios, yielding up to 100 A/D and 100 Signalling candidate ratios. The final feature set comprised the top 50 A/D and top 50 Signalling ratios (or fewer if insufficient candidates passed the threshold).

## Batch correction of ratios

Ratio-level batch correction was performed using limma::removeBatchEffect on the ratio matrix (not the gene-level matrix), with study as the batch variable and phenotype as the protected covariate. This ratio-level correction was preferred over gene-level correction for the classifier because ratios that are inherently batch-robust (i.e. both genes shift together across batches) require minimal adjustment, whereas gene-level correction may introduce artefacts into the ratio space.

## Dimensionality reduction and classification

Principal component analysis (PCA) was applied separately to the batch-corrected A/D and Signalling ratio matrices. The number of principal components retained was determined by the minimum number explaining ≥ 95% of variance in each space. A/D and Signalling PCs were concatenated into a joint feature matrix and used to train a linear discriminant analysis (LDA) classifier (MASS package) to predict the eight phenotypic states. Model performance was evaluated by leave-one-out cross-validation (LOOCV).

## Visualisation

PCA was performed on the ComBat-corrected expression matrix (all 19,208 genes and ECM gene subset) and on the ratio feature matrix. For each, the first three principal components were plotted pairwise (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3), with samples coloured by phenotypic condition, study of origin, sequencing platform, tissue of origin, cell type and library layout. Heatmaps of mean z-scored ratios per phenotype were generated using pheatmap with Ward's hierarchical clustering. A diagnostic heatmap was constructed by selecting the top three most phenotype-specific ratios per condition using a greedy isolation-score metric (minimum absolute z-score difference to the nearest other phenotype).

## Data and code availability

All raw sequencing data are available from NCBI SRA under the accessions listed in Supplementary Table 1. The alignment pipeline (`eleanor_align_script.sh`), tximport preprocessing (`tximport_preprocess.R`), normalisation and PCA (`ecm_normalise_pca.R`) and the full ECMExplainer pipeline (`new_preprocess_ECMExplainer.R`) are provided as supplementary code. The Salmon index was built from *Homo sapiens* GRCh38 Ensembl release 111 cDNA sequences. Batch-corrected expression matrices and trained model objects are available as RDS files in the accompanying data repository.

---

**Software versions:** R v4.5.1, Salmon v1.11.4, SRA Toolkit v3.2.1 (fasterq-dump), tximport v1.38.2, sva v3.58.0 (ComBat), limma, MASS, pheatmap.
