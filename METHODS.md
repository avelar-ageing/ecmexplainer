# Methods and provenance

This app implements the ECMExplainer method described in Sheekey, Avelar &
Fujiwara, *ECMExplainer: a functional matrisome framework for decoding cellular
stress states and aging trajectories*. Section numbers refer to that
manuscript's Methods.

## What the app computes

**Sign encoding (5.3).** Each feature is the within-sample rank order of a gene
pair, `x = sign(log2 X_A − log2 X_B) ∈ {−1, 0, +1}`. Any per-sample additive
offset — library-size rescaling, TMM — cancels inside the difference, and any
monotone within-sample transform preserves rank order. This is why log2-CPM and
log2-RPM inputs are interchangeable here.

**Length-normalised units are rejected by design.** RPKM, FPKM and TPM apply a
per-gene length factor that shifts the two members of a ratio unequally and does
*not* cancel. Do not feed them to this app.

**Deployment projection (5.5).** A ratio is computed only when both members are
present in your data; ratios with a missing member are neutralised to `0`,
contributing neither a positive nor a negative vote. Symbol matching is
case-insensitive (the manuscript's scripts matched case-sensitively; on every
validation cohort the two give identical results). Duplicate symbols are collapsed to the highest-variance row
within your cohort. Each sample's sign vector is passed through the **frozen**
training PCA (`predict.prcomp` with the stored centre, scale and eigenvectors,
retaining the stored number of components) and then the **frozen** LDA. Nothing
is re-estimated on your data, so your samples land in the same coordinate system
as the published figures and distances are directly comparable.

The app reports **ratio coverage** — the fraction of panel ratios with both genes
present. Every one of the 18 validation cohorts in the manuscript had at least 90%
for each panel; below 85% the app warns you.

**Confidence and out-of-distribution measures.** For each sample and panel the
app reports the largest posterior, the margin to the second largest, the Shannon
entropy of the posterior, the nearest training centroid and its distance, and the
Mahalanobis distance of the sample's PC scores from the training PC scores, with
a flag when it exceeds the 95th percentile of χ² on the number of PCs (SI Table 6,
sheet O). On the manuscript's validation cohorts that flag is raised for every
sample under every panel (median distance ≈ 900 against a threshold of 18.3), so
it does not separate typical from atypical samples and should not be read as one.

**Centroid geometry (5.6).** Training phenotype centroids are per-class means of
the training LDA scores. In this sphered discriminant space within-class
covariance is approximately isotropic, so plain Euclidean distance is equivalent
to the Mahalanobis distance the classifier uses.

**Trajectories (5.6–5.7).** Each condition is represented by its **medoid** — the
real sample minimising total distance to the others in that condition — computed
in the full LDA space, not in the 2-D distance projection. Using a real sample
guarantees the plotted point respects `D_A + D_B ≥ ‖μ_A − μ_B‖`, which
independent per-axis medians can violate. Protractor angles are
`θ = arccos((v·d)/(‖v‖‖d‖))` from a shared origin, with a 2,000-resample
percentile bootstrap CI and a minimum of 2 samples per group.

**ECM Aging Score (5.9).** Per-ratio age slopes are fitted by OLS on the
continuous log-ratio, retained at Benjamini–Hochberg FDR < 0.05, and applied to
the sign features. Because β is an age slope, the score increases with age by
construction — a lower score is oriented as younger.

**Retraining (5.2–5.4).** Expression filter CPM > 1 in at least as many samples as
the smallest class; gene-universe intersection; variance > 0.005; η² of the
continuous log-ratio with a floor of 0.10; TSP delta (one-versus-rest difference
in how often gene A exceeds gene B); one global phenotype-stratified greedy pass
in descending TSP delta, each gene usable at most once; PCA (centred and scaled,
≥95% cumulative variance, capped at 10 components); LDA; leave-one-out
cross-validation. Three designs share this pipeline and differ only in which
pairs are candidates:

| Design | Candidate pairs | Bundled panel |
|---|---|---|
| ADSIG | Assembly × Destruction, and Signalling × Signalling, from the gene list | ADSIG-180 (ECMExplainer) |
| POOL | Any two genes in the gene list | POOL-180 |
| ALLGENE | Any two of the 2,000 most variable expressed genes genome-wide | ALLGENE-180 |

Any gene list can be supplied; ADSIG additionally needs each gene's module.

**Nested leave-one-study-out (019_PanelSweep_Fast.r).** The gene universe is fixed
on all samples; then for each study in turn, η², TSP delta, selection, PCA and LDA
are recomputed on the other studies and the held-out study is classified. LOOCV
saturates near 100% because selection sees every sample; nested LOSO does not.
"Restricted" accuracy counts only phenotypes present in at least three training
studies, the ones a held-out study can test fairly.

## Provenance of the bundled models

`data/TrainedModel_ADSIG_180.rds` is the model from the manuscript, produced by
`003_FinalTrain_TwoPanels.R` from `salmon_counts_Combined_FINAL.rds`.
`data/TrainedModel_POOL_180.rds` and `data/TrainedModel_ALLGENE_180.rds` are built
by `scripts/build_models.R` with this repository's code and stop the build unless
they reproduce Supplementary Tables 3 and 5:

| Panel | LOOCV | Nested LOSO (all / restricted) | Ratios |
|---|---|---|---|
| ADSIG-180 | 92.71% | 47.92% / 56.79% | identical to SI Table 3 A |
| POOL-180 | 97.92% | 41.67% / 49.38% | identical to SI Table 5 A |
| ALLGENE-180 | 100.00% | 42.71% / 50.62% | identical to SI Table 5 I |

`data/validation_expr.rds` holds the 769 genes used by any of the three panels for
the 327 validation samples of Supplementary Table 6, each cohort on its own scale
(log2-CPM, or the batch-corrected VST of the collaborator cohorts).
`tests/test_validation.R` scores them and asserts that every prediction, nearest
centroid and confidence value, and the headline accuracies, equal SI Table 6.
Gene symbols for the recount3 cohorts come from `org.Hs.eg.db` 3.18 (as in the manuscript; 3.19 and 3.20 give the same result); later
releases rename a few non-matrisome genes in ALLGENE-180 and change 10 of its
predictions, while ADSIG-180 and POOL-180 are unaffected.

`tests/test_manuscript_numbers.R` retrains it from raw counts through this
repository's own `R/ecm_core.R` and asserts the published values. It passes,
reproducing the shipped panel **exactly — all 180 ratios in the same order**:

| Quantity | Manuscript | This code |
|---|---|---|
| Training samples / studies | 96 / 10 | 96 / 10 |
| Curated panel (A / D / SIG) | 1,139 (412 / 119 / 608) | identical |
| ECM genes detected in all 96 | 923 | 923 |
| Modelling universe (A / D / SIG) | 711 (303 / 88 / 320) | identical |
| ADSIG-180 composition | 76 A/D + 104 SIG | 76 + 104 |
| PCs retained | ≤10 | 10 |
| LOOCV accuracy | 92.7% | 92.71% |
| macro-F1 | 0.921 | 0.921 |

## Data handling

Uploaded files are read into memory for the duration of your session and are
never written to disk by the app. Nothing is logged, stored or transmitted
onward. The public instance additionally caps upload size, sample count and
retraining size; those caps are configuration, not method, and are absent when
you run the app yourself.
