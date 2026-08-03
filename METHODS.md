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
case-insensitive. Duplicate symbols are collapsed to the highest-variance row
within your cohort. Each sample's sign vector is passed through the **frozen**
training PCA (`predict.prcomp` with the stored centre, scale and eigenvectors,
retaining the stored number of components) and then the **frozen** LDA. Nothing
is re-estimated on your data, so your samples land in the same coordinate system
as the published figures and distances are directly comparable.

The app reports **ratio coverage** — the fraction of panel ratios with both genes
present. All 21 validation cohorts in the manuscript exceeded 85%; below that the
app warns you.

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
the smallest class; ECM-panel intersection; variance > 0.005; η² scoring of
candidate pairs with a floor of 0.10; phenotype-stratified greedy selection with
each gene usable at most once; PCA (centred and scaled, ≥95% cumulative variance,
capped at 10 components); LDA; leave-one-out cross-validation.

## Provenance of the bundled model

`data/TrainedModel_ADSIG_180.rds` is the model from the manuscript, produced by
`003_FinalTrain_TwoPanels.R` from `salmon_counts_Combined_FINAL.rds`.

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

## Known discrepancies between the manuscript text and the code

Recorded so anyone reading the paper alongside this repository can see where they
diverge. **These are documentation issues, not errors in the published panel** —
the model reproduces exactly — but the Methods text does not everywhere describe
what the code does.

1. **Panel selection ranks by TSP delta, not η² (5.3).** The text says the greedy
   selector adds ratios "in descending η² order". The code ranks by
   `TSP_Delta` — the absolute difference in the proportion of positive signs
   inside versus outside each phenotype — and uses η² only as a > 0.10 pre-filter.
   This is not cosmetic: at 6 ratios per phenotype the two criteria select
   **completely disjoint panels** (0 of 36 ratios shared). The app exposes both
   via *Selection ranking*; `tsp_delta` is the default because it is what
   produced the published model.

2. **The "at least five ratios per module" constraint (5.3) is not implemented.**
   No such rule exists in the selection code; the published 76/104 split satisfies
   it incidentally. Exposed here as *Minimum ratios per module*, default 0.

3. **ECM-Genes baseline size (5.4).** The text says 820 genes; the filter in
   `003_competing_models_TwoPanel.r` (curated list ∩ matrix, then variance >
   0.005) yields **895**, which is also the number 5.2 gives for the
   batch-corrected universe.

4. **"Common feature space of 18,957 genes" (5.1).** The common matrix holds
   **19,861** genes. 18,957 is the count *after* the variance > 0.005 filter — the
   All-Genes baseline figure from 5.4.

5. **LOSO results are not reproducible from the released code.** The 57.3% peak,
   the 45.8–57.3% sweep range and the reference-projection gains (All-Genes
   28→45%, ECM-Genes 40→49%) are read from `PanelSweep_LOOCV_LOSO.csv`, generated
   by a script (`002c-FASTv2`) not present in the released bundle. The §3.8
   benchmark figures (77.1 / 70.0 / 65.0 / 56.3%) likewise depend on validation
   scripts and cohort data that are not included.

6. **LOOCV cannot separate the four models.** POOL-120, All-Genes and ECM-Genes
   all reach 100% LOOCV against ADSIG-180's 92.7%. The model comparison therefore
   rests entirely on LOSO and the held-out validation set — see point 5.

7. **`final_version/` does not run as shipped.** `003_FinalTrain_TwoPanels.R`
   loads `salmon_counts_Combined_FINAL.rds` and
   `salmon_metadata_Combined_FINAL.csv`, which are absent from that folder.

8. **Dead code.** Lines 363–679 of `003_FinalTrain_TwoPanels.R` reference
   out-of-scope variables and files that do not exist, and error if the script is
   sourced whole.

## Data handling

Uploaded files are read into memory for the duration of your session and are
never written to disk by the app. Nothing is logged, stored or transmitted
onward. The public instance additionally caps upload size, sample count and
retraining size; those caps are configuration, not method, and are absent when
you run the app yourself.
