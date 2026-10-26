# ECMExplainer

**Eleanor Sheekey**<sup>1,3</sup>, **Roberto A. Avelar**<sup>2</sup>, **Hironobu Fujiwara**<sup>1</sup>

<sup>1</sup> RIKEN Center for Biosystems Dynamics Research, Kobe, Japan<br>
<sup>2</sup> Institute for Biostatistics and Informatics in Medicine and Ageing Research, Rostock University Medical Center, Germany<br>
<sup>3</sup> Corresponding author

*ECMExplainer: a functional matrisome framework for decoding cellular stress
states and aging trajectories.* The citation and link will be added on publication.

ECMExplainer places a fibroblast transcriptome in a phenotype space defined by six
reference states — Proliferating, Quiescence, Senescence, Fibrosis, Inflammation and
Hypoxia — using 180 matrisome gene-pair ratios. Each ratio records only *which of
its two genes is more highly expressed* in a sample, so it is insensitive to
sequencing depth and per-sample scaling, and a sample can be scored on its own.

**Web version:** <https://apps.avelar.science/ecmexplainer/> (capped for shared
use; run it locally for unrestricted analysis)

---

## Run it on your own computer

You need R ≥ 4.1.

```bash
git clone https://github.com/avelar-ageing/ecmexplainer.git
cd ecmexplainer
Rscript install_packages.R     # CRAN packages only
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

Or with Docker, which pins the whole environment:

```bash
docker build -t ecmexplainer .
docker run --rm -p 3838:3838 -e ECMX_MODE=local ecmexplainer
# open http://localhost:3838/ecmexplainer/
```

The local app has no upload, sample or training caps. The public instance sets
`ECMX_MODE=web`; nothing scientific differs between them.

## What you can do

| Tab | |
|---|---|
| **Classify** | Score samples with the three manuscript panels side by side — **ADSIG-180** (the ECMExplainer model), **POOL-180** and **ALLGENE-180** (controls) — and with any model you train. Phenotype space, posteriors, centroid distances, panel agreement, panel heatmap, a results table with confidence and out-of-distribution measures, and CSV/PNG/PDF downloads. Filter samples by any metadata column. With a known-phenotype column, the **Accuracy** view gives micro and balanced accuracy with Wilson CIs, recall per phenotype (Figure 7B), confusion matrices, accuracy per dataset and McNemar tests between panels. Opens on the manuscript's 327 validation samples. |
| **Trajectories** | Medoid trajectories between two training centroids and protractor angles with bootstrap CIs. |
| **Retrain** | Train a panel on the manuscript cohort or your own data, with **your own labels** and **your own gene list** (upload or paste), as ADSIG, POOL or ALLGENE; filter the training samples; estimate nested leave-one-study-out accuracy. The new model joins the Classify tab for testing, and can be downloaded and loaded again later. |
| **ECM Aging Score** | Fit per-ratio age slopes (OLS, BH FDR) and score samples. |
| **Manuscript** | The published validation recomputed live (SI Table 6), training cross-validation (SI Tables 3, 5), and Supplementary Tables 1–6 to download. |

## Input format

An expression table with HGNC symbols in the first column and one sample per
remaining column (CSV, optionally gzipped):

```
Gene,Sample1,Sample2,...
COL1A1,1420,980,...
MMP2,311,455,...
```

Raw counts, log2-CPM or VST all work, because each ratio compares two genes within
one sample. **TPM, FPKM and RPKM are not valid inputs**: length normalisation applies
a per-gene factor that does not cancel inside a ratio. Symbols are matched
case-insensitively; duplicates collapse to the highest-variance row.

Metadata is optional: a `SampleID` column plus any columns to group, filter or
score by (for example `Condition`, `Dataset`, `Truth`). Retraining needs a label
column and, for leave-one-study-out, a study column. Example files are linked in
the app under *Upload my own* and live in `data/examples/`.

A gene list for retraining is a CSV with a `Gene` column; add a `Function` column
(Assembly / Destruction / Signalling) to use the ADSIG design.

## Reproducing the manuscript

```bash
Rscript tests/test_manuscript_numbers.R   # retrains ADSIG-180 from raw counts
Rscript tests/test_validation.R           # reproduces SI Table 6 from the bundled cohorts
Rscript tests/test_app_server.R           # drives the app headlessly
Rscript tests/test_app_smoke.R
Rscript scripts/build_models.R            # rebuilds POOL-180 and ALLGENE-180, gated on SI Tables 3 and 5
```

| Check | Result |
|---|---|
| ADSIG-180 retrained from raw counts | all 180 ratios in the shipped order; LOOCV 92.71%, macro-F1 0.921 |
| POOL-180, ALLGENE-180 rebuilt | ratio sets identical to SI Table 5 sheets A and I; LOOCV 97.92% and 100% |
| Nested leave-one-study-out | 47.92 / 41.67 / 42.71%, as in SI Tables 3 (tab C) and 5 |
| 327 validation samples, 3 panels | every prediction, nearest centroid and confidence value equal to SI Table 6; ADSIG-180 balanced accuracy 72.7%, micro 83.6% |

`data/validation_expr.rds` holds, for the 327 validation samples, the 769 genes any
of the three panels uses, each cohort on its own scale. The cohorts were rebuilt
from recount3, GEO and the collaborator matrix exactly as the manuscript's
validation scripts built them, with gene symbols from `org.Hs.eg.db` 3.20.0.

## Layout

```
app.R                  Shiny UI and server
R/ecm_core.R           the method: encoding, projection, training, cross-validation, metrics
R/ecm_figures.R        manuscript-style figures
data/                  the three panels, curated gene list, training cohort,
                       validation cohorts, example upload files
supplementary/         Supplementary Tables 1–6
scripts/               panel builds and the manuscript's training-side analyses
tests/                 reproduction and app tests
legacy/v19_app/        an earlier app on a different model, kept for provenance
```

`legacy/v19_app/` was built on a different model (V19: 8 studies, 19 ratios,
`OxidativeStress` in place of `Hypoxia`) and does not reproduce the manuscript.
Nothing in the current app depends on it.

## Licence

MIT — see `LICENSE`.
