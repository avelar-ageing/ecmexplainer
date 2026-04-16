# ECMExplainer Shiny App

**By Eleanor Sheekey & Roberto A. Avelar**

> Version 0.2.0

A modular Shiny application for predicting fibroblast phenotypic states from bulk
RNA-seq data using ECM gene-expression ratios and a dual-panel PCA + LDA classifier.
Supports both **human** and **mouse** data, flexible batch correction, and fully
reproducible named model bundles.

---

## Current readiness status (April 2026)

### ✅ Ready to use

| Area | Detail |
|------|--------|
| **Backend pipeline** | Smoke test passes — normalise → ratios → PCA+LDA → predictions |
| **Training Wizard** | Upload counts + metadata → feature selection → train → save named bundle |
| **Test Wizard** | 3-way model selector: default model / upload bundle / use session-trained model |
| **Default model** | V19 built-in (83 samples, 6 phenotypes, LOOCV 90.4%) loads automatically |
| **Model bundles** | Save `.rds` + `.json` params; restore all training settings from either file |
| **Make Default button** | Saves a newly trained model to `training_data/` and registers it in the dropdown |
| **Unified model info panel** | Shows name, sample count, classes, LOOCV, studies for all 3 model sources |
| **Batch correction** | ComBat (gene-level) + limma `removeBatchEffect` (ratio-level) |
| **Species selector** | Human / Mouse / Other |
| **Packaging** | `DESCRIPTION` v0.2.0, `LICENSE` (MIT), `.gitignore`, `README.md` all in order |
| **All 11 R modules** | Parse cleanly; no syntax errors |

### ⚠️ Needs your input before it is fully complete

| Item | What to do |
|------|-----------|
| **Register your Salmon-aligned data as a default model** | You have `salmon_raw_counts.csv` + `salmon_metadata.csv` in `/home/ro/Downloads/el_test/`. Upload these in the Training Wizard, run training, then click **"Make this a default model"** — it will be saved to `training_data/` and appear in the Test Wizard dropdown for all future sessions. To make it permanent across app restarts, add one entry to `R/default_models.R`. |
| **Hypoxia samples decision** | The Salmon metadata includes `Hypoxia` (n=3) and `HypoxiaFibrosis` (n=3) samples from GSE139963 — phenotypes not in the V19 model. Decide whether to include them as new classes or filter them out before training. |
| **Combined PCA overlay (optional)** | The `plot_pca_combined()` function (test + training samples together in one PCA space) is a stub — it currently shows training-only PCA with a console message. The per-panel projections work correctly; this is a visualisation enhancement only. |
| **`renv.lock` refresh (optional)** | The lockfile reflects the original development environment. After confirming all packages are installed, run `renv::snapshot()` from inside `shiny_app/` to pin current versions. |

---

---

## Contents

- [What it does](#what-it-does)
- [Predicted phenotypes](#predicted-phenotypes)
- [Data format requirements](#data-format-requirements)
- [Installation](#installation)
- [Launching the app](#launching-the-app)
- [Training Wizard](#training-wizard)
- [Test Wizard](#test-wizard)
- [Saving and restoring models](#saving-and-restoring-models)
- [Prediction modes](#prediction-modes)
- [Folder structure](#folder-structure)
- [Troubleshooting](#troubleshooting)

---

## What it does

ECMExplainer learns a fibroblast phenotype classifier from a labelled bulk RNA-seq
training set. It computes log-ratio features from curated ECM gene pairs, applies
optional sample-level batch correction (ComBat / limma), and trains a two-stage
PCA + LDA model. The trained model can then be applied to new, unlabelled samples
to obtain per-sample phenotype predictions and posterior probabilities.

**Pipeline summary:**

```
Raw counts (genes × samples)
  → log2(CPM + 1) normalisation
  → optional gene-level ComBat batch correction
  → ECM ratio features  (log A − log B)
  → optional ratio-level limma batch correction
  → dual-panel PCA  (AD panel | SIG panel)
  → LDA classifier
  → Phenotype predictions + posterior probabilities
```

---

## Predicted phenotypes

The classifier distinguishes between the following fibroblast states
(derived automatically from the `Phenotype` column of your training metadata):

`Fibrosis` · `Inflammation` · `OxidativeStress` · `Proliferating` ·
`Quiescence` · `Senescence`

*(The exact classes depend on which phenotypes are present in the training metadata
you supply.)*

---

## Data format requirements

### Count matrix (CSV)

| Requirement | Detail |
|-------------|--------|
| Orientation | **Genes in rows, samples in columns** |
| First column | Gene symbols (HGNC for human, MGI for mouse). No header for this column, or use `"Gene"`. |
| Remaining columns | Sample IDs — these must exactly match `Sample_ID` in your metadata |
| Values | Raw integer counts (the app will normalise to log2(CPM + 1)) |
| Gene symbols | Human: `COL1A1`, `FN1`, etc. Mouse: `Col1a1`, `Fn1`, etc. |

### Metadata (CSV)

| Column | Required? | Description |
|--------|-----------|-------------|
| `Sample_ID` | **Yes** | Must match column headers of the count matrix |
| `Phenotype` | **Yes (training only)** | Class label for each sample |
| `Study_ID` | Recommended | Used for batch correction; samples from the same study share a value |

Any additional columns are allowed and ignored by the pipeline.

---

## Installation

### Requirements

- **R** ≥ 4.1 — <https://cran.r-project.org>
- **RStudio** (recommended) — <https://posit.co/download/rstudio-desktop>

### Option A — Reproducible install via renv (recommended)

A `renv.lock` file is included with exact package versions used during development.

```r
install.packages("renv")   # only needed once
renv::restore()            # run from inside the shiny_app/ folder
```

### Option B — Install latest versions

```r
source("shiny_app/install_packages.R")
```

> Run from the **project root**, or adjust the path accordingly.  
> This installs all CRAN and Bioconductor dependencies and prints a confirmation
> table. Only needed once.

---

## Launching the app

#### From RStudio

1. Open `shiny_app/app.R` in RStudio
2. Click the **Run App** button at the top of the editor

#### From an R console

```r
shiny::runApp("shiny_app")
```

#### From a terminal

```bash
Rscript -e 'shiny::runApp("shiny_app")'
```

---

## Training Wizard

The **Train** tab walks you through building a new ECMExplainer model.

### Step 1 — Feature selection options

| Setting | Description |
|---------|-------------|
| **Species** | Human (HGNC) or Mouse (MGI). Determines the ECM gene dictionary used. |
| **Min. ratio coverage (%)** | Ratios present in fewer than this fraction of samples are dropped. |
| **Batch correction** | Tick to enable ComBat gene-level batch correction before computing ratios. Requires `Study_ID` in metadata. |
| **Ratio batch correction** | Tick to enable limma `removeBatchEffect` on the ratio matrix after feature computation. |
| **Top N AD / SIG ratios** | Number of ratios to retain per panel based on variance. |
| **PCs for AD / SIG panel** | Number of principal components to pass to LDA. |
| **Save / restore params** | Upload a previous `.json` params log or `.rds` model bundle to auto-fill all settings. |

### Step 2 — Upload training data

| Upload | Format |
|--------|--------|
| **Count matrix** | CSV — genes × samples, gene symbols in first column |
| **Metadata** | CSV — must contain `Sample_ID` and `Phenotype` |

Click **Run Training** to execute the full pipeline. LOOCV accuracy is shown on completion.

### Step 3 — Save model

Give your model a name, then download:

| Button | Output |
|--------|--------|
| **Download Model (.rds)** | Full model bundle — pass to Test Wizard |
| **Params (.json)** | Human-readable log of all training settings and metrics |
| **Ratios (.csv)** | Table of ratio features used in this model |

---

## Test Wizard

The **Test** tab predicts phenotypes for new, unlabelled samples.

### Step 1 — Upload test data

| Upload | Format |
|--------|--------|
| **Count matrix** | CSV — same format as training |
| **Metadata** | CSV — must contain `Sample_ID`; `Study_ID` optional for batch QC |
| **Model bundle (.rds)** | Saved from the Training Wizard |

When a bundle is uploaded the app automatically displays the model name, species,
training LOOCV accuracy, and number of training samples.

### Step 2 — Run predictions

Click **Run Predictions**. The app projects test samples onto the training PCA
reference frame (batch-corrected) and applies the trained LDA classifier.

### Step 3 — Results

| Tab | Contents |
|-----|---------|
| **Predictions** | Per-sample predicted phenotype + posterior probabilities |
| **PCA plot** | Test samples projected onto training PCA space |
| **Downloads** | Prediction CSV and plot PNG |

---

## Saving and restoring models

### Model bundle (`.rds`)

The bundle contains everything needed to reproduce or re-apply the model:

- `train_obj` — trained PCA + LDA objects
- `features` — ratio feature table
- `train_meta` / `train_expr` — training metadata and expression matrix
- `ecm_dict` — ECM gene dictionary used
- `training_params` — all feature-selection settings, metrics, and dates
- `model_name`, `version`

### Params log (`.json`)

A human-readable JSON file capturing all Training Wizard settings. Upload it in
the **Feature Selection Options** box to auto-restore all numeric, slider, and
checkbox inputs for a new training run.

---

## Prediction modes

Three modes are available in the Test Wizard (Step 2 dropdown):

| Mode | Description |
|------|-------------|
| **Ratios (default)** | Compute ECM gene-pair ratios from test counts and classify |
| **ECM genes only** | Use raw log-CPM of ECM genes directly (no ratio computation) |
| **All genes** | Use full transcriptome log-CPM (high-dimensional; slower) |

---

## Folder structure

```
shiny_app/
├── app.R                    ← entry point (open this in RStudio)
├── install_packages.R       ← run once to install all dependencies
├── smoke_test.R             ← headless pipeline smoke test
├── DESCRIPTION              ← package metadata
├── LICENSE                  ← MIT licence
├── README.md                ← this file
├── renv.lock                ← pinned package versions for reproducibility
├── training_data/           ← V19 pre-trained ECMExplainer data
│   ├── ECMExplained_Master_List_Expanded_Fixed.csv   ← ECM gene dictionary
│   ├── ECMExplorer_v19_Foundation.rds                ← foundation object
│   ├── ECMExplorer_v19_TrainedModel.rds              ← default trained model
│   ├── ECMExplorer_v19_Metadata.csv                  ← training metadata
│   ├── ECMExplorer_v19_Final_Ratios.csv              ← ratio feature table
│   ├── ECMExplorer_v19_AllGenes_Merged.rds           ← full expression matrix
│   └── ECMExplorer_v19_CorrectedVis.rds              ← batch-corrected matrix
├── R/                       ← modular backend (no need to edit)
│   ├── mod_data_loader.R    ← load_training_data()
│   ├── mod_normalize.R      ← normalize_counts()
│   ├── mod_batch.R          ← apply_batch_correction(), apply_ratio_batch_correction()
│   ├── mod_ratios.R         ← compute_ratios()
│   ├── mod_train.R          ← train_ecmexplainer()
│   ├── mod_predict.R        ← predict_ecmexplainer()
│   ├── mod_plots.R          ← plot_* functions
│   ├── app_ui.R             ← UI definition
│   └── app_server.R         ← server logic
└── www/                     ← static assets (CSS, images)
```

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `could not find function "runApp"` | Run `install.packages("shiny")` first |
| `there is no package called '...'` | Re-run `source("install_packages.R")` or `renv::restore()` |
| No ratios found | Check that gene symbols are HGNC (human) or MGI (mouse) — not Ensembl IDs |
| `Sample_ID` mismatch | Column headers of count matrix must exactly match `Sample_ID` values in metadata |
| Plots blank after run | Check the R console for error messages; ensure training completed without warnings |
| LOOCV accuracy is NA | Usually means all samples belong to a single phenotype class — check your metadata |
| Param restore has no effect | Make sure you upload the `.json` **before** clicking Run Training |
| Mouse data not recognised | Select **Mouse** in the Species dropdown before uploading data |

---

## Citation

If you use ECMExplainer in your work, please cite:

> Sheekey E, Avelar RA, et al. *ECMExplainer: …* (manuscript in preparation)
