# ECMExplainer

**Eleanor Sheekey, Roberto A. Avelar & Hironobu Fujiwara**

A functional matrisome framework for decoding cellular stress states and aging
trajectories. ECMExplainer positions a fibroblast transcriptome in a phenotype
space defined by six reference states — Proliferating, Quiescence, Senescence,
Fibrosis, Inflammation, Hypoxia — using 180 sign-encoded matrisome gene-pair
ratios (76 Assembly/Destruction + 104 Signaling/Signaling).

Because each ratio records only *which of two genes is more highly expressed*, it
is insensitive to sequencing depth and per-sample scaling, so a sample can be
scored on its own without cohort-level normalisation.

**Live instance:** <https://apps.avelar.science/ecm>
(capped for shared use — run it locally for unrestricted analysis)

---

## Run it locally, with no limits

```bash
git clone https://github.com/avelar-ageing/ecmexplainer.git
cd ecmexplainer
Rscript install_packages.R     # CRAN packages only, no Bioconductor needed
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

Or with Docker, which pins the whole environment:

```bash
docker build -t ecmexplainer .
docker run --rm -p 3838:3838 ecmexplainer
```

The local build has no upload cap, no sample cap and no retraining cap. The
public instance sets `ECMX_MODE=web`; nothing scientific differs between them.

## What you can do

| Tab | |
|---|---|
| **Classify** | Project your samples through the frozen ADSIG-180 model: phenotype space, posterior probabilities, centroid distances, panel heatmap, downloadable results |
| **Trajectories** | Medoid-based condition trajectories in any two-centroid plane, plus protractor angles with a 2,000-resample bootstrap CI |
| **Retrain** | Train a new panel on your own data with **your own phenotype labels** — the six reference states are not required. Choose design, panel size, ranking criterion and module floor; download the model bundle |
| **ECM Aging Score** | Fit per-ratio age slopes by OLS with BH FDR control, then score samples |
| **Methods** | The implemented method, model provenance, and the manuscript/code discrepancy record |

## Input format

A CSV with genes in the first column and one sample per remaining column:

```
Gene,Sample1,Sample2,...
COL1A1,1420,980,...
MMP2,311,455,...
```

Raw counts or log2-CPM are both fine — the app detects which and transforms if
needed. **TPM, FPKM and RPKM are not valid inputs**: length normalisation applies
a per-gene factor that does not cancel inside a ratio. HGNC symbols, matched
case-insensitively; duplicates collapse to the highest-variance row.

Metadata is optional — `SampleID` plus any columns you want to group by.

## Reproducing the manuscript

```bash
Rscript tests/test_manuscript_numbers.R   # ~2 min
Rscript tests/test_app_smoke.R            # ~3 min
```

The first retrains ADSIG-180 from raw counts through this repository's own code
and asserts every published number. It reproduces the shipped model exactly —
all 180 ratios in the same order, LOOCV 92.71%, macro-F1 0.921.

`METHODS.md` records eight places where the manuscript text and the released code
disagree, including one that matters scientifically: the panel is selected by TSP
delta, not by the η² ordering the Methods describe, and the two criteria pick
disjoint panels.

## Layout

```
app.R                  Shiny UI + server
R/ecm_core.R           the method: encoding, projection, training, aging score
R/ecm_figures.R        manuscript-style figures
data/                  frozen model, curated gene list, training cohort, demo data
tests/                 manuscript-number and smoke tests
legacy/v19_app/        the previous V19 app, kept for reference (see below)
```

### On the previous app

`legacy/v19_app/` is the earlier application. It was built on a different model —
**V19: 8 studies, 19 ratios, LOOCV 90.4%, with `OxidativeStress` in place of
`Hypoxia`** — and therefore did not reproduce the manuscript. It is retained for
provenance only; nothing in the current app depends on it.

## Citation

Sheekey E, Avelar RA, Fujiwara H. *ECMExplainer: a functional matrisome framework
for decoding cellular stress states and aging trajectories.*

## Licence

MIT — see `LICENSE`.
