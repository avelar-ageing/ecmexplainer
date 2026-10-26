# ECMExplainer

**Eleanor Sheekey**<sup>1,3</sup>, **Roberto A. Avelar**<sup>2</sup>, **Hironobu Fujiwara**<sup>1</sup>

<sup>1</sup> RIKEN Center for Biosystems Dynamics Research, Kobe, Japan<br>
<sup>2</sup> Institute for Biostatistics and Informatics in Medicine and Ageing Research, Rostock University Medical Center, Germany<br>
<sup>3</sup> Corresponding author

*ECMExplainer: a functional matrisome framework for decoding cellular stress
states and aging trajectories.* Manuscript in preparation; the citation and link
will appear here on publication.

## What it does

ECMExplainer places a fibroblast transcriptome in a phenotype space defined by
six reference states — Proliferating, Quiescence, Senescence, Fibrosis,
Inflammation and Hypoxia — using 180 matrisome gene-pair ratios. Each ratio
records only which of its two genes is more highly expressed in a sample, so a
sample can be scored on its own, without normalising it against a cohort.

The app carries the manuscript's three panels:

- **ADSIG-180, the ECMExplainer model:** Assembly/Destruction and
  Signalling/Signalling ratios from a curated matrisome list.
- **POOL-180:** any pair within the same list, without the functional rule.
- **ALLGENE-180:** any pair within the 2,000 most variable genes genome-wide.

## Quick start

1. **Classify** — upload a gene × sample table (raw counts or log2-CPM; HGNC
   symbols in the first column) and, optionally, a metadata table with a
   `SampleID` column. Filter samples by any metadata column, tick the panels to
   run, and pick a known-phenotype column to score accuracy. Example files are
   linked under *Upload my own*.
2. **Retrain** — train a new panel on the manuscript cohort or your own data,
   with your own labels and your own gene list, then test it on the Classify tab.
3. **Manuscript** — the published validation, recomputed live, and the
   supplementary tables.

## Run it on your own computer

The public instance caps upload and training size. The local app has no caps:

```bash
git clone https://github.com/avelar-ageing/ecmexplainer.git
cd ecmexplainer
Rscript install_packages.R
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

or, with Docker, `docker build -t ecmexplainer . && docker run --rm -p 3838:3838 ecmexplainer`
and open <http://localhost:3838/ecmexplainer/>.

## Methods
