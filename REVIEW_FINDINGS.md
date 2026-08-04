# ECMExplainer — outstanding manuscript/code issues

**Date:** 4 August 2026 (second pass)
**Manuscript:** `Human_ECMExplainer_HF-1.pdf`
**Code/tables:** `final_version/`, `extra_for_claude/` (2026-08-04 13:22 drop),
Supplementary Tables 1–4 and **SuppTable5_Validation_MultiTab_v2.xlsx**

Everything below is still open. Items fixed in this manuscript revision, and in the
app and repository, are listed at the end rather than repeated as findings.

## Summary

The model is sound and reproduces exactly: re-running the training from raw counts
regenerates the published ADSIG-180 panel — all 180 ratios in the same order,
LOOCV 92.71%, macro-F1 0.921, 76 A/D + 104 SIG — and the nested LOSO sweep
reproduces Supplementary Table 5C row for row.

This revision fixed seven of the eleven results-level issues from the first pass.
Six remain, of which two are stale numbers left behind by the revision itself.

---

## A. Issues affecting reported numbers

### A1. The Discussion still quotes the old benchmark figure

§4.1 opens: "ECMExplainer achieved a mean balanced accuracy of **77%** across
independent validation cohorts". §3.8 now reports **75.7%**, which is what
Supplementary Table 5 v2 supports (mean of 96.8, 100.0, 61.1, 45.0 = 75.725).
The Discussion was not updated with the Results.

### A2. "100% recall (65/65)" no longer matches the sample count

§3.8 states ADSIG-180 "achieved 100% recall (65/65)" for Proliferating. In
Supplementary Table 5 v2 Proliferating is **N = 59**, not 65 — the scored set moved
from 226 samples/21 datasets to 211/19. The recall of 100% is right; the fraction
should read **59/59**. (POOL-120's 93.2% in the same sentence does match v2.)

### A3. §5.4 promises restricted LOSO figures that appear nowhere, and reuses the term

The paragraph states: "LOSO **balanced** accuracy was therefore reported both across
all six phenotypes and restricted to phenotypes represented by ≥2 and ≥3 source
studies." No such values appear in the manuscript or in Supplementary Table 5C,
which carries only LOOCV, LOSO and LOSO macro-F1. Two problems:

1. The sweep's LOSO column is **plain accuracy** (`mean(pred == truth)`), so calling
   it balanced accuracy is the first-pass A1 issue surviving in this sentence —
   four lines after the same section correctly says "peak LOSO **accuracy**".
2. The promised restricted values are missing. For the record, computed from the
   sweep's per-sample predictions:

| ADSIG-180 LOSO | Value |
|---|---|
| Plain accuracy (the reported 57.3%) | 57.29% |
| Balanced accuracy, all six phenotypes | **37.60%** |
| Balanced accuracy, ≥2 source studies | **37.60%** |
| Balanced accuracy, ≥3 source studies | **56.39%** |

Reproduce with `scripts/loso_balanced_accuracy.R`. Either add these numbers or drop
the sentence; as written it commits to a result the paper does not contain.

### A4. The η² admission filter is not computed as described

§3.1 and §5.3 define η² on the ratio's "binary sign vector", and
`002c_nested_loso_sweep.R` does exactly that. But
`003_FinalTrain_TwoPanels.R` — the script that produced the released model —
computes η² on the **continuous log-ratio** instead. Since η² gates which pairs
enter the candidate pool at all, the deployed panel was filtered by a different
statistic from the one described. Lower stakes now that the *ranking* is correctly
attributed to TSP delta, but the filter description still does not match the
released model.

### A5. LOOCV cannot separate the four models, and this is not stated

POOL-120, All-Genes and ECM-Genes all reach **100% LOOCV**; ADSIG-180 reaches
92.71%. On the in-sample metric the ratio model is the weakest of the four, so the
entire model comparison rests on LOSO and the 211-sample validation set. The
per-model LOOCV confusion matrices are in Supplementary Table 2 but the ceiling is
never noted in the text, which risks a reader assuming LOOCV supports the ranking.

### A6. The "52–57% vs 49%" comparison may mix statistics *(unverified)*

52–57% matches the sweep's plain LOSO accuracy range for ADSIG. The baseline
figures (All-Genes 28→45%, ECM-Genes 40→49%) come from the reference-projection
computation in `003_competing_models_TwoPanel.r`, which has not been run here for
lack of the validation cohort data. If those are balanced accuracies, the
comparison is not like-for-like.

---

## B. Code defects

### B1. `003_FinalTrain_TwoPanels.R` lines 363–679 are dead code that errors

They reference out-of-scope variables (`model_tag`, `PANEL`, `ratio_mat_sign`,
`cm`, `prec`, `rec`, `f1`) and files that do not exist
(`Final_Ratios_ADSIG_BALANCED*.csv`, `perSample_geneExpression_long*.csv`).
Sourcing the script whole fails after the summary table.

### B2. The "Top-3 ratio pairs per phenotype" figure does not rank by any score

In the same block: `slice_max(order_by = Ratio, n = 3)`. `Ratio` is the
`"GENEA / GENEB"` string, so this takes the alphabetically-last three names per
phenotype, not the top three by any statistic. The code comment even says "swap to
TSP_Delta to rank by score". Note that Figure 1G is captioned "top three
diagnostic ratios per phenotype" — if it came from this block, the caption is wrong.

### B3. `final_version/` cannot run as shipped

`003` loads `salmon_counts_Combined_FINAL.rds` and
`salmon_metadata_Combined_FINAL.csv`; the folder contains `salmon_raw_counts.rds`
and `salmon_metadata.csv` instead. The real inputs had to be taken from the older
bundle.

---

## C. Publication readiness

### C1. §6 Data Availability is still an unfinished placeholder

"ECMExplainer can be found on GitHub…" — no URL, no accession, no DOI. The app is
live at <https://apps.avelar.science/ecm> and the code is at
`github.com/avelar-ageing/ecmexplainer`.

### C2. The GitHub repository is private

Reviewers cannot see the code or the app source as things stand.

### C3. No archived snapshot

No Zenodo DOI and no published container image, so the artefact does not outlive the
server hosting it.

---

## Fixed in this revision

- §5.4 headline now reads "peak LOSO **accuracy** of 57.3%" (was "balanced accuracy")
- §5.3 now states the selector works "in descending **TSP delta** order", and that
  "η² acts as a variance-based admission filter, while TSP delta provides the
  ranking" — with the reason given
- The five-ratio module floor is now correctly scoped to the LOSO sweep, and noted
  as "never binding at any panel size tested"
- Figure 1E legend now states that per-fold panels differ from the fixed ADSIG-180
  model (76 A/D, 104 SIG), which was the first pass's panel-mismatch finding
- ECM-Genes universe corrected from 820 to **895** throughout
- §5.1 common feature space corrected from 18,957 to **19,861**
- Supplementary Table 5 v2 adds a **Metric_Definitions** sheet distinguishing
  balanced from micro and macro-per-dataset accuracy, and a `B_Model_Summary` sheet
- The aging-score ratio count is now recoverable from §3.6 (20+20+17+16+14+9 = 96)
- Verified correct in v2: 75.7 / 68.5 / 64.5 / 55.4% balanced accuracies, and the
  stress-phenotype sensitivities 78.9 / 76.3 / 53.9 / 42.1% (pooled recall)

App and repository issues from the first pass — the app serving the superseded V19
model, the container environment not reaching R so the public caps were inactive,
the stale `renv` lockfile, and a crash on empty candidate pools — are all fixed.

---

## Scope of verification

**Re-run from raw counts:** the training pipeline; the ADSIG-180 panel (exact
reproduction); LOOCV 92.71%; macro-F1 0.921; all gene-universe counts; the full
nested-LOSO sweep (all 14 rows); §3.8's balanced accuracies and stress
sensitivities against Supplementary Table 5 v2.

**Read but not executed:** the validation scripts (`005b`, `007b`, `009b`, `010b`,
`014`, SenNet), the figure scripts, `003_competing_models_TwoPanel.r`, and the
supplementary-table compilers — their cohort data is not available. A6 and any
number originating in those scripts remain unchecked.

## Reproducing these checks

```bash
Rscript tests/test_manuscript_numbers.R            # panel, LOOCV, gene counts
Rscript scripts/002c_nested_loso_sweep.R           # the LOSO sweep (~6 min)
Rscript scripts/loso_balanced_accuracy.R           # plain vs balanced accuracy
```
