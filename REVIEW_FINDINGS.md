# ECMExplainer — outstanding manuscript/code issues

**Date:** 4 August 2026
**Manuscript:** Sheekey E, Avelar RA, Fujiwara H. *ECMExplainer: a functional
matrisome framework for decoding cellular stress states and aging trajectories.*
**Reviewed:** `Rcode_ECMExplainer/` (app + k-TSP scripts), `final_version/`,
`extra_for_claude/` (nested-LOSO sweep, figure scripts, Supplementary Tables 1–5).

Everything below is **still open**. Items already fixed in the app and repository
have been removed.

## Summary

The model itself is sound. Re-running the training pipeline from raw counts
reproduces the published ADSIG-180 panel **exactly** — all 180 ratios in the same
order, LOOCV 92.71%, macro-F1 0.921, 76 A/D + 104 SIG, and every gene-universe
count (1,139 curated; 923 detected; 711 modelling universe = 303 A / 88 D /
320 SIG). The nested leave-one-study-out sweep also reproduces Supplementary
Table 5C exactly, all 14 rows.

The issues below are therefore about **what the paper says about the code**, not
about whether the model works. One of them (A1) changes a headline number.

---

## A. Issues that change reported results

### A1. §5.4 — "balanced accuracy" is actually plain accuracy

The sweep computes LOSO as `mean(loso_pred == phenotype, na.rm = TRUE)`, an
unweighted count of correct calls. The Methods report the result as "a peak LOSO
balanced accuracy of 57.3%".

| ADSIG-180 LOSO | Value |
|---|---|
| Plain accuracy (what the code computes) | **57.29%** |
| Balanced accuracy (mean per-class recall) | **37.60%** |
| Balanced accuracy, ≥3-study phenotypes only | 56.39% |

Cause: the cohort is 42% Proliferating (40/96) and Proliferating is the
best-recovered class, so the unweighted mean sits well above the mean of per-class
recalls. Reproduce with `scripts/loso_balanced_accuracy.R`.

**Options:** relabel 57.3% as overall accuracy; report 37.6% as the balanced
figure; or state the ≥3-study restriction explicitly, under which 56.39% is
defensible and close to what is written. §5.4 already mentions the ≥2/≥3-study
restriction, so the last is probably what was intended — but as written, 57.3% is
presented as the all-six-phenotype number.

### A2. Two phenotypes are never recovered cross-study

Per-class LOSO recall at ADSIG-180: Proliferating 92.5%, Quiescence 53.3%,
Fibrosis 44.4%, Senescence 35.3%, **Hypoxia 0.0%, Inflammation 0.0%**. Hypoxia and
Inflammation each have only two source studies. The paper does not state that two
of its six reference states have zero cross-study recall.

### A3. "Balanced accuracy" denotes two different quantities

§3.8's 77.1 / 70.0 / 65.0 / 56.3% **are** correct balanced accuracies — verified as
the mean of the four per-phenotype recalls in Supplementary Table 5B. §5.4's 57.3%
is not. The same term is used for both.

### A4. The LOSO figure describes a different panel from the deployed one

`002c_nested_loso_sweep.R` sets `MODULE_FLOOR <- 5L`;
`003_FinalTrain_TwoPanels.R` has no such constraint. Consequently:

- Sweep ADSIG-180 = **77 A/D + 103 SIG**
- Deployed ADSIG-180 = **76 A/D + 104 SIG**

The manuscript reports 76/104 and attributes the 57.3% LOSO to it. Nearly the same
panel, but not the same one.

### A5. §5.3 — selection is ranked by TSP delta, not η²

The text says the greedy selector adds ratios "in descending η² order". Both
scripts rank by `TSP_Delta` — the absolute difference in the proportion of positive
signs inside versus outside each phenotype — and use η² only as a > 0.10
pre-filter. This is not cosmetic: at 6 ratios per phenotype the two criteria select
**completely disjoint panels, 0 of 36 ratios shared**.

### A6. The two scripts define η² differently, and only one matches the text

§5.3 defines η² on the ratio's *sign vector*. `002c_nested_loso_sweep.R` does that
(between-group variance of `P(A > B)`). `003_FinalTrain_TwoPanels.R`, which
produced the released model, computes η² on the *continuous* log-ratio. The
description matches the sweep, not the deployed model.

### A7. §5.3 — the "at least five ratios per module" constraint is not in the final training

It exists only in the sweep. The published 76/104 split satisfies it incidentally,
not by design.

### A8. §5.4 — ECM-Genes baseline gene count is wrong

Text says 820 genes. The filter in `003_competing_models_TwoPanel.r` (curated list
∩ matrix, then variance > 0.005) yields **895**. 820 appears nowhere in
Supplementary Table 1 either. Note that 895 is also the number §5.2 gives for the
batch-corrected universe.

### A9. §5.1 — "a common feature space of 18,957 genes" is wrong

The common matrix holds **19,861** genes. 18,957 is the count *after* the
variance > 0.005 filter — the All-Genes baseline figure from §5.4, borrowed into
§5.1.

### A10. LOOCV cannot separate the four models

POOL-120 100%, All-Genes 100%, ECM-Genes 100%, ADSIG-180 92.71%. On LOOCV the ratio
model is the *worst* of the four; the model comparison rests entirely on LOSO and
the 226-sample validation set.

### A11. §5.4 — the "52–57% vs 49%" comparison may mix statistics *(unverified)*

52–57% matches the sweep's plain LOSO accuracy range for ADSIG. The baseline
figures (All-Genes 28→45%, ECM-Genes 40→49%) come from the separate
reference-projection computation in `003_competing_models_TwoPanel.r`, which was
not run here. If those are balanced accuracies, the comparison is
apples-to-oranges. Worth checking before submission.

---

## B. Code defects

### B1. `003_FinalTrain_TwoPanels.R` lines 363–679 are dead code that errors

They reference out-of-scope variables (`model_tag`, `PANEL`, `ratio_mat_sign`,
`cm`, `prec`, `rec`, `f1`) and files that do not exist
(`Final_Ratios_ADSIG_BALANCED*.csv`, `perSample_geneExpression_long*.csv`).
Sourcing the script whole fails after the summary table.

### B2. The "Top-3 ratio pairs per phenotype" figure does not rank by any score

In the same block: `slice_max(order_by = Ratio, n = 3)`. `Ratio` is the
`"GENEA / GENEB"` string, so this selects the alphabetically-last three names per
phenotype, not the top three by any statistic. The code comment even says "swap to
TSP_Delta to rank by score". If that figure reached the manuscript, its labelling
is wrong.

### B3. `final_version/` cannot run as shipped

`003` loads `salmon_counts_Combined_FINAL.rds` and
`salmon_metadata_Combined_FINAL.csv`; the folder contains `salmon_raw_counts.rds`
and `salmon_metadata.csv` instead. The real inputs had to be taken from the older
bundle.

---

## C. Reproducibility and publication readiness

### C1. §6 Data Availability is an unfinished placeholder

"ECMExplainer can be found on GitHub…" — no URL, no accession, no DOI.

### C2. The GitHub repository is private

Reviewers cannot currently see the code or the app source.

### C3. No archived snapshot

No Zenodo DOI and no published container image, so the artefact does not outlive
the server hosting it.

### C4. The aging-score ratio count is not stated

Supplementary Table 4 shows 96 of 180 ratios passed BH FDR < 0.05; §5.9 gives no
number.

---

## Scope of verification

**Verified by re-running from raw counts:** the training pipeline; the ADSIG-180
panel (reproduces exactly, all 180 ratios in the same order); LOOCV 92.71%;
macro-F1 0.921; every gene-universe count; the full nested-LOSO sweep (all 14
rows); and §3.8's benchmark balanced accuracies against Supplementary Table 5B.

**Read but not executed:** the validation scripts (`005b`, `007b`, `009b`, `010b`,
`014`, SenNet), the figure scripts, `003_competing_models_TwoPanel.r`, and the
supplementary-table compilers. Their cohort data was not available, so A11 and any
numbers originating in those scripts remain unchecked.

## Reproducing these checks

```bash
Rscript tests/test_manuscript_numbers.R            # panel, LOOCV, gene counts
Rscript scripts/002c_nested_loso_sweep.R           # the LOSO sweep (~6 min)
Rscript scripts/loso_balanced_accuracy.R           # plain vs balanced accuracy
```
