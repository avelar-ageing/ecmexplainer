# Analysis scripts

The original manuscript pipeline, kept here so every published number can be
regenerated. Run them from the repository root — they read from `data/`.

| Script | What it produces |
|---|---|
| `003_final_train_two_panels.R` | The deployed models: ADSIG-180 and POOL-120, with LOOCV, macro-F1, panel definitions and figures |
| `002c_nested_loso_sweep.R` | The nested leave-one-study-out sweep over panel sizes 60–180 — the source of Supplementary Table 5C and of the 57.3% figure |
| `002c_saturation.R` | The earlier saturation sweep |
| `loso_balanced_accuracy.R` | Recomputes the sweep metrics separating plain from balanced accuracy (see below) |

## Reproduced values

`002c_nested_loso_sweep.R` was re-run here and reproduces Supplementary Table 5C
exactly, all 14 rows — ADSIG peak 57.29% at 180 ratios, range 45.83–57.29%.

`003_final_train_two_panels.R` reproduces the shipped `TrainedModel_ADSIG_180.rds`
exactly. `tests/test_manuscript_numbers.R` verifies the same thing through the
app's own `R/ecm_core.R`.

## Which accuracy the sweep reports

The sweep writes its `LOSO` column as `mean(loso_pred == phenotype)` — **plain,
unweighted accuracy**. That is not the same as balanced accuracy, and on this
cohort the difference is large: 42% of the 96 samples are Proliferating, which is
also the best-recovered class, so the unweighted mean sits well above the mean of
per-class recalls.

```
ADSIG-180   plain accuracy 57.29%   balanced accuracy 37.60%
```

`loso_balanced_accuracy.R` prints both for every panel size, plus per-class
recalls and the number of source studies behind each phenotype. Be explicit about
which one you are quoting.

## Two notes on running them

`003_final_train_two_panels.R` errors after its summary table: lines 363 onward
are dead code referencing out-of-scope variables and absent files. Everything up
to the summary is sound; stop there.

The two scripts do not compute η² the same way: the sweep uses the ratio's sign
vector, the final training uses the continuous log-ratio. In both, η² is the
admission filter only — the greedy selection is ranked by TSP delta.
