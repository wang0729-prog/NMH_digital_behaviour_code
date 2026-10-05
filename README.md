# Digital behaviour profiles, generative AI use and mental health in US college students

Analysis code for the manuscript above. Two cohorts, four modelling layers:
latent class analysis, doubly robust estimation, causal forests, and
gradient-boosted prediction -- each paired with a conventional comparator run
on the same respondents, the same covariates and the same survey weights.

## Data availability

**No individual-level data is included in this repository, and none of it can
be.** Both cohorts are released under a data-use agreement that forbids
redistribution:

| Cohort | Waves used | How to obtain |
|---|---|---|
| Healthy Minds Study (HMS) | 2022-23, 2023-24, 2024-25 | https://healthymindsnetwork.org/hms/ -- institutional registration |
| China Family Panel Studies (CFPS) | 2020, 2022 | https://www.isss.pku.edu.cn/cfps/ -- application |

`00_config.R` therefore reads two directories that you must supply yourself.
Scripts that need them stop with an explicit message rather than failing deep
inside `read_dta()`.

The `results/` directory holds the aggregate numbers that appear in the paper --
odds ratios, risk differences, prevalences, test counts. Every table there is a
summary, not a respondent list: the largest is 120 rows.

## Layout

```
00_config.R          paths and environment; sourced first by every script
01_clean.R ... 31_   analysis pipeline, numbered in run order
99*.R                verification: recompute manuscript numbers from the outputs
probe_*.R, test_*.R  small compatibility probes
results/             aggregate output tables (CSV)
```

## Environment

R 4.4.3 on Windows. Package versions below are read off the analysis machine
when this README is generated:

| `dplyr` | 1.2.1 |
| `ggplot2` | 4.0.3 |
| `glmnet` | 4.1.10 |
| `grf` | 2.6.1 |
| `gtable` | 0.3.6 |
| `haven` | 2.5.5 |
| `lmtest` | 0.9.40 |
| `mice` | 3.19.0 |
| `openxlsx` | 4.2.8.1 |
| `pdftools` | 3.8.0 |
| `poLCA` | 1.6.0.2 |
| `sandwich` | 3.1.3 |
| `scales` | 1.4.0 |
| `survey` | 4.5 |
| `svglite` | 2.2.2 |
| `tidyr` | 1.3.2 |
| `xgboost` | 3.2.1.1 |
| `R` | 4.4.3 |

Two notes on the machine-learning stack. The TMLE in `06_tmle.R` uses a
**hand-coded discrete SuperLearner** (main-effects GLM, GLM with interactions,
XGBoost) with cluster-aware cross-fitting by school, not the `SuperLearner`
package -- so that package is *not* a dependency. `xgboost` is used directly
through its native interface.

## Running the code

```r
# 1. point the config at your own directories
Sys.setenv(HMS_DATA  = "path/to/HMS",        # HMS_*_PUBLIC_instchars.dta
           CFPS_DATA = "path/to/CFPS",       # cfps<year>person.dta
           HMS_OUT   = "path/for/outputs")
source("00_config.R")

# 2. run in numeric order; each script logs to <HMS_OUT>/<name>_log.txt
source("01_clean.R"); source("02_descriptive.R"); ...
```

Equivalently, export those three variables as environment variables before
starting R. `HMS_ANALYSIS` defaults to the working directory.

Run order and purpose:

| Script | Purpose |
|---|---|
| `01_clean.R` | Builds the analysis frames from the raw files: HMS 2023-24 and 2024-25 (screen time, online experiences, generative-AI uses, five outcomes, SDOH, covariates) and CFPS 2020 and 2022. |
| `02_descriptive.R` | Table 1: weighted sample characteristics of the HMS 2023-24 module sample, overall and by screen-time category. |
| `03_lca.R` | Latent class analysis of digital-behavior profiles, HMS 2023-24 main sample. Indicators (12): time5cat (ordinal 0-4, treated as polytomous), |
| `03b_lca_robust.R` | Robust LCA: fit k=1..8, save each model + stats IMMEDIATELY after fitting so a mid-run kill never loses completed models. |
| `03c_lca_noaiother.R` | Robust LCA: fit k=1..8, save each model + stats IMMEDIATELY after fitting so a mid-run kill never loses completed models. |
| `04_var_centered.R` | Variable-centered weighted associations: digital-behavior exposures -> mental-health outcomes HMS 2023-24 main analytic sample. |
| `04b_sens_missing.R` | Sensitivity to listwise deletion of covariates (main CC n=41,862 vs full 45,504). Variant A: missing-category (firstgen / race_cat / female / international get explicit |
| `05_lca_outcomes.R` | Correction: the row index used to re-attach LCA classes to the survey object was obtained with as.integer(rownames(dd)). That is only correct while `d` is a plain |
| `06_tmle.R` | Doubly-robust cross-fitted TMLE for the effect of digital-behavior exposures on mental-health outcomes, HMS 2023-24. |
| `07_cate.R` | Heterogeneous treatment effects via causal forests (grf). Treatment: A = time_gt3h (>3 h/d non-academic screen time); Outcomes: 5 binary MH outcomes. |
| `08_prediction.R` | Incremental predictive value of digital-behaviour features for mental-health outcomes. Design: outer 5-fold CV grouped by school (PSU) -> honest AUC / Brier / calibration; |
| `09_replication_cfps.R` | Correction: the CFPS analysis extract never carried a cluster identifier, so every svydesign below was silently built with ids = ~1 (no clustering) while the |
| `11_export_xlsx.R` | 11_export_figdata_xlsx.R Export figure source data + key result tables to a single Excel workbook |
| `12_sens_mi.R` | Sensitivity: multiple imputation of missing covariates (vs complete-case n=41,862 and missing-category n=44,679). m = 20 imputations; weighted svyglm (schools = PSU) |
| `13_evalue_bh.R` | (1) E-values for the variable-centred odds ratios (unmeasured-confounding sensitivity). RR approximated from OR and the outcome prevalence in the unexposed |
| `14_negcontrol.R` | Negative- and positive-control outcomes to probe residual confounding. Positive control : flourish (flourishing score) -- expected to DECREASE with risk |
| `15_trend_3waves.R` | Background outcome trends across the three HMS waves (2022-23 / 2023-24 / 2024-25). 2022-23 has NO digital-behaviour module -> descriptive only, no exposure analysis. |
| `16_method_benchmark.R` | Methodological benchmark for the "AI-driven analysis" claim. Question: does ML-based causal inference (TMLE + discrete SuperLearner, cross-fitted by |
| `16b_summarize.R` | Prints the relative divergence between the conventional g-computation and TMLE risk differences in table12, the number quoted as '11.05% versus 0.70% median divergence'. |
| `17_pred_benchmark.R` | Prediction layer of the methodological benchmark. Question: on this prediction task, does gradient boosting actually beat the conventional |
| `18_het_benchmark.R` | Heterogeneity layer of the methodological benchmark. Question: the conventional route to effect modification is a pre-specified product term in |
| `19_school_panel.R` | School-level two-period panel with school fixed effects (2023-24 -> 2024-25). Why: the negative-control outcome (institutional graduation rate) showed that the |
| `19a_check_compat.R` | Pre-flight for the school-level two-period fixed-effects panel. Verify that (a) the exposure item carries the same wording in both waves, |
| `19b_check_items.R` | Which exposure item is genuinely comparable across 2023-24 and 2024-25? Print full value labels for internet_1 (screen time) and all internet_2_* items. |
| `19c_check_time.R` | Is the screen-time item comparable across waves? Compare the raw internet_1 distribution in the SAME 132 overlapping schools. |
| `20_school_fe.R` | School fixed effects at the INDIVIDUAL level, 2023-24 only. Why this design instead of the two-period school panel (19_school_panel.R): |
| `22_final_figures.R` | Reorganised figure set for the manuscript. Problems with the previous set of 13 single-panel figures: |
| `23_export_png_preview.R` | Rasterises the PDF figures to 150 dpi PNG previews for quick review. Not part of the deliverable: the paper is submitted as SVG/PDF only. |
| `24_sens_weight_trim.R` | Sensitivity to extreme non-response weights (nrweight). Question addressed: HMS nrweight spans 0.21-76.7 (max/median = 67.5). If a handful of |
| `25_sens_school_type.R` | Effect modification by institutional context (stratified analysis + formal test). Why this matters: the school-level negative control (graduation rate) flagged residual |
| `26_cfps_panel.R` | Individual-level two-wave panel, CFPS 2020 -> 2022 (external-validation cohort). Purpose: every HMS estimate is cross-sectional, so time-invariant person-level |
| `27_central_illustration.R` | Central Illustration / graphical abstract. WHY THIS IS A SEPARATE FILE AND NOT A MAIN FIGURE: |
| `28_verify_manuscript_numbers.R` | Itemized verification of every number that appears in the journal manuscript (01_Manuscript.md) against the pipeline output files in analysis/outputs/. |
| `30_full_reaudit.R` | Independent verification: every number quoted in the paper is rebuilt from the raw files. Rule for this script: NOTHING is read from analysis/outputs except the manuscript text. |
| `31_lca_3step_globalbh.R` | Closes two methodological gaps: (a) LCA low entropy (0.603) -> does modal (hard) class assignment bias the class-outcome |
| `99_verify.R` | 99_verify.R -- independent recomputation of headline numbers (second code path) |
| `99b_check5.R` | Unit check of the model matrix used by the prediction layer: builds it on a 20,000-row subsample and confirms no row drops to NA. |
| `99c_verify2.R` | 99c_verify2.R -- independent recomputation of the new sensitivity / trend results |
| `99d_verify3.R` | Independent re-computation of the three methodological benchmarks (Tables 12-14). Reads only the written CSVs and re-derives every headline number from scratch. |
| `99e_verify4.R` | Verification for the school-confounding work: (1) the within transformation is arithmetically correct (school-demeaned variables |
| `99f_verify_figures.R` | Independent verification of the figure set produced by 22_final_figures.R (plus the Central Illustration from 27_central_illustration.R). |
| `99g_verify_sens.R` | Independent verification of 24 (weight trimming), 25 (school-type stratification) and 26 (CFPS individual-level panel). |
| `check_pkgs.R` | Reports which of the candidate packages are installed on this machine. |
| `probe_nc.R` | Probe: distributions of the negative-control candidate variables (height, BMI, graduation rate, flourishing) in HMS 2023-24. |
| `probe_vars.R` | Probe: which outcome and covariate names exist in the HMS 2022-23 file, used to decide what the trend analysis could carry. |
| `test_auc.R` | Unit test of the weighted AUC function: an informative predictor must score ~0.75 and a random one 0.50. |

## Verification

The `99*` and `3*` scripts are the self-check layer: they recompute every number
quoted in the manuscript from the saved outputs and print a PASS/FAIL line per
assertion. `30_full_reaudit.R` is the full pass over the main text;
`28_verify_manuscript_numbers.R` covers the earlier draft. Run them after the
pipeline:

```r
source("30_full_reaudit.R")     # prints N PASS / N FAIL
source("31_lca_3step_globalbh.R")
```

Random components (LCA starts, cross-fitting folds, imputation, causal forest)
are seeded at the top of each script.

## Deliberately not included

| Excluded | Reason |
|---|---|
| `*.rds` | fitted `poLCA` objects carry the 45,566 x 5 posterior matrix; class-assignment files carry a class per respondent |
| `hms_2324_module.csv`, `cfps_2022_young_hiedu.csv`, `cfps_panel_2020_2022.csv` | one row per respondent |
| `table15_school_panel_wide.csv` | one row per institution, keyed by the real institution id, carrying that institution's mental-health prevalence. The published estimates are in `table15_school_panel.csv` |
| `*.dta`, `*.sav` | restricted survey files, blocked by `.gitignore` |
| `*_log.txt`, `figures/` | regenerated on every run |

`.gitignore` blocks these patterns as well, so a stray `outputs/` copy cannot be
committed by accident.

## Ethics

Both surveys obtained informed consent and institutional review board approval
under their own protocols; this secondary analysis used de-identified public-use
files and was exempt from further review. No individual is identifiable from
anything in this repository.

## License

MIT -- see `LICENSE`. It covers the code only, and says nothing about the
survey data, which remains governed by its own agreements.
