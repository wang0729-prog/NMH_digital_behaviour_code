source("00_config.R")
## 11_export_figdata_xlsx.R
## Export figure source data + key result tables to a single Excel workbook
## (user requirement: raw data available for independent re-plotting).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(openxlsx) })
options(openxlsx.tempdir = OPENXLSX_TMP)   # ASCII temp dir: avoids the Chinese-path failure

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")

wb <- createWorkbook()

add <- function(sheet, obj) {
  if (is.null(obj) || nrow(obj) == 0) return(invisible(NULL))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, as.data.frame(obj))
  invisible(NULL)
}

## Sheet names follow the REORGANISED figure set (22_final_figures.R):
## F1A/F1B/F1C, F2D/F2E/F2F, F3G/F3H/F3I/F3J, FS1, FS2a/FS2b, FS3a/FS3b
add("F1A_LCA_condprob",    read.csv(file.path(out, "lca2", "condprob_k5.csv")))
add("F1A_profile_prev",    read.csv(file.path(out, "lca2_k5_prevalence.csv")))
add("F1B_outcomes_class",  read.csv(file.path(out, "lca2_k5_outcome_by_class.csv")))
add("F1B_adjusted_OR",     read.csv(file.path(out, "lca2_k5_adjusted_OR.csv")))
add("F1C_dose_response",   read.csv(file.path(out, "table2_var_centered_OR.csv")))
add("F2D_TMLE_RD",         read.csv(file.path(out, "table3_tmle_RD.csv")))
add("F2D_TMLE_aiheavy",    read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv")))
add("F2E_method_bench",    read.csv(file.path(out, "table12_method_benchmark.csv")))
add("F2F_CATE_ATE",        read.csv(file.path(out, "table4_cate_ate.csv")))
add("F2F_CATE_BLP",        read.csv(file.path(out, "table4_cate_blp.csv")))
add("F3G_het_bench",       read.csv(file.path(out, "table14_het_benchmark.csv")))
add("F3G_trad_inter",      read.csv(file.path(out, "table14b_trad_interactions.csv")))
add("F3H_pred_bench",      read.csv(file.path(out, "table13_pred_benchmark.csv")))
add("F3H_pred_bench_DCA",  read.csv(file.path(out, "table13_pred_dca.csv")))
add("F3I_control_outcomes", read.csv(file.path(out, "table9_control_outcomes.csv")))
add("F3J_school_FE",       read.csv(file.path(out, "table16_school_FE.csv")))
add("FS1_trend_prevalence", read.csv(file.path(out, "table8_trend_prevalence.csv")))
add("FS1_trend_tests",     read.csv(file.path(out, "table8_trend_tests.csv")))
add("FS1_trend_samples",   read.csv(file.path(out, "table8_trend_samples.csv")))
add("FS2a_SHAP",           read.csv(file.path(out, "table5_shap_importance.csv")))
add("FS2b_DCA",            read.csv(file.path(out, "table5_dca.csv")))
add("FS2b_prediction_AUC", read.csv(file.path(out, "table5_prediction.csv")))
add("FS3a_replication_2425", read.csv(file.path(out, "table6_replication_2425.csv")))
add("FS3b_CFPS_validation",  read.csv(file.path(out, "table7_cfps_validation.csv")))
add("Table1_sample",       read.csv(file.path(out, "table1_hms2324.csv")))
add("Sens_missing",        read.csv(file.path(out, "table_sens_missing.csv")))
add("Sens_MI_vs_CC",       read.csv(file.path(out, "table_sens_mi_vs_cc.csv")))
add("Evalues",             read.csv(file.path(out, "table10_evalues.csv")))
add("BH_corrected",        read.csv(file.path(out, "table11_bh_corrected.csv")))
add("Control_outcomes",    read.csv(file.path(out, "table9_control_outcomes.csv")))
add("Sens_lca2_k6_OR",     read.csv(file.path(out, "lca2_k6_adjusted_OR.csv")))
add("Sens_lca2_k7_OR",     read.csv(file.path(out, "lca2_k7_adjusted_OR.csv")))
add("Sens_lca2_k8_OR",     read.csv(file.path(out, "lca2_k8_adjusted_OR.csv")))
add("Sens_lca_k4_OR",      read.csv(file.path(out, "lca_k4_adjusted_OR.csv")))
add("Sens_lca_k6_OR",      read.csv(file.path(out, "lca_k6_adjusted_OR.csv")))
add("Arch_school_panel",   read.csv(file.path(out, "table15_school_panel.csv")))

## Robustness checks (scripts 24 / 25 / 26, verified by 99g):
add("Sens17_weight_diag",    read.csv(file.path(out, "table17_weight_diagnostics.csv")))
add("Sens17_weight_trim_OR", read.csv(file.path(out, "table17_weight_trim_OR.csv")))
add("Sens17_weight_summary", read.csv(file.path(out, "table17_weight_trim_summary.csv")))
add("Sens18_school_strata",  read.csv(file.path(out, "table18_school_type_stratified.csv")))
add("Sens18_school_inter",   read.csv(file.path(out, "table18_school_type_interaction.csv")))
add("Sens18_school_spread",  read.csv(file.path(out, "table18_school_type_spread.csv")))
add("Sens19_CFPS_panel_FE",  read.csv(file.path(out, "table19_cfps_panel_FE.csv")))
add("Sens19_CFPS_compare",   read.csv(file.path(out, "table19_cfps_panel_compare.csv")))
add("Sens19_CFPS_reconcile", read.csv(file.path(out, "table19_cfps_panel_reconcile.csv")))

saveWorkbook(wb, file.path(dir, "figure_source_data.xlsx"), overwrite = TRUE)
cat("XLSX DONE:", file.path(dir, "figure_source_data.xlsx"), "\n")
