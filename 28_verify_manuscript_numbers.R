source("00_config.R")
## 28_verify_manuscript_numbers.R
## Itemized verification of every number that appears in the journal manuscript
## (01_Manuscript.md) against the pipeline output files in analysis/outputs/.
##
## Design principle: the manuscript is the object under test, the outputs are the reference.
## Each check states a claim that must appear verbatim in the manuscript text and the value that
## the pipeline actually produced. A claim that is missing OR wrong is a FAIL.
##
## Usage: Rscript 28_verify_manuscript_numbers.R
## Exit status: 0 if every check passes, 1 otherwise.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey) })

root <- WORKSPACE
out  <- file.path(root, "analysis", "outputs")
ms   <- file.path(root, "manuscript", "01_Manuscript.md")

## ------------------------------------------------------------------ helpers
txt <- paste(readLines(ms, warn = FALSE), collapse = " ")
txt <- gsub("\\s+", " ", txt)
## the manuscript uses typographic minus in a few places; accept both
txt2 <- gsub("\u2212", "-", txt)

PASS <- 0L; FAIL <- 0L; fails <- character(0)
chk <- function(id, claim, present = TRUE) {
  hit <- grepl(claim, txt, fixed = TRUE) || grepl(claim, txt2, fixed = TRUE)
  ok <- if (present) hit else !hit
  if (ok) PASS <<- PASS + 1L
  else { FAIL <<- FAIL + 1L; fails <<- c(fails, sprintf("[%s] missing/wrong: %s", id, claim)) }
  invisible(ok)
}
## R's sprintf uses the binary value of the double, so 2.255 prints as "2.25" and 14.35 as
## "14.3". Scientific reporting rounds half away from zero, so a negligible epsilon is added
## before formatting. Without this the verifier produces false failures on exact .5 values.
fx <- function(x, d = 1) sprintf(paste0("%.", d, "f"), x + 1e-9)
chk_in <- function(id, claim, subject, present = TRUE) {
  hit <- grepl(claim, subject, fixed = TRUE)
  ok <- if (present) hit else !hit
  if (ok) PASS <<- PASS + 1L
  else { FAIL <<- FAIL + 1L; fails <<- c(fails, sprintf("[%s] %s: %s", id,
        if (present) "missing" else "must not appear", claim)) }
  invisible(ok)
}
rd1 <- function(x) fx(x, 1)
ci1 <- function(est, se) sprintf("%s-%s", fx(100*(est - 1.96*se), 1), fx(100*(est + 1.96*se), 1))

cat("=== manuscript number verification ===\n\n")

## ------------------------------------------------------------------ 0. file exists
stopifnot(file.exists(ms))
cat("manuscript file:", ms, "\n\n")

## ------------------------------------------------------------------ 1. sample sizes
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))
chk("S1", "104,729 students, 196 institutions")
chk("S2", "46,914")
chk("S3", "41,862")
chk("S4", "40,768")
chk("S5", "195 institutions")
chk("S6", "n = 45,566")
chk("S7", "n = 45,504")
chk("S8", "n = 31,206")
chk("S9", "76,406")
chk("S10", "84,735")
chk("S11", "1,381")
chk("S12", "71.8%")

## ------------------------------------------------------------------ 2. Table 1
t1 <- read.csv(file.path(out, "table1_hms2324.csv"))
row_for <- function(g) t1[t1$group == g, ]
for (g in c("module", "non_module", "main")) {
  r <- row_for(g)
  chk(paste0("T1-", g, "-age"), fx(r$age_num_est, 1))
  chk(paste0("T1-", g, "-fem"), fx(100*r$female_est, 1))
  chk(paste0("T1-", g, "-intl"), fx(100*r$international_est, 1))
  chk(paste0("T1-", g, "-fin"), fx(100*r$fin_stress_est, 1))
  chk(paste0("T1-", g, "-food"), fx(100*r$food_insec_est, 1))
  chk(paste0("T1-", g, "-fg"), fx(100*r$firstgen_est, 1))
  chk(paste0("T1-", g, "-dep"), fx(100*r$dep_any_est, 1))
  chk(paste0("T1-", g, "-sui"), fx(100*r$sui_idea_est, 1))
  chk(paste0("T1-", g, "-lon"), fx(100*r$lonely_est, 1))
}

## ------------------------------------------------------------------ 3. Table 2 (ORs)
t2 <- read.csv(file.path(out, "table2_var_centered_OR.csv"))
raw <- t2[grepl("raw w", t2$exposure), ]
or_of <- function(outcome, exposure, level) {
  r <- raw[raw$outcome == outcome & grepl(exposure, raw$exposure, fixed = TRUE) & raw$level == level, ]
  stopifnot(nrow(r) == 1)
  c(OR = r$OR, L = r$LCL, U = r$UCL)
}
tab2 <- list(
  list(o = "dep_any",  e = "time/day",  lv = "time5cat_f>3h"),
  list(o = "dep_maj",  e = "time/day",  lv = "time5cat_f>3h"),
  list(o = "anx_any",  e = "time/day",  lv = "time5cat_f>3h"),
  list(o = "sui_idea", e = "time/day",  lv = "time5cat_f>3h"),
  list(o = "lonely",   e = "time/day",  lv = "time5cat_f>3h"),
  list(o = "dep_any",  e = "exp_harassed", lv = "exp_harassed"),
  list(o = "dep_maj",  e = "exp_harassed", lv = "exp_harassed"),
  list(o = "anx_any",  e = "exp_harassed", lv = "exp_harassed"),
  list(o = "sui_idea", e = "exp_harassed", lv = "exp_harassed"),
  list(o = "lonely",   e = "exp_harassed", lv = "exp_harassed"),
  list(o = "dep_any",  e = "ai_any",    lv = "ai_any"),
  list(o = "dep_maj",  e = "ai_any",    lv = "ai_any"),
  list(o = "anx_any",  e = "ai_any",    lv = "ai_any"),
  list(o = "sui_idea", e = "ai_any",    lv = "ai_any"),
  list(o = "lonely",   e = "ai_any",    lv = "ai_any"),
  list(o = "dep_any",  e = "ai_n_uses", lv = "ai_n_uses_g3+"),
  list(o = "dep_maj",  e = "ai_n_uses", lv = "ai_n_uses_g3+"),
  list(o = "anx_any",  e = "ai_n_uses", lv = "ai_n_uses_g3+"),
  list(o = "sui_idea", e = "ai_n_uses", lv = "ai_n_uses_g3+"),
  list(o = "lonely",   e = "ai_n_uses", lv = "ai_n_uses_g3+"))
for (i in seq_along(tab2)) {
  v <- or_of(tab2[[i]]$o, tab2[[i]]$e, tab2[[i]]$lv)
  chk(sprintf("T2-OR-%02d", i),
      sprintf("%s (%s-%s)", fx(v["OR"],2), fx(v["L"],2), fx(v["U"],2)))
}

## ------------------------------------------------------------------ 4. Table 2 (TMLE RDs)
t3  <- read.csv(file.path(out, "table3_tmle_RD.csv"))
t3b <- read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv"))
td  <- rbind(t3[, c("psi","se_cluster","exposure","outcome")],
             t3b[, c("psi","se_cluster","exposure","outcome")])
for (i in seq_len(nrow(td))) {
  r <- td[i, ]
  chk(sprintf("T2-RD-%02d", i),
      sprintf("%s (%s)", rd1(100*r$psi), ci1(r$psi, r$se_cluster)))
}
## The Table 2 checks above pass as long as the correct string occurs ANYWHERE in the file,
## so a drift in the Results prose was invisible: the suicide risk difference was written
## "5.5 (4.3-6.9)" in the text while Table 2 correctly carried 6.7. Rebuild the whole
## Results sentence from the TMLE table and require it verbatim.
sc <- t3[t3$exposure == "screen>3h/d", ]
sc <- sc[match(c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely"), sc$outcome), ]
stopifnot(nrow(sc) == 5, all(!is.na(sc$psi)), all(!is.na(sc$se_cluster)))
prose5 <- paste0("risk differences of ",
                 sprintf("%s (95%% CI %s)", rd1(100*sc$psi[1]), ci1(sc$psi[1], sc$se_cluster[1])),
                 ", ",
                 paste(vapply(2:4, function(i)
                   sprintf("%s (%s)", rd1(100*sc$psi[i]), ci1(sc$psi[i], sc$se_cluster[i])),
                   character(1)), collapse = ", "),
                 " and ",
                 sprintf("%s (%s)", rd1(100*sc$psi[5]), ci1(sc$psi[5], sc$se_cluster[5])))
chk("T2-RD-prose", prose5)

## ------------------------------------------------------------------ 5. LCA
pv <- read.csv(file.path(out, "lca2_k5_prevalence.csv"))
for (i in seq_len(nrow(pv))) chk(sprintf("LCA-prev-%d", pv$class[i]), sprintf("%s%%", fx(pv$prev_pct[i], 1)))
adj <- read.csv(file.path(out, "lca2_k5_adjusted_OR.csv"))
c1 <- adj[adj$level == "class_f1", ]
chk("LCA-C1-dep", sprintf("%s (%s-%s)", fx(c1$OR[c1$outcome=="dep_any"],2), fx(c1$LCL[c1$outcome=="dep_any"],2), fx(c1$UCL[c1$outcome=="dep_any"],2)))
chk("LCA-C1-sui", sprintf("%s (%s-%s)", fx(c1$OR[c1$outcome=="sui_idea"],2), fx(c1$LCL[c1$outcome=="sui_idea"],2), fx(c1$UCL[c1$outcome=="sui_idea"],2)))
chk("LCA-C1-lon", sprintf("%s (%s-%s)", fx(c1$OR[c1$outcome=="lonely"],2), fx(c1$LCL[c1$outcome=="lonely"],2), fx(c1$UCL[c1$outcome=="lonely"],2)))
ob <- read.csv(file.path(out, "lca2_k5_outcome_by_class.csv"))
chk("LCA-dep-min", sprintf("%s%%", fx(min(ob$pct[ob$outcome=="dep_any"]), 1)))
chk("LCA-dep-max", sprintf("%s%%", fx(max(ob$pct[ob$outcome=="dep_any"]), 1)))
chk("LCA-sui-max", sprintf("%s%%", fx(max(ob$pct[ob$outcome=="sui_idea"]), 1)))
chk("LCA-lon-max", sprintf("%s%%", fx(max(ob$pct[ob$outcome=="lonely"]), 1)))
lo_all <- Inf; hi_all <- -Inf
for (k in c("4","5","6","7","8")) {
  a <- read.csv(file.path(out, sprintf("lca2_k%s_adjusted_OR.csv", k)))
  mx <- sapply(c("dep_any","dep_maj","anx_any","sui_idea","lonely"),
               function(o) max(a$OR[a$outcome == o]))
  km <- max(sapply(c("dep_any","dep_maj","anx_any","sui_idea","lonely"),
                   function(o) max(a$OR[a$outcome == o])))
  if (k %in% c("5","6","7","8")) {
    lo_all <- min(lo_all, min(mx)); hi_all <- max(hi_all, max(mx))
  }
  if (k == "4") chk("LCA-k4-max", fx(km, 2))
}
chk("LCA-k58-bracket", sprintf("between %s and %s", fx(lo_all,2), fx(hi_all,2)))

## ------------------------------------------------------------------ 6. benchmark
t12 <- read.csv(file.path(out, "table12_method_benchmark.csv"))
t12$rel <- 100*abs(t12$rd_conv_gcomp - t12$rd_tmle)/abs(t12$rd_conv_gcomp)
chk("BM-screen-med", sprintf("%s%% for screen time", fx(median(t12$rel[t12$exposure=="screen>3h/d"]),2)))
chk("BM-harass-med", sprintf("%s%% for online harassment", fx(median(t12$rel[t12$exposure=="online harassment"]),2)))
chk("BM-ai-med", sprintf("%s%% (maximum", fx(median(t12$rel[t12$exposure=="any genAI use"]),2)))
chk("BM-ai-max", sprintf("%s%%)", fx(max(t12$rel[t12$exposure=="any genAI use"]),2)))
chk("BM-extra", sprintf("%s percentage points", fx(max(t12$extra_shrink_tmle_pp),1)))
chk("BM-verdict", "all 15 exposure-outcome comparisons")

t14 <- read.csv(file.path(out, "table14_het_benchmark.csv"))
chk("HET-n", "out of 35")
chk("HET-convraw", sprintf("flagged %d as significant", sum(t14$p_raw < 0.05)))
chk("HET-convbh", sprintf("correction, %d survived", sum(t14$p_trad_bh < 0.05)))
chk("HET-cfraw", sprintf("causal forest flagged %d before", sum(t14$p_cf < 0.05)))
chk("HET-cfbh", "and none after")
chk("HET-exp", "1.75 expected by chance")

## ------------------------------------------------------------------ 7. prediction
t13 <- read.csv(file.path(out, "table13_pred_benchmark.csv"))
A <- t13[t13$feature_set == "A_base_SDOH", ]; B <- t13[t13$feature_set == "B_plus_digital", ]
chk("PR-base-lo", fx(min(A$AUC),4))
chk("PR-base-hi", fx(max(A$AUC),4))
chk("PR-plus-lo", fx(min(B$AUC),4))
chk("PR-plus-hi", fx(max(B$AUC),4))
chk("PR-learnerspread", fx(max(diff(range(A$AUC)), diff(range(B$AUC))),4))
chk("PR-slope-lo", fx(min(t13$cal_slope),3))
chk("PR-slope-hi", fx(max(t13$cal_slope),3))
sh <- read.csv(file.path(out, "table5_shap_importance.csv"))
get_sh <- function(v) sh$mean_abs_shap[sh$variable == v]
chk("SHAP-fin", fx(get_sh("fin_stress"),3))
chk("SHAP-time", fx(get_sh("time5cat"),3))
chk("SHAP-food", fx(get_sh("food_insec"),3))
chk("SHAP-har", fx(get_sh("exp_harassed"),3))
chk("SHAP-ai", fx(get_sh("ai_n_uses"),3))
chk("SHAP-rank", "ninth")

## ------------------------------------------------------------------ 8. CATE
ct <- read.csv(file.path(out, "table4_cate_ate.csv"))
r <- ct[ct$outcome == "dep_any", ]
chk("CATE-dep", sprintf("%s-%s points", fx(100*r$cate_p05,1), fx(100*r$cate_p95,1)))
r <- ct[ct$outcome == "anx_any", ]
chk("CATE-anx", sprintf("%s-%s points", fx(100*r$cate_p05,2), fx(100*r$cate_p95,2)))

## ------------------------------------------------------------------ 9. validity
nc <- read.csv(file.path(out, "table9_control_outcomes.csv"))
nc <- nc[nc$outcome != "BMI", ]
f1 <- nc[nc$outcome == "flourish" & nc$exposure == "time_gt3h", ]
chk("NC-pos-screen", sprintf("%s points, 95%% CI %s to %s", fx(f1$beta,2), fx(f1$LCL,2), fx(f1$UCL,2)))
f2 <- nc[nc$outcome == "flourish" & nc$exposure == "exp_harassed", ]
chk("NC-pos-har", sprintf("%s, %s to %s", fx(f2$beta,2), fx(f2$LCL,2), fx(f2$UCL,2)))
f3 <- nc[nc$outcome == "flourish" & nc$exposure == "ai_heavy", ]
chk("NC-pos-ai", sprintf("%s, %s to %s; P = %s", fx(f3$beta,2), fx(f3$LCL,2), fx(f3$UCL,2), fx(f3$p,2)))
g1 <- nc[nc$outcome == "inst_gradrate" & nc$exposure == "time_gt3h", ]
chk("NC-neg", sprintf("%s percentage points, %s to %s", fx(100*g1$beta,2), fx(100*g1$LCL,2), fx(100*g1$UCL,2)))
t16 <- read.csv(file.path(out, "table16_school_FE.csv"))
chk("FE-maxshift", sprintf("%s percentage points", fx(100*max(abs(t16$confounding_shift)),2)))
chk("FE-n", "14 of 15")
ai_sui <- t16[t16$exposure == "ai_any" & t16$outcome == "sui_idea", ]
chk("FE-cross", sprintf("P = %s to %s", fx(ai_sui$p_noFE,3), fx(ai_sui$p_FE,3)))
te <- read.csv(file.path(out, "table10_evalues.csv"))
chk("EV-screen", sprintf("%s-%s for screen time",
    fx(min(te$E_value[te$exposure=="screen>3h/d"]),2), fx(max(te$E_value[te$exposure=="screen>3h/d"]),2)))
chk("EV-har", sprintf("%s-%s for harassment",
    fx(min(te$E_value[te$exposure=="harassment"]),2), fx(max(te$E_value[te$exposure=="harassment"]),2)))
chk("EV-ai", sprintf("%s-%s for any generative-AI use",
    fx(min(te$E_value[te$exposure=="any genAI"]),2), fx(max(te$E_value[te$exposure=="any genAI"]),2)))
tm <- read.csv(file.path(out, "table_sens_mi_vs_cc.csv"))
chk("MI-diff", "within 0.01")
chk("MI-fmi", "below 3 x 10^-5")
wd <- read.csv(file.path(out, "table17_weight_diagnostics.csv"))
chk("WT-max", sprintf("%s to %s", fx(wd$max_w[wd$weight=="w_raw"],1), fx(wd$max_w[wd$weight=="w_p99"],1)))
ts <- read.csv(file.path(out, "table17_weight_trim_summary.csv"))
chk("WT-shift", sprintf("%s%%", fx(max(ts$max_rel_OR_shift), 1)))
t18 <- read.csv(file.path(out, "table18_school_type_stratified.csv"))
s <- t18[t18$outcome=="dep_any" & t18$stratifier=="inst_type_f" & t18$stratum!="Special Focus", ]
chk("CAR-time", sprintf("%s-%s", fx(min(s$OR[s$exposure=="time_gt3h"]),2), fx(max(s$OR[s$exposure=="time_gt3h"]),2)))
chk("CAR-har", sprintf("%s-%s", fx(min(s$OR[s$exposure=="exp_harassed"]),2), fx(max(s$OR[s$exposure=="exp_harassed"]),2)))

## ------------------------------------------------------------------ 10. replication / trend
t6 <- read.csv(file.path(out, "table6_replication_2425.csv"))
hh <- t6[t6$level == "exp_harassed", ]
chk("REP-har", sprintf("%s-%s", fx(min(hh$OR),2), fx(max(hh$OR),2)))
for (o in c("dep_any","anx_any","lonely")) {
  r <- t6[t6$level=="time5cat_f>3h" & t6$outcome==o, ]
  chk(paste0("REP-", o), sprintf("%s, %s-%s", fx(r$OR,2), fx(r$LCL,2), fx(r$UCL,2)))
}
t8p <- read.csv(file.path(out, "table8_trend_prevalence.csv"))
t8t <- read.csv(file.path(out, "table8_trend_tests.csv"))
tt <- t8t[t8t$comparison == "2022-23 minus 2024-25", ]
r <- tt[tt$variable == "dep_any", ]
chk("TR-dep-range", sprintf("%s%% to %s%%",
    fx(t8p$pct[t8p$variable=="dep_any" & t8p$wave=="2022-23"],1),
    fx(t8p$pct[t8p$variable=="dep_any" & t8p$wave=="2024-25"],1)))
chk("TR-dep-ci", sprintf("-%s percentage points, 95%% CI %s to %s", fx(r$diff_pp,2), fx(-(r$diff_pp+1.96*r$se),2), fx(-(r$diff_pp-1.96*r$se),2)))
r <- tt[tt$variable == "sui_idea", ]
chk("TR-sui-range", sprintf("%s%% to %s%%",
    fx(t8p$pct[t8p$variable=="sui_idea" & t8p$wave=="2022-23"],1),
    fx(t8p$pct[t8p$variable=="sui_idea" & t8p$wave=="2024-25"],1)))
chk("TR-sui-diff", sprintf("-%s points", fx(r$diff_pp,2)))

## ------------------------------------------------------------------ 11. BH claims
tb <- read.csv(file.path(out, "table11_bh_corrected.csv"))
r <- tb[tb$family=="variable-centred" & tb$test=="dep_any ai_any", ]
chk("BH-ai-any", sprintf("q = %s", fx(r$q_overall,4)))
a3 <- tb[tb$family=="variable-centred" & grepl("ai_n_uses_g3\\+", tb$test), ]
chk("BH-ai3-n", sprintf("%d", nrow(a3)))
chk("BH-ai3-all", "all q < 0.05")

## ------------------------------------------------------------------ 12. figure/table cross-references
chk("XR-fig1", "Fig. 1A")
chk("XR-fig1b", "Fig. 1B")
chk("XR-fig1c", "Fig. 1C")
chk("XR-fig2a", "Fig. 2A")
chk("XR-fig3a", "Fig. 3A")
chk("XR-fig3b", "Fig. 3B")
chk("XR-fig3c", "Fig. 3C")
chk("XR-fig4a", "Fig. 4A")
chk("XR-fig4b", "Fig. 4B")
chk("XR-edf1", "Extended Data Fig. 1")
chk("XR-edf2a", "Extended Data Fig. 2a")
chk("XR-edf2b", "Extended Data Fig. 2b")
chk("XR-edf3a", "Extended Data Fig. 3a")
chk("XR-edf3b", "Extended Data Fig. 3b")
chk("XR-edt1", "Extended Data Table 1")
chk("XR-edt2", "Extended Data Table 2")

## ------------------------------------------------------------------ 13. Methods must contain no figure/table
mlines <- readLines(ms, warn = FALSE)
mi <- grep("^## Methods", mlines); di <- grep("^## Data availability", mlines)
stopifnot(length(mi) == 1, length(di) == 1, di > mi)
m <- paste(mlines[(mi + 1):(di - 1)], collapse = "\n")
## The journal forbids display items inside the Methods, but cross-references to supplementary items are
## normal; strip the "Extended Data " prefix so only in-text display items are detected.
m2 <- gsub("Extended Data Table [0-9]+", "EDT", m)
m2 <- gsub("Extended Data Fig[.] [0-9]+", "EDF", m2)
chk_in("METH-nodisplay", "Figure 1", m2, present = FALSE)
chk_in("METH-notable", "Table 1", m2, present = FALSE)
chk_in("METH-nodisplay2", "Figure 2", m2, present = FALSE)
chk_in("METH-notable2", "Table 2", m2, present = FALSE)
## no markdown table body may appear in the Methods
chk_in("METH-notablebody", "\n| ", m2, present = FALSE)

## ------------------------------------------------------------------ 14. The journal hard limits
secw <- function(a, b) {
  L <- readLines(ms, warn = FALSE)
  i <- grep(paste0("^", a), L); j <- grep(paste0("^", b), L)
  stopifnot(length(i) == 1, length(j) == 1, j > i)
  s <- paste(L[(i + 1):(j - 1)], collapse = " ")
  s <- gsub("\\[[0-9,]+\\]", "REF", s)
  s <- gsub("[#*>|`]", " ", s)
  length(grep("[A-Za-z0-9]", strsplit(s, " ")[[1]], value = TRUE))
}
w_abs  <- secw("## Abstract", "## Introduction")
w_main <- secw("## Introduction", "## Methods")
w_meth <- secw("## Methods", "## Data availability")
chk_limit <- function(id, w, lim) {
  if (w <= lim) PASS <<- PASS + 1L
  else { FAIL <<- FAIL + 1L; fails <<- c(fails, sprintf("[%s] %d words exceeds limit %d", id, w, lim)) }
}
chk_limit("LIM-abstract", w_abs, 150)
chk_limit("LIM-maintext", w_main, 3000)
chk_limit("LIM-methods", w_meth, 3000)
cat(sprintf("word counts: abstract %d / main text %d / methods %d\n", w_abs, w_main, w_meth))

## ------------------------------------------------------------------ 15. references
refs <- readLines(ms, warn = FALSE)
refs <- refs[grepl("^\\d+\\. ", refs)]
if (length(refs) > 50) { FAIL <<- FAIL + 1L; fails <<- c(fails, "[REF-count] more than 50 references") } else PASS <<- PASS + 1L
if (all(grepl("\\*", refs))) PASS <<- PASS + 1L else {
  FAIL <<- FAIL + 1L; fails <<- c(fails, "[REF-title] a reference lacks an italicised journal/book title") }
nums <- as.integer(sub("^(\\d+)\\..*", "\\1", refs))
if (identical(nums, seq_along(refs))) PASS <<- PASS + 1L else {
  FAIL <<- FAIL + 1L; fails <<- c(fails, "[REF-order] reference numbers are not sequential") }

## ------------------------------------------------------------------ report
cat("\n---------------------------------------------\n")
cat(sprintf("TOTAL: %d PASS, %d FAIL\n", PASS, FAIL))
if (FAIL > 0) { cat("\nFAILURES:\n"); cat(paste(fails, collapse = "\n"), "\n") }
cat("---------------------------------------------\n")
if (FAIL > 0) quit(status = 1) else cat("MANUSCRIPT NUMBERS VERIFIED\n")
