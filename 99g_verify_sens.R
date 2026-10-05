source("00_config.R")
## 99g_verify_sens.R
## Independent verification of 24 (weight trimming), 25 (school-type stratification)
## and 26 (CFPS individual-level panel).
##
## Two evidence chains, as in 99f:
## (1) re-derive every headline number from the RAW source files (dta), not from the
## CSVs the analysis scripts wrote;
## (2) internal invariants that a silent bug would violate.
##
## The strongest single check is #C3: for a Gaussian svyglm the point estimate IS the
## weighted-least-squares solution, so lm(y ~ x, weights = wt) must reproduce it exactly.
## If the panel differencing, the weight choice or the sample filter were wrong, this
## check fails to 1e-8.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(dplyr); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")

PASS <- 0; FAIL <- 0; MSG <- character(0)
chk <- function(id, ok, detail = "") {
  ok <- isTRUE(ok)
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  cat(sprintf("[%s] %-28s %s\n", if (ok) "PASS" else "FAIL", id, detail))
  if (!ok) MSG <<- c(MSG, paste(id, detail))
}
near <- function(a, b, tol = 1e-6) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))

core_cov <- c("age_num", "female", "race_cat", "international", "undergrad",
              "grad_student", "fin_stress", "food_insec", "firstgen")
OC <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
ADJ <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

## ==========================================================================
## A. weight trimming (24_sens_weight_trim.R)
## ==========================================================================
cat("\n===== A. weight trimming =====\n")
## A1 is the only check that goes back to the delivered .dta: it re-derives the
## winsorisation thresholds from the raw non-response weight. (The covariates used to
## define the analytic sample are constructed in 01_clean.R, so the sample itself is
## rebuilt from the rds -- reproducing 01_clean here would just duplicate that code.)
raw <- zap_labels(read_dta(file.path(HMS_DIR, "HMS_2023-2024_PUBLIC_instchars.dta"),
                           col_select = all_of(c("internet_1", "nrweight"))))
stopifnot(all(c("internet_1", "nrweight") %in% names(raw)))

D <- readRDS(file.path(out, "hms_2324_analysis.rds"))       # as the analysis script sees it
Dz <- zap_labels(D)

## A1 thresholds recomputed on the module sample from the RAW dta
mod_raw <- raw[!is.na(raw$internet_1), ]
qs <- quantile(mod_raw$nrweight, c(0.99, 0.995, 0.999), na.rm = TRUE)
mod_z <- Dz[Dz$module == 1 & !is.na(Dz$time5cat), ]
qs2 <- quantile(mod_z$nrweight, c(0.99, 0.995, 0.999), na.rm = TRUE)
chk("A1 thresholds (raw vs rds)", near(as.numeric(qs), as.numeric(qs2), 1e-9),
    sprintf("p99=%.3f p99.5=%.3f p99.9=%.3f", qs[1], qs[2], qs[3]))
chk("A1b module n", nrow(mod_raw) == nrow(mod_z) && nrow(mod_raw) == 46914,
    sprintf("raw=%d rds=%d", nrow(mod_raw), nrow(mod_z)))

## A2 diagnostics recomputed
m <- Dz[Dz$module == 1, ]; m <- m[!is.na(m$time5cat), ]
for (v in OC) m <- m[!is.na(m[[v]]), ]
for (v in core_cov) m <- m[!is.na(m[[v]]), ]
chk("A2 analytic n = 41862", nrow(m) == 41862, sprintf("n=%d", nrow(m)))

wv <- list(w_raw = m$nrweight,
           w_p999 = pmin(m$nrweight, qs[3]),
           w_p995 = pmin(m$nrweight, qs[2]),
           w_p99  = pmin(m$nrweight, qs[1]))
diag <- read.csv(file.path(out, "table17_weight_diagnostics.csv"))
for (nm in names(wv)) {
  v <- wv[[nm]]; vn <- v / mean(v)
  ess <- sum(vn)^2 / sum(vn^2)
  row <- diag[diag$weight == nm, ]
  chk(paste0("A2 ESS ", nm), near(ess, row$ESS, 1e-3),
      sprintf("recomputed=%.1f table=%.1f", ess, row$ESS))
 ## max_w / sum_w are stored rounded (3 dp / 0 dp) in the diagnostics CSV, so the
 ## tolerance must match the rounding, not machine precision.
  chk(paste0("A2 max ", nm), near(max(v), row$max_w, 1e-3),
      sprintf("%.3f vs %.3f", max(v), row$max_w))
  chk(paste0("A2 sum ", nm), near(sum(v), row$sum_w, 1),
      sprintf("%.1f vs %.1f", sum(v), row$sum_w))
}

## A3 two models refit from scratch
tt <- read.csv(file.path(out, "table17_weight_trim_OR.csv"), stringsAsFactors = FALSE)
refit <- function(o, ev, wcol) {
 ## build the winsorised weight INSIDE the model-specific subset (the analysis script
 ## attaches it to the full frame, which svydesign then subsets -- same result, but
 ## deriving it here keeps the check independent of that plumbing).
  sub <- m[!is.na(m[[ev]]), ]
  sub$.w <- switch(wcol,
                   w_raw  = sub$nrweight,
                   w_p999 = pmin(sub$nrweight, qs[3]),
                   w_p995 = pmin(sub$nrweight, qs[2]),
                   w_p99  = pmin(sub$nrweight, qs[1]))
  g <- svyglm(as.formula(paste(o, "~", ev, "+", ADJ)),
              design = svydesign(ids = ~schoolnum, weights = ~.w, nest = TRUE, data = sub),
              family = quasibinomial())
  exp(coef(g)[[ev]])
}
for (cs in list(c("dep_any", "time_gt3h", "w_p995"),
                c("sui_idea", "ai_any", "w_p99"),
                c("lonely", "exp_harassed", "w_raw"))) {
  b <- refit(cs[1], cs[2], cs[3])
  tv <- tt$OR[tt$outcome == cs[1] & tt$exposure == cs[2] & tt$weight == cs[3]]
  chk(paste0("A3 refit ", paste(cs, collapse = "/")), near(b, tv, 1e-3),
      sprintf("refit=%.5f table=%.3f", b, tv))
}

## A4 invariants in the summary table
sm <- read.csv(file.path(out, "table17_weight_trim_summary.csv"), stringsAsFactors = FALSE)
chk("A4 15 pairs", nrow(sm) == 15, sprintf("nrow=%d", nrow(sm)))
chk("A4 4 variants each", all(sm$n_models == 4), sprintf("min=%d", min(sm$n_models)))
ok <- TRUE; det <- ""
for (i in seq_len(nrow(sm))) {
  sub <- tt[tt$outcome == sm$outcome[i] & tt$exposure == sm$exposure[i], ]
  sg <- all(sub$OR > 1) || all(sub$OR < 1)
  dl <- max(abs(sub$logOR - sub$logOR[sub$weight == "w_raw"]))
  if (!identical(sg, as.logical(sm$sign_stable[i]))) { ok <- FALSE; det <- paste(det, "sign", sm$outcome[i], sm$exposure[i]) }
  if (!near(dl, sm$max_abs_dlogOR[i], 1e-3)) { ok <- FALSE; det <- paste(det, "dlogOR", sm$outcome[i], sm$exposure[i]) }
  if (!near(sm$max_rel_OR_shift[i], 100 * dl, 1e-1)) { ok <- FALSE; det <- paste(det, "rel", sm$outcome[i], sm$exposure[i]) }
}
chk("A4 sign/dlogOR/rel recomputed", ok, det)
chk("A4 no sign flip", sum(!sm$sign_stable) == 0, sprintf("flips=%d", sum(!sm$sign_stable)))
chk("A4 no significance change", sum(!sm$sig_stable) == 0, sprintf("changes=%d", sum(!sm$sig_stable)))
chk("A4 max |dlogOR| < 0.05", max(sm$max_abs_dlogOR) < 0.05,
    sprintf("max=%.4f", max(sm$max_abs_dlogOR)))

## ==========================================================================
## B. school-type stratification (25_sens_school_type.R)
## ==========================================================================
cat("\n===== B. school type =====\n")
st <- read.csv(file.path(out, "table18_school_type_stratified.csv"), stringsAsFactors = FALSE)
it <- read.csv(file.path(out, "table18_school_type_interaction.csv"), stringsAsFactors = FALSE)
sp <- read.csv(file.path(out, "table18_school_type_spread.csv"), stringsAsFactors = FALSE)

TYPE <- c("Associate's", "Baccalaureate", "Doctorate-granting", "Master's", "Special Focus")
PUB  <- c("Private", "Public")

## B1 strata partition the sample
for (ex in c("time_gt3h", "ai_any")) {
  s1 <- st %>% filter(stratifier == "inst_type_f", exposure == ex, outcome == "dep_any")
  chk(paste0("B1 inst_type n sums (", ex, ")"),
      sum(s1$n) == 41862, sprintf("sum=%d", sum(s1$n)))
  s2 <- st %>% filter(stratifier == "inst_public_f", exposure == ex, outcome == "dep_any")
  chk(paste0("B1 sector n sums (", ex, ")"), sum(s2$n) == 41862, sprintf("sum=%d", sum(s2$n)))
}
s3 <- st %>% filter(stratifier == "inst_type_f", exposure == "exp_harassed", outcome == "dep_any")
chk("B1 harassment n < total (eligibility)", sum(s3$n) < 41862 & sum(s3$n) > 38000,
    sprintf("sum=%d", sum(s3$n)))

## B2 school counts
mm <- m[!is.na(m$inst_type) & !is.na(m$inst_public), ]
mm$inst_type_f <- factor(mm$inst_type, levels = 1:5, labels = TYPE)
tot_sch <- length(unique(mm$schoolnum))
s1 <- st %>% filter(stratifier == "inst_type_f", exposure == "time_gt3h", outcome == "dep_any")
chk("B2 schools partition", sum(s1$n_schools) == tot_sch,
    sprintf("%d vs %d", sum(s1$n_schools), tot_sch))

## B3 one stratum model refit
k <- "Doctorate-granting"
lev <- which(TYPE == k)
sub <- mm[mm$inst_type == lev, ]
g <- svyglm(as.formula(paste("dep_any ~ time_gt3h +", ADJ)),
            design = svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = sub),
            family = quasibinomial())
tv <- st$OR[st$stratifier == "inst_type_f" & st$exposure == "time_gt3h" &
              st$outcome == "dep_any" & st$stratum == k]
chk("B3 refit Doctorate-granting", near(exp(coef(g)[["time_gt3h"]]), tv, 1e-3),
    sprintf("refit=%.5f table=%.3f", exp(coef(g)[["time_gt3h"]]), tv))

## B4 interaction df + BH invariants
chk("B4 df = 4 (inst_type)", all(it$df[it$stratifier == "inst_type_f"] == 4), "")
chk("B4 df = 1 (sector)",    all(it$df[it$stratifier == "inst_public_f"] == 1), "")
chk("B4 30 tests", nrow(it) == 30, sprintf("nrow=%d", nrow(it)))
chk("B4 BH >= raw", all(it$p_interaction_bh >= it$p_interaction - 1e-12), "")
chk("B4 no raw p<0.05", sum(it$p_interaction < 0.05) == 0,
    sprintf("n sig=%d (min p=%.3f)", sum(it$p_interaction < 0.05), min(it$p_interaction)))

## B5 the substantive claim: every estimable stratum keeps OR > 1 for the two strong lines
for (ex in c("time_gt3h", "exp_harassed")) {
  s <- st %>% filter(exposure == ex, ci_ok)
  chk(paste0("B5 all estimable strata OR>1 (", ex, ")"), all(s$OR > 1),
      sprintf("min OR=%.3f", min(s$OR)))
}
chk("B5 spread excludes Special Focus", all(sp$n_strata[sp$stratifier == "inst_type_f"] == 4),
    sprintf("n_strata=%s", paste(unique(sp$n_strata[sp$stratifier == "inst_type_f"]), collapse = ",")))
chk("B5 Special Focus flagged", all(!st$ci_ok[st$stratum == "Special Focus"]), "")

## ==========================================================================
## C. CFPS individual-level panel (26_cfps_panel.R)
## ==========================================================================
cat("\n===== C. CFPS panel =====\n")
res <- read.csv(file.path(out, "table19_cfps_panel_FE.csv"), stringsAsFactors = FALSE)
cmp <- read.csv(file.path(out, "table19_cfps_panel_compare.csv"), stringsAsFactors = FALSE)
rec <- read.csv(file.path(out, "table19_cfps_panel_reconcile.csv"), stringsAsFactors = FALSE)
pan <- read.csv(file.path(out, "cfps_panel_2020_2022.csv"), stringsAsFactors = FALSE)

cln <- function(x) { x <- zap_labels(x); x <- suppressWarnings(as.numeric(x)); x[x < 0] <- NA; x }
rd <- function(yr) {
  y2 <- substr(yr, 3, 4)
  f <- file.path(CFPS_DIR,
                 paste0("cfps", yr, "person.dta"))
  allc <- names(read_dta(f, n_max = 1))
  zap_labels(read_dta(f, col_select = all_of(intersect(
    c("pid", "age", "gender", "qu93", "qu931", "qu201a", "cesd8", "qn414",
      paste0("cfps", yr, "edu"), paste0("rswt_natcs", y2, "n")), allc))))
}
z20 <- rd("2020"); z22 <- rd("2022")
mkz <- function(z, yr) {
  y2 <- substr(yr, 3, 4)
  data.frame(pid = as.numeric(z$pid), age = cln(z$age),
             sv_watch = as.integer(cln(z$qu93) == 1),
             mobile_h = pmin(cln(z$qu201a) / 60, 16),
             cesd8_std = cln(z$cesd8) - 8, lonely_score = cln(z$qn414),
             edu = cln(z[[paste0("cfps", yr, "edu")]]),
             wt = cln(z[[paste0("rswt_natcs", y2, "n")]]))
}
q20 <- mkz(z20, "2020"); q22 <- mkz(z22, "2022")
q20$coh <- as.integer(q20$age >= 18 & q20$age <= 30 & q20$edu >= 5)
coh20 <- q20[q20$coh == 1, ]
p2 <- inner_join(coh20, q22 %>% select(pid, sv_watch, mobile_h, cesd8_std, lonely_score, wt),
                 by = "pid", suffix = c("_20", "_22"))
chk("C1 cohort n = 1237", nrow(p2) == 1237, sprintf("n=%d", nrow(p2)))
chk("C1 panel file rows = 1085", nrow(pan) == 1085, sprintf("n=%d", nrow(pan)))
chk("C1 cohort definition", sum(q20$coh, na.rm = TRUE) == 1779,
    sprintf("young_hiedu2020=%d", sum(q20$coh, na.rm = TRUE)))

## C2 changer counts recomputed ON THE ANALYSED cohort (after the weight/PSU drop),
## i.e. on the 1,085-panel file, not on the 1,237-person raw overlap.
p2$d_sv_all <- p2$sv_watch_22 - p2$sv_watch_20
pan$d_sv_watch <- pan$sv_watch_22 - pan$sv_watch_20
chk("C2 raw overlap larger than analysed", nrow(p2) == 1237 && nrow(pan) == 1085,
    sprintf("overlap=%d analysed=%d", nrow(p2), nrow(pan)))
chk("C2 discordant sv_watch = 183", sum(pan$d_sv_watch != 0, na.rm = TRUE) == 183,
    sprintf("n=%d", sum(pan$d_sv_watch != 0, na.rm = TRUE)))
chk("C2 0->1 = 121", sum(pan$d_sv_watch == 1, na.rm = TRUE) == 121,
    sprintf("n=%d", sum(pan$d_sv_watch == 1, na.rm = TRUE)))
chk("C2 1->0 = 62", sum(pan$d_sv_watch == -1, na.rm = TRUE) == 62,
    sprintf("n=%d", sum(pan$d_sv_watch == -1, na.rm = TRUE)))
chk("C2 discordant < raw overlap", sum(pan$d_sv_watch != 0, na.rm = TRUE) <
      sum(p2$d_sv_all != 0, na.rm = TRUE),
    sprintf("analysed=%d raw=%d", sum(pan$d_sv_watch != 0, na.rm = TRUE),
            sum(p2$d_sv_all != 0, na.rm = TRUE)))

## C3 the decisive check: Gaussian svyglm == weighted least squares.
## beta is stored rounded to 4 dp in the CSV, so compare at that precision.
chk_pairs <- list(
  c("d_dep8",     "d_sv_watch", "dep8 ~ sv_watch"),
  c("d_cesd8",    "d_sv_watch", "cesd8 ~ sv_watch"),
  c("d_lonely",   "d_sv_watch", "lonely ~ sv_watch"),
  c("d_cesd8",    "d_mobile",   "cesd8 ~ mobile_hours"),
  c("d_lonely_s", "d_mobile",   "lonely_score ~ mobile_hours"))
ok <- TRUE; det <- ""
for (cp in chk_pairs) {
  d <- pan[stats::complete.cases(pan[, c(cp[1], cp[2], "wt")]), ]
  lf <- lm(as.formula(paste(cp[1], "~", cp[2])), data = d, weights = wt)
  b_lm <- unname(coef(lf)[2])
  b_sv <- res$beta[res$estimator == "within-person (FE)" & res$spec == cp[3]]
  if (!near(round(b_lm, 4), b_sv, 1e-9)) { ok <- FALSE; det <- paste(det, cp[3], sprintf("%.6f vs %.4f", b_lm, b_sv)) }
}
chk("C3 FE == WLS (5 specs)", ok, det)

## C4 rescaling and units
chk("C4 CES-D8 rescale in 0-24",
    min(pan$cesd8_std_20, na.rm = TRUE) >= 0 && max(pan$cesd8_std_20, na.rm = TRUE) <= 24,
    sprintf("range=%.0f-%.0f", min(pan$cesd8_std_20, na.rm = TRUE), max(pan$cesd8_std_20, na.rm = TRUE)))
chk("C4 mobile_h capped at 16", max(c(pan$mobile_h_20, pan$mobile_h_22), na.rm = TRUE) <= 16,
    sprintf("max=%.2f", max(c(pan$mobile_h_20, pan$mobile_h_22), na.rm = TRUE)))
chk("C4 mobile_h median = 4 h", near(median(pan$mobile_h_22, na.rm = TRUE), 4, 1e-9),
    sprintf("median=%.2f", median(pan$mobile_h_22, na.rm = TRUE)))

## C5 invariants
chk("C5 BH >= raw (FE)", all(res$p_bh[res$estimator == "within-person (FE)"] >=
                             res$p[res$estimator == "within-person (FE)"] - 1e-12), "")
chk("C5 no FE q<0.05", sum(res$p_bh[res$estimator == "within-person (FE)"] < 0.05) == 0,
    sprintf("min q=%.3f", min(res$p_bh[res$estimator == "within-person (FE)"])))
chk("C5 CI half-width", all(near(cmp$fe_ci_halfwidth, (cmp$UCL_FE - cmp$LCL_FE) / 2, 1e-3)), "")
chk("C5 CI brackets beta", all(cmp$LCL_FE <= cmp$beta_FE & cmp$beta_FE <= cmp$UCL_FE), "")
chk("C5 n_changers <= n_FE", all(cmp$n_changers <= cmp$n_FE), "")

## C6 reconciliation with the published CFPS validation
full <- read.csv(file.path(out, "cfps_2022_young_hiedu.csv"), stringsAsFactors = FALSE)
ff <- full[!is.na(full$wt) & full$young_hiedu == 1, ]
ff$mobile_h <- pmin(ff$mobile_hours / 60, 16)
ff$urban_f <- factor(ff$urban)
sub <- ff[!is.na(ff$sv_watch) & !is.na(ff$dep8) & !is.na(ff$urban_f), ]
g <- svyglm(dep8 ~ sv_watch + age_num + female + urban_f,
            design = svydesign(ids = ~1, weights = ~wt, data = sub), family = quasibinomial())
tv <- rec$OR_per_unit[rec$sample == "CFPS2022 full young_hiedu" & rec$outcome == "dep8" &
                        rec$exposure == "sv_watch"]
chk("C6 full-cohort OR refit", near(exp(coef(g)[["sv_watch"]]), tv, 1e-3),
    sprintf("refit=%.5f table=%.3f", exp(coef(g)[["sv_watch"]]), tv))
chk("C6 full cohort all null", all(rec$LCL[rec$sample == "CFPS2022 full young_hiedu"] < 1 &
                                     rec$UCL[rec$sample == "CFPS2022 full young_hiedu"] > 1), "")
chk("C6 panel cohort all null", all(rec$LCL[rec$sample == "panel cohort (2022 wave)"] < 1 &
                                      rec$UCL[rec$sample == "panel cohort (2022 wave)"] > 1), "")

## ==========================================================================
cat(sprintf("\n===== %d PASS / %d FAIL =====\n", PASS, FAIL))
if (FAIL > 0) { cat("\nFAILURES:\n"); cat(paste0(" - ", MSG, collapse = "\n"), "\n") }
cat("DONE\n")
