source("00_config.R")
## 26_cfps_panel.R
## Individual-level two-wave panel, CFPS 2020 -> 2022 (external-validation cohort).
##
## Purpose: every HMS estimate is cross-sectional, so time-invariant person-level
## confounding (stable traits, family background, personality) is unaddressed by design.
## CFPS is the only dataset in this project that follows the SAME person across two
## waves, so it is the only place where a within-person (fixed-effects) estimator can be
## computed. This script estimates it and puts it side by side with the cross-sectional
## estimator computed on the identical people -- the difference is the part of the
## cross-sectional association attributable to stable between-person differences.
##
## Cohort: 18-30 years old with tertiary education (cfps2020edu >= 5) in 2020 AND
## present in 2022. Defining the cohort on the 2020 wave (rather than requiring the
## criteria in both waves) avoids selecting on survival into the 2022 age window,
## which would drop everyone who aged past 30.
##
## Estimator: first-difference (within-person) linear model on the change score,
## (y_2022 - y_2020) ~ (x_2022 - x_2020)
## weighted by the 2022 national cross-sectional weight, PSUs = CFPS psu.
## Time-invariant confounders difference out exactly. Two hard caveats:
## (i) with two waves, the estimator cannot be separated from a common period
## effect (everyone ages two years between waves);
## (ii) identification comes only from CHANGERS, so power is far below the
## cross-sectional comparison -- a null here is inconclusive, not evidence of
## no effect. Changer counts are reported alongside every estimate.
##
## Secondary: ANCOVA-style lagged model y_2022 ~ x_2020 + y_2020 (adjusts baseline
## outcome; still not immune to time-varying confounding).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(dplyr); library(haven) })

dir  <- ANALYSIS_DIR
out  <- file.path(dir, "outputs")
root <- CFPS_DIR

cln <- function(x) { x <- zap_labels(x); x <- suppressWarnings(as.numeric(x)); x[x < 0] <- NA; x }

read_wave <- function(yr) {
 ## NB: CFPS suffixes use the TWO-DIGIT year (provcd20, rswt_natcs20n), while the
 ## education variable uses the FOUR-DIGIT year (cfps2020edu). Getting this wrong
 ## makes the variable silently disappear from col_select.
  y2 <- substr(yr, 3, 4)
  f <- file.path(root, paste0("cfps", yr, "person.dta"))
  allc <- names(read_dta(f, n_max = 1))
  want <- c("pid", "psu", paste0("provcd", y2), "age", "gender",
            "qu93", "qu931", "qu201", "qu201a", "qu11",
            "cesd8", "qn414", paste0("cfps", yr, "edu"),
            "qc1", "qc3", paste0("urban", y2), paste0("rswt_natcs", y2, "n"))
  sel <- intersect(want, allc)
  stopifnot(all(c("pid", paste0("rswt_natcs", y2, "n")) %in% sel))
  d <- zap_labels(read_dta(f, col_select = all_of(sel)))
  d$.yr <- yr
  d
}

w20 <- read_wave("2020"); w22 <- read_wave("2022")
cat("CFPS 2020 rows:", nrow(w20), "| 2022 rows:", nrow(w22), "\n")

mk <- function(d, yr) {
  y2 <- substr(yr, 3, 4)
  edu <- intersect(paste0("cfps", yr, "edu"), names(d))
  tibble(
    pid = as.numeric(d$pid),
    psu  = suppressWarnings(as.numeric(d$psu)),
    prov = suppressWarnings(as.numeric(d[[paste0("provcd", y2)]])),
    age  = cln(d$age),
    female = as.integer(cln(d$gender) == 1),
    sv_watch = as.integer(cln(d$qu93) == 1),
    sv_daily = as.integer(cln(d$qu931) == 1),
 ## qu201a is recorded in MINUTES per day (median 120, max 1440). Convert to hours
 ## and cap at 16 h/day: 24 h/day is impossible and would dominate a change score.
    mobile_h = pmin(cln(d$qu201a) / 60, 16),
    wechat = as.integer(cln(d$qu11) == 1),
    cesd8_std = cln(d$cesd8) - 8,                 # CFPS CES-D8 is 8-32; rescale to 0-24
    lonely_score = cln(d$qn414),
    edu_max = if (length(edu) == 1) cln(d[[edu]]) else NA_real_,
    urban = cln(d[[paste0("urban", y2)]]),
    wt = cln(d[[paste0("rswt_natcs", y2, "n")]]))
}

a <- mk(w20, "2020"); b <- mk(w22, "2022")
a <- a %>% mutate(dep8 = as.integer(cesd8_std >= 8), lonely = as.integer(lonely_score >= 3),
                  young_hiedu = as.integer(age >= 18 & age <= 30 & edu_max >= 5))
b <- b %>% mutate(dep8 = as.integer(cesd8_std >= 8), lonely = as.integer(lonely_score >= 3),
                  young_hiedu = as.integer(age >= 18 & age <= 30 & edu_max >= 5))
cat("young_hiedu: 2020 =", sum(a$young_hiedu, na.rm = TRUE),
    "| 2022 =", sum(b$young_hiedu, na.rm = TRUE), "\n")

## ---- cohort: 2020-defined, present in both waves -----------------------------------
p <- inner_join(
  a %>% filter(young_hiedu == 1) %>% select(-young_hiedu),
  b %>% select(pid, psu, prov, age, sv_watch, sv_daily, mobile_h, wechat,
               cesd8_std, dep8, lonely_score, lonely, urban, wt),
  by = "pid", suffix = c("_20", "_22"))
cat("\npanel cohort n =", nrow(p), "\n")

## weights / PSUs: keep 2022 where available, fall back to 2020; drop if neither.
## (CFPS assigns no national cross-sectional weight to some respondents -- about
## 2,900 per wave -- so this is a real loss, reported explicitly.)
p <- p %>% mutate(
  psu = coalesce(psu_22, psu_20),
  prov = coalesce(prov_22, prov_20),
  wt = coalesce(wt_22, wt_20))
n_pre <- nrow(p)
p <- p %>% filter(!is.na(wt), !is.na(psu), wt > 0)
cat("dropped for missing weight/PSU:", n_pre - nrow(p), "-> analysis cohort n =", nrow(p), "\n")

p <- p %>% mutate(
  d_sv_watch = sv_watch_22 - sv_watch_20,
  d_sv_daily = sv_daily_22 - sv_daily_20,
  d_mobile   = mobile_h_22 - mobile_h_20,
  d_cesd8    = cesd8_std_22 - cesd8_std_20,
  d_dep8     = dep8_22 - dep8_20,
  d_lonely   = lonely_22 - lonely_20,
  d_lonely_s = lonely_score_22 - lonely_score_20)

cat("\n--- changer counts (the only source of identification) ---\n")
cat("sv_watch both valid :", sum(!is.na(p$sv_watch_20) & !is.na(p$sv_watch_22)),
    "| discordant:", sum(!is.na(p$d_sv_watch) & p$d_sv_watch != 0),
    " (0->1:", sum(p$d_sv_watch == 1, na.rm = TRUE), ", 1->0:", sum(p$d_sv_watch == -1, na.rm = TRUE), ")\n")
cat("mobile_h both       :", sum(!is.na(p$mobile_h_20) & !is.na(p$mobile_h_22)),
    "| changed:", sum(!is.na(p$d_mobile) & p$d_mobile != 0),
    "| |delta| > 8 h:", sum(abs(p$d_mobile) > 8, na.rm = TRUE), "\n")
cat("dep8 both valid     :", sum(!is.na(p$dep8_20) & !is.na(p$dep8_22)),
    "| discordant:", sum(!is.na(p$d_dep8) & p$d_dep8 != 0), "\n")
cat("cesd8 both valid    :", sum(!is.na(p$cesd8_std_20) & !is.na(p$cesd8_std_22)), "\n")
cat("lonely both valid   :", sum(!is.na(p$lonely_20) & !is.na(p$lonely_22)),
    "| discordant:", sum(!is.na(p$d_lonely) & p$d_lonely != 0), "\n")
cat("\nmobile hours/day: 2020 median", median(p$mobile_h_20, na.rm = TRUE),
    "| 2022 median", median(p$mobile_h_22, na.rm = TRUE),
    "| range", paste(range(c(p$mobile_h_20, p$mobile_h_22), na.rm = TRUE), collapse = "-"), "\n")

## ---------- models -------------------------------------------------------------------
des_of <- function(dat, w = "wt") svydesign(ids = ~psu, weights = as.formula(paste0("~", w)),
                                            nest = TRUE, data = dat)
fit <- function(dat, fo) {
  g <- svyglm(as.formula(fo), design = des_of(dat), family = gaussian())
  nm <- setdiff(names(coef(g)), "(Intercept)")
  ci <- confint(g)
  b <- coef(g)[[nm[1]]]; se <- sqrt(diag(vcov(g)))[[nm[1]]]
  data.frame(term = nm[1], beta = b, se = se,
             LCL = ci[nm[1], 1], UCL = ci[nm[1], 2],
             p = 2 * pt(-abs(b / se), df = degf(des_of(dat))), n = nrow(dat))
}

## (A) within-person (first difference) and (B) cross-sectional 2022, same people.
## The cross-sectional model carries the SAME adjustment as the published CFPS
## external validation (09_replication_cfps.R: age + female + urban) so that the
## two are directly comparable and any discrepancy is attributable to the sample
## (panel cohort vs full 2022 cohort), not to the specification.
SPECS <- list(
  list(id = "dep8 ~ sv_watch",
       fe = "d_dep8 ~ d_sv_watch", cs = "dep8_22 ~ sv_watch_22",
       ylab = "any depression (CES-D8>=8)", xlab = "short-video use"),
  list(id = "cesd8 ~ sv_watch",
       fe = "d_cesd8 ~ d_sv_watch", cs = "cesd8_std_22 ~ sv_watch_22",
       ylab = "CES-D8 score (0-24)", xlab = "short-video use"),
  list(id = "lonely ~ sv_watch",
       fe = "d_lonely ~ d_sv_watch", cs = "lonely_22 ~ sv_watch_22",
       ylab = "loneliness", xlab = "short-video use"),
  list(id = "cesd8 ~ mobile_hours",
       fe = "d_cesd8 ~ d_mobile", cs = "cesd8_std_22 ~ mobile_h_22",
       ylab = "CES-D8 score (0-24)", xlab = "mobile hours per day"),
  list(id = "lonely_score ~ mobile_hours",
       fe = "d_lonely_s ~ d_mobile", cs = "lonely_score_22 ~ mobile_h_22",
       ylab = "loneliness score", xlab = "mobile hours per day")
)
CS_ADJ <- "+ female + age_22 + urban_22"

cc <- function(dat, fo) {
  v <- all.vars(as.formula(fo))
  stopifnot(all(v %in% names(dat)))
  dat[stats::complete.cases(dat[, v, drop = FALSE]), , drop = FALSE]
}
## number of people who actually changed the exposure = the identifying information
n_changers <- function(dat, dvar) sum(!is.na(dat[[dvar]]) & dat[[dvar]] != 0, na.rm = TRUE)

rows <- list()
for (s in SPECS) {
  dfe <- cc(p, s$fe); dcs <- cc(p, paste(s$cs, CS_ADJ))
  rfe <- fit(dfe, s$fe); rcs <- fit(dcs, paste(s$cs, CS_ADJ))
  nch <- n_changers(dfe, all.vars(as.formula(s$fe))[2])
  rows[[length(rows) + 1]] <- data.frame(
    spec = s$id, outcome_lab = s$ylab, exposure_lab = s$xlab,
    estimator = "within-person (FE)", beta = rfe$beta, se = rfe$se,
    LCL = rfe$LCL, UCL = rfe$UCL, p = rfe$p, n = rfe$n, n_changers = nch)
  rows[[length(rows) + 1]] <- data.frame(
    spec = s$id, outcome_lab = s$ylab, exposure_lab = s$xlab,
    estimator = "cross-sectional 2022", beta = rcs$beta, se = rcs$se,
    LCL = rcs$LCL, UCL = rcs$UCL, p = rcs$p, n = rcs$n, n_changers = NA)
  cat(sprintf("%-30s FE  b=%8.4f [%7.4f,%7.4f] p=%.4f n=%d changers=%d\n",
              s$id, rfe$beta, rfe$LCL, rfe$UCL, rfe$p, rfe$n, nch))
  cat(sprintf("%-30s CS  b=%8.4f [%7.4f,%7.4f] p=%.4f n=%d\n",
              "", rcs$beta, rcs$LCL, rcs$UCL, rcs$p, rcs$n))
  flush.console()
}
res <- bind_rows(rows)

## (C) lagged / ANCOVA: outcome 2022 ~ exposure 2020 + outcome 2020
lag_rows <- list()
for (s in SPECS) {
  y22 <- all.vars(as.formula(s$cs))[1]; x20 <- sub("_22", "_20", all.vars(as.formula(s$cs))[2])
  y20 <- sub("_22", "_20", y22)
  if (!(y20 %in% names(p)) || !(x20 %in% names(p))) next
  fo <- paste(y22, "~", x20, "+", y20)
  dl <- cc(p, fo)
  if (nrow(dl) < 100) next
  rl <- fit(dl, fo)
  lag_rows[[length(lag_rows) + 1]] <- data.frame(
    spec = s$id, outcome_lab = s$ylab, exposure_lab = s$xlab,
    estimator = "lagged (2020 exposure, adj. 2020 outcome)",
    beta = rl$beta, se = rl$se, LCL = rl$LCL, UCL = rl$UCL, p = rl$p, n = rl$n)
}
lag <- bind_rows(lag_rows)
res <- bind_rows(res, lag)

## BH within each estimator family (5 tests per family: 3 outcomes x 2 exposures + 2)
res <- res %>% group_by(estimator) %>% mutate(p_bh = p.adjust(p, method = "BH")) %>% ungroup()
res <- res %>% mutate(beta = round(beta, 4), se = round(se, 4),
                      LCL = round(LCL, 4), UCL = round(UCL, 4),
                      p = signif(p, 4), p_bh = signif(p_bh, 4))

## summary: does the cross-sectional association survive the within-person test?
## An attenuation percentage is meaningless when the two signs disagree, so we report
## sign agreement and whether the FE interval EXCLUDES the cross-sectional point
## estimate (the informative case) instead of a ratio.
fe <- res %>% filter(estimator == "within-person (FE)") %>%
  select(spec, beta_FE = beta, LCL_FE = LCL, UCL_FE = UCL, p_FE = p, q_FE = p_bh,
         n_FE = n, n_changers)
cs <- res %>% filter(estimator == "cross-sectional 2022") %>%
  select(spec, beta_CS = beta, p_CS = p, q_CS = p_bh, n_CS = n)
cmp <- inner_join(fe, cs, by = "spec") %>%
  mutate(sign_agrees = sign(beta_FE) == sign(beta_CS),
         CS_outside_FE_CI = beta_CS < LCL_FE | beta_CS > UCL_FE,
 ## half-width of the FE interval = the smallest effect this design can still
 ## distinguish from zero at 95% confidence. Reported because a null FE result
 ## is only interpretable against it.
         fe_ci_halfwidth = round((UCL_FE - LCL_FE) / 2, 4))
cat("\n=== within-person vs cross-sectional (same people) ===\n")
print(as.data.frame(cmp), row.names = FALSE)

## ---- reconciliation with the published CFPS validation (09 / Figure S3b) ----------
## The published validation used the FULL 2022 young_hiedu cohort and LOGISTIC models
## (binary outcomes). Re-run it here on the full cohort to confirm the panel cohort is
## not producing a contradictory cross-sectional signal.
c22 <- zap_labels(readRDS(file.path(out, "cfps_2022_analysis.rds")))
full <- c22[c22$young_hiedu == 1 & !is.na(c22$wt), ]
full$mobile_h <- pmin(full$mobile_hours / 60, 16)
full$urban_f <- factor(full$urban)
rec <- list()
for (e in c("sv_watch", "mobile_h")) {
  for (o in c("dep8", "lonely")) {
    sub <- full[!is.na(full[[e]]) & !is.na(full[[o]]) & !is.na(full$urban_f), ]
    d2 <- svydesign(ids = ~1, weights = ~wt, data = sub)
    g <- svyglm(as.formula(paste(o, "~", e, "+ age_num + female + urban_f")),
                design = d2, family = quasibinomial())
    ci <- confint(g); b <- coef(g)[[e]]
    rec[[length(rec) + 1]] <- data.frame(sample = "CFPS2022 full young_hiedu", outcome = o,
      exposure = e, OR_per_unit = round(exp(b), 3), LCL = round(exp(ci[e, 1]), 3),
      UCL = round(exp(ci[e, 2]), 3), n = nrow(sub))
  }
}
## same logistic specification, restricted to the panel cohort, for a like-for-like check
for (e in c("sv_watch_22", "mobile_h_22")) {
  for (o in c("dep8_22", "lonely_22")) {
    sub <- p[stats::complete.cases(p[, c(e, o, "female", "age_22", "urban_22")]), ]
    sub$urban_f <- factor(sub$urban_22)
    d2 <- svydesign(ids = ~psu, weights = ~wt, data = sub)
    g <- svyglm(as.formula(paste(o, "~", e, "+ female + age_22 + urban_f")),
                design = d2, family = quasibinomial())
    ci <- confint(g); b <- coef(g)[[e]]
    rec[[length(rec) + 1]] <- data.frame(sample = "panel cohort (2022 wave)", outcome = o,
      exposure = e, OR_per_unit = round(exp(b), 3), LCL = round(exp(ci[e, 1]), 3),
      UCL = round(exp(ci[e, 2]), 3), n = nrow(sub))
  }
}
rec <- bind_rows(rec)
cat("\n=== reconciliation with 09/Figure S3b (logistic, adjusted) ===\n")
print(as.data.frame(rec), row.names = FALSE)

write.csv(res, file.path(out, "table19_cfps_panel_FE.csv"), row.names = FALSE)
write.csv(cmp, file.path(out, "table19_cfps_panel_compare.csv"), row.names = FALSE)
write.csv(rec, file.path(out, "table19_cfps_panel_reconcile.csv"), row.names = FALSE)
write.csv(p,   file.path(out, "cfps_panel_2020_2022.csv"), row.names = FALSE)
cat("\nDONE\n")
