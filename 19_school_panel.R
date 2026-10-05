source("00_config.R")
## 19_school_panel.R
## School-level two-period panel with school fixed effects (2023-24 -> 2024-25).
##
## Why: the negative-control outcome (institutional graduation rate) showed that the
## screen-time association carries school-level residual confounding. Individual-level
## adjustment cannot remove confounding from time-invariant school characteristics.
## A within-school (first-difference) design absorbs ALL time-invariant school confounding.
##
## Cross-wave comparability (verified in 19a/19b, not assumed):
## internet_1 scale differs in the number of bins (7 vs 8) but the COLLAPSED categories
## are identical: <1h / 1-2h / 2-3h / >3h. So time5cat and time_gt3h are comparable.
## internet_2_2 (harassment) wording CHANGED in 2024-25: "and/or been called names" was
## added, broadening the item. A broader item should RAISE prevalence, yet observed
## prevalence fell (12.81% -> 9.62%), so the true decline is at least as large as
## observed; the resulting measurement error attenuates the harassment coefficient
## towards zero, i.e. the harassment estimate is conservative.
##
## Estimators reported side by side (this is itself a benchmark):
## pooled OLS on 264 school-waves -- ignores school confounding
## first difference (= two-way FE with T=2) -- absorbs time-invariant school confounding

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

p <- OUT
a <- readRDS(file.path(p, "hms_2324_analysis.rds"))
b <- readRDS(file.path(p, "hms_2425_analysis.rds"))

outs <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
exps <- c("time_gt3h", "exp_harassed")

## ---- school-level weighted prevalence for one wave ----
agg <- function(d, wave) {
  d <- d[d$module == 1, ]
  for (v in c(exps, outs, "nrweight", "schoolnum", "female", "age_num",
              "fin_stress", "firstgen")) d <- d[!is.na(d[[v]]), ]
  sp <- split(d, d$schoolnum)
  do.call(rbind, lapply(names(sp), function(s) {
    x <- sp[[s]]
    w <- x$nrweight
    wp <- function(v) sum(w * x[[v]]) / sum(w)
    ig <- suppressWarnings(mean(x$inst_gradrate, na.rm = TRUE))
    if (is.nan(ig)) ig <- NA_real_
    data.frame(schoolnum = s, wave = wave,
               n_eff = sum(!is.na(w)),
               wsum = sum(w),
               inst_gradrate = ig,
               p_time_gt3h  = wp("time_gt3h"),
               p_exp_harassed = wp("exp_harassed"),
               dep_any = wp("dep_any"), dep_maj = wp("dep_maj"),
               anx_any = wp("anx_any"), sui_idea = wp("sui_idea"),
               lonely = wp("lonely"),
               female = wp("female"), age_num = sum(w * x$age_num) / sum(w),
               fin_stress = wp("fin_stress"), firstgen = wp("firstgen"),
               stringsAsFactors = FALSE)
  }))
}

A <- agg(a, "2324"); B <- agg(b, "2425")
cat("school-level cells: 2324 =", nrow(A), " 2425 =", nrow(B), "\n")

sch <- intersect(A$schoolnum, B$schoolnum)
cat("overlapping schools:", length(sch), "\n")
A <- A[A$schoolnum %in% sch, ]; B <- B[B$schoolnum %in% sch, ]

## ---- wide (one row per school, .1 = 2023-24, .2 = 2024-25) ----
W <- merge(A, B, by = "schoolnum", suffixes = c(".1", ".2"))
stopifnot(nrow(W) == length(sch))

## harmonic-mean effective sample size as the school-level weight
W$w <- 2 / (1 / W$n_eff.1 + 1 / W$n_eff.2)

## outcomes are stored bare (dep_any.1 / dep_any.2); exposures carry the p_ prefix
for (o in outs) W[[paste0("d_", o)]] <- W[[paste0(o, ".2")]] - W[[paste0(o, ".1")]]
for (e in exps) W[[paste0("d_", e)]] <- W[[paste0("p_", e, ".2")]] - W[[paste0("p_", e, ".1")]]
W$d_gradrate <- W$inst_gradrate.2 - W$inst_gradrate.1

cat("\n===== school-level change between waves (weighted by harmonic-mean n) =====\n")
for (v in c(paste0("p_", exps), outs)) {
  d1 <- W[[paste0(v, ".1")]]; d2 <- W[[paste0(v, ".2")]]
  cat(sprintf("%-14s 2324 = %.4f  2425 = %.4f  mean delta = %+.4f  (SD %.4f)\n",
              v, weighted.mean(d1, W$w), weighted.mean(d2, W$w),
              weighted.mean(d2 - d1, W$w), sd(d2 - d1)))
}

## ---- estimators ----
rows <- list()
for (e in exps) {
  xe <- paste0("p_", e)
  for (o in outs) {
 ## pooled OLS across 264 school-wave observations, no school FE
    pl <- rbind(
      data.frame(y = W[[paste0(o, ".1")]], x = W[[paste0(xe, ".1")]],
                 w = W$n_eff.1, school = W$schoolnum),
      data.frame(y = W[[paste0(o, ".2")]], x = W[[paste0(xe, ".2")]],
                 w = W$n_eff.2, school = W$schoolnum))
    fp <- lm(y ~ x, data = pl, weights = w)
    cp <- coeftest(fp, vcov = vcovCL(fp, cluster = ~school))[2, ]

 ## first difference: dY ~ dX, intercept = common time trend
    fd <- lm(W[[paste0("d_", o)]] ~ W[[paste0("d_", e)]], weights = W$w)
    cf <- coeftest(fd, vcov = vcovHC(fd, type = "HC1"))[2, ]

 ## adjusted first difference: also difference out changing school composition
    fd2 <- lm(W[[paste0("d_", o)]] ~ W[[paste0("d_", e)]] +
                (W$female.2 - W$female.1) + (W$age_num.2 - W$age_num.1) +
                (W$fin_stress.2 - W$fin_stress.1) + (W$firstgen.2 - W$firstgen.1),
              weights = W$w)
    cf2 <- coeftest(fd2, vcov = vcovHC(fd2, type = "HC1"))[2, ]

    rows[[length(rows) + 1]] <- data.frame(
      exposure = e, outcome = o, n_schools = nrow(W),
      beta_pooled = cp[1], se_pooled = cp[2], p_pooled = cp[4],
      beta_FE = cf[1], se_FE = cf[2], p_FE = cf[4],
      beta_FE_adj = cf2[1], se_FE_adj = cf2[2], p_FE_adj = cf2[4],
      confounding_shift = cp[1] - cf[1],
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
res[, 4:ncol(res)] <- round(res[, 4:ncol(res)], 4)
write.csv(res, file.path(p, "table15_school_panel.csv"), row.names = FALSE)

## ---- negative control: institutional graduation rate ----
## inst_gradrate is a school characteristic; its change should be unrelated to the
## change in student digital behaviour once school FE are absorbed.
cat(sprintf("\ndelta inst_gradrate: mean %+.4f  SD %.4f  n valid %d\n",
            mean(W$d_gradrate, na.rm = TRUE), sd(W$d_gradrate, na.rm = TRUE),
            sum(!is.na(W$d_gradrate))))
if (sd(W$d_gradrate, na.rm = TRUE) > 0) {
  cat("\n===== negative control: delta graduation rate ~ delta exposure =====\n")
  for (e in exps) {
    f <- lm(d_gradrate ~ W[[paste0("d_", e)]], data = W, weights = W$w)
    c1 <- coeftest(f, vcov = vcovHC(f, type = "HC1"))[2, ]
    cat(sprintf("%-14s beta = %+.4f  SE %.4f  p = %.4f\n", e, c1[1], c1[2], c1[4]))
  }
} else {
  cat("\nNOTE: inst_gradrate does not vary over time; negative control not informative.\n")
}

cat("\n===== school panel results =====\n")
print(res[, c("exposure", "outcome", "beta_pooled", "p_pooled",
              "beta_FE", "p_FE", "beta_FE_adj", "p_FE_adj", "confounding_shift")],
      row.names = FALSE)
write.csv(W, file.path(p, "table15_school_panel_wide.csv"), row.names = FALSE)
cat("DONE\n")
