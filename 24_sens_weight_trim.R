source("00_config.R")
## 24_sens_weight_trim.R
## Sensitivity to extreme non-response weights (nrweight).
##
## Question addressed: HMS nrweight spans 0.21-76.7 (max/median = 67.5). If a handful of
## very large weights drive the three headline associations (screen >3 h/d, online
## harassment, any generative-AI use), the paper's conclusions are an artefact of the
## weighting scheme rather than a population signal.
##
## Design: identical model, identical covariates, identical clustering (schools = PSU);
## only the weight vector changes. Four variants:
## w_raw : nrweight as delivered
## w_p999 : winsorised at the 99.9th percentile of the module sample
## w_p995 : winsorised at the 99.5th percentile (pre-registered variant in 04_)
## w_p99 : winsorised at the 99th percentile
## Winsorisation thresholds are computed ONCE on the full module sample (n = 46,914),
## i.e. the sample the non-response weights were constructed for -- NOT on each
## outcome-specific complete-case subset (that would make the threshold itself vary
## with the model and confound the comparison).
##
## Metrics reported per variant: Kish effective sample size, design effect, CV of weights.
## Per model: OR, 95% CI, p; plus stability flags (sign flip, p crossing 0.05).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(dplyr); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

## CRITICAL: outcomes/covariates are still haven_labelled (vctrs_vctr). With dplyr
## loaded, `binomial$initialize` does `y < 0` on the response and vctrs refuses to
## cast labelled-vs-double -> "Execution halted" for sui_idea (the other four outcomes
## survive only because of GLORP in how svyglm builds the model frame). Zap once here.
d <- zap_labels(d)
stopifnot(all(vapply(d, function(x) !inherits(x, "haven_labelled"), logical(1))))

core_cov <- c("age_num", "female", "race_cat", "international", "undergrad",
              "grad_student", "fin_stress", "food_insec", "firstgen")
OC <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")

## --- main analytic sample (identical construction to 04_var_centered.R) -------------
m <- d[d$module == 1, ]
m <- m[!is.na(m$time5cat), ]
for (v in OC) m <- m[!is.na(m[[v]]), ]
for (v in core_cov) m <- m[!is.na(m[[v]]), ]
cat("main analytic sample n =", nrow(m), "\n")

## --- weight variants ---------------------------------------------------------------
## thresholds from the FULL module sample (n = 46,914), not the regression subset
mod <- d[d$module == 1 & !is.na(d$time5cat), ]
qs <- quantile(mod$nrweight, c(0.99, 0.995, 0.999), na.rm = TRUE)
cat("winsorisation thresholds (module sample n =", nrow(mod), "):\n"); print(round(qs, 3))

m$w_raw  <- m$nrweight
m$w_p99  <- pmin(m$nrweight, qs[["99%"]])
m$w_p995 <- pmin(m$nrweight, qs[["99.5%"]])
m$w_p999 <- pmin(m$nrweight, qs[["99.9%"]])

WV <- c("w_raw", "w_p999", "w_p995", "w_p99")

## --- weight diagnostics ------------------------------------------------------------
kish <- function(w) { w <- w / mean(w); sum(w)^2 / sum(w^2) }
diag_rows <- list()
for (w in WV) {
  v <- m[[w]]
  ess <- kish(v)
  diag_rows[[w]] <- data.frame(weight = w, n = nrow(m), sum_w = round(sum(v)),
                               mean_w = round(mean(v), 3), max_w = round(max(v), 3),
                               cv_w = round(sd(v) / mean(v), 3),
                               ESS = round(ess, 1), DEFF = round(nrow(m) / ess, 2),
                               pct_winsorised = round(100 * mean(m$nrweight > max(v)), 3))
}
diag <- do.call(rbind, diag_rows)
cat("\nweight diagnostics:\n"); print(diag, row.names = FALSE)

## --- models ------------------------------------------------------------------------
adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

EXPS <- list(list(v = "time_gt3h",    lab = "screen>3h/d"),
             list(v = "exp_harassed", lab = "online harassment"),
             list(v = "ai_any",       lab = "any genAI use"))

one <- function(o, ev, w) {
  sub <- m[!is.na(m[[ev]]), ]
  des <- svydesign(ids = ~schoolnum, weights = as.formula(paste0("~", w)),
                   nest = TRUE, data = sub)
  g <- svyglm(as.formula(paste(o, "~", ev, "+", adj)), design = des, family = quasibinomial())
  ci <- confint(g); b <- coef(g)[[ev]]; se <- sqrt(diag(vcov(g)))[[ev]]
  data.frame(outcome = o, exposure = ev, weight = w,
             n = nrow(sub),
             OR = exp(b), LCL = exp(ci[ev, 1]), UCL = exp(ci[ev, 2]),
             logOR = b, se_logOR = se, p = 2 * pnorm(-abs(b / se)))
}

res <- list()
for (e in EXPS) for (o in OC) for (w in WV) {
  r <- try(one(o, e$v, w), silent = TRUE)
  if (inherits(r, "try-error")) { cat("FAILED:", o, e$v, w, "\n"); next }
  r$exposure_lab <- e$lab
  res[[length(res) + 1]] <- r
}
all <- bind_rows(res)
all <- all %>% mutate(OR = round(OR, 3), LCL = round(LCL, 3), UCL = round(UCL, 3),
                      logOR = round(logOR, 4), p = signif(p, 4)) %>%
  select(outcome, exposure, exposure_lab, weight, n, OR, LCL, UCL, logOR, se_logOR, p)

## --- stability summary -------------------------------------------------------------
sm <- all %>% group_by(outcome, exposure, exposure_lab) %>% summarise(
  n_models = n(),
  OR_raw   = OR[weight == "w_raw"],
  OR_min   = min(OR), OR_max = max(OR),
  max_abs_dlogOR   = round(max(abs(logOR - logOR[weight == "w_raw"])), 4),
  max_rel_OR_shift = round(100 * max(abs(logOR - logOR[weight == "w_raw"])), 2),
  sign_stable  = all(OR > 1) || all(OR < 1),
  sig_raw      = p[weight == "w_raw"] < 0.05,
  sig_stable   = all((p < 0.05) == (p[weight == "w_raw"] < 0.05)),
  .groups = "drop")
cat("\nstability summary (15 exposure x outcome pairs):\n")
print(as.data.frame(sm), row.names = FALSE)
cat("\nmax |delta log OR| across all variants:", max(sm$max_abs_dlogOR), "\n")
cat("sign flips:", sum(!sm$sign_stable), "| significance changes:", sum(!sm$sig_stable), "\n")

write.csv(diag, file.path(out, "table17_weight_diagnostics.csv"), row.names = FALSE)
write.csv(all,  file.path(out, "table17_weight_trim_OR.csv"),   row.names = FALSE)
write.csv(sm,   file.path(out, "table17_weight_trim_summary.csv"), row.names = FALSE)
cat("DONE\n")
