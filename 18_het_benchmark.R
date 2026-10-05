source("00_config.R")
## 18_het_benchmark.R
## Heterogeneity layer of the methodological benchmark.
##
## Question: the conventional route to effect modification is a pre-specified product term in
## a weighted logistic model, one test per candidate moderator. The ML route is a causal
## forest with best-linear projection onto the same moderators. Do they agree?
##
## Part A: head-to-head on the exposure the forest was fitted to (time_gt3h), 7 moderators
## x 5 outcomes = 35 tests per method, identical sample, identical moderators.
## Part B: the conventional product-term route applied to the two other exposures
## (exp_harassed, ai_any) -- 70 further tests -- to see how it behaves on a weak
## association, where the forest was never run.
##
## Reported: raw and BH-adjusted significance counts, and the share of conventional p-values
## below 0.05 (expected 5% under a correctly specified null).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

vars <- c("age_num", "female", "race_cat", "international", "undergrad", "grad_student",
          "fin_stress", "food_insec", "firstgen")
mods <- c("fin_stress", "food_insec", "firstgen", "female", "undergrad", "age_num",
          "international")
outs <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")

## ---- identical sample to 07_cate.R so the comparison is exact ----
m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c(outs, "time_gt3h", "exp_harassed", "ai_any", "nrweight", "schoolnum"))
  m <- m[!is.na(m[[v]]), ]
for (v in vars) m <- m[!is.na(m[[v]]), ]
rownames(m) <- NULL
cat("heterogeneity benchmark n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = m)

## conventional product-term test for one (outcome, exposure, moderator)
trad <- function(o, a, mod) {
  rest <- setdiff(vars, mod)
  fo <- as.formula(paste(o, "~", a, "*", mod, "+", paste(rest, collapse = " + ")))
  g <- svyglm(fo, design = des, family = quasibinomial())
  cf <- summary(g)$coefficients
  inm <- grep(":", rownames(cf), value = TRUE)
  inm <- inm[grepl(a, inm)]
  if (length(inm) != 1) return(NULL)
  data.frame(coef = unname(cf[inm, "Estimate"]),
             se = unname(cf[inm, "Std. Error"]),
             p = unname(cf[inm, "Pr(>|t|)"]),
             stringsAsFactors = FALSE)
}

rows <- list()
for (a in c("time_gt3h", "exp_harassed", "ai_any")) {
  for (o in outs) {
    for (mod in mods) {
      r <- trad(o, a, mod)
      if (is.null(r)) next
      rows[[length(rows) + 1]] <- data.frame(
        exposure = a, outcome = o, moderator = mod,
        coef = r$coef, se = r$se, p_raw = r$p, stringsAsFactors = FALSE)
    }
  }
  cat("  conventional product terms done:", a, "\n"); flush.console()
}
trad_all <- do.call(rbind, rows)

## ---- merge with the causal-forest BLP (35 tests on time_gt3h) ----
blp <- read.csv(file.path(out, "table4_cate_blp.csv"), stringsAsFactors = FALSE)
blp <- blp[blp$moderator != "(Intercept)", ]
blp$exposure <- "time_gt3h"
names(blp)[names(blp) == "p"] <- "p_cf"
names(blp)[names(blp) == "coef"] <- "coef_cf"
names(blp)[names(blp) == "se"] <- "se_cf"

tt <- trad_all[trad_all$exposure == "time_gt3h", ]
cmp <- merge(tt, blp[, c("outcome", "moderator", "coef_cf", "se_cf", "p_cf")],
             by = c("outcome", "moderator"), all.x = TRUE)
stopifnot(nrow(cmp) == 35)

cmp$p_trad_bh <- p.adjust(cmp$p_raw, "BH")
cmp$p_cf_bh <- p.adjust(cmp$p_cf, "BH")
cmp$sig_trad_raw <- cmp$p_raw < 0.05
cmp$sig_trad_bh <- cmp$p_trad_bh < 0.05
cmp$sig_cf_raw <- cmp$p_cf < 0.05
cmp$sig_cf_bh <- cmp$p_cf_bh < 0.05
cmp$sign_agree <- sign(cmp$coef) == sign(cmp$coef_cf)
cmp[, 4:ncol(cmp)] <- round(cmp[, 4:ncol(cmp)], 5)
write.csv(cmp, file.path(out, "table14_het_benchmark.csv"), row.names = FALSE)

## BH within exposure. NOT split+unlist: split reorders groups alphabetically and
## silently misaligns with the row order (this produced p_bh < p_raw on a previous run).
trad_all$p_bh <- NA_real_
for (a in unique(trad_all$exposure)) {
  idx <- trad_all$exposure == a
  trad_all$p_bh[idx] <- p.adjust(trad_all$p_raw[idx], "BH")
}
## invariant: BH-adjusted p can never be smaller than the raw p
stopifnot(all(trad_all$p_bh >= trad_all$p_raw - 1e-12),
          all(cmp$p_trad_bh >= cmp$p_raw - 1e-12),
          all(cmp$p_cf_bh >= cmp$p_cf - 1e-12))
trad_all[, c("coef", "se", "p_raw", "p_bh")] <-
  round(trad_all[, c("coef", "se", "p_raw", "p_bh")], 5)
write.csv(trad_all, file.path(out, "table14b_trad_interactions.csv"), row.names = FALSE)

summ <- data.frame(
  method = c("conventional product term", "conventional product term",
             "causal forest BLP", "causal forest BLP"),
  correction = c("raw", "BH", "raw", "BH"),
  n_tests = c(35, 35, 35, 35),
  n_significant = c(sum(cmp$sig_trad_raw), sum(cmp$sig_trad_bh),
                    sum(cmp$sig_cf_raw), sum(cmp$sig_cf_bh)),
  stringsAsFactors = FALSE)
summ$expected_under_null <- 0.05 * summ$n_tests
summ$pct_significant <- round(100 * summ$n_significant / summ$n_tests, 1)

extra <- data.frame(
  exposure = c("exp_harassed", "ai_any"),
  n_tests = 35,
  sig_raw = c(sum(trad_all$exposure == "exp_harassed" & trad_all$p_raw < 0.05),
              sum(trad_all$exposure == "ai_any" & trad_all$p_raw < 0.05)),
  sig_bh = c(sum(trad_all$exposure == "exp_harassed" & trad_all$p_bh < 0.05),
             sum(trad_all$exposure == "ai_any" & trad_all$p_bh < 0.05)),
  stringsAsFactors = FALSE)

cat("\n===== heterogeneity benchmark: 35 tests per method (exposure = time_gt3h) =====\n")
print(summ, row.names = FALSE)
cat("\n===== conventional route on the other two exposures (no forest fitted) =====\n")
print(extra, row.names = FALSE)
cat("\n===== discordant items: conventional significant, forest not =====\n")
print(cmp[cmp$sig_trad_raw & !cmp$sig_cf_raw,
          c("outcome", "moderator", "coef", "p_raw", "coef_cf", "p_cf")], row.names = FALSE)
cat("sign agreement across 35 items:", sum(cmp$sign_agree), "/ 35\n")
cat("DONE\n")
