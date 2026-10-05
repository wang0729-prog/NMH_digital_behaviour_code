source("00_config.R")
## 99d_verify3.R
## Independent re-computation of the three methodological benchmarks (Tables 12-14).
## Reads only the written CSVs and re-derives every headline number from scratch.

dir <- OUT
chk <- function(label, got, want, tol = 1e-6) {
  ok <- isTRUE(all.equal(got, want, tolerance = tol))
  cat(sprintf("[%s] %-52s got=%s want=%s\n", if (ok) "PASS" else "FAIL",
              label, format(got, digits = 6), format(want, digits = 6)))
  ok
}
res <- c()

## ---------- Table 12: causal estimation layer ----------
t12 <- read.csv(file.path(dir, "table12_method_benchmark.csv"))
res <- c(res, chk("T12 sample-identity self-check all PASS",
                  all(t12$sample_check == "PASS"), TRUE))
res <- c(res, chk("T12 n rows", nrow(t12), 15))
t12$rel <- 100 * (t12$rd_conv_gcomp - t12$rd_tmle) / t12$rd_conv_gcomp
for (e in unique(t12$exposure)) {
  s <- t12[t12$exposure == e, ]
  want <- switch(e,
                 "screen>3h/d"       = c(0.26, 2.05),
                 "online harassment" = c(0.24, 6.66),
                 "any genAI use"     = c(6.44, 30.98))
  res <- c(res, chk(paste("T12 rel-diff min", e), round(min(abs(s$rel)), 2), want[1], 0.02))
  res <- c(res, chk(paste("T12 rel-diff max", e), round(max(abs(s$rel)), 2), want[2], 0.02))
}
ai <- t12[t12$exposure == "any genAI use", ]
res <- c(res, chk("T12 AI dep_maj RD conv -> TMLE",
                  c(round(ai$rd_conv_gcomp[ai$outcome == "dep_maj"], 4),
                    round(ai$rd_tmle[ai$outcome == "dep_maj"], 4)),
                  c(0.0126, 0.0087), 1e-3))
res <- c(res, chk("T12 AI: 4 of 5 outcomes non-significant under TMLE",
                  sum(ai$p_tmle_cluster >= 0.05), 4))

## ---------- Table 13: prediction layer ----------
t13 <- read.csv(file.path(dir, "table13_pred_benchmark.csv"))
res <- c(res, chk("T13 n rows", nrow(t13), 6))
b <- t13[t13$feature_set == "B_plus_digital", ]
res <- c(res, chk("T13 B: XGB - logistic delta AUC",
                  round(b$AUC[b$learner == "xgboost"] - b$AUC[b$learner == "logistic"], 4),
                  0.0013, 1e-4))
a <- t13[t13$feature_set == "A_base_SDOH", ]
res <- c(res, chk("T13 A: XGB - logistic delta AUC",
                  round(a$AUC[a$learner == "xgboost"] - a$AUC[a$learner == "logistic"], 4),
                  -0.0002, 1e-4))
res <- c(res, chk("T13 digital-feature gain (B logistic - A logistic)",
                  round(b$AUC[b$learner == "logistic"] - a$AUC[a$learner == "logistic"], 4),
                  0.0318, 1e-3))
## cross-check against the earlier standalone run (08_prediction.R)
t5 <- read.csv(file.path(dir, "table5_prediction.csv"))
res <- c(res, chk("T13 XGB AUC reproduces 08_prediction.R (set A)",
                  round(a$AUC[a$learner == "xgboost"], 4),
                  round(t5$AUC[t5$model == "A_base_SDOH"], 4), 1e-4))
res <- c(res, chk("T13 XGB AUC reproduces 08_prediction.R (set B)",
                  round(b$AUC[b$learner == "xgboost"], 4),
                  round(t5$AUC[t5$model == "B_plus_digital"], 4), 1e-4))

## ---------- Table 14: heterogeneity layer ----------
t14 <- read.csv(file.path(dir, "table14_het_benchmark.csv"))
res <- c(res, chk("T14 n tests per method", nrow(t14), 35))
res <- c(res, chk("T14 conventional raw significant", sum(t14$p_raw < 0.05), 8))
res <- c(res, chk("T14 conventional BH significant", sum(t14$p_trad_bh < 0.05), 2))
res <- c(res, chk("T14 forest raw significant", sum(t14$p_cf < 0.05), 4))
res <- c(res, chk("T14 forest BH significant", sum(t14$p_cf_bh < 0.05), 0))
res <- c(res, chk("T14 discordant (conv sig, forest ns)",
                  sum(t14$p_raw < 0.05 & t14$p_cf >= 0.05), 6))
## invariants
res <- c(res, chk("T14 BH-adjusted p >= raw p (conventional)",
                  all(t14$p_trad_bh >= t14$p_raw - 1e-12), TRUE))
res <- c(res, chk("T14 BH-adjusted p >= raw p (forest)",
                  all(t14$p_cf_bh >= t14$p_cf - 1e-12), TRUE))
## sign reversal on dep_maj
dm <- t14[t14$outcome == "dep_maj", ]
res <- c(res, chk("T14 sign reversals on dep_maj",
                  sum(sign(dm$coef) != sign(dm$coef_cf)), 2))
res <- c(res, chk("T14 sign agreement across 35 items", sum(t14$sign_agree == TRUE), 28))

t14b <- read.csv(file.path(dir, "table14b_trad_interactions.csv"))
res <- c(res, chk("T14b n conventional tests", nrow(t14b), 105))
res <- c(res, chk("T14b exp_harassed raw significant",
                  sum(t14b$exposure == "exp_harassed" & t14b$p_raw < 0.05), 0))
res <- c(res, chk("T14b ai_any raw significant",
                  sum(t14b$exposure == "ai_any" & t14b$p_raw < 0.05), 2))
res <- c(res, chk("T14b ai_any BH significant",
                  sum(t14b$exposure == "ai_any" & t14b$p_bh < 0.05), 1))
res <- c(res, chk("T14b BH-adjusted p >= raw p (all 105)",
                  all(t14b$p_bh >= t14b$p_raw - 1e-12), TRUE))

cat(sprintf("\n===== %d checks: %d PASS, %d FAIL =====\n",
            length(res), sum(res), sum(!res)))
