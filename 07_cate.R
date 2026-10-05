source("00_config.R")
## 07_cate.R
## Heterogeneous treatment effects via causal forests (grf).
## Treatment: A = time_gt3h (>3 h/d non-academic screen time); Outcomes: 5 binary MH outcomes.
## Moderators examined via best linear projection: fin_stress, food_insec, firstgen,
## female, undergrad, age_num, international.
## School-level clustering and survey weights (nrweight) passed to the forest.

suppressPackageStartupMessages({ library(grf) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

vars <- c("age_num","female","race_cat","international","undergrad","grad_student",
          "fin_stress","food_insec","firstgen")
m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely","time_gt3h",
            "nrweight","schoolnum")) m <- m[!is.na(m[[v]]), ]
for (v in vars) m <- m[!is.na(m[[v]]), ]
cat("CATE analytic n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

X <- model.matrix(~ age_num + female + race_cat + international + undergrad +
                    grad_student + fin_stress + food_insec + firstgen, data = m)[, -1]
X <- X[, apply(X, 2, sd) > 0]
W <- as.numeric(m$time_gt3h)
cl <- as.factor(as.character(m$schoolnum))
levels(cl) <- seq_along(levels(cl))
cl <- as.integer(cl)
sw <- m$nrweight

set.seed(20260924)
res <- list(); blp <- list(); vi <- list()
di <- m[, c("fin_stress","food_insec","firstgen","female","undergrad",
            "age_num","international")]
Dm <- as.matrix(cbind(1, di[, c("fin_stress","food_insec","firstgen","female",
                                "undergrad","age_num","international")]))
colnames(Dm)[1] <- "intercept"

for (o in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) {
  Y <- as.numeric(m[[o]])
  cf <- causal_forest(X, Y, W, clusters = cl, sample.weights = sw,
                      tune.parameters = "all", num.trees = 2000, seed = 20260924)
  ate <- average_treatment_effect(cf, target.sample = "all")
  tau <- predict(cf)$predictions
  b <- best_linear_projection(cf, Dm[, -1, drop = FALSE])
  vimp <- variable_importance(cf)
  res[[o]] <- data.frame(outcome = o,
                         ate = ate[1], se = ate[2],
                         cate_p05 = quantile(tau, .05), cate_p25 = quantile(tau, .25),
                         cate_p50 = quantile(tau, .50), cate_p75 = quantile(tau, .75),
                         cate_p95 = quantile(tau, .95))
 ## grf returns a coeftest matrix (Estimate / Std. Error / t value / Pr(>|t|))
  blp[[o]] <- data.frame(outcome = o, moderator = rownames(b),
                         coef = round(as.numeric(b[, "Estimate"]), 4),
                         se = round(as.numeric(b[, "Std. Error"]), 4),
                         p = round(as.numeric(b[, "Pr(>|t|)"]), 4))
  vi[[o]] <- data.frame(outcome = o, variable = colnames(X),
                        importance = round(as.numeric(vimp), 4))
  cat(sprintf("%-9s ATE=%+.4f (SE %.4f)  CATE IQR=[%.3f, %.3f]\n",
              o, ate[1], ate[2], quantile(tau, .25), quantile(tau, .75)))
  flush.console()
  saveRDS(cf, file.path(out, sprintf("grf_cf_%s.rds", o)))
}

all <- do.call(rbind, res); allb <- do.call(rbind, blp); allv <- do.call(rbind, vi)
write.csv(all,  file.path(out, "table4_cate_ate.csv"), row.names = FALSE)
write.csv(allb, file.path(out, "table4_cate_blp.csv"), row.names = FALSE)
write.csv(allv, file.path(out, "table4_cate_vimp.csv"), row.names = FALSE)
cat("DONE\n")
