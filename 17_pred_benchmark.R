source("00_config.R")
## 17_pred_benchmark.R
## Prediction layer of the methodological benchmark.
##
## Question: on this prediction task, does gradient boosting actually beat the conventional
## learners (weighted main-effects logistic regression, LASSO)? Or is the ML machinery
## decorative?
##
## Design: identical data, identical feature matrix, identical outer 5-fold CV grouped by
## school, identical inner 3-fold for hyper-parameter selection, identical weighted metrics.
## Only the learner differs.
##
## Metrics: weighted AUC, Brier score, calibration slope/intercept, decision-curve net benefit.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(xgboost); library(glmnet) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

base <- c("age_num", "female", "race_cat", "international", "undergrad", "grad_student")
sdoh <- c("fin_stress", "food_insec", "firstgen")
digi <- c("time5cat", "exp_relation", "exp_harassed", "exp_create", "exp_activism",
          "ai_academic", "ai_forbidden", "ai_work", "ai_comm", "ai_health", "ai_fun",
          "ai_n_uses")

oc <- "dep_any"
m <- d[d$module == 1, ]
for (v in c(oc, base, sdoh, digi, "nrweight", "schoolnum")) m <- m[!is.na(m[[v]]), ]
rownames(m) <- NULL
cat("prediction benchmark n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

mkX <- function(dat, vars) {
  X <- model.matrix(as.formula(paste("~", paste(vars, collapse = " + "))), data = dat)
  X <- X[, -1, drop = FALSE]
  X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  colnames(X) <- make.names(colnames(X))
  X
}

Y <- as.numeric(m[[oc]])
w <- m$nrweight

set.seed(20260924)
sch <- unique(m$schoolnum)
fs <- sample(rep(1:5, length.out = length(sch))); names(fs) <- sch
fold <- fs[as.character(m$schoolnum)]
inner <- sample(rep(1:3, length.out = nrow(m)))     ## shared inner folds for all learners

auc_w <- function(p, y, w) {           ## weighted AUC = P(score_pos > score_neg)
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  n1 <- sum(w * y); n0 <- sum(w * (1 - y))
  if (n1 == 0 || n0 == 0) return(NA_real_)
  cn <- cumsum(w * (1 - y))
  (sum((y * w) * cn)) / (n1 * n0)
}
brier_w <- function(p, y, w) sum(w * (y - p)^2) / sum(w)
nb <- function(p, y, w, pt) {
  fl <- as.integer(p >= pt)
  tp <- sum(w * y * fl); fp <- sum(w * (1 - y) * fl); N <- sum(w)
  tp / N - (fp / N) * (pt / (1 - pt))
}

xgb_fit <- function(X, y, w, depth, nr) {
  xgboost::xgboost(x = X, y = factor(y, levels = c(0, 1)), weight = w,
                   nrounds = nr, max_depth = depth, learning_rate = 0.05,
                   subsample = 0.8, colsample_bytree = 0.8,
                   objective = "binary:logistic", eval_metric = "auc", nthread = 2)
}

## ---- three learners, identical fit/predict interface: (Xtr, ytr, wtr, inn) ----
L_logistic <- list(
  name = "logistic",
  fit = function(Xtr, ytr, wtr, inn) {
    df <- data.frame(.y = ytr, Xtr)
    suppressWarnings(glm(.y ~ ., data = df, weights = wtr, family = binomial()))
  },
  pred = function(obj, Xte) {
    df <- data.frame(Xte)
    as.numeric(predict(obj, newdata = df, type = "response"))
  })

L_lasso <- list(
  name = "lasso",
  fit = function(Xtr, ytr, wtr, inn) {
    cv <- suppressWarnings(cv.glmnet(Xtr, ytr, weights = wtr, family = "binomial",
                                     alpha = 1, foldid = inn, type.measure = "deviance"))
    list(cv = cv, lam = cv$lambda.min)
  },
  pred = function(obj, Xte) {
    as.numeric(predict(obj$cv, Xte, s = obj$lam, type = "response"))
  })

L_xgb <- list(
  name = "xgboost",
  fit = function(Xtr, ytr, wtr, inn) {
    best <- -Inf; bd <- 3
    for (dd in c(3, 5)) {
      a <- sapply(1:3, function(j) {
        itr <- inn != j
        mo <- xgb_fit(Xtr[itr, , drop = FALSE], ytr[itr], wtr[itr], dd, nr = 150)
        pr <- as.numeric(predict(mo, Xtr[!itr, , drop = FALSE]))
        -brier_w(pr, ytr[!itr], wtr[!itr])
      })
      sc <- mean(a)
      if (sc > best) { best <- sc; bd <- dd }
    }
    xgb_fit(Xtr, ytr, wtr, bd, nr = 300)
  },
  pred = function(obj, Xte) as.numeric(predict(obj, Xte)))

run_learner <- function(X, L, label) {
  p <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    tr <- fold != k; te <- fold == k
    obj <- L$fit(X[tr, , drop = FALSE], Y[tr], w[tr], inner[tr])
    p[te] <- L$pred(obj, X[te, , drop = FALSE])
    cat(sprintf("  %s / %s fold %d done\n", label, L$name, k)); flush.console()
  }
  lp <- log(p / (1 - p))
  cal <- suppressWarnings(glm(Y ~ lp, weights = w, family = binomial()))
  list(p = p, auc = auc_w(p, Y, w), brier = brier_w(p, Y, w),
       cal_slope = unname(coef(cal)[2]), cal_int = unname(coef(cal)[1]),
       nb = sapply(seq(0.05, 0.60, 0.05), function(pt) nb(p, Y, w, pt)))
}

sets <- list(A_base_SDOH = c(base, sdoh), B_plus_digital = c(base, sdoh, digi))
learners <- list(L_logistic, L_lasso, L_xgb)

rows <- list(); dca <- list()
for (sn in names(sets)) {
  X <- mkX(m, sets[[sn]])
  cat("=== feature set", sn, " p =", ncol(X), "===\n")
  for (L in learners) {
    r <- run_learner(X, L, sn)
    rows[[length(rows) + 1]] <- data.frame(
      feature_set = sn, learner = L$name, n_features = ncol(X),
      AUC = r$auc, Brier = r$brier,
      cal_slope = r$cal_slope, cal_intercept = r$cal_int,
      stringsAsFactors = FALSE)
    dca[[length(dca) + 1]] <- data.frame(
      feature_set = sn, learner = L$name,
      threshold = seq(0.05, 0.60, 0.05), net_benefit = r$nb,
      stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, rows)
## delta vs the conventional logistic learner, within each feature set
res$delta_AUC_vs_logistic <- NA_real_
for (sn in unique(res$feature_set)) {
  idx <- res$feature_set == sn
  ref <- res$AUC[idx & res$learner == "logistic"]
  res$delta_AUC_vs_logistic[idx] <- round(res$AUC[idx] - ref, 4)
}
res[, c("AUC", "Brier")] <- round(res[, c("AUC", "Brier")], 4)
res[, c("cal_slope", "cal_intercept")] <- round(res[, c("cal_slope", "cal_intercept")], 3)
write.csv(res, file.path(out, "table13_pred_benchmark.csv"), row.names = FALSE)

dca_all <- do.call(rbind, dca)
dca_all$net_benefit <- round(dca_all$net_benefit, 5)
write.csv(dca_all, file.path(out, "table13_pred_dca.csv"), row.names = FALSE)

cat("\n===== prediction benchmark =====\n")
print(res, row.names = FALSE)
cat("DONE\n")
