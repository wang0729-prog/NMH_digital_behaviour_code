source("00_config.R")
## 08_prediction.R
## Incremental predictive value of digital-behaviour features for mental-health outcomes.
## Design: outer 5-fold CV grouped by school (PSU) -> honest AUC / Brier / calibration;
## inner 3-fold for xgboost depth selection.
## Model A: demographics + SDOH only.
## Model B: A + digital behaviour (time, online experiences, genAI uses).
## Output: AUC delta, calibration slope/intercept, decision-curve net benefit,
## and built-in TreeSHAP importance (predcontrib).

suppressPackageStartupMessages({ library(xgboost) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

base <- c("age_num","female","race_cat","international","undergrad","grad_student")
sdoh <- c("fin_stress","food_insec","firstgen")
digi <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
          "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun",
          "ai_n_uses")

oc <- "dep_any"
m <- d[d$module == 1, ]
for (v in c(oc, base, sdoh, digi, "nrweight","schoolnum")) m <- m[!is.na(m[[v]]), ]
cat("prediction n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

mkX <- function(dat, vars) {
  X <- model.matrix(as.formula(paste("~", paste(vars, collapse = " + "))), data = dat)
  X <- X[, -1, drop = FALSE]
  X[, apply(X, 2, sd) > 0, drop = FALSE]
}
XA <- mkX(m, c(base, sdoh))
XB <- mkX(m, c(base, sdoh, digi))
Y  <- as.numeric(m[[oc]])
w  <- m$nrweight

set.seed(20260924)
sch <- unique(m$schoolnum)
fs <- sample(rep(1:5, length.out = length(sch))); names(fs) <- sch
fold <- fs[as.character(m$schoolnum)]

fit_xgb <- function(X, y, w, depth, nr = 300) {
  xgboost::xgboost(x = X, y = factor(y, levels = c(0, 1)), weight = w,
                   nrounds = nr, max_depth = depth, learning_rate = 0.05,
                   subsample = 0.8, colsample_bytree = 0.8,
                   objective = "binary:logistic", eval_metric = "auc",
                   nthread = 2)
}

auc_w <- function(p, y, w) {           # weighted AUC = P(score_pos > score_neg)
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  n1 <- sum(w * y); n0 <- sum(w * (1 - y))
  if (n1 == 0 || n0 == 0) return(NA)
  cn <- cumsum(w * (1 - y))            # weight of negatives scored below i
  pos <- y * w
  sum(pos * cn) / (n1 * n0)
}
brier_w <- function(p, y, w) sum(w * (y - p)^2) / sum(w)

nb <- function(p, y, w, pt) {          # net benefit at threshold pt
  fl <- as.integer(p >= pt)
  tp <- sum(w * y * fl); fp <- sum(w * (1 - y) * fl); N <- sum(w)
  tp / N - (fp / N) * (pt / (1 - pt))
}

run_model <- function(X, label) {
  p <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    tr <- fold != k; te <- fold == k
 ## inner selection of depth
    inn <- sample(rep(1:3, length.out = sum(tr)))
    best <- -Inf; bd <- 3
    for (dd in c(3, 5)) {
      a <- sapply(1:3, function(j) {
        itr <- tr; itr[tr] <- inn != j
        mo <- fit_xgb(X[itr, , drop = FALSE], Y[itr], w[itr], dd, nr = 150)
        pr <- as.numeric(predict(mo, X[tr & !itr, , drop = FALSE]))
        -brier_w(pr, Y[tr & !itr], w[tr & !itr])
      })
      sc <- mean(a)
      if (sc > best) { best <- sc; bd <- dd }
    }
    mo <- fit_xgb(X[tr, , drop = FALSE], Y[tr], w[tr], bd, nr = 300)
    p[te] <- as.numeric(predict(mo, X[te, , drop = FALSE]))
    cat(sprintf("  %s fold %d: depth=%d\n", label, k, bd)); flush.console()
  }
 ## calibration: weighted logistic of Y on logit(p)
  lp <- log(p / (1 - p))
  cal <- suppressWarnings(glm(Y ~ lp, weights = w, family = binomial()))
  list(p = p, auc = auc_w(p, Y, w), brier = brier_w(p, Y, w),
       cal_slope = coef(cal)[2], cal_int = coef(cal)[1],
       nb = sapply(seq(0.05, 0.60, 0.05), function(pt) nb(p, Y, w, pt)),
       xgb_model = NULL)
}

cat("=== Model A: demographics + SDOH ===\n"); rA <- run_model(XA, "A")
cat("=== Model B: + digital behaviour ===\n"); rB <- run_model(XB, "B")

## TreeSHAP importance for model B.
## NOTE: xgboost >= 3.x silently ignores predcontrib for models fit with xgboost;
## TreeSHAP requires the xgb.train + xgb.DMatrix route.
mb <- fit_xgb(XB, Y, w, 5, nr = 300)
bst <- xgboost::xgb.train(
  params = list(objective = "binary:logistic", max_depth = 5, eta = 0.05,
                subsample = 0.8, colsample_bytree = 0.8),
  data = xgboost::xgb.DMatrix(XB, label = Y, weight = w),
  nrounds = 300, nthread = 2)
sh <- predict(bst, xgboost::xgb.DMatrix(XB), predcontrib = TRUE)
sh <- sh[, seq_len(ncol(XB)), drop = FALSE]
imp <- data.frame(variable = colnames(XB),
                  mean_abs_shap = round(colMeans(abs(sh)), 5))
imp <- imp[order(-imp$mean_abs_shap), ]
write.csv(imp, file.path(out, "table5_shap_importance.csv"), row.names = FALSE)
saveRDS(mb, file.path(out, "xgb_modelB.rds"))

res <- data.frame(
  model = c("A_base_SDOH","B_plus_digital"),
  AUC = round(c(rA$auc, rB$auc), 4),
  Brier = round(c(rA$brier, rB$brier), 4),
  cal_slope = round(c(rA$cal_slope, rB$cal_slope), 3),
  cal_intercept = round(c(rA$cal_int, rB$cal_int), 3))
res$delta_AUC <- c(NA, round(rB$auc - rA$auc, 4))
write.csv(res, file.path(out, "table5_prediction.csv"), row.names = FALSE)
print(res, row.names = FALSE)

dca <- data.frame(threshold = seq(0.05, 0.60, 0.05),
                  nb_treat_all = sapply(seq(0.05, 0.60, 0.05),
                                        function(pt) {
                                          pr <- mean(Y)
                                          sum(w * Y) / sum(w) -
                                            (sum(w * (1 - Y)) / sum(w)) * (pt / (1 - pt))
                                        }),
                  nb_A = rA$nb, nb_B = rB$nb)
write.csv(dca, file.path(out, "table5_dca.csv"), row.names = FALSE)
cat("Top SHAP variables:\n"); print(head(imp, 12), row.names = FALSE)
cat("DONE\n")
