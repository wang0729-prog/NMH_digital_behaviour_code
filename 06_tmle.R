source("00_config.R")
## 06_tmle.R
## Doubly-robust cross-fitted TMLE for the effect of digital-behavior exposures on
## mental-health outcomes, HMS 2023-24.
## Hand-coded TMLE (binary point treatment) with discrete SuperLearner nuisance estimation
## (candidates: main-effects GLM, GLM + interactions, XGBoost), cluster-aware cross-fitting
## by school (schools are the PSU), and survey-weighted fluctuation / aggregation.
## Exposures analysed: time_gt3h (>3 h/d non-academic screen time), exp_harassed, ai_any.
## Outcomes: dep_any, dep_maj, anx_any, sui_idea, lonely.

suppressPackageStartupMessages({ library(xgboost) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")

m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely","ai_any","ai_n_uses",
            "exp_harassed","time_gt3h","nrweight","schoolnum"))
  m <- m[!is.na(m[[v]]), ]
for (v in W) m <- m[!is.na(m[[v]]), ]
cat("TMLE analytic n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

## model matrix shared by all learners
mm <- model.matrix(~ age_num + female + race_cat + international + undergrad +
                     grad_student + fin_stress + food_insec + firstgen, data = m)
mm_int <- model.matrix(~ (age_num + female + race_cat + international + undergrad +
                            grad_student + fin_stress + food_insec + firstgen)^2,
                       data = m)
mm_int <- mm_int[, apply(mm_int, 2, sd) > 0]

set.seed(20260924)
schools <- unique(m$schoolnum)
set.seed(11)
fold <- rep(NA_integer_, nrow(m))
fs <- sample(rep(1:5, length.out = length(schools)))
names(fs) <- schools
fold <- fs[as.character(m$schoolnum)]
cat("fold sizes:", table(fold), "\n")

## ---------- learners ----------
fit_glm <- function(X, y, w) {
  df <- as.data.frame(X); names(df) <- make.names(names(df))
  df$.y <- y
  suppressWarnings(glm(.y ~ ., data = df, weights = w, family = binomial()))
}
pred_glm <- function(obj, X) {
  df <- as.data.frame(X); names(df) <- make.names(names(df))
  p <- predict(obj, newdata = df, type = "response")
  p <- as.numeric(p); p[is.na(p)] <- mean(obj$y, na.rm = TRUE)
  pmin(pmax(p, 0.01), 0.99)
}
fit_xgb <- function(X, y, w) {
 ## xgboost >= 3.x requires a factor response for binary:logistic in xgboost
  xgboost::xgboost(x = X, y = factor(y, levels = c(0, 1)), weight = w,
                   nrounds = 150, max_depth = 4,
                   learning_rate = 0.08, subsample = 0.8, colsample_bytree = 0.8,
                   objective = "binary:logistic", nthread = 2)
}
pred_xgb <- function(obj, X) {
  p <- predict(obj, newdata = X); pmin(pmax(p, 0.01), 0.99)
}

cv_risk <- function(Xlist, y, w, folds) {
 ## negative log-likelihood risk, weighted
  n <- length(y)
  R <- matrix(NA_real_, n, length(Xlist))
  for (j in seq_along(Xlist)) {
    X <- Xlist[[j]]
    for (k in unique(folds)) {
      tr <- folds != k; te <- folds == k
      if (sum(te) == 0) next
      obj <- if (j == 3) fit_xgb(X[tr, , drop = FALSE], y[tr], w[tr])
             else fit_glm(X[tr, , drop = FALSE], y[tr], w[tr])
      p <- if (j == 3) pred_xgb(obj, X[te, , drop = FALSE])
           else pred_glm(obj, X[te, , drop = FALSE])
      R[te, j] <- -(y[te] * log(p) + (1 - y[te]) * log(1 - p))
    }
  }
  colMeans(R * w) / mean(w)
}

nuisance <- function(A, Y, w) {
 ## discrete SL for g = P(A=1|W) and Qbar(A,W) = E[Y|A,W]; cross-fitted
  g <- rep(NA_real_, length(A)); Q1 <- g; Q0 <- g
  Xlist <- list(mm, mm_int, mm)
  for (k in unique(fold)) {
    tr <- fold != k; te <- fold == k
    wt <- w[tr]
    rg <- cv_risk(lapply(Xlist, function(X) X[tr, , drop = FALSE]), A[tr], wt, fold[tr])
    jg <- which.min(rg)
 ## propensity
    gmod <- if (jg == 3) fit_xgb(Xlist[[jg]][tr, , drop = FALSE], A[tr], wt)
            else fit_glm(Xlist[[jg]][tr, , drop = FALSE], A[tr], wt)
    g[te] <- if (jg == 3) pred_xgb(gmod, Xlist[[jg]][te, , drop = FALSE])
             else pred_glm(gmod, Xlist[[jg]][te, , drop = FALSE])
 ## outcome: fit separately within treated / untreated (saturated in A)
    for (a in 0:1) {
      idx <- tr & A == a
      rQ <- cv_risk(lapply(Xlist, function(X) X[idx, , drop = FALSE]),
                    Y[idx], w[idx], fold[idx])
      jQ <- which.min(rQ)
      qm <- if (jQ == 3) fit_xgb(Xlist[[jQ]][idx, , drop = FALSE], Y[idx], w[idx])
            else fit_glm(Xlist[[jQ]][idx, , drop = FALSE], Y[idx], w[idx])
      p <- if (jQ == 3) pred_xgb(qm, Xlist[[jQ]][te, , drop = FALSE])
           else pred_glm(qm, Xlist[[jQ]][te, , drop = FALSE])
      if (a == 1) Q1[te] <- p else Q0[te] <- p
    }
  }
  list(g = g, Q1 = Q1, Q0 = Q0)
}

tmle <- function(A, Y, w, cluster) {
  n <- length(Y)
  cf <- nuisance(A, Y, w)
  g <- pmin(pmax(cf$g, 0.025), 0.975)
  Q1 <- cf$Q1; Q0 <- cf$Q0
  Q <- ifelse(A == 1, Q1, Q0)
  H1 <- 1 / g; H0 <- -1 / (1 - g)
  H <- ifelse(A == 1, H1, H0)
 ## weighted logistic fluctuation on clever covariate, offset = logit(Q)
  eps <- suppressWarnings(coef(glm(Y ~ -1 + H + offset(qlogis(Q)),
                                   family = binomial(), weights = w)))
  eps <- if (is.na(eps)) 0 else eps
  Q1s <- plogis(qlogis(Q1) + eps / g)
  Q0s <- plogis(qlogis(Q0) - eps / (1 - g))
  Qs <- ifelse(A == 1, Q1s, Q0s)
  psi <- sum(w * (Q1s - Q0s)) / sum(w)
 ## influence curve for the weighted (population) parameter
  D <- (w / mean(w)) * (H * (Y - Qs) + (Q1s - Q0s) - psi)
  se_naive <- sqrt(var(D) / n)
 ## cluster-robust (schools = PSU): Var = G * var(cluster sums) / n^2
  S <- tapply(D, cluster, sum)
  G <- length(S)
  se_cluster <- sqrt(G * var(S)) / n
  p1 <- sum(w[A == 1] * Y[A == 1]) / sum(w[A == 1])
  p0 <- sum(w[A == 0] * Y[A == 0]) / sum(w[A == 0])
  data.frame(psi = psi, eps = eps, rd_unadj = p1 - p0,
             se_naive = se_naive, se_cluster = se_cluster,
             p_naive = 2 * pnorm(-abs(psi / se_naive)),
             p_cluster = 2 * pnorm(-abs(psi / se_cluster)))
}

exps <- list(list(v = "time_gt3h", lab = "screen>3h/d"),
             list(v = "exp_harassed", lab = "online harassment"),
             list(v = "ai_any", lab = "any genAI use"))
ocs <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")

if (Sys.getenv("AIHEAVY") == "1") {       # heavy multi-purpose genAI use instead
  m$ai_heavy <- as.integer(m$ai_n_uses >= 3)
  exps <- list(list(v = "ai_heavy", lab = "genAI 3+ uses"))
  cat("AI-HEAVY MODE; prevalence =", round(mean(m$ai_heavy), 3), "\n")
}

if (Sys.getenv("SMOKE") == "1") {         # fast smoke test
  set.seed(1); keep <- sample(seq_len(nrow(m)), 6000)
  m <- m[keep, ]
  mm <- mm[keep, , drop = FALSE]; mm_int <- mm_int[keep, , drop = FALSE]
  fold <- fold[keep]
  ocs <- "dep_any"; exps <- exps[1]
  cat("SMOKE MODE n =", nrow(m), "\n")
}

res <- list()
for (e in exps) {
  A <- as.integer(m[[e$v]])
  for (o in ocs) {
    Y <- as.integer(m[[o]])
    r <- tmle(A, Y, m$nrweight, m$schoolnum)
    r$exposure <- e$lab; r$outcome <- o
    res[[length(res) + 1]] <- r
    cat(sprintf("%-14s -> %-9s RD=%+.4f (naive SE %.4f, cluster SE %.4f) p_cluster=%.2e\n",
                e$lab, o, r$psi, r$se_naive, r$se_cluster, r$p_cluster))
    flush.console()
  }
}

all <- do.call(rbind, res)
all$OR_tmle <- NA
fn <- if (Sys.getenv("AIHEAVY") == "1") "table3b_tmle_aiheavy_RD.csv" else "table3_tmle_RD.csv"
write.csv(all, file.path(out, fn), row.names = FALSE)
cat("DONE\n")
