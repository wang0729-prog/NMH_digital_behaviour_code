source("00_config.R")
## 12_sens_mi.R
## Sensitivity: multiple imputation of missing covariates (vs complete-case n=41,862
## and missing-category n=44,679). m = 20 imputations; weighted svyglm (schools = PSU)
## within each imputed data set; Rubin's rules for pooling (Barnard-Rubin df).
## Focus: 3 exposures x 3 outcomes (keeps the run tractable).

suppressPackageStartupMessages({ library(mice); library(survey); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

oc <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")
W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")

## base: module + exposure + all outcomes observed (covariates allowed to be missing)
d <- zap_labels(d)          # haven_labelled -> plain numeric (mice cannot coerce them)
base <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in oc) base <- base[!is.na(base[[v]]), ]
base$ai_heavy <- as.integer(base$ai_n_uses >= 3)
cat("MI base n:", nrow(base), "\n")
print(colSums(is.na(base[, W])))

X <- base[, c(W, "dep_any","sui_idea","lonely","time5cat","exp_harassed",
              "ai_n_uses","nrweight")]
## race_cat coded numerically + pmm: only 86 values missing (0.2%); polyreg with 9
## levels on 44k rows would dominate run time for negligible gain.
X$race_code <- as.integer(factor(base$race_cat))
X <- X[, -which(names(X) == "race_cat")]
for (v in c("female","international","undergrad","grad_student","fin_stress",
            "food_insec","firstgen")) X[[v]] <- factor(X[[v]])
X$time5cat <- factor(X$time5cat)

meth <- c(age_num = "pmm", female = "logreg", race_code = "pmm",
          international = "logreg", undergrad = "logreg", grad_student = "logreg",
          fin_stress = "logreg", food_insec = "logreg", firstgen = "logreg",
          dep_any = "", sui_idea = "", lonely = "", time5cat = "",
          exp_harassed = "", ai_n_uses = "", nrweight = "")

set.seed(20260924)
t0 <- Sys.time()
imp <- mice(X, m = 20, maxit = 8, method = meth, printFlag = FALSE)
cat("imputation time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(imp, file.path(out, "mi_imp.rds"))

## ---- pooled estimation ----
adj <- "age_num + female + race_code + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"
## Exposure coding is matched EXACTLY to Table 2 (04_var_centered.R) so that the
## imputed and complete-case estimates are directly comparable:
## time5cat_f : 5-level factor, reference = "none" -> term time5cat_f>3h
## exp_harassed: binary -> term exp_harassed
## ai_n_uses_g : 0 / 1 / 2 / 3+, reference = 0 -> term ai_n_uses_g3+
exps <- list(list(v = "time5cat_f", lab = "screen>3h/d", term = "time5cat_f>3h"),
             list(v = "exp_harassed", lab = "harassment", term = "exp_harassed"),
             list(v = "ai_n_uses_g", lab = "genAI 3+ uses", term = "ai_n_uses_g3+"))
ocs <- c("dep_any","sui_idea","lonely")

pool_rubin <- function(bs, Vs, n_clust, n_par) {
  m <- length(bs)
  bm <- do.call(rbind, bs)
  Q <- colMeans(bm)
  U <- Reduce("+", Vs) / m
  Udiag <- mean(sapply(Vs, function(v) diag(v)))   # mean within-imputation variance (diag)
  B <- apply(bm, 2, var)
  Td <- Udiag + (1 + 1 / m) * B
  se <- sqrt(Td)
  lam <- (1 + 1 / m) * B / Td                      # fraction of missing information
  nu <- (m - 1) / lam^2                            # Barnard-Rubin
  nu_com <- max(n_clust - n_par, 5)                # complete-data df (clusters - params)
  nu_adj <- nu * nu_com / (nu + nu_com)
  p <- 2 * pt(-abs(Q / se), df = nu_adj)
  data.frame(term = names(Q), est = Q, se = se, p = p, fmi = lam)
}

res <- list()
for (e in exps) {
  for (o in ocs) {
    bs <- list(); Vs <- list(); nc <- NA
    for (j in 1:imp$m) {
      dat <- complete(imp, j)
      dat$race_code <- factor(as.integer(round(dat$race_code)))
      dat$time5cat_f <- factor(base$time5cat, levels = 0:4,
                               labels = c("none","<1h","1-2h","2-3h",">3h"))
      dat$time5cat_f <- relevel(dat$time5cat_f, ref = "none")
      dat$ai_n_uses_g <- cut(base$ai_n_uses, breaks = c(-0.5, 0.5, 1.5, 2.5, 100),
                             labels = c("0","1","2","3+"))
      dat$ai_n_uses_g <- relevel(dat$ai_n_uses_g, ref = "0")
      dat$schoolnum <- base$schoolnum
      dat$nrweight <- base$nrweight
      dd <- dat[!is.na(dat[[e$v]]) & !is.na(dat[[o]]), ]
      des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dd)
      g <- svyglm(as.formula(paste(o, "~", e$v, "+", adj)), design = des,
                  family = quasibinomial())
      cc <- coef(g); keep <- which(names(cc) == e$term)
      bs[[j]] <- cc[keep]; Vs[[j]] <- vcov(g)[keep, keep, drop = FALSE]
      nc <- length(unique(dd$schoolnum))
    }
    pr <- pool_rubin(bs, Vs, nc, n_par = 10)
    pr$exposure <- e$lab; pr$outcome <- o
    pr$OR <- exp(pr$est); pr$LCL <- exp(pr$est - 1.96 * pr$se)
    pr$UCL <- exp(pr$est + 1.96 * pr$se)
    res[[length(res) + 1]] <- pr[, c("outcome","exposure","term","OR","LCL","UCL","p","fmi")]
    cat(sprintf("%-14s -> %-9s OR=%.3f [%.3f-%.3f] p=%.4g  fmi=%.3f\n",
                e$lab, o, pr$OR, pr$LCL, pr$UCL, pr$p, pr$fmi)); flush.console()
  }
}
all <- do.call(rbind, res)

## ---- like-for-like complete-case estimates on the SAME base sample ----
ccb <- base
complete_cov <- complete.cases(ccb[, W])
ccb <- ccb[complete_cov, ]
ccb$race_code <- factor(as.integer(factor(ccb$race_cat)))  # same coding as imputed sets
ccb$time5cat_f <- relevel(factor(ccb$time5cat, levels = 0:4,
                                 labels = c("none","<1h","1-2h","2-3h",">3h")), ref = "none")
ccb$ai_n_uses_g <- relevel(cut(ccb$ai_n_uses, breaks = c(-0.5, 0.5, 1.5, 2.5, 100),
                               labels = c("0","1","2","3+")), ref = "0")
ccres <- list()
for (e in exps) {
  for (o in ocs) {
    dd <- ccb[!is.na(ccb[[e$v]]) & !is.na(ccb[[o]]), ]
    des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dd)
    g <- svyglm(as.formula(paste(o, "~", e$v, "+", adj)), design = des,
                family = quasibinomial())
    i <- which(names(coef(g)) == e$term)
    ci <- confint(g)
    ccres[[length(ccres) + 1]] <- data.frame(outcome = o, exposure = e$lab,
                                             OR = round(exp(coef(g)[i]), 3),
                                             LCL = round(exp(ci[i, 1]), 3),
                                             UCL = round(exp(ci[i, 2]), 3), n = nrow(dd))
  }
}
ccres <- do.call(rbind, ccres)
cmp <- merge(all[, c("outcome","exposure","OR","LCL","UCL","p","fmi")], ccres,
             by = c("outcome","exposure"), suffixes = c("_MI","_CC"))
cat("\n=== MI (m=20, pooled) vs complete-case, identical specification ===\n")
print(cmp[, c("outcome","exposure","OR_MI","LCL_MI","UCL_MI","OR_CC","LCL_CC","UCL_CC","p","fmi","n")],
      row.names = FALSE)
write.csv(cmp, file.path(out, "table_sens_mi_vs_cc.csv"), row.names = FALSE)
write.csv(all, file.path(out, "table_sens_mi_pooled.csv"), row.names = FALSE)
cat("DONE\n")
