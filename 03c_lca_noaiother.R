source("00_config.R")
## 03c_lca_noaiother.R
## Robust LCA: fit k=1..8, save each model + stats IMMEDIATELY after fitting
## so a mid-run kill never loses completed models.
## nrep reduced to 5 (k<=6 previously fitted with nrep=10; results in 03_lca.log
## serve as a convergence cross-check).

suppressPackageStartupMessages({ library(poLCA) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
lcaout <- file.path(out, "lca2")
dir.create(lcaout, showWarnings = FALSE)

d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

ind <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
         "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun" )

dd <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in ind) dd <- dd[!is.na(dd[[v]]), ]
cat("LCA analytic n:", nrow(dd), "\n")
saveRDS(dd$schoolnum, file.path(lcaout, "id_schoolnum.rds"))

X <- as.data.frame(dd[, ind])
X$time5cat <- as.integer(X$time5cat) + 1
for (v in setdiff(ind, "time5cat")) X[[v]] <- as.integer(X[[v]]) + 1
f <- as.formula(paste("cbind(", paste(ind, collapse = ","), ") ~ 1"))

set.seed(20260924)

for (k in 1:8) {
  fk <- file.path(lcaout, sprintf("fit_k%d.rds", k))
  if (file.exists(fk)) { cat("skip k =", k, "(exists)\n"); next }
  cat("---- fitting k =", k, "----\n"); flush.console()
  m <- try(poLCA(f, X, nclass = k, nrep = 5, maxiter = 5000,
                 calc.se = FALSE, verbose = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) {
    cat("k =", k, "FAILED\n"); flush.console(); next
  }
  saveRDS(m, fk)

  P <- m$posterior
  ent <- 1 + sum(P * log(P + 1e-12)) / (nrow(P) * log(k))
  pr <- colMeans(P)
  st <- data.frame(k = k, logLik = m$llik, BIC = m$bic, AIC = m$aic,
                   Gsq = m$Gsq, Chisq = m$Chisq, entropy = round(ent, 4),
                   min_class_pct = round(100 * min(pr), 2),
                   classes_pct = paste(round(100 * pr, 1), collapse = "/"))
  write.csv(st, file.path(lcaout, sprintf("stats_k%d.csv", k)), row.names = FALSE)
  cat(sprintf("k=%d BIC=%.0f AIC=%.0f entropy=%.3f min%%=%.2f classes=%s\n",
              k, m$bic, m$aic, ent, 100 * min(pr), st$classes_pct)); flush.console()

  rows <- list()
  for (v in ind) {
    cp <- m$probs[[v]]
    for (lv in 1:ncol(cp)) {
      rows[[length(rows) + 1]] <- data.frame(k = k, item = v, level = lv - 1,
                                             t(data.frame(round(cp[, lv], 4))))
    }
  }
  tab <- do.call(rbind, rows)
  names(tab)[4:ncol(tab)] <- paste0("class", 1:k)
  write.csv(tab, file.path(lcaout, sprintf("condprob_k%d.csv", k)), row.names = FALSE)
}

cat("DONE\n")
