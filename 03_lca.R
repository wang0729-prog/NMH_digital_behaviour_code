source("00_config.R")
## 03_lca.R
## Latent class analysis of digital-behavior profiles, HMS 2023-24 main sample.
## Indicators (12): time5cat (ordinal 0-4, treated as polytomous),
## exp_relation, exp_harassed, exp_create, exp_activism (binary),
## ai_academic, ai_forbidden, ai_work, ai_comm, ai_health, ai_fun, ai_other (binary)
## exp_none excluded: it is the logical complement of the other 4 experience items
## and would violate local independence by construction.
## Unweighted LCA (poLCA has no native weights); weighted prevalence computed post-hoc.

suppressPackageStartupMessages({ library(poLCA); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")

d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

ind <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
         "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun","ai_other")

dd <- d[d$module==1 & !is.na(d$time5cat), ]
for (v in ind) dd <- dd[!is.na(dd[[v]]), ]
cat("LCA analytic n:", nrow(dd), "\n")

X <- as.data.frame(dd[, ind])
# poLCA requires consecutive integers starting at 1
X$time5cat <- X$time5cat + 1
for (v in setdiff(ind, "time5cat")) X[[v]] <- X[[v]] + 1
# ensure numeric not labelled
for (v in ind) X[[v]] <- as.integer(X[[v]])

f <- as.formula(paste("cbind(", paste(ind, collapse=","), ") ~ 1"))

set.seed(20260923)
fits <- list()
stats <- data.frame()

for (k in 1:8) {
  cat("---- fitting k =", k, "----\n")
  m <- poLCA(f, X, nclass=k, nrep=10, maxiter=5000, verbose=FALSE)
  fits[[k]] <- m
 # entropy from posterior matrix
  P <- m$posterior
  ent <- 1 + sum(P * log(P + 1e-12)) / (nrow(P) * log(k))
  pr <- colMeans(P)
  stats <- rbind(stats, data.frame(
    k=k, logLik=m$llik, BIC=m$bic, AIC=m$aic, Gsq=m$Gsq, Chisq=m$Chisq,
    entropy=round(ent,4), min_class_pct=round(100*min(pr),2),
    classes_pct=paste(round(100*pr,1), collapse="/")
  ))
  cat(sprintf("k=%d  BIC=%.0f  AIC=%.0f  entropy=%.3f  min%%=%.2f\n",
              k, m$bic, m$aic, ent, 100*min(pr)))
}

write.csv(stats, file.path(out, "lca_fit_stats.csv"), row.names=FALSE)

## Save conditional-probability tables for every k >= 2
for (k in 2:8) {
  m <- fits[[k]]
  pr <- colMeans(m$posterior)
  rows <- list()
  for (v in ind) {
    cp <- m$probs[[v]]  # matrix nclass x nlevels
    for (lv in 1:ncol(cp)) {
      rows[[length(rows)+1]] <- data.frame(
        k=k, item=v, level=lv - 1,  # back to original coding (0/1; time 0-4)
        t(data.frame(round(cp[,lv], 4)))
      )
    }
  }
  tab <- do.call(rbind, rows)
  names(tab)[4:ncol(tab)] <- paste0("class", 1:k)
  write.csv(tab, file.path(out, sprintf("lca_condprob_k%d.csv", k)), row.names=FALSE)
}

saveRDS(fits, file.path(out, "lca_fits_k1_8.rds"))
cat("DONE\n")
