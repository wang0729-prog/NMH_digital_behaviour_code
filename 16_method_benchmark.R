source("00_config.R")
## 16_method_benchmark.R
## Methodological benchmark for the "AI-driven analysis" claim.
##
## Question: does ML-based causal inference (TMLE + discrete SuperLearner, cross-fitted by
## school, doubly robust) give a DIFFERENT answer from the conventional route used throughout
## this literature (survey-weighted logistic regression + plug-in g-computation, main effects
## only, no cross-fitting)?
##
## Design: identical complete-case sample, identical covariate set, identical survey weights.
## Only the estimation machinery differs, so any divergence is attributable to the method.
##
## Self-check: the crude (unadjusted) weighted risk difference recomputed here must reproduce
## rd_unadj stored by 06_tmle.R bit-for-bit. If it does, the samples are provably identical.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")

B <- as.integer(Sys.getenv("BOOT", "0"))   ## set BOOT=200 for cluster bootstrap SEs

d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

W <- c("age_num", "female", "race_cat", "international", "undergrad", "grad_student",
       "fin_stress", "food_insec", "firstgen")

## ---- complete-case sample: byte-identical construction to 06_tmle.R ----
m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely", "ai_any", "ai_n_uses",
            "exp_harassed", "time_gt3h", "nrweight", "schoolnum"))
  m <- m[!is.na(m[[v]]), ]
for (v in W) m <- m[!is.na(m[[v]]), ]
rownames(m) <- NULL                       ## guard against the model.matrix row-name trap
cat("benchmark analytic n:", nrow(m), " schools:", length(unique(m$schoolnum)), "\n")

for (v in c("time_gt3h", "exp_harassed", "ai_any")) {
  stopifnot(all(m[[v]] %in% c(0, 1)))
  m[[v]] <- as.integer(m[[v]])
}

tm <- read.csv(file.path(out, "table3_tmle_RD.csv"), stringsAsFactors = FALSE)
adj <- paste(W, collapse = " + ")

## ---- conventional route: main-effects weighted logistic + plug-in g-computation ----
gcomp <- function(dat, outcome, exposure) {
  fo <- as.formula(paste(outcome, "~", exposure, "+", adj))
  g  <- glm(fo, data = dat, weights = nrweight, family = binomial())
  d1 <- dat; d1[[exposure]] <- 1L
  d0 <- dat; d0[[exposure]] <- 0L
  w  <- dat$nrweight
  p1 <- sum(w * predict(g, newdata = d1, type = "response")) / sum(w)
  p0 <- sum(w * predict(g, newdata = d0, type = "response")) / sum(w)
  c(p1 = p1, p0 = p0, rd = p1 - p0, or_cond = exp(unname(coef(g)[exposure])))
}

crude_rd <- function(dat, outcome, exposure) {
  w <- dat$nrweight; a <- dat[[exposure]]; y <- dat[[outcome]]
  sum(w[a == 1] * y[a == 1]) / sum(w[a == 1]) - sum(w[a == 0] * y[a == 0]) / sum(w[a == 0])
}

exps <- c(time_gt3h = "screen>3h/d", exp_harassed = "online harassment", ai_any = "any genAI use")
outs <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")

## ---- fast g-computation: avoid predict.glm, which rebuilds model.matrix twice per call ----
## Because the exposure enters linearly, eta(A=1) = eta + beta_a * (1 - a) and
## eta(A=0) = eta - beta_a * a. No newdata construction needed.
gcomp_fast <- function(dat, outcome, exposure) {
  fo <- as.formula(paste(outcome, "~", exposure, "+", adj))
  g  <- glm(fo, data = dat, weights = nrweight, family = binomial())
  X  <- model.matrix(g)
  stopifnot(nrow(X) == nrow(dat))
  b  <- coef(g)
  a  <- dat[[exposure]]
  eta <- as.vector(X %*% b)
  ba  <- b[[exposure]]
  w   <- dat$nrweight
  p1 <- sum(w * plogis(eta + ba * (1 - a))) / sum(w)
  p0 <- sum(w * plogis(eta - ba * a)) / sum(w)
  c(rd = p1 - p0, or_cond = exp(ba))
}

## ---- school-cluster bootstrap for the conventional estimator (optional) ----
set.seed(20260924)
boot_rd <- function(outcome, exposure, B) {
  idx_by_school <- split(seq_len(nrow(m)), factor(m$schoolnum))
  ns <- length(idx_by_school)
  r <- numeric(B)
  for (b in seq_len(B)) {
    s  <- sample(idx_by_school, ns, replace = TRUE)
    ii <- unlist(s, use.names = FALSE)
    r[b] <- suppressWarnings(gcomp_fast(m[ii, , drop = FALSE], outcome, exposure)["rd"])
  }
  r
}

## verify the fast estimator reproduces the reference implementation exactly
chk  <- gcomp(m, "dep_any", "time_gt3h")
chkf <- gcomp_fast(m, "dep_any", "time_gt3h")
stopifnot(abs(chk["rd"] - chkf["rd"]) < 1e-10,
          abs(chk["or_cond"] - chkf["or_cond"]) < 1e-10)
cat("fast estimator check: PASS\n")

rows <- list()
for (en in names(exps)) {
  for (o in outs) {
    gc_  <- gcomp(m, o, en)
    cr   <- crude_rd(m, o, en)
    tr   <- tm[tm$exposure == exps[en] & tm$outcome == o, ]
    stopifnot(nrow(tr) == 1)

 ## sample-identity self-check against the TMLE run
    dev <- abs(cr - tr$rd_unadj)
    flag <- if (dev < 1e-8) "PASS" else "FAIL"

    rd_tmle <- tr$psi
    p0      <- gc_["p0"]
    or_of   <- function(rd) {
      pa <- p0 + rd
      (pa / (1 - pa)) / (p0 / (1 - p0))
    }

    se_b <- p_b <- NA_real_
    if (B > 0) {
      bs <- boot_rd(o, en, B)
      se_b <- sd(bs)
      p_b  <- 2 * pnorm(-abs(gc_["rd"]) / se_b)
    }

    rows[[length(rows) + 1]] <- data.frame(
      exposure = exps[en], outcome = o,
      rd_unadj = cr, rd_unadj_tmle_run = tr$rd_unadj, sample_check = flag,
      rd_conv_gcomp = gc_["rd"], rd_tmle = rd_tmle,
      or_cond_conv  = gc_["or_cond"], or_marg_conv = or_of(gc_["rd"]),
      or_marg_tmle  = or_of(rd_tmle),
      shrink_conv_pct = 100 * (cr - gc_["rd"]) / cr,
      shrink_tmle_pct = 100 * (cr - rd_tmle) / cr,
      extra_shrink_tmle_pp = 100 * (gc_["rd"] - rd_tmle) / cr,
      se_conv_boot = se_b, p_conv_boot = p_b,
      p_tmle_cluster = tr$p_cluster,
      verdict = if (is.na(p_b)) NA_character_ else
        if (p_b < 0.05 & tr$p_cluster >= 0.05) "FLIP: conv sig -> TMLE ns" else
        if (p_b >= 0.05 & tr$p_cluster < 0.05) "FLIP: conv ns -> TMLE sig" else
        if (p_b < 0.05 & tr$p_cluster < 0.05) "both significant" else "both ns",
      stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, rows)
res[, sapply(res, is.numeric)] <- round(res[, sapply(res, is.numeric)], 5)
write.csv(res, file.path(out, "table12_method_benchmark.csv"), row.names = FALSE)

cat("\n===== sample-identity self-check =====\n")
print(table(res$sample_check))
cat("\n===== method benchmark (risk difference, percentage points) =====\n")
print(res[, c("exposure", "outcome", "rd_unadj", "rd_conv_gcomp", "rd_tmle",
              "shrink_conv_pct", "shrink_tmle_pct")], row.names = FALSE)
cat("\n===== odds ratios =====\n")
print(res[, c("exposure", "outcome", "or_cond_conv", "or_marg_conv", "or_marg_tmle")],
      row.names = FALSE)
