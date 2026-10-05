source("00_config.R")
## 20_school_fe.R
## School fixed effects at the INDIVIDUAL level, 2023-24 only.
##
## Why this design instead of the two-period school panel (19_school_panel.R):
## the screen-time item is NOT measurement-equivalent across waves (see 19c): within the
## same 132 schools, >3h/day rose from 30.64% to 54.44% while the 1-2h and 2-3h bins both
## FELL -- a pattern that cannot be a real behaviour change and reflects response-category
## anchoring (2024-25 splits >3h into 4 bins). The harassment item also changed wording.
## A cross-wave panel is therefore not defensible for either exposure.
##
## What this does instead: absorb ALL school-level confounding -- observed and unobserved,
## time-invariant -- by including school fixed effects in the individual-level model,
## using one wave only. Identification comes from within-school contrasts.
##
## Estimator: weighted linear probability model with the within (school-demeaned)
## transformation, school-cluster robust SEs. Reported against the same model WITHOUT
## school FE, so the confounding shift is explicit.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

p <- OUT
d <- readRDS(file.path(p, "hms_2324_analysis.rds"))

exps <- c("time_gt3h", "exp_harassed", "ai_any")
outs <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
covs <- c("age_num", "female", "international", "undergrad", "grad_student",
          "fin_stress", "food_insec", "firstgen")
races <- c("Multiracial", "APIDA", "Hispanic", "Black", "AI/AN", "MENA", "Other")

m <- d[d$module == 1, ]
for (v in c(exps, outs, covs, "race_cat", "nrweight", "schoolnum")) m <- m[!is.na(m[[v]]), ]
rownames(m) <- NULL
m$race_cat <- factor(m$race_cat)
m$race_cat <- relevel(m$race_cat, ref = "White")
m$schoolnum <- factor(m$schoolnum)

## drop schools with fewer than 10 respondents: they contribute nothing to within-school
## identification and destabilise the demeaned fit
tab <- table(m$schoolnum)
keep <- names(tab)[tab >= 10]
cat("schools before:", length(tab), " >=10 respondents:", length(keep), "\n")
m <- m[m$schoolnum %in% keep, ]
m$schoolnum <- droplevels(m$schoolnum)
cat("individual n:", nrow(m), " schools:", nlevels(m$schoolnum), "\n")

Xextra <- model.matrix(~ race_cat, data = m)[, -1, drop = FALSE]
cov_mat <- cbind(as.matrix(m[, covs]), Xextra)

## within transformation, weighted by nrweight
demean <- function(v, g, w) {
  if (is.matrix(v)) {
    out <- v
    for (j in seq_len(ncol(v))) out[, j] <- demean(v[, j], g, w)
    return(out)
  }
  wm <- tapply(w * v, g, sum) / tapply(w, g, sum)
  v - wm[as.character(g)]
}

rows <- list()
for (e in exps) {
  for (o in outs) {
    X <- cbind(m[[e]], cov_mat)
    colnames(X)[1] <- e
    y <- as.numeric(m[[o]])
    g <- m$schoolnum
    w <- m$nrweight

 ## (A) no school FE
    fa <- lm(y ~ X, weights = w)
    ca <- coeftest(fa, vcov = vcovCL(fa, cluster = ~g))[2, ]

 ## (B) school FE via within transformation
    Xt <- cbind(demean(X, g, w), 1)          # intercept absorbs the common level
    yt <- demean(y, g, w)
    fb <- lm(yt ~ Xt - 1, weights = w)
    cb <- coeftest(fb, vcov = vcovCL(fb, cluster = ~g))[1, ]

    rows[[length(rows) + 1]] <- data.frame(
      exposure = e, outcome = o,
      n = nrow(m), n_schools = nlevels(g),
      rd_noFE = ca[1], se_noFE = ca[2], p_noFE = ca[4],
      rd_FE = cb[1], se_FE = cb[2], p_FE = cb[4],
      confounding_shift = ca[1] - cb[1],
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
## round the estimates first, then recompute the shift from the rounded columns so the
## written CSV is internally consistent (otherwise shift differs from rd_noFE - rd_FE by
## up to 1e-4 and fails arithmetic checks).
res[, 4:(ncol(res) - 1)] <- round(res[, 4:(ncol(res) - 1)], 4)
res$confounding_shift <- round(res$rd_noFE - res$rd_FE, 4)
write.csv(res, file.path(p, "table16_school_FE.csv"), row.names = FALSE)

cat("\n===== individual-level school fixed effects (2023-24) =====\n")
print(res[, c("exposure", "outcome", "rd_noFE", "p_noFE", "rd_FE", "p_FE",
              "confounding_shift")], row.names = FALSE)
cat("\nDONE\n")
