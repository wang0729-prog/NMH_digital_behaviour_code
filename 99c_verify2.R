source("00_config.R")
## 99c_verify2.R -- independent recomputation of the new sensitivity / trend results
suppressPackageStartupMessages({ library(survey); library(haven) })
dir <- ANALYSIS_DIR; out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

## (1) control outcomes: UNADJUSTED weighted difference as a second code path
m <- d[d$module == 1 & !is.na(d$time5cat), ]
m$hi <- as.integer(m$time5cat == 4)
des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = m)
for (y in c("flourish","inst_gradrate")) {
  s <- m[!is.na(m[[y]]) & !is.na(m$hi), ]
  ds <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
  tb <- svyby(as.formula(paste0("~", y)), by = ~hi, design = ds, FUN = svymean)
  cat(sprintf("CHECK %-14s unadjusted diff (%s): %.3f   [adjusted beta in table9]\n",
              y, ifelse(y == "inst_gradrate", "pp x100", "score"),
              ifelse(y == "inst_gradrate",
                     100 * (tb[[2]][2] - tb[[2]][1]), tb[[2]][2] - tb[[2]][1])))
}

## (2) E-value spot check by hand for screen>3h / dep_any (P0 = 29.4%, OR = 2.479)
p0 <- 0.294; or <- 2.479
rr <- or / (1 - p0 + p0 * or)
cat("CHECK E-value hand calc: RR =", round(rr, 3),
    " E =", round(rr + sqrt(rr * (rr - 1)), 2), " (table10: 1.728 / 2.85)\n")

## (3) trend: unweighted prevalences as a directional sanity check
f <- file.path(HMS_DIR, "HMS_2022-2023_PUBLIC_instchars.dta")
h1 <- zap_labels(read_dta(f, col_select = all_of(c("dep_any","sui_idea","lonely"))))
h2 <- readRDS(file.path(out, "hms_2324_analysis.rds"))
h3 <- readRDS(file.path(out, "hms_2425_analysis.rds"))
for (v in c("dep_any","sui_idea","lonely")) {
  u1 <- mean(h1[[v]], na.rm = TRUE); u2 <- mean(h2[[v]], na.rm = TRUE)
  u3 <- mean(h3[[v]], na.rm = TRUE)
  cat(sprintf("CHECK %-9s unweighted %%: 22-23 %.1f | 23-24 %.1f | 24-25 %.1f\n",
              v, 100 * u1, 100 * u2, 100 * u3))
}

## (4) BH boundary check: recompute q for the two tests straddling q = .05
tb <- read.csv(file.path(out, "table11_bh_corrected.csv"))
tb <- tb[order(tb$p), ]
n <- nrow(tb)
tb$q_manual <- pmin(1, cummin(rev(pmin(1, rev(tb$p * n / seq(n, 1))))))
i <- which(diff(tb$sig_BH_overall) != 0)
if (length(i)) {
  j <- max(i)
  cat("CHECK BH boundary rows", j, "-", j + 1, ": p =", signif(tb$p[j], 3),
      "/", signif(tb$p[j + 1], 3), " q =", round(tb$q_overall[j], 4),
      "/", round(tb$q_overall[j + 1], 4), " manual q =",
      round(tb$q_manual[j], 4), "/", round(tb$q_manual[j + 1], 4), "\n")
}
cat("CHECK n tests:", n, "| raw<0.05:", sum(tb$p < 0.05),
    "| BH<0.05:", sum(tb$q_overall < 0.05), "\n")

## (5) MI vs CC: verify the two pooled ORs differ by <2%
mi <- read.csv(file.path(out, "table_sens_mi_vs_cc.csv"))
cat("CHECK max |log OR_MI - log OR_CC| =",
    round(max(abs(log(mi$OR_MI) - log(mi$OR_CC))), 4), "\n")
cat("VERIFY2 DONE\n")
