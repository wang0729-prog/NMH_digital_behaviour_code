source("00_config.R")
## 99e_verify4.R
## Verification for the school-confounding work:
## (1) the within transformation is arithmetically correct (school-demeaned variables
## must have zero weighted sum within every school);
## (2) Table 16 numbers re-derive from the written CSV;
## (3) the cross-wave screen-time incomparability evidence in 19c is reproducible.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(haven) })

dir <- OUT
chk <- function(label, got, want, tol = 1e-4) {
  ok <- isTRUE(all.equal(got, want, tolerance = tol))
  cat(sprintf("[%s] %-54s got=%s want=%s\n", if (ok) "PASS" else "FAIL",
              label, format(got, digits = 5), format(want, digits = 5)))
  ok
}
res <- c()

## ---------- (1) within-transformation identity ----------
d <- readRDS(file.path(dir, "hms_2324_analysis.rds"))
m <- d[d$module == 1 & !is.na(d$time_gt3h) & !is.na(d$dep_any) &
         !is.na(d$nrweight) & !is.na(d$schoolnum), ]
rownames(m) <- NULL
tab <- table(m$schoolnum); m <- m[m$schoolnum %in% names(tab)[tab >= 10], ]
m$schoolnum <- droplevels(factor(m$schoolnum))

demean <- function(v, g, w) {
  wm <- tapply(w * v, g, sum) / tapply(w, g, sum)
  v - wm[as.character(g)]
}
w <- m$nrweight; g <- m$schoolnum
for (v in c("time_gt3h", "dep_any")) {
  dv <- demean(m[[v]], g, w)
  s <- tapply(w * dv, g, sum)
  res <- c(res, chk(paste("within-transform: max |weighted school sum| for", v),
                    max(abs(s)), 0, 1e-6))
}

## ---------- (2) Table 16 ----------
t16 <- read.csv(file.path(dir, "table16_school_FE.csv"))
res <- c(res, chk("T16 n rows", nrow(t16), 15))
res <- c(res, chk("T16 individual n", unique(t16$n), 40732))
res <- c(res, chk("T16 n schools", unique(t16$n_schools), 190))
res <- c(res, chk("T16 confounding_shift == rd_noFE - rd_FE",
                  max(abs(t16$confounding_shift - (t16$rd_noFE - t16$rd_FE))), 0, 1e-9))
res <- c(res, chk("T16 max |confounding_shift| across 15 models",
                  max(abs(t16$confounding_shift)), 0.0099, 1e-4))
res <- c(res, chk("T16 all school-FE estimates keep the same sign as no-FE",
                  all(sign(t16$rd_FE) == sign(t16$rd_noFE)), TRUE))
tg <- t16[t16$exposure == "time_gt3h", ]
res <- c(res, chk("T16 time_gt3h dep_any rd_noFE -> rd_FE",
                  c(tg$rd_noFE[tg$outcome == "dep_any"], tg$rd_FE[tg$outcome == "dep_any"]),
                  c(0.1470, 0.1403), 1e-4))
eh <- t16[t16$exposure == "exp_harassed", ]
res <- c(res, chk("T16 exp_harassed dep_any rd_noFE -> rd_FE",
                  c(eh$rd_noFE[eh$outcome == "dep_any"], eh$rd_FE[eh$outcome == "dep_any"]),
                  c(0.1823, 0.1724), 1e-4))
ai <- t16[t16$exposure == "ai_any", ]
res <- c(res, chk("T16 ai_any: dep_any + sui_idea significant under school FE",
                  sum(ai$p_FE < 0.05), 2))
res <- c(res, chk("T16 ai_any sui_idea crosses significance under FE",
                  (ai$p_noFE[ai$outcome == "sui_idea"] >= 0.05) &&
                    (ai$p_FE[ai$outcome == "sui_idea"] < 0.05), TRUE))

## ---------- (3) cross-wave screen-time incomparability ----------
D <- HMS_DIR
ra <- read_dta(file.path(D, "HMS_2023-2024_PUBLIC_instchars.dta"),
               col_select = all_of(c("internet_1", "schoolnum")))
rb <- read_dta(file.path(D, "HMS_2024-2025_PUBLIC_instchars.dta"),
               col_select = all_of(c("internet_1", "schoolnum")))
a <- readRDS(file.path(dir, "hms_2324_analysis.rds"))
b <- readRDS(file.path(dir, "hms_2425_analysis.rds"))
sch <- intersect(unique(a$schoolnum[!is.na(a$time_gt3h)]),
                 unique(b$schoolnum[!is.na(b$time_gt3h)]))
ra <- ra[ra$schoolnum %in% sch & !is.na(ra$internet_1), ]
rb <- rb[rb$schoolnum %in% sch & !is.na(rb$internet_1), ]
va <- as.numeric(ra$internet_1); vb <- as.numeric(rb$internet_1)
res <- c(res, chk("19c overlapping schools", length(sch), 134))
res <- c(res, chk("19c 2324 share >3h (bin 7)", round(100 * mean(va == 7), 2), 30.64, 0.01))
res <- c(res, chk("19c 2425 share >3h (bins 5-8)", round(100 * mean(vb >= 5), 2), 54.44, 0.01))
res <- c(res, chk("19c 2324 bin 5 (1-2h)", round(100 * mean(va == 5), 2), 23.06, 0.01))
res <- c(res, chk("19c 2425 bin 3 (1-2h)", round(100 * mean(vb == 3), 2), 16.63, 0.01))
res <- c(res, chk("19c 2324 bin 6 (2-3h)", round(100 * mean(va == 6), 2), 22.80, 0.01))
res <- c(res, chk("19c 2425 bin 4 (2-3h)", round(100 * mean(vb == 4), 2), 19.48, 0.01))
## the internal contradiction: >3h rises while both adjacent bins fall
res <- c(res, chk("19c contradiction: >3h rises AND 1-2h falls AND 2-3h falls",
                  (mean(vb >= 5) > mean(va == 7)) &&
                    (mean(vb == 3) < mean(va == 5)) &&
                    (mean(vb == 4) < mean(va == 6)), TRUE))

cat(sprintf("\n===== %d checks: %d PASS, %d FAIL =====\n",
            length(res), sum(res), sum(!res)))
