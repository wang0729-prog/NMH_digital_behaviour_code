source("00_config.R")
## 15_trend_3waves.R
## Background outcome trends across the three HMS waves (2022-23 / 2023-24 / 2024-25).
## 2022-23 has NO digital-behaviour module -> descriptive only, no exposure analysis.
## Each wave is independently weight-standardised (schools equally weighted within wave),
## so prevalences are comparable across waves.

suppressPackageStartupMessages({ library(haven); library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
D <- HMS_DIR

bin <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")
con <- c("anx_score","flourish")

wprev <- function(dat, vars, w, psu) {
  des <- svydesign(ids = as.formula(paste0("~", psu)), weights = as.formula(paste0("~", w)),
                   nest = TRUE, data = dat)
  rows <- list()
  for (v in vars) {
    ok <- !is.na(dat[[v]])
    dv <- svydesign(ids = as.formula(paste0("~", psu)),
                    weights = as.formula(paste0("~", w)), nest = TRUE, data = dat[ok, ])
    mn <- svymean(as.formula(paste0("~", v)), design = dv)
    rows[[length(rows) + 1]] <- data.frame(variable = v,
                                           pct = round(100 * as.numeric(mn[1]), 2),
                                           se = round(100 * as.numeric(SE(mn)), 2),
                                           n_valid = sum(ok))
  }
  do.call(rbind, rows)
}

wmean <- function(dat, vars, w, psu) {
  rows <- list()
  for (v in vars) {
    ok <- !is.na(dat[[v]])
    dv <- svydesign(ids = as.formula(paste0("~", psu)),
                    weights = as.formula(paste0("~", w)), nest = TRUE, data = dat[ok, ])
    mn <- svymean(as.formula(paste0("~", v)), design = dv)
    rows[[length(rows) + 1]] <- data.frame(variable = v,
                                           mean = round(as.numeric(mn[1]), 3),
                                           se = round(as.numeric(SE(mn)), 3),
                                           n_valid = sum(ok))
  }
  do.call(rbind, rows)
}

## ---------- 2022-23 (build from raw dta) ----------
f <- file.path(D, "HMS_2022-2023_PUBLIC_instchars.dta")
av <- names(read_dta(f, n_max = 1))
sel <- intersect(c(bin, con, "nrweight","schoolnum","sui_plan","sui_att"), av)
h1 <- zap_labels(read_dta(f, col_select = all_of(sel)))
cat("2022-23 loaded:", nrow(h1), "x", ncol(h1), "\n")
t1 <- wprev(h1, bin, "nrweight", "schoolnum"); t1$wave <- "2022-23"
c1 <- wmean(h1, con, "nrweight", "schoolnum"); c1$wave <- "2022-23"
cat("2022-23 sui_plan valid:", sum(!is.na(h1$sui_plan)),
    "| sui_att valid:", sum(!is.na(h1$sui_att)), "(conditional sub-sample)\n")

## ---------- 2023-24 / 2024-25 ----------
h2 <- readRDS(file.path(out, "hms_2324_analysis.rds"))
h3 <- readRDS(file.path(out, "hms_2425_analysis.rds"))
t2 <- wprev(h2, bin, "nrweight", "schoolnum"); t2$wave <- "2023-24"
t3 <- wprev(h3, bin, "nrweight", "schoolnum"); t3$wave <- "2024-25"
c2 <- wmean(h2, con, "nrweight", "schoolnum"); c2$wave <- "2023-24"
c3 <- wmean(h3, con, "nrweight", "schoolnum"); c3$wave <- "2024-25"

tab <- rbind(t1, t2, t3)
cat("\n=== weighted prevalence (%) by wave ===\n"); print(tab, row.names = FALSE)
cat("\n=== weighted means by wave ===\n"); print(rbind(c1, c2, c3), row.names = FALSE)

nsch <- data.frame(wave = c("2022-23","2023-24","2024-25"),
                   n = c(nrow(h1), nrow(h2), nrow(h3)),
                   schools = c(length(unique(h1$schoolnum)),
                               length(unique(h2$schoolnum)),
                               length(unique(h3$schoolnum))),
                   module_n = c(NA, sum(h2$module == 1), sum(h3$module == 1)))
cat("\n=== sample sizes ===\n"); print(nsch, row.names = FALSE)

## ---- pairwise Wald tests (waves are independent samples of students) ----
cmp <- function(tab, v, w1, w2) {
  a <- tab[tab$variable == v & tab$wave == w1, ]
  b <- tab[tab$variable == v & tab$wave == w2, ]
  diff <- a$pct - b$pct; se <- sqrt(a$se^2 + b$se^2)
  data.frame(variable = v, comparison = paste0(w1, " minus ", w2),
             diff_pp = round(diff, 2), se = round(se, 2),
             z = round(diff / se, 2), p = signif(2 * pnorm(-abs(diff / se)), 3))
}
tests <- list()
for (v in bin) {
  tests[[length(tests) + 1]] <- cmp(tab, v, "2022-23", "2024-25")
  tests[[length(tests) + 1]] <- cmp(tab, v, "2022-23", "2023-24")
  tests[[length(tests) + 1]] <- cmp(tab, v, "2023-24", "2024-25")
}
tests <- do.call(rbind, tests)
cat("\n=== pairwise differences (percentage points) ===\n"); print(tests, row.names = FALSE)

cmp_m <- function(cm, v, w1, w2) {
  a <- cm[cm$variable == v & cm$wave == w1, ]
  b <- cm[cm$variable == v & cm$wave == w2, ]
  diff <- a$mean - b$mean; se <- sqrt(a$se^2 + b$se^2)
  data.frame(variable = v, comparison = paste0(w1, " minus ", w2),
             diff_pp = round(diff, 3), se = round(se, 3),
             z = round(diff / se, 2), p = signif(2 * pnorm(-abs(diff / se)), 3))
}
cm <- rbind(c1, c2, c3)
tm <- do.call(rbind, lapply(con, function(v) cmp_m(cm, v, "2022-23", "2024-25")))
cat("\n=== continuous outcomes, 2022-23 vs 2024-25 ===\n"); print(tm, row.names = FALSE)

write.csv(rbind(tests, tm), file.path(out, "table8_trend_tests.csv"), row.names = FALSE)
write.csv(tab, file.path(out, "table8_trend_prevalence.csv"), row.names = FALSE)
write.csv(rbind(c1, c2, c3), file.path(out, "table8_trend_means.csv"), row.names = FALSE)
write.csv(nsch, file.path(out, "table8_trend_samples.csv"), row.names = FALSE)

saveRDS(h1[, intersect(c(bin, con, "nrweight","schoolnum"), names(h1))],
        file.path(out, "hms_2223_outcomes.rds"))
cat("DONE\n")
