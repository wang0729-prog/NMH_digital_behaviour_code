source("00_config.R")
## 19a_check_compat.R
## Pre-flight for the school-level two-period fixed-effects panel.
## Verify that (a) the exposure item carries the same wording in both waves,
## (b) the outcomes are defined identically, (c) enough schools overlap.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(haven) })

D <- HMS_DIR

lab <- function(year, v) {
  f <- file.path(D, paste0("HMS_", year, "_PUBLIC_instchars.dta"))
  allc <- names(read_dta(f, n_max = 1))
  v <- intersect(v, allc)
  if (length(v) == 0) return(setNames(rep("<absent>", length(v)), v))
  d <- read_dta(f, n_max = 1, col_select = all_of(v))
  out <- vapply(d, function(x) {
    l <- attr(x, "label"); if (is.null(l) || is.na(l)) "<none>" else as.character(l)
  }, character(1))
  unname(out)
}

vs <- c(paste0("internet_2_", 1:5), "internet_1",
        "dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")

cat("========== variable labels: 2023-2024 vs 2024-2025 ==========\n")
for (v in vs) {
  l1 <- lab("2023-2024", v); l2 <- lab("2024-2025", v)
  cat("\n---", v, "---\n")
  cat("  2324:", if (!length(l1)) "<absent>" else l1, "\n")
  cat("  2425:", if (!length(l2)) "<absent>" else l2, "\n")
  cat("  IDENTICAL WORDING:", identical(l1, l2), "\n")
}

## ---- value labels for internet_2_2 (the harassment item) ----
cat("\n========== value labels of internet_2_2 ==========\n")
for (y in c("2023-2024", "2024-2025")) {
  f <- file.path(D, paste0("HMS_", y, "_PUBLIC_instchars.dta"))
  d <- read_dta(f, n_max = 1, col_select = all_of("internet_2_2"))
  cat(y, ":\n"); print(attr(d$internet_2_2, "labels"))
}

## ---- overlap and cell sizes ----
p <- OUT
a <- readRDS(file.path(p, "hms_2324_analysis.rds"))
b <- readRDS(file.path(p, "hms_2425_analysis.rds"))

cat("\n========== wave-level exposure availability ==========\n")
for (nm in c("2324", "2425")) {
  d <- if (nm == "2324") a else b
  dd <- d[!is.na(d$exp_harassed), ]
  cat(sprintf("%s: exp_harassed valid n=%d  schools=%d  mean=%.4f\n",
              nm, nrow(dd), length(unique(dd$schoolnum)),
              mean(dd$exp_harassed, na.rm = TRUE)))
}

sa <- unique(a$schoolnum[!is.na(a$exp_harassed)])
sb <- unique(b$schoolnum[!is.na(b$exp_harassed)])
cat("\n====================== overlap ======================\n")
cat("schools with exposure, 2324:", length(sa), "\n")
cat("schools with exposure, 2425:", length(sb), "\n")
cat("OVERLAP:", length(intersect(sa, sb)), "\n")
