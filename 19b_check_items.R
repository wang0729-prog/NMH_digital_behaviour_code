source("00_config.R")
## 19b_check_items.R
## Which exposure item is genuinely comparable across 2023-24 and 2024-25?
## Print full value labels for internet_1 (screen time) and all internet_2_* items.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(haven) })

D <- HMS_DIR

vl <- function(year, v) {
  f <- file.path(D, paste0("HMS_", year, "_PUBLIC_instchars.dta"))
  allc <- names(read_dta(f, n_max = 1))
  if (!v %in% allc) return(NULL)
  d <- read_dta(f, n_max = 1, col_select = all_of(v))
  attr(d[[v]], "labels")
}

show <- function(year, v) {
  L <- vl(year, v)
  cat("\n[", year, "] ", v, ":\n", sep = "")
  if (is.null(L)) { cat("   <absent>\n"); return(invisible(NULL)) }
  for (i in seq_along(L)) cat(sprintf("   %3s = %s\n", names(L)[i], as.character(L[[i]])))
}

for (y in c("2023-2024", "2024-2025")) show(y, "internet_1")

cat("\n\n############ internet_2_* items ############")
for (y in c("2023-2024", "2024-2025")) {
  cat("\n=================", y, "=================\n")
  for (i in 1:14) show(y, paste0("internet_2_", i))
}

## ---- outcome distributions as a cross-wave sanity check ----
p <- OUT
a <- readRDS(file.path(p, "hms_2324_analysis.rds"))
b <- readRDS(file.path(p, "hms_2425_analysis.rds"))
cat("\n\n############ outcome prevalence by wave (unweighted) ############\n")
for (o in c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")) {
  cat(sprintf("%-9s  2324 = %.4f (n=%d)   2425 = %.4f (n=%d)\n", o,
              mean(a[[o]], na.rm = TRUE), sum(!is.na(a[[o]])),
              mean(b[[o]], na.rm = TRUE), sum(!is.na(b[[o]]))))
}
