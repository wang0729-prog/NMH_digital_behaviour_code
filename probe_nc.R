source("00_config.R")
suppressPackageStartupMessages({ library(haven) })
f23 <- file.path(HMS_DIR, "HMS_2023-2024_PUBLIC_instchars.dta")
v <- c("weight","height_ft","height_in","height_ft_clean","height_in_clean","bmi",
       "inst_gradrate","flourish")
sub <- read_dta(f23, col_select = all_of(intersect(v, names(read_dta(f23, n_max = 1)))))
for (nm in names(sub)) {
  x <- suppressWarnings(as.numeric(sub[[nm]]))
  cat(sprintf("%-18s valid=%7d  min=%8.2f  med=%8.2f  max=%8.2f\n",
              nm, sum(!is.na(x)), min(x, na.rm = TRUE),
              median(x, na.rm = TRUE), max(x, na.rm = TRUE)))
}
