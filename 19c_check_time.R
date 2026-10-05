source("00_config.R")
## 19c_check_time.R
## Is the screen-time item comparable across waves? Compare the raw internet_1
## distribution in the SAME 132 overlapping schools.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(haven) })

D <- HMS_DIR
p <- OUT

ra <- read_dta(file.path(D, "HMS_2023-2024_PUBLIC_instchars.dta"),
               col_select = all_of(c("internet_1", "nrweight", "schoolnum")))
rb <- read_dta(file.path(D, "HMS_2024-2025_PUBLIC_instchars.dta"),
               col_select = all_of(c("internet_1", "nrweight", "schoolnum")))

a <- readRDS(file.path(p, "hms_2324_analysis.rds"))
b <- readRDS(file.path(p, "hms_2425_analysis.rds"))
sa <- unique(a$schoolnum[!is.na(a$time_gt3h)])
sb <- unique(b$schoolnum[!is.na(b$time_gt3h)])
sch <- intersect(sa, sb)

ra <- ra[ra$schoolnum %in% sch & !is.na(ra$internet_1), ]
rb <- rb[rb$schoolnum %in% sch & !is.na(rb$internet_1), ]

cat("overlapping schools:", length(sch), "\n")
cat("raw internet_1 valid n: 2324 =", nrow(ra), " 2425 =", nrow(rb), "\n\n")

fa <- table(as.numeric(ra$internet_1)); fb <- table(as.numeric(rb$internet_1))
cat("--- 2023-2024 internet_1 distribution (7 bins) ---\n")
print(round(100 * fa / sum(fa), 2))
cat("\n--- 2024-2025 internet_1 distribution (8 bins) ---\n")
print(round(100 * fb / sum(fb), 2))

cat("\n--- share above 3 h/day under each scale definition ---\n")
cat("2324 (internet_1 == 7):", round(100 * mean(as.numeric(ra$internet_1) == 7), 2), "%\n")
cat("2425 (internet_1 >= 5):", round(100 * mean(as.numeric(rb$internet_1) >= 5), 2), "%\n")

cat("\n--- 2024-2025 bins 5-8 broken out (3-4h / 4-6h / 6-9h / 10+h) ---\n")
print(round(100 * fb[as.character(5:8)] / sum(fb), 2))
