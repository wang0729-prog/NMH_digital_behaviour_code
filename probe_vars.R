source("00_config.R")
suppressPackageStartupMessages({ library(haven) })
f22 <- file.path(HMS_DIR, "HMS_2022-2023_PUBLIC_instchars.dta")
allv <- names(read_dta(f22, n_max = 1))
cat("2022-23 total columns:", length(allv), "\n")
want <- c("dep_any","dep_maj","anx_any","anx_score","sui_idea","sui_plan","sui_att",
          "lonely","flourish","nrweight","schoolnum","age","sex_birth","weight",
          "fincur","finpast","food_worry","educ_par1","educ_par2","race","international",
          "undergrad","grad_student","year_in_school","sex")
for (v in want) cat(sprintf("%-16s %s\n", v, v %in% allv))
cat("--- candidates matching weight/height ---\n")
print(grep("weight|height|bmi", allv, value = TRUE, ignore.case = TRUE)[1:20])
cat("--- 2023-24 weight distribution (from raw dta) ---\n")
f23 <- file.path(HMS_DIR, "HMS_2023-2024_PUBLIC_instchars.dta")
if ("weight" %in% names(read_dta(f23, n_max = 1))) {
  w <- read_dta(f23, col_select = "weight")[[1]]
  print(summary(w)); print(table(w >= 999))
} else cat("no weight var in 2023-24\n")
