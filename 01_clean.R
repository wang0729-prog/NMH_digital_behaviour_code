source("00_config.R")
if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({library(haven); library(dplyr)})

## OUT is set in 00_config.R
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(OUT, "01_clean_log.txt"), open = "w", encoding = "UTF-8")
W <- function(...) writeLines(paste0(...), log_con)

## ============================================================
## HMS helper
## ============================================================
build_hms <- function(year) {
  f <- file.path(HMS_DIR, paste0("HMS_", year, "_PUBLIC_instchars.dta"))
  allc <- names(read_dta(f, n_max = 1))
  want <- c("internet_1", paste0("internet_2_", 1:14), paste0("internet_4_", 1:7),
            "dep_any", "dep_maj", "anx_any", "anx_score", "sui_idea", "sui_plan", "sui_att",
            "lonely", "flourish", "fincur", "finpast", "food_worry", "food_notlast",
            "educ_par1", "educ_par2", "age", "sex_birth",
            "gender_male", "gender_female", "gender_trans", "gender_nonbin",
            paste0("race_", c("white", "black", "asian", "his", "pi", "ainaan", "mides", "other")),
            "international", "degree_bach", "degree_ma", "degree_phd", "transfer",
            "nrweight", "schoolnum", "inst_type", "inst_size", "inst_public", "inst_geo",
            "inst_gradrate", "inst_msi")
  sel <- intersect(want, allc)
  d <- read_dta(f, col_select = all_of(sel))
  W("\n===== HMS ", year, " loaded: ", nrow(d), " x ", ncol(d), " =====")

 ## --- exposure ---
  d$module <- as.integer(!is.na(d$internet_1))
  mp <- if (year == "2023-2024") c("1"=0,"2"=1,"3"=1,"4"=1,"5"=2,"6"=3,"7"=4) else c("1"=0,"2"=1,"3"=2,"4"=3,"5"=4,"6"=4,"7"=4,"8"=4)
  d$time5cat <- as.integer(mp[as.character(d$internet_1)])
  d$time_gt3h <- as.integer(!is.na(d$internet_1) & d$internet_1 == if (year == "2023-2024") 7 else 5)
  if (year == "2024-2025") d$time_gt3h <- as.integer(!is.na(d$internet_1) & d$internet_1 >= 5)

 ## --- online experiences (2023-24: 5 items; multiselect: selected=1, unselected=NA -> 0, ONLY among internet_1 != 1) ---
  elig_exp <- !is.na(d$internet_1) & d$internet_1 != 1
  if (year == "2023-2024") {
    d$exp_relation  <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_1)))
    d$exp_harassed  <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_2)))
    d$exp_create    <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_3)))
    d$exp_activism  <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_4)))
    d$exp_none      <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_5)))
  } else {
    d$exp_harassed  <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_2)))
    d$exp_neg_none  <- ifelse(!elig_exp, NA, as.integer(!is.na(d$internet_2_1)))
  }

 ## --- generative AI (2023-24 only) ---
  if (year == "2023-2024") {
    ai_cols <- paste0("internet_4_", 1:7)
    d$ai_academic  <- as.integer(!is.na(d$internet_4_1))
    d$ai_forbidden <- as.integer(!is.na(d$internet_4_2))
    d$ai_work      <- as.integer(!is.na(d$internet_4_3))
    d$ai_comm      <- as.integer(!is.na(d$internet_4_4))
    d$ai_health    <- as.integer(!is.na(d$internet_4_5))
    d$ai_fun       <- as.integer(!is.na(d$internet_4_6))
    d$ai_other     <- as.integer(!is.na(d$internet_4_7))
    M <- as.matrix(d[, ai_cols])
    nsel <- rowSums(!is.na(M))
    d$ai_any    <- ifelse(d$module == 0, NA, as.integer(nsel > 0))
    d$ai_n_uses <- ifelse(d$module == 0, NA, nsel)
    W("AI: any = ", sum(d$ai_any == 1, na.rm = TRUE), " | mean n_uses among users = ",
      round(mean(d$ai_n_uses[d$ai_any == 1], na.rm = TRUE), 2))
  }

 ## --- SDOH ---
  d$fin_stress  <- as.integer(!is.na(d$fincur) & d$fincur <= 2)
  d$food_insec  <- as.integer(!is.na(d$food_worry) & d$food_worry >= 2)
  e1 <- suppressWarnings(as.numeric(d$educ_par1)); e2 <- suppressWarnings(as.numeric(d$educ_par2))
  e1[e1 == 8] <- NA; e2[e2 == 8] <- NA
  both_known <- !is.na(e1) & !is.na(e2)
  d$firstgen <- ifelse(!both_known, NA, as.integer(pmax(e1, e2, na.rm = TRUE) < 6))
  d$parent_edu_max <- pmax(e1, e2, na.rm = FALSE)
  d$parent_edu_max[is.infinite(d$parent_edu_max)] <- NA

 ## --- covariates ---
  d$female <- as.integer(d$sex_birth == 1)
  d$male   <- as.integer(d$sex_birth == 2)
  d$age_num <- suppressWarnings(as.numeric(d$age))
  d$age_gt40 <- as.integer(!is.na(d$age_num) & d$age_num >= 40)

  rcols <- paste0("race_", c("white", "black", "asian", "his", "pi", "ainaan", "mides", "other"))
  R <- as.matrix(d[, rcols])
  nsel <- rowSums(!is.na(R))
  d$race_cat <- ifelse(nsel == 0, NA,
                ifelse(nsel > 1, "Multiracial",
                ifelse(!is.na(d$race_white), "White",
                ifelse(!is.na(d$race_asian), "APIDA",
                ifelse(!is.na(d$race_his), "Hispanic",
                ifelse(!is.na(d$race_black), "Black",
                ifelse(!is.na(d$race_ainaan), "AI/AN",
                ifelse(!is.na(d$race_mides), "MENA", "Other"))))))))
  d$international <- as.integer(d$international == 1)
  d$undergrad <- as.integer(!is.na(d$degree_bach) & is.na(d$degree_ma) & is.na(d$degree_phd))
  d$grad_student <- as.integer(!is.na(d$degree_ma) | !is.na(d$degree_phd))

  d$school_wave <- paste(year, d$schoolnum, sep = "_")

  keep <- c("module", "time5cat", "time_gt3h",
            grep("^exp_", names(d), value = TRUE), grep("^ai_", names(d), value = TRUE),
            "dep_any", "dep_maj", "anx_any", "anx_score", "sui_idea", "sui_plan", "sui_att",
            "lonely", "flourish", "fin_stress", "food_insec", "firstgen", "parent_edu_max",
            "fincur", "finpast", "food_worry", "educ_par1", "educ_par2",
            "age_num", "age_gt40", "female", "male", "race_cat", "international",
            "undergrad", "grad_student", "nrweight", "schoolnum", "school_wave",
            "inst_type", "inst_size", "inst_public", "inst_geo", "inst_gradrate", "inst_msi")
  out <- d[, intersect(keep, names(d))]
  out$wave <- year
  W("analysis frame: ", nrow(out), " x ", ncol(out), " | module n = ", sum(out$module))
  out
}

h2324 <- build_hms("2023-2024")
h2425 <- build_hms("2024-2025")
saveRDS(h2324, file.path(OUT, "hms_2324_analysis.rds"))
saveRDS(h2425, file.path(OUT, "hms_2425_analysis.rds"))
write.csv(h2324[h2324$module == 1, ], file.path(OUT, "hms_2324_module.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ============================================================
## CFPS
## ============================================================
cln <- function(x) { x <- suppressWarnings(as.numeric(x)); x[x < 0] <- NA; x }
build_cfps <- function(yr) {
  f <- file.path(CFPS_DIR, paste0("cfps", yr, "person.dta"))
  allc <- names(read_dta(f, n_max = 1))
  edu <- intersect(paste0("cfps", yr, "edu"), allc)
  want <- c("pid", "age", "gender", "qc1", "qc3", "qu93", "qu931", "qu201", "qu201a", "qu11", "qu111",
            "qu91", "cesd8", "qn414", "qm2016", "urban20", "urban22", "rswt_natcs20n", "rswt_natcs22n", edu)
  sel <- intersect(want, allc)
  d <- read_dta(f, col_select = all_of(sel))
  W("\n===== CFPS ", yr, " loaded: ", nrow(d), " x ", ncol(d), " =====")

  d$age_num  <- cln(d$age)
  d$female   <- as.integer(cln(d$gender) == 1)
  d$urban    <- cln(d[[intersect(c("urban20", "urban22"), sel)[1]]])
  d$sv_watch <- as.integer(cln(d$qu93) == 1)
  d$sv_watch[is.na(cln(d$qu93))] <- NA
  d$sv_daily <- as.integer(cln(d$qu931) == 1)
  d$sv_daily[is.na(cln(d$qu931))] <- NA
  d$mobile_hours <- cln(d$qu201a)
  d$wechat    <- as.integer(cln(d$qu11) == 1)
  d$cesd8_std <- cln(d$cesd8) - 8
  d$dep8      <- as.integer(d$cesd8_std >= 8)
  d$dep8[is.na(d$cesd8_std)] <- NA
  d$dep10     <- as.integer(d$cesd8_std >= 10)
  d$dep10[is.na(d$cesd8_std)] <- NA
  d$lonely_score <- cln(d$qn414)
  d$lonely    <- as.integer(d$lonely_score >= 3)
  d$lonely[is.na(d$lonely_score)] <- NA
  wv <- intersect(c("rswt_natcs20n", "rswt_natcs22n"), sel)
  d$wt <- cln(d[[wv[1]]])
  if (length(edu) == 1) d$edu_max <- cln(d[[edu]]) else d$edu_max <- NA_real_

  d$college_strict <- as.integer(!is.na(d$qc1) & d$qc1 == 1 & !is.na(d$qc3) & d$qc3 %in% 6:9)
  d$young_hiedu <- as.integer(!is.na(d$age_num) & d$age_num >= 18 & d$age_num <= 30 &
                              !is.na(d$edu_max) & d$edu_max >= 5)
  W("college_strict = ", sum(d$college_strict), " | young_hiedu = ", sum(d$young_hiedu))
  W("young_hiedu: sv_watch valid = ", sum(!is.na(d$sv_watch[d$young_hiedu == 1])),
    " | cesd8 valid = ", sum(!is.na(d$cesd8_std[d$young_hiedu == 1])))

  keep <- c("pid", "age_num", "female", "urban", "sv_watch", "sv_daily", "mobile_hours",
            "wechat", "cesd8_std", "dep8", "dep10", "lonely", "lonely_score", "wt",
            "edu_max", "college_strict", "young_hiedu")
  out <- d[, keep]
  out$wave <- yr
  out
}

c20 <- build_cfps("2020")
c22 <- build_cfps("2022")
saveRDS(c20, file.path(OUT, "cfps_2020_analysis.rds"))
saveRDS(c22, file.path(OUT, "cfps_2022_analysis.rds"))
write.csv(c22[c22$young_hiedu == 1, ], file.path(OUT, "cfps_2022_young_hiedu.csv"), row.names = FALSE, fileEncoding = "UTF-8")

close(log_con)
cat("DONE\n")
