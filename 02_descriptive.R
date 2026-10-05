source("00_config.R")
if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({library(haven); library(dplyr); library(survey)})

## OUT is set in 00_config.R
log_con <- file(file.path(OUT, "02_desc_log.txt"), open = "w", encoding = "UTF-8")
W <- function(...) writeLines(paste0(...), log_con)

d <- readRDS(file.path(OUT, "hms_2324_analysis.rds"))

## ---------- helpers ----------
wm <- function(des, v) {  # weighted mean for binary/continuous
  m <- svymean(as.formula(paste0("~", v)), design = des, na.rm = TRUE)
  round(c(est = as.numeric(coef(m)[1]), se = as.numeric(SE(m)[1])), 3)
}
catw <- function(des, v) {  # weighted proportions for factor
  m <- svymean(as.formula(paste0("~factor(", v, ")")), design = des, na.rm = TRUE)
  out <- round(as.numeric(coef(m)) * 100, 1)
  paste(names(coef(m)), out, sep = "=", collapse = " | ")
}

report <- function(dd, tag) {
  dd <- dd[!is.na(dd$nrweight), ]
  des <- svydesign(ids = ~schoolnum, weights = ~nrweight, data = dd, nest = TRUE)
  W("\n--- ", tag, " | n = ", nrow(dd), " | weighted N = ", round(sum(weights(des)), 0), " ---")
  for (v in c("age_num", "female", "international", "undergrad", "grad_student",
              "fin_stress", "food_insec", "firstgen")) {
    if (v %in% names(dd)) W("  ", v, ": est=", wm(des, v)[1], " se=", wm(des, v)[2])
  }
  W("  race_cat: ", catw(des, "race_cat"))
  for (v in c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")) {
    if (v %in% names(dd)) W("  ", v, ": est=", wm(des, v)[1], " se=", wm(des, v)[2])
  }
  for (v in c("anx_score", "flourish")) {
    if (v %in% names(dd)) W("  ", v, ": mean=", wm(des, v)[1], " se=", wm(des, v)[2])
  }
  invisible(des)
}

## ---------- 1) module vs non-module (randomization check) ----------
W("========== RANDOMIZATION CHECK: digital-behavior module vs rest ==========")
mod  <- d[d$module == 1, ]
rest <- d[d$module == 0, ]
report(mod,  "MODULE RESPONDENTS (n)")
report(rest, "NON-MODULE (rest of sample)")

## ---------- 2) main analytic sample ----------
main <- d[d$module == 1 & !is.na(d$dep_any) & !is.na(d$sui_idea) &
          !is.na(d$fin_stress) & !is.na(d$age_num) & !is.na(d$nrweight), ]
W("\n========== MAIN ANALYTIC SAMPLE ==========")
report(main, "MAIN (module + outcomes + core covariates)")

## ---------- 3) by >3h/day online ----------
W("\n========== MAIN by time_gt3h ==========")
report(main[main$time_gt3h == 0 & !is.na(main$time_gt3h), ], "<=3h/day")
report(main[main$time_gt3h == 1 & !is.na(main$time_gt3h), ], ">3h/day")

## ---------- 4) exposure profile (main) ----------
W("\n========== EXPOSURE PROFILE (main, weighted %) ==========")
des_main <- svydesign(ids = ~schoolnum, weights = ~nrweight, data = main[!is.na(main$nrweight), ], nest = TRUE)
W("  time5cat: ", catw(des_main, "time5cat"))
for (v in c("exp_relation", "exp_harassed", "exp_create", "exp_activism", "exp_none",
            "ai_academic", "ai_forbidden", "ai_work", "ai_comm", "ai_health", "ai_fun", "ai_other", "ai_any")) {
  if (v %in% names(main)) W("  ", v, ": est=", wm(des_main, v)[1], " se=", wm(des_main, v)[2])
}
W("  ai_n_uses: mean=", wm(des_main, "ai_n_uses")[1], " se=", wm(des_main, "ai_n_uses")[2])

## ---------- 5) save Table 1 csv ----------
mk <- function(dd, tag) {
  dd <- dd[!is.na(dd$nrweight), ]
  des <- svydesign(ids = ~schoolnum, weights = ~nrweight, data = dd, nest = TRUE)
  out <- data.frame(group = tag, n = nrow(dd), stringsAsFactors = FALSE)
  vars <- c("age_num", "female", "international", "undergrad", "fin_stress", "food_insec",
            "firstgen", "time_gt3h", "ai_any", "dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
  for (v in vars) {
    if (v %in% names(dd)) {
      m <- wm(des, v)
      out[[paste0(v, "_est")]] <- m[1]; out[[paste0(v, "_se")]] <- m[2]
    }
  }
  out
}
tab <- rbind(mk(mod, "module"), mk(rest, "non_module"), mk(main, "main"),
             mk(main[main$time_gt3h == 0, ], "main_le3h"), mk(main[main$time_gt3h == 1, ], "main_gt3h"))
write.csv(tab, file.path(OUT, "table1_hms2324.csv"), row.names = FALSE, fileEncoding = "UTF-8")

close(log_con)
cat("DONE\n")
