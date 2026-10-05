source("00_config.R")
## 09_replication_cfps.R
## Correction: the CFPS analysis extract never carried a cluster identifier, so
## every svydesign below was silently built with ids = ~1 (no clustering) while the
## manuscript claimed province-level clustering. The province code is now merged in from the
## raw CFPS person file by pid (100% of the analytic rows match), and the design declares it.
## (1) Internal replication in HMS 2024-25: duration (aligned 5-level) + harassment.
## NOTE: the two waves use different duration scales (7 vs 8 response options), so
## levels are NOT comparable; only the SHAPE of the association is replicated.
## (2) External validation in CFPS 2022: short-video use / mobile hours -> depressive
## symptoms (CES-D8 >= 8 and >= 10) and loneliness, in 18-30-year-olds with
## tertiary education (the closest available analogue of US college students).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(haven); library(dplyr) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")

## ---------------- (1) HMS 2024-25 replication ----------------
h2 <- readRDS(file.path(out, "hms_2425_analysis.rds"))
W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")
oc <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")

s <- h2[h2$module == 1 & !is.na(h2$time5cat), ]
for (v in oc) s <- s[!is.na(s[[v]]), ]
for (v in W) s <- s[!is.na(s[[v]]), ]
s$time5cat_f <- relevel(factor(s$time5cat, levels = 0:4,
                               labels = c("none","<1h","1-2h","2-3h",">3h")), ref = "none")
cat("HMS 2024-25 replication n:", nrow(s), "\n")
des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

rep_res <- list()
for (e in c("time5cat_f","exp_harassed")) {
  for (o in oc) {
    g <- svyglm(as.formula(paste(o, "~", e, "+", adj)), design = des, family = quasibinomial())
    ci <- confint(g)
    tt <- data.frame(wave = "2024-25", outcome = o, level = rownames(ci),
                     OR = round(exp(coef(g)), 3), LCL = round(exp(ci[, 1]), 3),
                     UCL = round(exp(ci[, 2]), 3))
    tt <- tt[grep(e, tt$level), ]
    rep_res[[length(rep_res) + 1]] <- tt
    cat(sprintf("%-14s %-9s: %s\n", e, o,
                paste(sprintf("%s=%.2f[%.2f-%.2f]", tt$level, tt$OR, tt$LCL, tt$UCL),
                      collapse = " ")))
  }
}
rep_res <- do.call(rbind, rep_res)
write.csv(rep_res, file.path(out, "table6_replication_2425.csv"), row.names = FALSE)

## ---------------- (2) CFPS 2022 external validation ----------------
c22 <- readRDS(file.path(out, "cfps_2022_analysis.rds"))
cf <- c22[c22$young_hiedu == 1, ]
cf <- cf[!is.na(cf$wt), ]

## ---- cluster identifier -------------------------------------------------------
## The extract keeps only pid, so the sampling design has to be re-attached from the raw
## file. provcd22 is the province code (32 provinces present in this subsample); psu is the
## true primary sampling unit and is carried too so the choice can be audited.
rawf <- file.path(CFPS_DIR, "cfps2022person.dta")
geo <- read_dta(rawf, col_select = c("pid", "psu", "provcd22"))
geo <- distinct(geo, pid, .keep_all = TRUE)
cf <- left_join(cf, geo, by = "pid")
stopifnot(sum(is.na(cf$provcd22)) == 0)
cat("\nCFPS 2022 young_hiedu n:", nrow(cf), "\n")
cat("clusters: province", length(unique(cf$provcd22)), "| psu", length(unique(cf$psu)), "\n")
cat("sv_watch valid:", sum(!is.na(cf$sv_watch)),
    "| mobile_hours valid:", sum(!is.na(cf$mobile_hours)),
    "| dep8 valid:", sum(!is.na(cf$dep8)),
    "| lonely valid:", sum(!is.na(cf$lonely)), "\n")

cf$urban_f <- factor(cf$urban)
## NOTE: qu201a is recorded in MINUTES per day (median 240); convert to hours
cf$mobile_h <- cf$mobile_hours / 60
cat("mobile hours/day: median", median(cf$mobile_h, na.rm = TRUE),
    "| IQR", paste(round(quantile(cf$mobile_h, c(.25,.75), na.rm = TRUE), 1), collapse = "-"), "\n")
cf$mobile_g <- cut(cf$mobile_h, breaks = c(-0.01, 2, 4, 6, Inf),
                   labels = c("<=2h","2-4h","4-6h",">6h"))
des2 <- svydesign(ids = ~provcd22, weights = ~wt, data = cf)

cf_res <- list()
for (e in c("sv_watch","mobile_g","mobile_h")) {
  for (o in c("dep8","dep10","lonely")) {
    sub <- cf[!is.na(cf[[e]]) & !is.na(cf[[o]]), ]
    d2 <- svydesign(ids = ~provcd22, weights = ~wt, data = sub)
    g <- svyglm(as.formula(paste(o, "~", e, "+ age_num + female + urban_f")),
                design = d2, family = quasibinomial())
    ci <- confint(g)
    tt <- data.frame(outcome = o, level = rownames(ci),
                     OR = round(exp(coef(g)), 3), LCL = round(exp(ci[, 1]), 3),
                     UCL = round(exp(ci[, 2]), 3), n = nrow(sub))
    tt <- tt[grep(e, tt$level), ]
    cf_res[[length(cf_res) + 1]] <- tt
    cat(sprintf("CFPS %-10s %-7s n=%5d: %s\n", e, o, nrow(sub),
                paste(sprintf("%s=%.2f[%.2f-%.2f]", tt$level, tt$OR, tt$LCL, tt$UCL),
                      collapse = " ")))
  }
}
cf_res <- do.call(rbind, cf_res)
write.csv(cf_res, file.path(out, "table7_cfps_validation.csv"), row.names = FALSE)
cat("DONE\n")
