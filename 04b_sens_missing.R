source("00_config.R")
## 04b_sens_missing.R
## Sensitivity to listwise deletion of covariates (main CC n=41,862 vs full 45,504).
## Variant A: missing-category (firstgen / race_cat / female / international get explicit
## "Unknown" level) -> retains near-full sample.
## Variant B: complete-case restricted to the same rows (reference, reproduced here for
## a like-for-like comparison on identical covariate handling).
## Reports ORs for the two headline exposures: time5cat (>3h vs none) and exp_harassed.

suppressPackageStartupMessages({ library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) m <- m[!is.na(m[[v]]), ]
for (v in c("age_num","fin_stress","food_insec","undergrad","grad_student"))
  m <- m[!is.na(m[[v]]), ]
cat("base sample (outcomes + fully-observed covariates):", nrow(m), "\n")

## --- Variant A: explicit unknown categories ---
a <- m
for (v in c("female","race_cat","international","firstgen")) {
  a[[v]] <- as.character(a[[v]])
  a[[v]][is.na(a[[v]])] <- "Unknown"
  a[[v]] <- factor(a[[v]])
}
a$time5cat_f <- relevel(factor(a$time5cat, levels = 0:4,
                               labels = c("none","<1h","1-2h","2-3h",">3h")), ref = "none")

## --- Variant B: complete cases ---
b <- m[!is.na(m$female) & !is.na(m$race_cat) & !is.na(m$international) &
         !is.na(m$firstgen), ]
b$time5cat_f <- relevel(factor(b$time5cat, levels = 0:4,
                               labels = c("none","<1h","1-2h","2-3h",">3h")), ref = "none")
cat("variant A n:", nrow(a), " variant B n:", nrow(b), "\n")

adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

one <- function(dat, outcome, exposure) {
  des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dat)
  g <- svyglm(as.formula(paste(outcome, "~", exposure, "+", adj)),
              design = des, family = quasibinomial())
  ci <- confint(g)
  tt <- data.frame(OR = round(exp(coef(g)), 3), LCL = round(exp(ci[, 1]), 3),
                   UCL = round(exp(ci[, 2]), 3))
  tt$level <- rownames(tt)
  tt[grep(exposure, tt$level), ]
}

res <- list()
for (o in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) {
  for (e in c("time5cat_f","exp_harassed","ai_any")) {
    ra <- one(a, o, e); rb <- one(b, o, e)
    res[[length(res) + 1]] <- data.frame(outcome = o, exposure = e,
                                         level = ra$level,
                                         OR_A = ra$OR, LCL_A = ra$LCL, UCL_A = ra$UCL,
                                         OR_B = rb$OR, LCL_B = rb$LCL, UCL_B = rb$UCL)
  }
}
all <- do.call(rbind, res)
write.csv(all, file.path(out, "table_sens_missing.csv"), row.names = FALSE)
print(all, row.names = FALSE)
cat("DONE\n")
