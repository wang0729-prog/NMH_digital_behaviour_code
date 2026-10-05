source("00_config.R")
## 14_negcontrol.R
## Negative- and positive-control outcomes to probe residual confounding.
## Positive control : flourish (flourishing score) -- expected to DECREASE with risk
## exposures; verifies the pipeline can detect a known association.
## Negative control : inst_gradrate (institution 6-year graduation rate) -- an
## individual student's digital behaviour cannot causally affect a
## school-level graduation rate; any association signals residual
## school-level / unmeasured confounding.
## Supplementary NC : BMI (lb/in^2) in the height-weight sub-sample (n ~ 3,165, 3%).

suppressPackageStartupMessages({ library(survey); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")
adj <- paste(W, collapse = " + ")

m <- d[d$module == 1 & !is.na(d$time5cat), ]
m$ai_heavy <- as.integer(m$ai_n_uses >= 3)
exps <- c("time_gt3h","exp_harassed","ai_heavy")

## ---- main controls on the module sub-sample ----
run <- function(yvar) {
  rows <- list()
  for (e in exps) {
    s <- m[!is.na(m[[e]]) & !is.na(m[[yvar]]), ]
    for (v in W) s <- s[!is.na(s[[v]]), ]
    des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
    g <- svyglm(as.formula(paste(yvar, "~", e, "+", adj)), design = des)
    ci <- confint(g)
    b <- coef(g)[e]; se <- sqrt(diag(vcov(g)))[e]
    rows[[length(rows) + 1]] <- data.frame(outcome = yvar, exposure = e, n = nrow(s),
                                           beta = round(b, 4), se = round(se, 4),
                                           LCL = round(ci[e, 1], 4), UCL = round(ci[e, 2], 4),
                                           p = signif(2 * pnorm(-abs(b / se)), 3))
  }
  do.call(rbind, rows)
}
res <- rbind(run("flourish"), run("inst_gradrate"))
res$scale <- ifelse(res$outcome == "inst_gradrate",
                    "percentage points of graduation rate", "flourishing score points")
cat("=== positive control (flourish) and negative control (inst_gradrate) ===\n")
print(res, row.names = FALSE)

## ---- supplementary: BMI in the height/weight sub-sample ----
## The analysis RDS keeps the raw dta row order, so columns can be bound by position.
f <- file.path(HMS_DIR, "HMS_2023-2024_PUBLIC_instchars.dta")
hv <- c("weight","height_ft","height_in")
if (all(hv %in% names(read_dta(f, n_max = 1)))) {
  hw <- zap_labels(read_dta(f, col_select = all_of(hv)))
  stopifnot(nrow(hw) == nrow(d))
 ## height/weight are stored as character after zap_labels -> coerce explicitly
  w_lb <- suppressWarnings(as.numeric(as.character(hw$weight)))
  ft   <- suppressWarnings(as.numeric(as.character(hw$height_ft)))
  inch <- suppressWarnings(as.numeric(as.character(hw$height_in)))
  h_in <- 12 * ft + inch
  bmi <- ifelse(!is.na(w_lb) & !is.na(h_in) & h_in > 48 & h_in < 84 &
                  w_lb > 70 & w_lb < 500, 703 * w_lb / h_in^2, NA)
  db <- d
  db$bmi <- bmi
  mb <- db[db$module == 1 & !is.na(db$time5cat) & !is.na(db$bmi), ]
  cat("\nBMI sub-sample n:", nrow(mb), "(",
      round(100 * nrow(mb) / sum(d$module == 1), 1), "% of module sub-sample )\n")
  mb$ai_heavy <- as.integer(mb$ai_n_uses >= 3)
  rows <- list()
  for (e in exps) {
    s <- mb[!is.na(mb[[e]]), ]
    for (v in W) s <- s[!is.na(s[[v]]), ]
    des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
    g <- svyglm(as.formula(paste("bmi ~", e, "+", adj)), design = des)
    ci <- confint(g); b <- coef(g)[e]; se <- sqrt(diag(vcov(g)))[e]
    rows[[length(rows) + 1]] <- data.frame(outcome = "BMI", exposure = e, n = nrow(s),
                                           beta = round(b, 4), se = round(se, 4),
                                           LCL = round(ci[e, 1], 4), UCL = round(ci[e, 2], 4),
                                           p = signif(2 * pnorm(-abs(b / se)), 3))
  }
  bm <- do.call(rbind, rows)
  bm$scale <- "kg/m^2"
  print(bm, row.names = FALSE)
  res <- rbind(res, bm)
}
write.csv(res, file.path(out, "table9_control_outcomes.csv"), row.names = FALSE)
cat("DONE\n")
