source("00_config.R")
## 04_var_centered.R
## Variable-centered weighted associations: digital-behavior exposures -> mental-health outcomes
## HMS 2023-24 main analytic sample.
## Survey design: ids = schoolnum (schools = PSU), weights = nrweight, nest = TRUE.
## Outcomes (5): dep_any, dep_maj, anx_any, sui_idea, lonely
## Exposures: time5cat (0..4 ordinal-ish factor, ref = 0), exp_* (binary), ai_any, ai_n_uses
## Covariates: age_num, female, race_cat, international, undergrad, grad_student,
## fin_stress, food_insec, firstgen

suppressPackageStartupMessages({ library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

## --- reconstruct main analytic sample exactly as in 02_descriptive.R ---
core_cov <- c("age_num","female","race_cat","international","undergrad",
              "grad_student","fin_stress","food_insec","firstgen")
m <- d[d$module == 1, ]
m <- m[!is.na(m$time5cat), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) m <- m[!is.na(m[[v]]), ]
for (v in core_cov) m <- m[!is.na(m[[v]]), ]
cat("main sample n:", nrow(m), "\n")

m$time5cat_f <- factor(m$time5cat, levels = 0:4,
                       labels = c("none","<1h","1-2h","2-3h",">3h"))
m$ai_n_uses_g <- cut(m$ai_n_uses, breaks = c(-0.5, 0.5, 1.5, 2.5, 100),
                     labels = c("0","1","2","3+"))
m$time5cat_f <- relevel(m$time5cat_f, ref = "none")
m$ai_n_uses_g <- relevel(m$ai_n_uses_g, ref = "0")

## weights: raw + winsorized at 99.5 pct (sensitivity for extreme nrweight)
q995 <- quantile(m$nrweight, 0.995, na.rm = TRUE)
m$nrweight_w <- pmin(m$nrweight, q995)
cat("nrweight: median", median(m$nrweight), "p99.5", q995, "max", max(m$nrweight), "\n")

design <- function(w) svydesign(ids = ~schoolnum, weights = w, nest = TRUE, data = m)

adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

## missingness audit on covariates (explains listwise loss)
cat("covariate NA counts (module+outcome-complete set):\n")
print(colSums(is.na(d[d$module == 1 & !is.na(d$time5cat) & !is.na(d$dep_any), core_cov])))

res <- list()
run <- function(outcome, exposure, des, label) {
  fo <- as.formula(paste(outcome, "~", exposure, "+", adj))
  g <- svyglm(fo, design = des, family = quasibinomial())
  cc <- coef(g); se <- sqrt(diag(vcov(g)))
  ci <- confint(g)
  or <- exp(cc); lo <- exp(ci[, 1]); hi <- exp(ci[, 2])
  tt <- data.frame(outcome = outcome, exposure = label, level = names(cc),
                   OR = round(or, 3), LCL = round(lo, 3), UCL = round(hi, 3),
                   p = round(2 * pnorm(-abs(cc / se)), 5))
  tt <- tt[grep(exposure, tt$level, fixed = FALSE), ]
  tt
}

exps <- list(
  list(v = "time5cat_f", lab = "time/day"),
  list(v = "exp_relation", lab = "exp_relation"),
  list(v = "exp_harassed", lab = "exp_harassed"),
  list(v = "exp_create",   lab = "exp_create"),
  list(v = "exp_activism", lab = "exp_activism"),
  list(v = "ai_any",       lab = "ai_any"),
  list(v = "ai_n_uses_g",  lab = "ai_n_uses")
)
oc <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")

d1 <- design(~nrweight)
d2 <- design(~nrweight_w)

for (e in exps) {
  for (o in oc) {
    r1 <- run(o, e$v, d1, paste0(e$lab, " [raw w]"))
    r2 <- run(o, e$v, d2, paste0(e$lab, " [winsor w]"))
    res[[length(res) + 1]] <- rbind(r1, r2)
    cat(sprintf("%-12s x %-14s : %s\n", o, e$lab,
                paste(sprintf("%s=%.2f[%.2f-%.2f]", r1$level, r1$OR, r1$LCL, r1$UCL),
                      collapse = " ")))
    flush.console()
  }
}

all <- do.call(rbind, res)
write.csv(all, file.path(out, "table2_var_centered_OR.csv"), row.names = FALSE)
cat("DONE\n")
