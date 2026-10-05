source("00_config.R")
## 05_lca_outcomes.R
## Correction: the row index used to re-attach LCA classes to the survey object
## was obtained with as.integer(rownames(dd)). That is only correct while `d` is a plain
## data.frame. The extract is a tibble, so as soon as tibble is on the search path (any
## library(dplyr) earlier in the session) rownames(dd) returns 1..n instead of the original
## row numbers, every class label shifts, and the run SILENTLY produces a different answer
## (n 40,768 -> 37,325, 195 -> 196 schools, and the highest-risk class disappears). The
## index is now built from the selection rule itself, which is class-agnostic.
## Person-centred step: assign LCA classes and estimate weighted class prevalence
## plus adjusted associations with the 5 outcomes (weighted svyglm, schools = PSU).
## Usage: Rscript 05_lca_outcomes.R K (K = number of classes chosen from LCA fit stats)

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(poLCA) })

args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else 5
set <- if (length(args) >= 2) args[2] else "lca"      # "lca" (12 items) or "lca2" (11, no ai_other)

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs"); lcaout <- file.path(out, set)
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

ind <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
         "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun",
         "ai_other")
if (set == "lca2") ind <- setdiff(ind, "ai_other")
cat("indicator set:", set, "(", length(ind), "items )\n")

## rebuild the LCA data set in exactly the same order as 03b
## keep is a logical mask on the FULL extract, so idx is independent of data-frame class
keep <- d$module == 1 & !is.na(d$time5cat)
for (v in ind) keep <- keep & !is.na(d[[v]])
keep[is.na(keep)] <- FALSE
dd <- d[keep, ]
idx <- which(keep)
cat("LCA n:", nrow(dd), "\n")
stopifnot(nrow(dd) == length(idx), max(idx) <= nrow(d))
## guards against the silent mislabelling this script used to produce
cat("idx range:", min(idx), "-", max(idx), "| n:", length(idx), "\n")

fitf <- file.path(lcaout, sprintf("fit_k%d.rds", K))
stopifnot(file.exists(fitf))
m <- readRDS(fitf)
cls <- apply(m$posterior, 1, which.max)
cat("class sizes (unweighted):", table(cls), "\n")
cat("modal-class assignment entropy check: mean max posterior =",
    round(mean(apply(m$posterior, 1, max)), 3), "\n")

## attach class to the full data frame by row position
d$class <- NA_integer_
d$class[idx] <- cls

## sample for regression: class assigned + outcomes + covariates
W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")
s <- d[!is.na(d$class), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) s <- s[!is.na(s[[v]]), ]
for (v in W) s <- s[!is.na(s[[v]]), ]
s$class_f <- factor(s$class)
cat("regression n:", nrow(s), "\n")

des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)

## weighted class prevalence
prev <- svymean(~factor(class_f), design = des)
pv <- data.frame(class = 1:K,
                 prev_pct = round(100 * as.numeric(coef(prev)), 2),
                 se_pct = round(100 * as.numeric(SE(prev)), 2))
cat("\nweighted class prevalence (%):\n"); print(pv, row.names = FALSE)

## weighted outcome prevalence by class
oc <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")
byc <- list()
for (o in oc) {
  f <- as.formula(paste0("~factor(class_f):", o))
  tb <- svyby(as.formula(paste0("~", o)), by = ~factor(class_f), design = des,
              FUN = svymean, keep.names = FALSE)
  x <- data.frame(outcome = o, class = 1:K,
                  pct = round(100 * as.numeric(tb[[2]]), 2),
                  se = round(100 * as.numeric(tb[[3]]), 2))
  byc[[o]] <- x
}
byc <- do.call(rbind, byc)
cat("\nweighted outcome prevalence (%) by class:\n"); print(byc, row.names = FALSE)

## adjusted associations: reference = largest class
ref <- as.integer(names(sort(table(s$class_f), decreasing = TRUE))[1])
cat("\nreference class:", ref, "\n")
s$class_f <- relevel(s$class_f, ref = as.character(ref))
des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"

res <- list()
for (o in oc) {
  g <- svyglm(as.formula(paste(o, "~ class_f +", adj)), design = des,
              family = quasibinomial())
  ci <- confint(g)
  tt <- data.frame(outcome = o, level = rownames(ci),
                   OR = round(exp(coef(g)), 3), LCL = round(exp(ci[, 1]), 3),
                   UCL = round(exp(ci[, 2]), 3))
  tt <- tt[grep("class_f", tt$level), ]
  res[[o]] <- tt
  cat(sprintf("%-9s: %s\n", o,
              paste(sprintf("%s=%.2f[%.2f-%.2f]", tt$level, tt$OR, tt$LCL, tt$UCL),
                    collapse = "  ")))
}
all <- do.call(rbind, res)

write.csv(pv, file.path(out, sprintf("%s_k%d_prevalence.csv", set, K)), row.names = FALSE)
write.csv(byc, file.path(out, sprintf("%s_k%d_outcome_by_class.csv", set, K)), row.names = FALSE)
write.csv(all, file.path(out, sprintf("%s_k%d_adjusted_OR.csv", set, K)), row.names = FALSE)
saveRDS(d$class, file.path(lcaout, sprintf("class_assignment_k%d.rds", K)))
cat("DONE\n")
