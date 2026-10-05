source("00_config.R")
## 25_sens_school_type.R
## Effect modification by institutional context (stratified analysis + formal test).
##
## Why this matters: the school-level negative control (graduation rate) flagged residual
## school-context confounding on the DURATION line only (report section 9.3). If a single
## institution type drives the association, the pooled estimate is not describing a
## general student-level phenomenon. Stratifying by institutional type / sector answers
## that directly, and a formal interaction test guards against over-reading stratum
## noise (the same discipline the het-benchmark applies to the causal forest).
##
## Stratifiers (HMS institutional characteristics, IPEDS-derived):
## inst_type : 1=Associate's, 2=Baccalaureate, 3=Doctorate-granting,
## 4=Master's, 5=Special Focus
## inst_public : 0=Private, 1=Public
##
## Note on inference: inst_type / inst_public are SCHOOL-level variables. Stratum-specific
## fits therefore use disjoint school sets, so their SEs are valid; the pooled interaction
## model reuses all schools. We treat the interaction test as the inferential statement
## and the stratum-specific ORs as descriptive.
##
## Interaction test: survey::regTermTest (F test with design denominator df).
## Fallback: manual Wald chi-square on the interaction block (used if regTermTest errors).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(dplyr); library(haven) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- zap_labels(readRDS(file.path(out, "hms_2324_analysis.rds")))
stopifnot(all(vapply(d, function(x) !inherits(x, "haven_labelled"), logical(1))))

core_cov <- c("age_num", "female", "race_cat", "international", "undergrad",
              "grad_student", "fin_stress", "food_insec", "firstgen")
OC <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")

m <- d[d$module == 1, ]
m <- m[!is.na(m$time5cat), ]
for (v in OC) m <- m[!is.na(m[[v]]), ]
for (v in core_cov) m <- m[!is.na(m[[v]]), ]
m <- m[!is.na(m$inst_type) & !is.na(m$inst_public), ]
cat("analytic sample with school characteristics n =", nrow(m), "\n")

TYPE_LAB <- c("Associate's", "Baccalaureate", "Doctorate-granting", "Master's", "Special Focus")
PUB_LAB  <- c("Private", "Public")
m$inst_type_f   <- factor(m$inst_type,   levels = 1:5, labels = TYPE_LAB)
m$inst_public_f <- factor(m$inst_public, levels = 0:1, labels = PUB_LAB)
cat("\nunweighted n by institution type:\n"); print(table(m$inst_type_f))
cat("\nunweighted n by sector:\n");           print(table(m$inst_public_f))

adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"
EXPS <- list(list(v = "time_gt3h",    lab = "screen>3h/d"),
             list(v = "exp_harassed", lab = "online harassment"),
             list(v = "ai_any",       lab = "any genAI use"))

des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = m)

## ---------- stratum-specific weighted model (fresh design on the subset) -------------
strat_or <- function(o, ev, strat, lev) {
  dat <- m[m[[strat]] == lev & !is.na(m[[ev]]), ]
  if (nrow(dat) < 50) return(NULL)
  sub <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dat)
  g <- svyglm(as.formula(paste(o, "~", ev, "+", adj)), design = sub, family = quasibinomial())
  ci <- confint(g); b <- coef(g)[[ev]]; se <- sqrt(diag(vcov(g)))[[ev]]
  ns <- length(unique(dat$schoolnum))
  lo <- suppressWarnings(try(ci[ev, 1], silent = TRUE))
  hi <- suppressWarnings(try(ci[ev, 2], silent = TRUE))
  lo <- if (inherits(lo, "try-error")) NA_real_ else as.numeric(lo)
  hi <- if (inherits(hi, "try-error")) NA_real_ else as.numeric(hi)
 ## confint returns limits on the LOG-ODDS scale; exponentiate so OR/LCL/UCL share the
 ## odds scale (previously LCL/UCL were written un-exponentiated while OR was exp'd).
  data.frame(outcome = o, exposure = ev, stratifier = strat, stratum = as.character(lev),
             n = nrow(dat), n_schools = ns,
             OR = round(exp(b), 3), LCL = round(exp(lo), 3), UCL = round(exp(hi), 3),
             logOR = b, p = 2 * pnorm(-abs(b / se)), failed = FALSE,
 ## with < 10 PSUs the design denominator df is exhausted and the survey
 ## CI is not estimable; the point estimate is retained but flagged.
             ci_ok = ns >= 10)
}

interaction_test <- function(g, ev, strat) {
  tt <- try(regTermTest(g, paste(ev, strat, sep = ":")), silent = TRUE)
  if (!inherits(tt, "try-error")) {
    return(list(stat = as.numeric(tt$Ftest[1, 1]), df = as.integer(tt$df),
                ddf = as.numeric(tt$ddf), p = as.numeric(tt$p[1, 1]),
                method = "regTermTest F"))
  }
  nm <- grep(paste0("^", ev, ":", strat), names(coef(g)), value = TRUE)
  V <- vcov(g)[nm, nm, drop = FALSE]; b <- coef(g)[nm]
  W <- as.numeric(t(b) %*% qr.solve(V) %*% b)
  list(stat = W, df = length(nm), ddf = NA_real_,
       p = pchisq(W, length(nm), lower.tail = FALSE), method = "manual Wald chi2")
}

rows <- list(); inter <- list()
for (strat in c("inst_type_f", "inst_public_f")) {
  for (e in EXPS) {
    for (o in OC) {
 ## ---- stratum-specific ----
      for (k in levels(m[[strat]])) {
        r <- try(strat_or(o, e$v, strat, k), silent = TRUE)
        if (inherits(r, "try-error") || is.null(r)) {
          cat("  stratum FAIL:", strat, "|", k, "|", o, "|", e$v, "\n")
          rows[[length(rows) + 1]] <- data.frame(
            outcome = o, exposure = e$v, exposure_lab = e$lab, stratifier = strat,
            stratum = as.character(k), n = sum(m[[strat]] == k, na.rm = TRUE),
            n_schools = NA, OR = NA, LCL = NA, UCL = NA, logOR = NA, p = NA,
            failed = TRUE, ci_ok = FALSE)
        } else {
          r$exposure_lab <- e$lab
          rows[[length(rows) + 1]] <- r
        }
      }
 ## ---- formal interaction test ----
      g <- try(svyglm(as.formula(paste(o, "~", e$v, "*", strat, "+", adj)),
                      design = des, family = quasibinomial()), silent = TRUE)
      if (inherits(g, "try-error")) { cat("  interaction FAIL:", o, e$v, strat, "\n"); next }
      it <- interaction_test(g, e$v, strat)
      inter[[length(inter) + 1]] <- data.frame(
        outcome = o, exposure = e$v, exposure_lab = e$lab, stratifier = strat,
        stat = it$stat, df = it$df, ddf = it$ddf,
        p_interaction = it$p, method = it$method)
      cat(sprintf("%-9s x %-18s | %-13s : %s  p = %.4f\n",
                  o, e$lab, strat, it$method, it$p))
      flush.console()
    }
  }
}

strat_df <- bind_rows(rows)
inter_df <- bind_rows(inter)
inter_df <- inter_df %>% group_by(stratifier) %>%
  mutate(p_interaction_bh = p.adjust(p_interaction, method = "BH")) %>% ungroup()

## per exposure x outcome x stratifier: spread of stratum estimates.
## Only strata with an estimable CI enter the spread (Special Focus has 6 schools ->
## design df exhausted; its point estimate is reported but cannot carry inference).
spread <- strat_df %>% filter(!failed, ci_ok) %>%
  group_by(outcome, exposure, exposure_lab, stratifier) %>% summarise(
    n_strata = n(),
    OR_min = min(OR), OR_max = max(OR),
    spread_logOR = round(max(logOR) - min(logOR), 4),
    all_same_sign = all(logOR > 0) || all(logOR < 0),
    .groups = "drop")

write.csv(strat_df, file.path(out, "table18_school_type_stratified.csv"), row.names = FALSE)
write.csv(inter_df, file.path(out, "table18_school_type_interaction.csv"), row.names = FALSE)
write.csv(spread,   file.path(out, "table18_school_type_spread.csv"),     row.names = FALSE)

cat("\n=== interaction tests (BH within stratifier) ===\n")
print(as.data.frame(inter_df %>% select(outcome, exposure_lab, stratifier, df, ddf,
                                        p_interaction, p_interaction_bh)), row.names = FALSE)
cat("\n=== stratum spread ===\n")
print(as.data.frame(spread), row.names = FALSE)
cat("\ninteraction p<0.05 raw:", sum(inter_df$p_interaction < 0.05),
    "| BH q<0.05:", sum(inter_df$p_interaction_bh < 0.05),
    "| stratum failures:", sum(strat_df$failed == TRUE), "/", nrow(strat_df), "\n")
cat("DONE\n")
