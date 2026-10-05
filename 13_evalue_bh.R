source("00_config.R")
## 13_evalue_bh.R
## (1) E-values for the variable-centred odds ratios (unmeasured-confounding sensitivity).
## RR approximated from OR and the outcome prevalence in the unexposed
## (VanderWeele & Ding): RR ~ OR / (1 - P0 + P0 * OR)
## E-value = RR + sqrt(RR * (RR - 1)); reported for the point estimate and for the
## confidence limit closest to the null (conservative).
## (2) Benjamini-Hochberg correction across all primary (discovery-set) tests.

suppressPackageStartupMessages({ library(survey) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

## ---------------- (1) E-values ----------------
cc <- c("age_num","female","race_cat","international","undergrad","grad_student",
        "fin_stress","food_insec","firstgen")
m <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) m <- m[!is.na(m[[v]]), ]
for (v in cc) m <- m[!is.na(m[[v]]), ]

cat("nrweight NA in analytic set:", sum(is.na(m$nrweight)), "\n")
prev0 <- function(idx, o) {                     # weighted prevalence in reference group
  s <- m[idx & !is.na(m[[o]]) & !is.na(m$nrweight), ]
  des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
  as.numeric(svymean(as.formula(paste0("~", o)), design = des)[1])
}
## NOTE (recurring R trap): a logical index containing NA makes data.frame[idx, ] return
## all-NA rows. Every comparison below is guarded with !is.na.
ref <- list(
  list(lab = "screen>3h/d",  idx = !is.na(m$time5cat) & m$time5cat == 0,
       term = "time5cat_f>3h"),
  list(lab = "harassment",   idx = !is.na(m$exp_harassed) & m$exp_harassed == 0,
       term = "exp_harassed"),
  list(lab = "any genAI",    idx = !is.na(m$ai_any) & m$ai_any == 0,
       term = "ai_any"),
  list(lab = "genAI 3+ uses", idx = !is.na(m$ai_n_uses) & m$ai_n_uses < 3,
       term = "ai_n_uses_g3+"))

t2 <- read.csv(file.path(out, "table2_var_centered_OR.csv"))
t2 <- t2[grepl("raw w", t2$exposure, fixed = TRUE), ]

ev <- function(or) {                            # E-value for RR >= 1 (OR -> RR approx done outside)
  rr <- or
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
rr_from_or <- function(or, p0) or / (1 - p0 + p0 * or)

rows <- list()
for (r in ref) {
  for (o in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) {
    line <- t2[t2$outcome == o & t2$level == r$term, ]
    if (nrow(line) == 0) next
    p0 <- prev0(r$idx, o)
    or <- line$OR; lo <- line$LCL; hi <- line$UCL
    rr <- rr_from_or(or, p0)
    ci_or <- if (abs(lo - 1) < abs(hi - 1)) lo else hi   # limit closest to null
    rr_ci <- rr_from_or(ci_or, p0)
    rows[[length(rows) + 1]] <- data.frame(
      exposure = r$lab, outcome = o, prev_unexposed_pct = round(100 * p0, 1),
      OR = or, LCL = lo, UCL = hi, RR_approx = round(rr, 3),
      E_value = round(ev(rr), 2), E_value_CI = round(ev(rr_ci), 2))
  }
}
ev_tab <- do.call(rbind, rows)
cat("=== E-values (risk-ratio scale, approximation for common outcomes) ===\n")
print(ev_tab, row.names = FALSE)
write.csv(ev_tab, file.path(out, "table10_evalues.csv"), row.names = FALSE)

## ---------------- (2) Benjamini-Hochberg ----------------
p_from_ci <- function(or, lo, hi) {
  se <- (log(hi) - log(lo)) / (2 * 1.96)
  z <- log(or) / se
  2 * pnorm(-abs(z))
}
pool <- list()

## (a) variable-centred weighted models
a <- t2[, c("outcome","exposure","level","OR","LCL","UCL","p")]
a$family <- "variable-centred"; a$test <- paste(a$outcome, a$level)
pool[["a"]] <- a[, c("family","test","outcome","OR","p")]

## (b) TMLE (cluster-robust p)
t3 <- read.csv(file.path(out, "table3_tmle_RD.csv"))
t3b <- read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv"))
b <- rbind(t3, t3b)[, c("outcome","exposure","psi","se_cluster","p_cluster")]
names(b) <- c("outcome","test","OR","se","p"); b$family <- "TMLE"
pool[["b"]] <- b[, c("family","test","outcome","OR","p")]

## (c) LCA profile contrasts (p back-calculated from CI)
lk <- read.csv(file.path(out, "lca2_k5_adjusted_OR.csv"))
c_ <- data.frame(family = "LCA profile", test = paste(lk$outcome, lk$level),
                 outcome = lk$outcome, OR = lk$OR,
                 p = p_from_ci(lk$OR, lk$LCL, lk$UCL))
pool[["c"]] <- c_

all <- do.call(rbind, pool)
all$q_overall <- p.adjust(all$p, method = "BH")
all <- all[order(all$p), ]
byfam <- split(seq_len(nrow(all)), all$outcome)
all$q_by_outcome <- NA
for (i in seq_along(byfam)) {
  idx <- byfam[[i]]
  all$q_by_outcome[idx] <- p.adjust(all$p[idx], method = "BH")
}
all$sig_raw <- all$p < 0.05
all$sig_BH_overall <- all$q_overall < 0.05
all$sig_BH_by_outcome <- all$q_by_outcome < 0.05

cat("\n=== all primary tests with BH q-values (sorted by p) ===\n")
print(all[, c("family","test","OR","p","q_overall","q_by_outcome",
              "sig_BH_overall")], row.names = FALSE)
cat("\ntests:", nrow(all), "| raw p<.05:", sum(all$sig_raw),
    "| BH-overall q<.05:", sum(all$sig_BH_overall),
    "| BH-by-outcome q<.05:", sum(all$sig_BH_by_outcome), "\n")
write.csv(all, file.path(out, "table11_bh_corrected.csv"), row.names = FALSE)
cat("DONE\n")
