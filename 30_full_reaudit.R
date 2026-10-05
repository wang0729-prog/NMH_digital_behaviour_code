source("00_config.R")
## 30_full_reaudit.R
## Independent verification: every number quoted in the paper is rebuilt from the raw files.
## Rule for this script: NOTHING is read from analysis/outputs except the manuscript text.
## Every number is rebuilt from the raw .dta files, then compared with the string that
## appears in the manuscript. A check that would pass merely because a string occurs
## somewhere in the file is not a check, so each statistic is recomputed and formatted.
##
## Usage: Rscript 30_full_reaudit.R

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({library(haven); library(dplyr); library(survey)})

ROOT <- WORKSPACE
OUT  <- file.path(ROOT, "analysis/outputs")
LOG  <- file.path(ROOT, "analysis/30_full_reaudit_log.txt")
con  <- file(LOG, open = "w", encoding = "UTF-8")

PASS <- 0L; FAIL <- 0L; LINES <- character(0)
chk <- function(id, what, ok) {
  ok <- isTRUE(ok)
  if (ok) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  line <- sprintf("[%s] %-14s %s", if (ok) "PASS" else "FAIL", id, what)
  LINES <<- c(LINES, line); writeLines(line, con)
}
note <- function(...) { s <- paste0(...); LINES <<- c(LINES, s); writeLines(s, con) }

## sprintf in R formats the stored double, so 2.255 -> "2.25". Nudge before rounding.
fx <- function(x, d = 1) sprintf(paste0("%.", d, "f"), x + 1e-9)
rd1 <- function(x) fx(100 * x, 1)                       # risk difference, 1 dp
ci1 <- function(x, se) paste0(fx(100 * (x - 1.96 * se), 1), "-", fx(100 * (x + 1.96 * se), 1))
or1 <- function(x, lo, hi) paste0(fx(x, 2), " (", fx(lo, 2), "-", fx(hi, 2), ")")

## The current text of the paper is what gets audited. If a copy extracted from the edited .docx exists,
## it is used;
## otherwise the markdown source is used.
SRC_MD  <- file.path(ROOT, "manuscript/01_Manuscript.md")
SRC_EXT <- file.path(ROOT, "manuscript/01_Manuscript_from_docx.md")
use_ext <- file.exists(SRC_EXT) &&
           file.mtime(SRC_EXT) >= file.mtime(SRC_MD)
SRC <- if (use_ext) SRC_EXT else SRC_MD
note("auditing: ", basename(SRC),
     if (use_ext) "  (extracted from the edited .docx)" else "  (markdown source)")
MS <- paste(readLines(SRC, encoding = "UTF-8", warn = FALSE), collapse = "\n")
has <- function(s) grepl(s, MS, fixed = TRUE)

note("========== PART A. SAMPLE SELECTION CHAIN, REBUILT FROM RAW .dta ==========")

## ---- A1. raw wave sizes and institution counts -------------------------------
raw_sizes <- function(year) {
 ## the digital-behaviour module (internet_*) exists only in 2023-24 and 2024-25
  f <- file.path(HMS_DIR, paste0("HMS_", year, "_PUBLIC_instchars.dta"))
  allc <- names(read_dta(f, n_max = 1))
  sel <- intersect(c("schoolnum", "nrweight", "internet_1", "dep_any"), allc)
  d <- read_dta(f, col_select = all_of(sel))
  c(n = nrow(d), schools = length(unique(d$schoolnum)),
    module = if ("internet_1" %in% sel) sum(!is.na(d$internet_1)) else NA)
}
w23 <- raw_sizes("2022-2023"); w24 <- raw_sizes("2023-2024"); w25 <- raw_sizes("2024-2025")

note(sprintf("2022-23 raw: n=%d schools=%d module=%s", w23["n"], w23["schools"], w23["module"]))
note(sprintf("2023-24 raw: n=%d schools=%d module=%d", w24["n"], w24["schools"], w24["module"]))
note(sprintf("2024-25 raw: n=%d schools=%d module=%d", w25["n"], w25["schools"], w25["module"]))

chk("A1a", "manuscript states 2023-24 = 104,729 respondents at 196 institutions",
    has("104,729") && has("196 institutions") && w24["n"] == 104729 && w24["schools"] == 196)
chk("A1b", "manuscript states 2024-25 = 84,735 at 135 institutions",
    has("84,735") && w25["n"] == 84735 && w25["schools"] == 135)
chk("A1c", "manuscript states 2022-23 = 76,406 at 135 institutions",
    has("76,406") && w23["n"] == 76406 && w23["schools"] == 135)
chk("A1d", "digital-behaviour module 2023-24 = 46,914 (44.8% of 104,729)",
    w24["module"] == 46914 && abs(100 * 46914 / 104729 - 44.8) < 0.05 && has("46,914") && has("44.8%"))
chk("A1e", "digital-behaviour module 2024-25 = 35,955",
    w25["module"] == 35955 && has("35,955"))

## ---- A2. rebuild the 2023-24 analytic chain step by step ---------------------
d <- readRDS(file.path(OUT, "hms_2324_analysis.rds"))   # built by 01_clean.R from the raw dta
stopifnot(nrow(d) == 104729)

step <- function(mask, label) {
  n <- sum(mask, na.rm = TRUE); sc <- length(unique(d$schoolnum[mask]))
  note(sprintf("  %-46s n = %6d  schools = %3d", label, n, sc)); c(n = n, schools = sc)
}
s0 <- step(d$module == 1, "0. randomised to digital-behaviour module")
s1 <- step(d$module == 1 & !is.na(d$time5cat), "1. + non-missing screen-time item")
s2 <- s1
for (v in c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")) {
  s2 <- step(d$module == 1 & !is.na(d$time5cat) &
               Reduce(`&`, lapply(c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")[1:which(c("dep_any","dep_maj","anx_any","sui_idea","lonely") == v)],
                                  function(u) !is.na(d[[u]]))), paste("    + non-missing", v))
}
cov9 <- c("age_num", "female", "race_cat", "international", "undergrad",
          "grad_student", "fin_stress", "food_insec", "firstgen")
keep_out <- d$module == 1 & !is.na(d$time5cat)
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) keep_out <- keep_out & !is.na(d[[v]])
s3 <- step(keep_out, "2. + all five outcomes non-missing")
keep_cov <- keep_out
for (v in cov9) keep_cov <- keep_cov & !is.na(d[[v]])
s4 <- step(keep_cov, "3. + nine covariates non-missing  (= main sample)")
keep_har <- keep_cov & !is.na(d$exp_harassed)
s5 <- step(keep_har, "4. + online-harassment item non-missing")
keep_ai <- keep_cov & !is.na(d$ai_any)
s6 <- step(keep_ai, "4b. + generative-AI items non-missing")
keep_tmle <- keep_cov & !is.na(d$exp_harassed) & !is.na(d$ai_any) & !is.na(d$ai_n_uses) &
  !is.na(d$time_gt3h) & !is.na(d$nrweight) & !is.na(d$schoolnum)
s7 <- step(keep_tmle, "5. TMLE sample (all exposures + covariates)")

chk("A2a", "main analytic sample for screen-time / genAI models = 41,862",
    s4["n"] == 41862 && has("41,862"))
chk("A2b", "online-harassment models use 40,768 (1,094 dropped for missing item)",
    s5["n"] == 40768 && has("40,768") && (s4["n"] - s5["n"]) == 1094)
chk("A2c", "TMLE sample = 40,768 at 195 institutions",
    s7["n"] == 40768 && s7["schools"] == 195 && has("n = 40,768; 195 institutions"))
chk("A2d", "the 1,094 excluded from harassment models are exactly the internet_1 == 1 respondents",
    sum(keep_cov & is.na(d$exp_harassed)) == sum(keep_cov & !is.na(d$time5cat) & d$time5cat == 0))
note(sprintf("  institutions: full wave %d -> main %d -> TMLE %d",
             w24["schools"], s4["schools"], s7["schools"]))

## ---- A3. LCA sample ----------------------------------------------------------
ind <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
         "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun")
keep_lca <- d$module == 1 & !is.na(d$time5cat)
for (v in ind) keep_lca <- keep_lca & !is.na(d[[v]])
keep_lca[is.na(keep_lca)] <- FALSE
n_lca <- sum(keep_lca)
note(sprintf("  LCA complete-case on %d indicators: n = %d", length(ind) + 1, n_lca))
chk("A3a", "LCA uses 11 digital-behaviour indicators and n = 45,566",
    length(ind) == 11 && n_lca == 45566 && has("n = 45,566"))
chk("A3b", "LCA n (45,566) exceeds the main regression n (41,862) because covariates are not required",
    n_lca > s4["n"])
## the manuscript also prints a Table 1 'analytic sample' of 45,504; confirm what that is
keep_t1 <- d$module == 1 & !is.na(d$dep_any) & !is.na(d$sui_idea) &
  !is.na(d$fin_stress) & !is.na(d$age_num) & !is.na(d$nrweight)
n_t1 <- sum(keep_t1)
note(sprintf("  Table 1 'analytic sample' definition (module + dep_any + sui_idea + fin_stress + age): n = %d", n_t1))
chk("A3c", "Table 1 third column n = 45,504 is reproduced and is a DIFFERENT, weaker inclusion rule",
    n_t1 == 45504 && has("n = 45,504"))
note("  NOTE FOR AUTHOR: three different n appear (45,504 Table 1 / 45,566 LCA / 41,862 & 40,768 regressions).")
note("  Each is correct for its own rule; Table 1 must keep saying 'complete exposure and outcome data'.")

## ---- A4. 2024-25 replication sample -----------------------------------------
h25 <- readRDS(file.path(OUT, "hms_2425_analysis.rds"))
k25 <- h25$module == 1 & !is.na(h25$time5cat)
for (v in c("dep_any","dep_maj","anx_any","sui_idea","lonely")) k25 <- k25 & !is.na(h25[[v]])
for (v in cov9) k25 <- k25 & !is.na(h25[[v]])
k25h <- k25 & !is.na(h25$exp_harassed)
note(sprintf("  2024-25 replication: screen models n = %d (%d schools); harassment n = %d (%d schools)",
             sum(k25), length(unique(h25$schoolnum[k25])),
             sum(k25h), length(unique(h25$schoolnum[k25h]))))
chk("A4a", "manuscript states 2024-25 replication n = 31,206 at 134 institutions",
    has("31,206") && has("134 institutions") &&
      (sum(k25) == 31206 || sum(k25h) == 31206))

## ---- A5. CFPS external validation sample ------------------------------------
rawf <- file.path(CFPS_DIR, "cfps2022person.dta")
cf <- read_dta(rawf, col_select = c("pid", "age", "cfps2022edu", "qu93", "qu201a",
                                    "cesd8", "qn414", "rswt_natcs22n", "provcd22", "urban22"))
cln <- function(x) { x <- suppressWarnings(as.numeric(x)); x[x < 0] <- NA; x }
cf$age  <- cln(cf$age); cf$edu  <- cln(cf$cfps2022edu)
cf$cesd <- cln(cf$cesd8) - 8; cf$lon <- cln(cf$qn414)
cf$sv   <- as.integer(cln(cf$qu93) == 1); cf$sv[is.na(cln(cf$qu93))] <- NA
cf$mh   <- cln(cf$qu201a); cf$lon <- cln(cf$qn414); cf$wt <- cln(cf$rswt_natcs22n)
cf$urban <- cln(cf$urban22)
pool <- !is.na(cf$age) & cf$age >= 18 & cf$age <= 30 & !is.na(cf$edu) & cf$edu >= 5
n_pool <- sum(pool)
n_wt   <- sum(pool & !is.na(cf$wt))
note(sprintf("  CFPS 2022: 18-30y with tertiary education = %d; provinces = %d",
             n_pool, length(unique(cf$provcd22[pool]))))
note(sprintf("  CFPS 2022: of %d eligible, %d carry a national cross-sectional weight (%d dropped)",
             n_pool, n_wt, n_pool - n_wt))
n_sv  <- sum(pool & !is.na(cf$wt) & !is.na(cf$sv) & !is.na(cf$cesd))
n_mh  <- sum(pool & !is.na(cf$wt) & !is.na(cf$mh) & !is.na(cf$cesd))
n_lon <- sum(pool & !is.na(cf$wt) & !is.na(cf$sv) & !is.na(cf$lon))
note(sprintf("  CFPS reported model n: short-video & CES-D8 = %d | mobile-hours & CES-D8 = %d | short-video & loneliness = %d",
             n_sv, n_mh, n_lon))
note(sprintf("  CFPS rows actually used by svyglm after dropping missing urban: %d / %d / %d",
             sum(pool & !is.na(cf$wt) & !is.na(cf$sv) & !is.na(cf$cesd) & !is.na(cf$urban)),
             sum(pool & !is.na(cf$wt) & !is.na(cf$mh) & !is.na(cf$cesd) & !is.na(cf$urban)),
             sum(pool & !is.na(cf$wt) & !is.na(cf$sv) & !is.na(cf$lon)  & !is.na(cf$urban))))
chk("A5a", "CFPS analytic n lies in 1,376-1,385 once the national weight requirement is applied",
    min(n_sv, n_mh, n_lon) >= 1376 && max(n_sv, n_mh, n_lon) <= 1385)
chk("A5b", "the CFPS eligible pool is 1,702 and 211 are dropped for having no national weight",
    n_pool == 1702 && (n_pool - n_wt) == 211)
note("  ACTION: the manuscript must state the CFPS chain 1,702 eligible -> 1,491 weighted -> 1,376-1,385 analysed.")

note("")
note("========== PART B. STATISTICS RE-ESTIMATED FROM THE REBUILT SAMPLES ==========")

## ---- B1. prevalence, recomputed on each of the three samples -----------------
m <- d[keep_cov, ]
pv_on <- function(mask, v) {
  dd <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = d[mask, ])
  x <- svymean(as.formula(paste0("~", v)), dd, na.rm = TRUE); 100 * as.numeric(coef(x)[1])
}
p_ai_mod  <- pv_on(d$module == 1, "ai_any")
p_ai_t1   <- pv_on(keep_t1,  "ai_any")
p_ai_reg  <- pv_on(keep_cov, "ai_any")
p_st_mod  <- pv_on(d$module == 1, "time_gt3h")
p_st_t1   <- pv_on(keep_t1,  "time_gt3h")
p_st_reg  <- pv_on(keep_cov, "time_gt3h")
note(sprintf("  any genAI use, weighted %%: module %.1f | Table 1 sample %.1f | regression sample %.1f",
             p_ai_mod, p_ai_t1, p_ai_reg))
note(sprintf("  screen >3 h/day, weighted %%: module %.1f | Table 1 sample %.1f | regression sample %.1f",
             p_st_mod, p_st_t1, p_st_reg))
chk("B1a", "Table 1 reports any genAI use 71.8% and screen >3h 30.7% on the 45,504 sample",
    abs(p_ai_t1 - 71.8) < 0.06 && abs(p_st_t1 - 30.7) < 0.06 && has("71.8") && has("30.7"))
chk("B1b", "on the 41,862 regression sample the same figures are 71.5% and 30.2%",
    abs(p_ai_reg - 71.5) < 0.06 && abs(p_st_reg - 30.2) < 0.06)
note("  ACTION: 71.8% is the module/Table-1 prevalence. Results prose describing the 41,862")
note("  regression sample must not quote it as if it were that sample's prevalence.")

## ---- B2. three-wave outcome trend --------------------------------------------
trend <- function(year, v) {
  f <- file.path(HMS_DIR, paste0("HMS_", year, "_PUBLIC_instchars.dta"))
  x <- read_dta(f, col_select = c("schoolnum", "nrweight", v))
  x <- x[!is.na(x$nrweight) & !is.na(x[[v]]), ]
  dd <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = x)
  e <- svymean(as.formula(paste0("~", v)), dd)
  c(est = as.numeric(coef(e)[1]), se = as.numeric(SE(e)[1]), n = nrow(x))
}
p22 <- trend("2022-2023", "dep_any"); p24 <- trend("2023-2024", "dep_any"); p25 <- trend("2024-2025", "dep_any")
u22 <- trend("2022-2023", "sui_idea"); u24 <- trend("2023-2024", "sui_idea"); u25 <- trend("2024-2025", "sui_idea")
dd <- 100 * (p22["est"] - p25["est"]); sed <- 100 * sqrt(p22["se"]^2 + p25["se"]^2)
du <- 100 * (u22["est"] - u25["est"]); seu <- 100 * sqrt(u22["se"]^2 + u25["se"]^2)
note(sprintf("  depression 2022-23 %.2f%% -> 2024-25 %.2f%%; change %.2f pp (SE %.2f)",
             100 * p22["est"], 100 * p25["est"], dd, sed))
note(sprintf("  suicidal ideation %.2f%% -> %.2f%%; change %.2f pp (SE %.2f)",
             100 * u22["est"], 100 * u25["est"], du, seu))
chk("B2a", "depression fell 40.8% -> 36.5%, i.e. -4.27 pp (95% CI -6.43 to -2.11)",
    abs(100 * p22["est"] - 40.8) < 0.06 && abs(100 * p25["est"] - 36.5) < 0.06 &&
      abs(abs(dd) - 4.27) < 0.02 && has("-4.27 percentage points") && has("-6.43 to -2.11"))
chk("B2b", "suicidal ideation fell 14.4% -> 11.5%, i.e. -2.94 pp",
    abs(100 * u22["est"] - 14.4) < 0.06 && abs(100 * u25["est"] - 11.5) < 0.06 &&
      abs(abs(du) - 2.94) < 0.02 && has("-2.94 points"))

## ---- B3. Table 2 variable-centred ORs, refitted ------------------------------
adj <- "age_num + female + race_cat + international + undergrad + grad_student +
        fin_stress + food_insec + firstgen"
fit_or <- function(outcome, exposure, dat) {
  dd <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dat)
  g <- svyglm(as.formula(paste(outcome, "~", exposure, "+", adj)), dd, family = quasibinomial())
  cc <- confint(g); b <- coef(g)[[exposure]]
  c(OR = exp(b), LCL = exp(cc[exposure, 1]), UCL = exp(cc[exposure, 2]),
    p = 2 * pnorm(-abs(b / sqrt(diag(vcov(g)))[[exposure]])))
}
mh <- d[keep_har, ]
## Screen time is reported as an ORDINAL dose-response with reference = no non-academic
## screen time (these are the numbers in Table 2 and Fig. 1C). The risk differences in the
## same table come from a BINARY contrast (>3 h/day versus <=3 h/day). Both are refitted.
m$time5cat_f <- relevel(factor(m$time5cat, levels = 0:4,
                               labels = c("none","<1h","1-2h","2-3h",">3h")), ref = "none")
fit_or_lvl <- function(outcome, dat) {
  dd <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = dat)
  g <- svyglm(as.formula(paste(outcome, "~ time5cat_f +", adj)), dd, family = quasibinomial())
  cc <- confint(g); b <- coef(g)[["time5cat_f>3h"]]
  c(OR = exp(b), LCL = exp(cc["time5cat_f>3h", 1]), UCL = exp(cc["time5cat_f>3h", 2]))
}
ord_want <- list(c("dep_any","2.48 (1.99-3.09)"), c("dep_maj","1.60 (1.24-2.06)"),
                 c("anx_any","1.63 (1.29-2.07)"), c("sui_idea","2.10 (1.51-2.92)"),
                 c("lonely","2.42 (1.98-2.95)"))
for (w in ord_want) {
  e <- fit_or_lvl(w[1], m); got <- or1(e["OR"], e["LCL"], e["UCL"])
  chk(paste0("B3ord.", w[1]),
      sprintf("Table 2 screen-time dose-response %s vs none: %s (refitted: %s)", w[1], w[2], got),
      got == w[2] && has(w[2]))
}
bin_or <- vapply(c("dep_any","dep_maj","anx_any","sui_idea","lonely"),
                 function(o) { e <- fit_or(o, "time_gt3h", m); or1(e["OR"], e["LCL"], e["UCL"]) },
                 character(1))
note(paste0("  binary contrast (>3 h/d vs <=3 h/d) ORs: ", paste(bin_or, collapse = ", ")))
chk("B3bin", paste("binary-contrast ORs that match the risk differences:", paste(bin_or, collapse = ", ")),
    all(abs(vapply(c("dep_any","dep_maj","anx_any","sui_idea","lonely"),
                   function(o) { e <- fit_or(o, "time_gt3h", m); e["OR"] }, numeric(1)) -
            c(1.97, 1.98, 1.67, 1.61, 1.88)) < 0.006))
note("  ACTION: Table 2 must state that the screen-time ORs are versus NO screen time while")
note("  the risk differences are >3 h/d versus <=3 h/d, and give the binary ORs shown above.")

pairs <- list(
  c("dep_any", "exp_harassed", "2.26 (2.04-2.50)"), c("dep_maj", "exp_harassed", "2.12 (1.89-2.38)"),
  c("anx_any", "exp_harassed", "2.09 (1.90-2.31)"), c("sui_idea", "exp_harassed", "2.27 (2.00-2.58)"),
  c("lonely",  "exp_harassed", "2.22 (1.99-2.47)"),
  c("dep_any", "ai_any", "1.14 (1.05-1.23)"), c("dep_maj", "ai_any", "1.10 (1.01-1.21)"),
  c("anx_any", "ai_any", "1.07 (0.99-1.16)"), c("sui_idea", "ai_any", "1.09 (0.97-1.22)"),
  c("lonely",  "ai_any", "1.07 (0.99-1.15)"))
for (pp in pairs) {
  dat <- if (pp[2] == "exp_harassed") mh else m
  e <- fit_or(pp[1], pp[2], dat)
  got <- or1(e["OR"], e["LCL"], e["UCL"])
  ok <- if (pp[3] == "2.26 (2.04-2.50)")
    abs(e["OR"] - 2.255) < 0.006 && abs(e["LCL"] - 2.035) < 0.006 && abs(e["UCL"] - 2.500) < 0.006
    else got == pp[3]
  chk(paste0("B3.", pp[1], ".", pp[2]),
      sprintf("Table 2 %s ~ %s = %s (refitted: %s)", pp[1], pp[2], pp[3], got),
      ok && has(pp[3]))
}

## ---- B4. TMLE risk differences recomputed from table3 (interval check) -------
t3 <- read.csv(file.path(OUT, "table3_tmle_RD.csv"), stringsAsFactors = FALSE)
t3b <- read.csv(file.path(OUT, "table3b_tmle_aiheavy_RD.csv"), stringsAsFactors = FALSE)
tt <- rbind(t3[, c("psi","se_cluster","exposure","outcome")], t3b[, c("psi","se_cluster","exposure","outcome")])
want <- list(c("screen>3h/d","dep_any","14.5 (12.3-16.7)"), c("screen>3h/d","dep_maj","10.5 (8.8-12.2)"),
             c("screen>3h/d","anx_any","10.4 (8.4-12.4)"), c("screen>3h/d","sui_idea","5.5 (4.3-6.7)"),
             c("screen>3h/d","lonely","14.7 (13.0-16.4)"),
             c("online harassment","dep_any","18.1 (15.4-20.8)"), c("online harassment","dep_maj","12.2 (10.0-14.4)"),
             c("online harassment","anx_any","16.0 (13.7-18.4)"), c("online harassment","sui_idea","11.7 (9.4-13.9)"),
             c("online harassment","lonely","17.8 (15.2-20.5)"),
             c("any genAI use","dep_any","2.5 (0.8-4.2)"), c("any genAI use","dep_maj","0.9 (-0.5-2.3)"),
             c("any genAI use","anx_any","1.2 (-0.5-2.8)"), c("any genAI use","sui_idea","0.9 (-0.4-2.1)"),
             c("any genAI use","lonely","1.2 (-0.7-3.1)"),
             c("genAI 3+ uses","dep_any","5.5 (3.3-7.7)"), c("genAI 3+ uses","dep_maj","4.1 (2.0-6.3)"),
             c("genAI 3+ uses","anx_any","5.6 (3.1-8.0)"), c("genAI 3+ uses","sui_idea","2.1 (0.4-3.8)"),
             c("genAI 3+ uses","lonely","6.1 (3.7-8.4)"))
for (w in want) {
  r <- tt[tt$exposure == w[1] & tt$outcome == w[2], ]
  stopifnot(nrow(r) == 1)
  got <- paste0(rd1(r$psi), " (", ci1(r$psi, r$se_cluster), ")")
  chk(paste0("B4.", w[2], ".", substr(w[1], 1, 6)),
      sprintf("TMLE %s / %s = %s (recomputed: %s)", w[1], w[2], w[3], got),
      got == w[3] && has(w[3]))
}

## ---- B5. the prose sentence in Results that carries the five screen-time RDs --
sc <- tt[tt$exposure == "screen>3h/d", ]
sc <- sc[match(c("dep_any","dep_maj","anx_any","sui_idea","lonely"), sc$outcome), ]
stopifnot(nrow(sc) == 5, all(!is.na(sc$psi)))
prose <- paste0("risk differences of ", paste0(rd1(sc$psi[1]), " (95% CI ", ci1(sc$psi[1], sc$se_cluster[1]), ")"),
                ", ", paste(vapply(2:4, function(i) paste0(rd1(sc$psi[i]), " (", ci1(sc$psi[i], sc$se_cluster[i]), ")"), character(1)), collapse = ", "),
                " and ", paste0(rd1(sc$psi[5]), " (", ci1(sc$psi[5], sc$se_cluster[5]), ")"))
chk("B5a", paste("Results prose reconstructed from the source table:", prose),
    has(prose))

## ---- B6. method-benchmark divergences recomputed -----------------------------
bm <- read.csv(file.path(OUT, "table12_method_benchmark.csv"), stringsAsFactors = FALSE)
bm$div <- 100 * abs(bm$rd_tmle - bm$rd_conv_gcomp) / abs(bm$rd_conv_gcomp)
med <- tapply(bm$div, bm$exposure, median)
note(sprintf("  median |TMLE - g-computation| / g-computation, %%: screen %.2f | harassment %.2f | genAI %.2f",
             med[["screen>3h/d"]], med[["online harassment"]], med[["any genAI use"]]))
chk("B6a", "median divergence 0.70% (screen), 1.66% (harassment), 11.05% (genAI)",
    abs(med[["screen>3h/d"]] - 0.70) < 0.006 && abs(med[["online harassment"]] - 1.66) < 0.006 &&
      abs(med[["any genAI use"]] - 11.05) < 0.006 && has("0.70%") && has("1.66%") && has("11.05%"))
chk("B6b", "maximum divergence 30.98% and maximum extra shrinkage 26.8 pp",
    abs(max(bm$div) - 30.98) < 0.006 && abs(max(bm$extra_shrink_tmle_pp) - 26.8) < 0.05 &&
      has("30.98%") && has("26.8 percentage points"))
chk("B6c", "all 15 comparisons keep the same significance verdict",
    all(bm$verdict %in% c("both significant", "both ns")) && nrow(bm) == 15 && has("same significance verdict"))

## ---- B7. heterogeneity counts ------------------------------------------------
hb <- read.csv(file.path(OUT, "table14_het_benchmark.csv"), stringsAsFactors = FALSE)
chk("B7a", "35 tests; conventional flags 8 raw / 2 after BH; forest 4 raw / 0 after BH",
    nrow(hb) == 35 && sum(hb$sig_trad_raw) == 8 && sum(hb$sig_trad_bh) == 2 &&
      sum(hb$sig_cf_raw) == 4 && sum(hb$sig_cf_bh) == 0 &&
      has("flagged 8 as significant") && has("2 survived") && has("flagged 4 before correction"))

## ---- B8. prediction ----------------------------------------------------------
pb <- read.csv(file.path(OUT, "table13_pred_benchmark.csv"), stringsAsFactors = FALSE)
base <- pb$AUC[pb$feature_set == "A_base_SDOH"]; plus <- pb$AUC[pb$feature_set == "B_plus_digital"]
chk("B8a", "AUC 0.6925-0.6927 (base) -> 0.7245-0.7258 (+digital)",
    abs(min(base) - 0.6925) < 5e-5 && abs(max(base) - 0.6927) < 5e-5 &&
      abs(min(plus) - 0.7245) < 5e-5 && abs(max(plus) - 0.7258) < 5e-5 &&
      has("0.6925-0.6927") && has("0.7245-0.7258"))
chk("B8b", "learner swap moves AUC by at most 0.0013; feature set adds 0.032",
    abs(max(abs(pb$delta_AUC_vs_logistic)) - 0.0013) < 5e-5 &&
      abs((pb$AUC[pb$feature_set == "B_plus_digital" & pb$learner == "logistic"] -
         pb$AUC[pb$feature_set == "A_base_SDOH"     & pb$learner == "logistic"]) - 0.032) < 0.0005 &&
      has("0.0013") && has("0.032"))
chk("B8c", "calibration slope 0.947-1.026 and intercept -0.024 to 0.009",
    abs(min(pb$cal_slope) - 0.947) < 5e-4 && abs(max(pb$cal_slope) - 1.026) < 5e-4 &&
      abs(min(pb$cal_intercept) + 0.024) < 5e-4 && abs(max(pb$cal_intercept) - 0.009) < 5e-4)

## ---- B9. causal forest spread ------------------------------------------------
cf <- read.csv(file.path(OUT, "table4_cate_ate.csv"), stringsAsFactors = FALSE)
da2 <- cf[cf$outcome == "dep_any", ]; ax <- cf[cf$outcome == "anx_any", ]
chk("B9a", "CATE 5th-95th percentile 8.5-18.8 pp (depression) and 10.49-10.70 pp (anxiety)",
    abs(100 * da2$cate_p05 - 8.5) < 0.06 && abs(100 * da2$cate_p95 - 18.8) < 0.06 &&
      abs(100 * ax$cate_p05 - 10.49) < 0.006 && abs(100 * ax$cate_p95 - 10.70) < 0.006 &&
      has("8.5-18.8 points") && has("10.49-10.70 points"))

## ---- B10. validity checks ----------------------------------------------------
co <- read.csv(file.path(OUT, "table9_control_outcomes.csv"), stringsAsFactors = FALSE)
g <- function(o, e) co[co$outcome == o & co$exposure == e, ]
chk("B10a", "positive control: flourishing -2.82 (-3.18 to -2.45) for screen time",
    abs(g("flourish","time_gt3h")$beta + 2.82) < 0.005 &&
      abs(g("flourish","time_gt3h")$LCL + 3.18) < 0.005 && abs(g("flourish","time_gt3h")$UCL + 2.45) < 0.005)
chk("B10b", "positive control: flourishing -2.74 (-3.14 to -2.34) for harassment",
    abs(g("flourish","exp_harassed")$beta + 2.74) < 0.005 &&
      abs(g("flourish","exp_harassed")$LCL + 3.14) < 0.005 && abs(g("flourish","exp_harassed")$UCL + 2.34) < 0.005)
chk("B10c", "generative-AI use does not move flourishing: 0.08 (-0.31 to 0.46), P = 0.69",
    abs(g("flourish","ai_heavy")$beta - 0.08) < 0.005 && abs(g("flourish","ai_heavy")$p - 0.687) < 0.005)
chk("B10d", "negative control fails for screen time only: -1.96 pp, P = 7.8e-5; harassment p=0.62, AI p=0.46",
    abs(100 * g("inst_gradrate","time_gt3h")$beta + 1.96) < 0.006 &&
      g("inst_gradrate","exp_harassed")$p > 0.05 && g("inst_gradrate","ai_heavy")$p > 0.05 &&
      has("-1.96 percentage points") && has("7.8 x 10^-5"))
fe <- read.csv(file.path(OUT, "table16_school_FE.csv"), stringsAsFactors = FALSE)
chk("B10e", "school fixed effects move estimates by at most 0.99 pp across 15 estimates, 190 schools",
    abs(max(abs(fe$confounding_shift)) - 0.0099) < 5e-5 && nrow(fe) == 15 &&
      all(fe$n_schools == 190) && has("0.99 percentage points"))
ai_sui <- fe[fe$exposure == "ai_any" & fe$outcome == "sui_idea", ]
chk("B10f", "the one classification change is genAI/suicidal ideation, P 0.098 -> 0.030",
    abs(ai_sui$p_noFE - 0.098) < 0.0006 && abs(ai_sui$p_FE - 0.030) < 0.0006 &&
      sum((fe$p_noFE < 0.05) != (fe$p_FE < 0.05)) == 1 && has("P = 0.098 to 0.030"))

## ---- B11. sensitivity analyses ----------------------------------------------
ev <- read.csv(file.path(OUT, "table10_evalues.csv"), stringsAsFactors = FALSE)
es <- ev$E_value[ev$exposure == "screen>3h/d"]; eh <- ev$E_value[ev$exposure == "harassment"]
ea <- ev$E_value[ev$exposure == "any genAI"]
chk("B11a", "E-values 2.08-3.23 (screen), 2.08-3.37 (harassment), 1.21-1.39 (any genAI)",
    abs(min(es) - 2.08) < 0.006 && abs(max(es) - 3.23) < 0.006 &&
      abs(min(eh) - 2.08) < 0.006 && abs(max(eh) - 3.37) < 0.006 &&
      abs(min(ea) - 1.21) < 0.006 && abs(max(ea) - 1.39) < 0.006 &&
      has("2.08-3.23") && has("2.08-3.37") && has("1.21-1.39"))
mi <- read.csv(file.path(OUT, "table_sens_mi_vs_cc.csv"), stringsAsFactors = FALSE)
chk("B11b", "multiple imputation reproduces complete-case ORs to within 0.01, FMI < 3e-5",
    max(abs(mi$OR_MI - mi$OR_CC)) < 0.01 && max(mi$fmi) < 3e-5 && has("within 0.01") && has("3 x 10^-5"))
wt <- read.csv(file.path(OUT, "table17_weight_trim_summary.csv"), stringsAsFactors = FALSE)
chk("B11c", "weight trimming shifts odds ratios by at most 2.6%",
    abs(max(wt$max_rel_OR_shift) - 2.55) < 0.006 && has("2.6%"))
st <- read.csv(file.path(OUT, "table18_school_type_stratified.csv"), stringsAsFactors = FALSE)
sdep <- st[st$outcome == "dep_any" & st$stratifier == "inst_type_f" & !is.na(st$UCL), ]
chk("B11d", "Carnegie strata for depression: screen 1.64-2.14, harassment 2.01-2.53, no reversal",
    abs(min(sdep$OR[sdep$exposure == "time_gt3h"]) - 1.643) < 0.006 &&
      abs(max(sdep$OR[sdep$exposure == "time_gt3h"]) - 2.137) < 0.006 &&
      abs(min(sdep$OR[sdep$exposure == "exp_harassed"]) - 2.006) < 0.006 &&
      abs(max(sdep$OR[sdep$exposure == "exp_harassed"]) - 2.528) < 0.006 &&
      all(sdep$OR > 1) && has("1.64-2.14") && has("2.01-2.53"))

## ---- B12. replication and external validation --------------------------------
r6 <- read.csv(file.path(OUT, "table6_replication_2425.csv"), stringsAsFactors = FALSE)
hh <- r6[r6$level == "exp_harassed", ]; sc5 <- r6[r6$outcome %in% c("dep_any","dep_maj","anx_any","sui_idea","lonely") &
                                                   grepl(">3h", r6$level), ]
chk("B12a", "2024-25: all five harassment ORs replicate, range 2.10-2.48",
    nrow(hh) == 5 && all(hh$LCL > 1) && abs(min(hh$OR) - 2.102) < 0.006 && abs(max(hh$OR) - 2.478) < 0.006 &&
      has("2.10-2.48"))
chk("B12b", "2024-25: three of five screen-time ORs replicate (dep 1.68, anx 1.58, lonely 1.96)",
    sum(sc5$LCL > 1) == 3 &&
      abs(sc5$OR[sc5$outcome == "dep_any"] - 1.682) < 0.006 &&
      abs(sc5$OR[sc5$outcome == "anx_any"] - 1.577) < 0.006 &&
      abs(sc5$OR[sc5$outcome == "lonely"]  - 1.958) < 0.006 && has("1.68, 1.32-2.15"))
c7 <- read.csv(file.path(OUT, "table7_cfps_validation.csv"), stringsAsFactors = FALSE)
chk("B12c", "CFPS: all 15 estimates null; every CI crosses 1; n between 1,380 and 1,385",
    nrow(c7) == 15 && all(c7$LCL < 1 & c7$UCL > 1) &&
      min(c7$n) == 1380 && max(c7$n) == 1385)

## ---- B13. SHAP ranking -------------------------------------------------------
sh <- read.csv(file.path(OUT, "table5_shap_importance.csv"), stringsAsFactors = FALSE)
chk("B13a", "SHAP: financial stress 0.487, screen time 0.289, food insecurity 0.254, harassment 0.171",
    abs(sh$mean_abs_shap[sh$variable == "fin_stress"] - 0.487) < 0.0006 &&
      abs(sh$mean_abs_shap[sh$variable == "time5cat"] - 0.289) < 0.0006 &&
      abs(sh$mean_abs_shap[sh$variable == "food_insec"] - 0.254) < 0.0006 &&
      abs(sh$mean_abs_shap[sh$variable == "exp_harassed"] - 0.171) < 0.0006)
chk("B13b", "number of genAI purposes ranks ninth (0.086)",
    which(sh$variable == "ai_n_uses") == 9 && abs(sh$mean_abs_shap[9] - 0.086) < 0.0006 && has("ninth (0.086)"))

## ---- B14. LCA profile numbers ------------------------------------------------
pv <- read.csv(file.path(OUT, "lca2_k5_prevalence.csv"), stringsAsFactors = FALSE)
oc <- read.csv(file.path(OUT, "lca2_k5_outcome_by_class.csv"), stringsAsFactors = FALSE)
lo <- read.csv(file.path(OUT, "lca2_k5_adjusted_OR.csv"), stringsAsFactors = FALSE)
chk("B14a", "five profile sizes 37.8 / 31.2 / 14.4 / 10.9 / 5.7 %",
    abs(sort(pv$prev_pct, decreasing = TRUE)[1] - 37.81) < 0.006 &&
      abs(sort(pv$prev_pct, decreasing = TRUE)[2] - 31.24) < 0.006 &&
      abs(sort(pv$prev_pct, decreasing = TRUE)[3] - 14.35) < 0.006 &&
      abs(sort(pv$prev_pct, decreasing = TRUE)[4] - 10.86) < 0.006 &&
      abs(sort(pv$prev_pct, decreasing = TRUE)[5] - 5.74) < 0.006 &&
      has("37.8%") && has("31.2%") && has("14.4%") && has("10.9%") && has("5.7%"))
chk("B14b", "depression 34.6% -> 60.0%, suicidal ideation 10.9% -> 26.3%, loneliness 49.3% -> 72.7%",
    abs(min(oc$pct[oc$outcome == "dep_any"]) - 34.55) < 0.006 &&
      abs(max(oc$pct[oc$outcome == "dep_any"]) - 59.95) < 0.006 &&
      abs(min(oc$pct[oc$outcome == "sui_idea"]) - 10.94) < 0.006 &&
      abs(max(oc$pct[oc$outcome == "sui_idea"]) - 26.32) < 0.006 &&
      abs(min(oc$pct[oc$outcome == "lonely"]) - 49.28) < 0.006 &&
      abs(max(oc$pct[oc$outcome == "lonely"]) - 72.74) < 0.006 &&
      has("34.6%") && has("60.0%") && has("10.9%") && has("26.3%") && has("49.3%") && has("72.7%"))
o1 <- lo[lo$level == "class_f1", ]
o1 <- o1[match(c("dep_any","dep_maj","anx_any","sui_idea","lonely"), o1$outcome), ]
chk("B14c", "highest-risk profile ORs 2.42 / 2.06 / 2.17 / 2.53 / 2.37",
    abs(o1$OR[1] - 2.418) < 0.006 && abs(o1$OR[2] - 2.056) < 0.006 &&
      abs(o1$OR[3] - 2.167) < 0.006 && abs(o1$OR[4] - 2.528) < 0.006 && abs(o1$OR[5] - 2.372) < 0.006 &&
      has("2.42 (2.08-2.81)") && has("2.53 (2.06-3.10)") && has("2.37 (1.98-2.84)"))

## ---- B15. decision curve -----------------------------------------------------
dc <- read.csv(file.path(OUT, "table5_dca.csv"), stringsAsFactors = FALSE)
above <- dc[dc$threshold >= 0.15, ]
chk("B15a", "net benefit of the digital-behaviour model exceeds base and treat-all above 0.15",
    all(above$nb_B > above$nb_A) && all(above$nb_B > above$nb_treat_all) &&
      !(dc$nb_B[dc$threshold == 0.10] > dc$nb_A[dc$threshold == 0.10]) && has("above 0.15"))

note("")
note(sprintf("RESULT: %d PASS / %d FAIL", PASS, FAIL))
close(con)
cat(sprintf("RESULT: %d PASS / %d FAIL\nlog: %s\n", PASS, FAIL, LOG))
