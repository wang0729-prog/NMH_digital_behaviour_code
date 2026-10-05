source("00_config.R")
## 31_lca_3step_globalbh.R
## Closes two methodological gaps:
## (a) LCA low entropy (0.603) -> does modal (hard) class assignment bias the class-outcome
## odds ratios? Quantified with a 3-step (R3STEP-style) regression that replaces the hard
## class dummies with the classification posterior probabilities. This is the
## BCH/R3STEP correction of Asparouhov & Muthen (2014) / Vermunt (2010) for
## classification error, and it also yields the correct clustered SEs.
## (b) Multiplicity was described in Methods as "the pooled family of variable-centred, TMLE
## and latent-class tests" but the manuscript reports within-outcome q values. This
## script recomputes BH over the genuinely pooled family so the claim matches the number.
## Neither script reads stored estimates for (a); (b) reads the stored p-value tables only,
## which is legitimate because BH is a deterministic function of the p-value vector.
##
## Usage: Rscript 31_lca_3step_globalbh.R

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(survey); library(dplyr) })

ROOT <- ANALYSIS_DIR
out  <- file.path(ROOT, "outputs")
LOG  <- file.path(ROOT, "31_lca_3step_globalbh_log.txt")
con <- file(LOG, "w")
say <- function(...) { cat(...); cat(..., file = con, append = TRUE) }

K <- 5
ind <- c("time5cat","exp_relation","exp_harassed","exp_create","exp_activism",
         "ai_academic","ai_forbidden","ai_work","ai_comm","ai_health","ai_fun")

d <- readRDS(file.path(out, "hms_2324_analysis.rds"))
m <- readRDS(file.path(out, "lca2", sprintf("fit_k%d.rds", K)))

## ---- rebuild the LCA analysis set in the fitted order ------------------------
keep <- d$module == 1 & !is.na(d$time5cat)
for (v in ind) keep <- keep & !is.na(d[[v]])
keep[is.na(keep)] <- FALSE
idx <- which(keep)
dd  <- d[keep, ]
post <- m$posterior
stopifnot(nrow(dd) == nrow(post))
say("== LCA sample rebuilt ==")
say("n = ", nrow(dd), " | posterior dim = ", paste(dim(post), collapse = " x "))

## ---- classification quality -------------------------------------------------
ent <- function(p) -sum(p * log(p))
H   <- mean(apply(post, 1, ent))                 # mean posterior entropy (nats)
R2  <- exp(-H)                                    # relative entropy, the reported "entropy"
mmp <- mean(apply(post, 1, max))                 # mean max posterior = modal confidence
cls <- apply(post, 1, which.max)
say("\n== classification quality (k = 5) ==")
say("mean posterior entropy H (nats) = ", round(H, 4))
say("relative entropy exp(-H)        = ", round(R2, 4), "   (manuscript quotes 0.603)")
say("mean maximum posterior          = ", round(mmp, 4))
say("share assigned with max post > 0.80 = ", round(100 * mean(apply(post, 1, max) > 0.80), 1), "%")
say("share assigned with max post > 0.90 = ", round(100 * mean(apply(post, 1, max) > 0.90), 1), "%")

## ---- attach both hard classes and posterior probabilities -------------------
d$class <- NA_integer_
d$class[idx] <- cls
for (j in seq_len(K)) d[[paste0("pp", j)]] <- NA_real_
for (j in seq_len(K)) d[[paste0("pp", j)]][idx] <- post[, j]

W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")
oc <- c("dep_any","dep_maj","anx_any","sui_idea","lonely")
s <- d[!is.na(d$class), ]
for (v in c(oc, W)) s <- s[!is.na(s[[v]]), ]
## sui_idea and friends arrive as haven_labelled, which quasibinomial cannot initialise
for (v in oc) s[[v]] <- as.numeric(as.vector(s[[v]]))
s$class_f <- factor(s$class)
say("\nregression n = ", nrow(s))

## largest class is the reference (class 5, 37.8%)
prevtab <- table(s$class_f)
ref <- as.integer(names(sort(prevtab, decreasing = TRUE))[1])
say("reference class = ", ref)
s$class_f <- relevel(s$class_f, ref = as.character(ref))
des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)

adj <- paste(W, collapse = " + ")

## ---- (a) hard assignment vs 3-step -----------------------------------------
say("\n== modal assignment vs 3-step (posterior-probability) regression ==")
say(sprintf("%-9s %-22s %8s %15s %8s %15s %7s",
            "outcome", "contrast", "hard OR", "hard 95% CI", "3step OR", "3step 95% CI", "%atten"))

cmp <- list()
z95 <- 1.959964
mkci <- function(g) {
  b <- coef(g); se <- SE(g); b <- as.numeric(b); se <- as.numeric(se)
  data.frame(term = names(coef(g)), OR = exp(b), lo = exp(b - z95 * se),
             hi = exp(b + z95 * se), p = 2 * pnorm(abs(b / se), lower.tail = FALSE),
             stringsAsFactors = FALSE)
}
for (o in oc) {
  g_hard <- svyglm(as.formula(paste(o, "~ class_f +", adj)), design = des,
                   family = quasibinomial())
  tab_h  <- mkci(g_hard)
  tab_h  <- tab_h[grepl("^class_f", tab_h$term), ]
  tab_h$class <- as.integer(sub("^class_f", "", tab_h$term))
  names(tab_h)[c(2, 3, 4, 5)] <- c("hard_OR", "hard_lo", "hard_hi", "hard_p")
  tab_h <- tab_h[, c("class", "hard_OR", "hard_lo", "hard_hi", "hard_p")]

 ## 3-step: reference-class probability omitted, the other K-1 enter as predictors
  rhs <- paste(sprintf("pp%d", setdiff(seq_len(K), ref)), collapse = " + ")
  g_3   <- svyglm(as.formula(paste(o, "~", rhs, "+", adj)), design = des,
                  family = quasibinomial())
  tab_3 <- mkci(g_3)
  tab_3 <- tab_3[grepl("^pp", tab_3$term), ]
  tab_3$class <- as.integer(sub("^pp", "", tab_3$term))
  names(tab_3)[c(2, 3, 4, 5)] <- c("step3_OR", "step3_lo", "step3_hi", "step3_p")
  tab_3 <- tab_3[, c("class", "step3_OR", "step3_lo", "step3_hi", "step3_p")]

  mm <- merge(tab_h, tab_3, by = "class")
  mm$outcome <- o
  mm$pct_change <- round(100 * (mm$step3_OR - mm$hard_OR) / mm$hard_OR, 1)
  cmp[[o]] <- mm
  for (i in seq_len(nrow(mm))) {
    say(sprintf("%-9s class %d vs %d      %8.3f %6.3f-%6.3f %8.3f %6.3f-%6.3f %+6.1f%%",
                o, mm$class[i], ref, mm$hard_OR[i], mm$hard_lo[i], mm$hard_hi[i],
                mm$step3_OR[i], mm$step3_lo[i], mm$step3_hi[i], mm$pct_change[i]))
  }
}
cmp <- do.call(rbind, cmp)
cmp <- cmp[order(cmp$outcome, cmp$class), ]

## the high-risk class is the one with the largest dep_any OR
hr <- cmp$class[which.max(cmp$hard_OR[cmp$outcome == "dep_any"])]
say("\nhigh-risk class (largest dep_any OR) = class ", hr)
say("3-step vs hard OR change for that class, across the five outcomes: ",
    paste(sprintf("%+.1f%%", cmp$pct_change[cmp$class == hr]), collapse = ", "))
say("all 20 contrasts move by at most ", round(max(abs(cmp$pct_change)), 1), "%")
sig_hard <- sum(cmp$hard_lo > 1 | cmp$hard_hi < 1)
sig_3    <- sum(cmp$step3_lo > 1 | cmp$step3_hi < 1)
say("contrasts excluding the null: hard ", sig_hard, " of 20 | 3-step ", sig_3, " of 20")

write.csv(cmp, file.path(out, "lca2_k5_3step_vs_hard_OR.csv"), row.names = FALSE)

## ---- (b) pooled-family multiplicity ----------------------------------------
say("\n== pooled-family Benjamini-Hochberg ==")
## The pooled family is every primary inferential test the paper reports, counted once:
## 35 variable-centred = 7 exposure contrasts x 5 outcomes
## (4 screen-time dose-response categories vs no use, online
## harassment, any generative-AI use, AI 3+ purposes)
## 20 TMLE = 4 exposures x 5 outcomes
## 20 latent class = 4 non-reference profiles x 5 outcomes
## Winsorised-weight duplicates are robustness copies of the same hypotheses and are
## deliberately NOT counted again.
vc <- read.csv(file.path(out, "table2_var_centered_OR.csv"), stringsAsFactors = FALSE)
vc <- vc[grepl("\\[raw w\\]", vc$exposure), ]
lv <- c("time5cat_f<1h", "time5cat_f1-2h", "time5cat_f2-3h", "time5cat_f>3h",
        "exp_harassed", "ai_any", "ai_n_uses_g3+")
vc2 <- vc[vc$level %in% lv, ]
stopifnot(nrow(vc2) == 35)
fam <- list()
fam[["variable-centred"]] <- data.frame(
  test = paste(vc2$outcome, vc2$exposure, vc2$level), p = vc2$p)
tm1 <- read.csv(file.path(out, "table3_tmle_RD.csv"), stringsAsFactors = FALSE)
tm2 <- read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv"), stringsAsFactors = FALSE)
tm  <- rbind(tm1, tm2)
fam[["TMLE"]] <- data.frame(test = paste(tm$exposure, tm$outcome), p = tm$p_cluster)
## latent-class p-values: the hard-assignment ORs are what the manuscript reports, so the
## pooled family must use those. Recovered from the stored OR and 95% CI.
lca_or <- read.csv(file.path(out, "lca2_k5_adjusted_OR.csv"), stringsAsFactors = FALSE)
lca_se <- (log(lca_or$UCL) - log(lca_or$LCL)) / (2 * 1.959964)
lca_p  <- 2 * pnorm(abs(log(lca_or$OR) / lca_se), lower.tail = FALSE)
fam[["latent class"]] <- data.frame(
  test = paste(lca_or$outcome, lca_or$level), p = lca_p)

pooled <- do.call(rbind, lapply(names(fam), function(n) {
  z <- fam[[n]]; z$family <- n; z$test <- z$test; z[, c("family", "test", "p")]
}))
pooled <- pooled[is.finite(pooled$p), ]
pooled$q_pooled <- p.adjust(pooled$p, "BH")
pooled$q_within <- ave(pooled$p, pooled$family, FUN = function(x) p.adjust(x, "BH"))
pooled$sig_raw     <- pooled$p < 0.05
pooled$sig_pooled  <- pooled$q_pooled < 0.05
pooled$sig_within  <- pooled$q_within  < 0.05
pooled <- pooled[order(pooled$p), ]
say("tests in pooled family: ", nrow(pooled))
for (f in unique(pooled$family))
  say(sprintf("  %-18s n=%2d  raw p<0.05: %2d  BH within family: %2d  BH pooled: %2d",
              f, sum(pooled$family == f), sum(pooled$sig_raw[pooled$family == f]),
              sum(pooled$sig_within[pooled$family == f]),
              sum(pooled$sig_pooled[pooled$family == f])))
say(sprintf("  %-18s n=%2d  raw p<0.05: %2d  BH within family: %2d  BH pooled: %2d",
            "TOTAL", nrow(pooled), sum(pooled$sig_raw), sum(pooled$sig_within),
            sum(pooled$sig_pooled)))
flip <- pooled[pooled$sig_within != pooled$sig_pooled, ]
say("tests that change verdict when the family is pooled rather than split by family: ",
    nrow(flip))
if (nrow(flip)) for (i in seq_len(nrow(flip)))
  say("   ", flip$family[i], "|", flip$test[i], "| q_within ", signif(flip$q_within[i], 3),
      "| q_pooled ", signif(flip$q_pooled[i], 3))
write.csv(pooled, file.path(out, "table19_pooled_family_bh.csv"), row.names = FALSE)

close(con)
cat("log:", LOG, "\n")
