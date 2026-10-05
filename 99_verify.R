source("00_config.R")
## 99_verify.R -- independent recomputation of headline numbers (second code path)
suppressPackageStartupMessages({ library(survey) })
dir <- ANALYSIS_DIR; out <- file.path(dir, "outputs")
d <- readRDS(file.path(out, "hms_2324_analysis.rds"))

## (1) crude weighted prevalence of dep_any by time_gt3h [compare: 33.7 / 50.4 in Table 1]
m <- d[d$module == 1 & !is.na(d$time5cat) & !is.na(d$dep_any), ]
des <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = m)
tb <- svyby(~dep_any, ~time_gt3h, design = des, FUN = svymean)
cat("CHECK1 crude dep_any by time_gt3h (%):",
    paste(round(100 * as.numeric(tb[[2]]), 1), collapse = " / "), "\n")
cat("CHECK1 crude RD (pp):",
    round(100 * (as.numeric(tb[[2]][2]) - as.numeric(tb[[2]][1])), 2), "\n")

## (2) adjusted OR for >3h on dep_any via a DIFFERENT construction (design subset)
W <- c("age_num","female","race_cat","international","undergrad","grad_student",
       "fin_stress","food_insec","firstgen")
mm <- d[d$module == 1 & !is.na(d$time5cat), ]
for (v in c("dep_any", W)) mm <- mm[!is.na(mm[[v]]), ]
mm$hi <- as.integer(mm$time5cat == 4)          # >3h only, binary, no factor levels
des2 <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = mm)
g <- svyglm(dep_any ~ hi + age_num + female + race_cat + international + undergrad +
              grad_student + fin_stress + food_insec + firstgen,
            design = subset(des2, time5cat %in% c(0, 4)), family = quasibinomial())
ci <- confint(g)
cat("CHECK2 OR(>3h vs none) =", round(exp(coef(g)["hi"]), 3),
    "[", round(exp(ci["hi", 1]), 3), "-", round(exp(ci["hi", 2]), 3), "]", "\n")

## (3) LCA class 1 (lca2, k=5) weighted dep_any prevalence [compare: 59.95]
cl <- readRDS(file.path(out, "lca2", "class_assignment_k5.rds"))
d$cl <- cl
s <- d[!is.na(d$cl) & !is.na(d$dep_any) & !is.na(d$lonely) & !is.na(d$sui_idea), ]
s <- s[!is.na(s$age_num) & !is.na(s$race_cat) & !is.na(s$firstgen), ]
des3 <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = s)
p <- svyby(~dep_any, ~factor(cl), design = des3, FUN = svymean)
cat("CHECK3 class dep_any %:", paste(round(100 * as.numeric(p[[2]]), 2), collapse = " / "), "\n")
q <- svyby(~sui_idea, ~factor(cl), design = des3, FUN = svymean)
cat("CHECK3 class sui_idea %:", paste(round(100 * as.numeric(q[[2]]), 2), collapse = " / "), "\n")

## (4) crude RD for harassment [compare with TMLE rd_unadj = 0.2381]
h <- d[d$module == 1 & !is.na(d$exp_harassed) & !is.na(d$dep_any), ]
des4 <- svydesign(ids = ~schoolnum, weights = ~nrweight, nest = TRUE, data = h)
tb4 <- svyby(~dep_any, ~exp_harassed, design = des4, FUN = svymean)
cat("CHECK4 crude RD harassment (pp):",
    round(100 * (as.numeric(tb4[[2]][2]) - as.numeric(tb4[[2]][1])), 2), "\n")

## (5) unweighted AUC of model B out-of-fold predictions (rank formula, independent path)
## refit quickly on a 50% sample with a single split, unweighted, as a sanity bound
set.seed(7)
ii <- sample(seq_len(nrow(d)), 20000)
sub <- d[ii, ]
sub <- sub[sub$module == 1 & !is.na(sub$time5cat) & !is.na(sub$dep_any) &
             !is.na(sub$ai_n_uses), ]
## NOTE: model.matrix silently drops rows with NA -> build the frame first and
## subset the outcome by the SAME retained rows.
Xs <- model.matrix(~ age_num + female + international + undergrad + grad_student +
                     fin_stress + food_insec + time5cat + exp_harassed + exp_relation +
                     exp_create + exp_activism + ai_n_uses, data = sub)[, -1]
Xs <- Xs[, apply(Xs, 2, sd) > 0, drop = FALSE]
yy <- sub$dep_any[as.integer(rownames(Xs))]
cat("CHECK5 rows kept:", nrow(Xs), "of", nrow(sub), "\n")
tr <- sample(c(TRUE, FALSE), nrow(Xs), replace = TRUE, prob = c(.7, .3))
mo <- suppressWarnings(glm(yy[tr] ~ ., data = as.data.frame(Xs[tr, ]), family = binomial()))
pr <- predict(mo, newdata = as.data.frame(Xs[!tr, ]), type = "response")
y <- yy[!tr]
r <- rank(pr); n1 <- sum(y); n0 <- sum(1 - y)
auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
cat("CHECK5 unweighted GLM AUC (sanity bound, expect 0.65-0.75):", round(auc, 3), "\n")

## (6) CFPS mobile hours: unweighted GLM as sanity check on the weighted OR 1.03
c22 <- readRDS(file.path(out, "cfps_2022_analysis.rds"))
cf <- c22[c22$young_hiedu == 1 & !is.na(c22$mobile_hours) & !is.na(c22$dep8), ]
cf$mh <- cf$mobile_hours / 60
gf <- glm(dep8 ~ mh + age_num + female, data = cf, family = binomial())
cat("CHECK6 CFPS unweighted OR per +1 h:", round(exp(coef(gf)["mh"]), 3),
    " (weighted est 1.03)\n")
cat("CHECK6 CFPS n:", nrow(cf), " | dep8 prevalence:", round(mean(cf$dep8), 3), "\n")
cat("VERIFY DONE\n")
