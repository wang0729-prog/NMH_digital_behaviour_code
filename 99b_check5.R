source("00_config.R")
d <- readRDS(file.path(OUT, "hms_2324_analysis.rds"))
set.seed(7); ii <- sample(seq_len(nrow(d)), 20000); sub <- d[ii, ]
sub <- sub[sub$module == 1 & !is.na(sub$time5cat) & !is.na(sub$dep_any) &
             !is.na(sub$ai_n_uses), ]
rownames(sub) <- NULL   # MUST precede model.matrix: it keeps original (large) row names
Xs <- model.matrix(~ age_num + female + international + undergrad + grad_student +
                     fin_stress + food_insec + time5cat + exp_harassed + exp_relation +
                     exp_create + exp_activism + ai_n_uses, data = sub)[, -1]
Xs <- Xs[, apply(Xs, 2, sd) > 0, drop = FALSE]
yy <- sub$dep_any[as.integer(rownames(Xs))]
cat("rows kept:", nrow(Xs), "| NA in yy:", sum(is.na(yy)), "\n")
tr <- sample(c(TRUE, FALSE), nrow(Xs), replace = TRUE, prob = c(.7, .3))
mo <- suppressWarnings(glm(yy[tr] ~ ., data = as.data.frame(Xs[tr, ]), family = binomial()))
pr <- predict(mo, newdata = as.data.frame(Xs[!tr, ]), type = "response")
y <- yy[!tr]
cat("NA in pr:", sum(is.na(pr)), "| n1:", sum(y), "| n0:", sum(1 - y), "\n")
r <- rank(pr); n1 <- sum(y); n0 <- sum(1 - y)
cat("AUC (rank formula):", round((sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 3), "\n")
