source("00_config.R")
auc_w <- function(p, y, w) {
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  n1 <- sum(w * y); n0 <- sum(w * (1 - y))
  cw <- cumsum(w * y)
  neg <- (1 - y) * w
  sum(neg * cw) / (n1 * n0)
}
set.seed(1); n <- 5000
x <- rnorm(n); pr <- plogis(x); y <- rbinom(n, 1, pr); w <- runif(n, 0.5, 3)
cat("AUC informative (expect ~0.75):", round(auc_w(pr, y, w), 3), "\n")
cat("AUC random      (expect 0.50):", round(auc_w(runif(n), y, w), 3), "\n")
cat("AUC flipped     (expect 0.25):", round(auc_w(-pr, y, w), 3), "\n")
cat("AUC unweighted check:", round(auc_w(pr, y, rep(1, n)), 3), "\n")
