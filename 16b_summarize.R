source("00_config.R")
x <- read.csv(file.path(OUT, "table12_method_benchmark.csv"))
x$rel <- 100 * (x$rd_conv_gcomp - x$rd_tmle) / x$rd_conv_gcomp
for (e in unique(x$exposure)) {
  s <- x[x$exposure == e, ]
  cat(sprintf("%-20s rel-diff range %6.2f%% to %6.2f%%   median %5.2f%%\n",
              e, min(abs(s$rel)), max(abs(s$rel)), median(abs(s$rel))))
}
cat("\nper-row:\n")
print(x[, c("exposure", "outcome", "rd_conv_gcomp", "rd_tmle", "rel")], row.names = FALSE)
