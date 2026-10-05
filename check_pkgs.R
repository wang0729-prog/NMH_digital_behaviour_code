source("00_config.R")
if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
pkgs <- c("survey", "poLCA", "mclust", "dplyr", "tidyr", "ggplot2", "svglite", "tidyLPA", "xgboost", "tidymodels", "grf", "mice")
avail <- sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))
cat(paste(names(avail), avail, sep = "=", collapse = "\n"), "\n")
need <- names(avail)[!avail]
if (length(need) > 0) cat("MISSING:", paste(need, collapse = ", "), "\n") else cat("ALL_PRESENT\n")
