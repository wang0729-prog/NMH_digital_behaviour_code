source("00_config.R")
if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages(library(pdftools))
fig <- FIG_DIR
for (f in list.files(fig, "^(Figure|CentralIllustration).*\\.pdf$")) {
  pdf_convert(file.path(fig, f), dpi = 150,
              filenames = file.path(fig, sub("\\.pdf$", ".png", f)))
  cat("rendered", f, "\n")
}
