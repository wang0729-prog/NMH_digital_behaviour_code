source("00_config.R")
## 27_central_illustration.R
## Central Illustration / graphical abstract.
##
## WHY THIS IS A SEPARATE FILE AND NOT A MAIN FIGURE:
## The journal's content-type specification states, verbatim, "we do not use
## schemes". A flow/schematic figure is therefore not admissible as a display item, so
## this is produced as a GRAPHICAL ABSTRACT and does not consume any of the 6 display
## items the journal allows. It is also the one artefact that states the study's thesis
## in a single view -- none of the 15 result panels does that.
##
## It is deliberately NOT a pure schematic: every number printed here is taken from the
## analysis outputs, so the figure can be checked rather than merely admired.
##
## Rendered with base grid (no DiagrammeR / graphviz in this library). All text is ASCII
## by assertion: svglite/cairo corrupt non-ASCII glyphs.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(grid) })

dir <- ANALYSIS_DIR
fig <- file.path(dir, "figures")
dir.create(fig, showWarnings = FALSE)

W_MM <- 180; H_MM <- 122
SZ_BAND <- 8; SZ_HEAD <- 7; SZ_BODY <- 6; SZ_NOTE <- 6

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
           "#D55E00", "#CC79A7", "#999999")
COL_DATA <- "#DDDDDD"; COL_EXP <- okabe[2]; COL_PIPE <- okabe[6]
COL_STRONG <- okabe[6]; COL_MID <- okabe[1]; COL_WEAK <- okabe[8]

## ---- everything printed must survive svglite/cairo --------------------------------
assert_ascii <- function(txt, where) {
  for (s in txt)
    if (!is.null(s) && any(as.integer(charToRaw(as.character(s))) > 127L))
      stop("non-ASCII text in ", where, ": ", s)
  invisible(NULL)
}

## ---- low-level helpers (all coordinates in mm, origin = bottom-left) ---------------
mkbox <- function(x, y, w, h, head, body, fill, col_txt = "black", head_pt = SZ_HEAD) {
  list(x = x, y = y, w = w, h = h, head = head, body = body,
       fill = fill, col_txt = col_txt, head_pt = head_pt)
}
draw_box <- function(b) {
  pushViewport(viewport(x = unit(b$x, "mm"), y = unit(b$y, "mm"),
                        just = c("left", "bottom"),
                        width = unit(b$w, "mm"), height = unit(b$h, "mm")))
  grid.rect(gp = gpar(fill = b$fill, col = "black", lwd = 0.4))
  grid.text(b$head, x = unit(1.6, "mm"), y = unit(1, "npc") - unit(1.6, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = b$head_pt, fontface = "bold",
                      fontfamily = "Arial", col = b$col_txt, lineheight = 1.05))
  grid.text(paste(b$body, collapse = "\n"),
            x = unit(1.6, "mm"), y = unit(1, "npc") - unit(1.6 + b$head_pt * 0.42, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = SZ_BODY, fontfamily = "Arial",
                      col = b$col_txt, lineheight = 1.15))
  popViewport()
}
band_label <- function(y, h, txt) {
  pushViewport(viewport(x = unit(0.5, "mm"), y = unit(y, "mm"), just = c("left", "bottom"),
                        width = unit(25, "mm"), height = unit(h, "mm")))
  grid.text(txt, x = unit(0, "mm"), y = unit(0.5, "npc"), just = c("left", "centre"),
            gp = gpar(fontsize = SZ_BAND, fontface = "bold", fontfamily = "Arial"))
  popViewport()
}
arrow_down <- function(x, y1, y2) {
  grid.lines(x = unit(c(x, x), "mm"), y = unit(c(y1, y2), "mm"),
             gp = gpar(lwd = 0.6, col = "grey35"),
             arrow = arrow(length = unit(2.2, "mm"), type = "closed"))
}

## =============================== content ============================================
## Every number below is copied from outputs/ (see the source noted per panel).
boxes <- list()

## Band geometry. The panel grid only has H_MM - 13 = 109 mm of height, so the four
## bands plus three arrow gaps must fit inside that; the first render put the data band
## at y = 96 with h = 20 (top = 116 mm), which silently clipped the box headings off the
## top of the page. Box heights are now sized to their content (1 heading + up to 3
## body lines ~= 13 mm) rather than to a round number.
XD <- 27
H_BOX <- 15
YB <- c(data = 82, exposure = 58, pipeline = 34, findings = 10)

## --- Band 1: data -------------------------------------------------------------------
boxes[[length(boxes) + 1]] <- mkbox(XD, YB[["data"]], 36, H_BOX, "HMS 2022-23",
  c("n = 76,406", "135 schools", "outcomes only"), COL_DATA)
boxes[[length(boxes) + 1]] <- mkbox(XD + 38, YB[["data"]], 36, H_BOX, "HMS 2023-24",
  c("n = 104,729", "digital module 46,914", "analysis 41,862"), COL_DATA)
boxes[[length(boxes) + 1]] <- mkbox(XD + 76, YB[["data"]], 36, H_BOX, "HMS 2024-25",
  c("n = 84,735", "replication 31,206"), COL_DATA)
boxes[[length(boxes) + 1]] <- mkbox(XD + 114, YB[["data"]], 36, H_BOX, "CFPS 2020-22",
  c("panel n = 1,085", "external validation"), COL_DATA)

## --- Band 2: exposure domains --------------------------------------------------------
boxes[[length(boxes) + 1]] <- mkbox(XD, YB[["exposure"]], 48, H_BOX, "Screen time",
  c("non-academic, 7 levels", "internet_1"), COL_EXP)
boxes[[length(boxes) + 1]] <- mkbox(XD + 51, YB[["exposure"]], 48, H_BOX, "Online experiences",
  c("5 items: relation,", "harassment, creation,", "activism, none"), COL_EXP)
boxes[[length(boxes) + 1]] <- mkbox(XD + 102, YB[["exposure"]], 48, H_BOX, "Generative-AI use",
  c("7 purposes", "internet_4_1-7", "72% of module users"), COL_EXP)

## --- Band 3: the AI-driven pipeline (four layers) -------------------------------------
boxes[[length(boxes) + 1]] <- mkbox(XD, YB[["pipeline"]], 36, H_BOX, "1. Construct",
  c("LCA, 11 items, k = 5", "entropy 0.603"), COL_PIPE, "white")
boxes[[length(boxes) + 1]] <- mkbox(XD + 38, YB[["pipeline"]], 36, H_BOX, "2. Causal",
  c("TMLE + SuperLearner", "cross-fitted, school-", "clustered SE"), COL_PIPE, "white")
boxes[[length(boxes) + 1]] <- mkbox(XD + 76, YB[["pipeline"]], 36, H_BOX, "3. Prediction",
  c("XGBoost, SHAP, DCA", "nested CV by school"), COL_PIPE, "white")
boxes[[length(boxes) + 1]] <- mkbox(XD + 114, YB[["pipeline"]], 36, H_BOX, "4. Heterogeneity",
  c("causal forest CATE", "best-linear-predictor", "test + BH"), COL_PIPE, "white")

## --- Band 4: findings, ordered by evidence strength ------------------------------------
boxes[[length(boxes) + 1]] <- mkbox(XD, YB[["findings"]], 48, H_BOX, "Online harassment",
  c("OR 2.1 - 2.3, 5/5 outcomes", "replicated 5/5 in 2024-25"), COL_STRONG, "white")
boxes[[length(boxes) + 1]] <- mkbox(XD + 51, YB[["findings"]], 48, H_BOX, "Screen time > 3 h/d",
  c("OR 1.6 - 2.0, dose-dependent", "3/5 replicated in 2024-25"), COL_MID, "white")
boxes[[length(boxes) + 1]] <- mkbox(XD + 102, YB[["findings"]], 48, H_BOX, "Generative-AI use",
  c("OR 1.07 - 1.14, 2/5 significant", "E-value 1.5 - 1.8"), COL_WEAK, "white")

BANDS <- list(
  list(y = YB[["data"]],     h = H_BOX, lab = "Data"),
  list(y = YB[["exposure"]], h = H_BOX, lab = "Exposure"),
  list(y = YB[["pipeline"]], h = H_BOX, lab = "AI pipeline"),
  list(y = YB[["findings"]], h = H_BOX, lab = "Findings"))

NOTE <- paste(
  "Machine-learning estimators diverge from conventional weighted regression by up to 11% where associations are weak,",
  "but add no discriminative accuracy (maximum cross-validated AUC change 0.0013). Estimates are survey-weighted with",
  "schools as primary sampling units and are adjusted for age, sex/gender, race/ethnicity, international status, degree",
  "level, financial stress, food insecurity and first-generation status.")

## ---- assertions before any device is opened ----------------------------------------
assert_ascii(c(unlist(lapply(boxes, function(b) c(b$head, b$body))),
               unlist(lapply(BANDS, function(b) b$lab)), NOTE), "central illustration")
for (b in boxes) {
  if (nchar(b$head) > 22) stop("box heading too long for its width: ", b$head)
  if (b$x + b$w > W_MM - 2) stop("box overflows canvas: ", b$head)
  if (b$y + b$h > H_MM - 2) stop("box overflows canvas height: ", b$head)
}

## ---- render -------------------------------------------------------------------------
draw <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0, y = unit(13, "mm"), just = c("left", "bottom"),
                        width = unit(W_MM, "mm"),
                        height = unit(H_MM - 13, "mm")))
  for (bd in BANDS) band_label(bd$y, bd$h, bd$lab)
  for (b in boxes) draw_box(b)
 ## flow arrows between bands, centred on the canvas (from the foot of one band's boxes
 ## to the head of the next)
  arrow_down(13.5, YB[["data"]],     YB[["exposure"]] + H_BOX + 1)
  arrow_down(13.5, YB[["exposure"]], YB[["pipeline"]] + H_BOX + 1)
  arrow_down(13.5, YB[["pipeline"]], YB[["findings"]] + H_BOX + 1)
  popViewport()
  pushViewport(viewport(x = unit(0.5, "mm"), y = 0, just = c("left", "bottom"),
                        width = unit(W_MM - 1, "mm"), height = unit(13, "mm")))
  grid.text(paste(strwrap(NOTE, width = 168), collapse = "\n"),
            x = unit(0, "mm"), y = unit(1, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = SZ_NOTE, fontfamily = "Arial",
                      col = "grey25", lineheight = 1.15))
  popViewport()
}

grDevices::cairo_pdf(file.path(fig, "CentralIllustration.pdf"),
                     width = W_MM / 25.4, height = H_MM / 25.4, family = "Arial")
suppressWarnings(draw()); dev.off()
svglite::svglite(file.path(fig, "CentralIllustration.svg"),
                 width = W_MM / 25.4, height = H_MM / 25.4,
                 system_fonts = list(sans = "Arial"))
suppressWarnings(draw()); dev.off()
cat("wrote CentralIllustration\n")
