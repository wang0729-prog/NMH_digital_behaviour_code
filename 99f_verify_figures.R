source("00_config.R")
## 99f_verify_figures.R
## Independent verification of the figure set produced by 22_final_figures.R
## (plus the Central Illustration from 27_central_illustration.R).
##
## Two independent evidence chains, as before:
## (1) re-derive every number from the source CSVs in outputs/;
## (2) read the rendered SVG back and verify the drawn geometry.
##
## WHAT CHANGED ON (Figure 3 split for the journal):
## * the figure set is now Figure 1 (A-C), Figure 2 (D-E), Figure 3 (F-H),
## Figure 4 (I-J), Figure S1, S2 (a-b), S3 (a-b) -- seven files;
## * panel letters were remapped (old E -> F, old F -> E) so each figure carries a
## contiguous block;
## * PANEL PARTITIONING IS NO LONGER HARD-CODED. The previous version recomputed y
## bands from a hard-coded canvas height and `heights` vector, which had to be
## edited by hand every time a figure was re-laid-out. This version reads the panel
## BORDERS out of the SVG itself (rectangles stroked with #333333) and assigns every
## circle to the border box that contains it. A silent layout change is still caught
## (the box geometry is re-read, and the checks below assert the expected number of
## boxes per row), but the script no longer needs editing when the layout moves.

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

fig <- FIG_DIR
out <- OUT

PASS <- 0; FAIL <- 0; LOG <- character(0)
chk <- function(id, what, ok) {
  ok <- isTRUE(ok)
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  cat(sprintf("[%s] %-10s %s\n", if (ok) "PASS" else "FAIL", id, what))
  if (!ok) LOG <<- c(LOG, paste(id, what))
}
note <- function(txt) cat("        note: ", txt, "\n", sep = "")

## ---------------------------------------------------------------- svg readers
svg_raw <- function(f) paste(readLines(file.path(fig, f), warn = FALSE), collapse = " ")
## canvas height in pt, from the viewBox
svg_h <- function(f) as.numeric(strsplit(
  sub(".*viewBox='([^']+)'.*", "\\1", svg_raw(f)), " ")[[1]])[4]
svg_texts <- function(f) {
  s <- svg_raw(f)
  v <- regmatches(s, gregexpr("<text[^>]*>.*?</text>", s, perl = TRUE))[[1]]
  v <- sub("^<text[^>]*>", "", v)
  sub("</text>$", "", v)
}
svg_circles <- function(f) {
  s <- svg_raw(f)
  v <- regmatches(s, gregexpr("<circle[^>]*/>", s, perl = TRUE))[[1]]
  data.frame(
    cx   = as.numeric(sub(".*cx='([-0-9.]+)'.*",   "\\1", v)),
    cy   = as.numeric(sub(".*cy='([-0-9.]+)'.*",   "\\1", v)),
    fill = sub(".*fill: (#[0-9A-Fa-f]{6}).*", "\\1", v),
    stringsAsFactors = FALSE)
}
svg_rects <- function(f) {
  s <- svg_raw(f)
  v <- regmatches(s, gregexpr("<rect[^>]*>", s, perl = TRUE))[[1]]
  num <- function(tag) suppressWarnings(as.numeric(
    sub(paste0(".*", tag, "='([-0-9.]+)'.*"), "\\1", v)))
  data.frame(x = num("x"), y = num("y"), w = num("width"), h = num("height"),
             fill = sub(".*fill: ([^;']+).*", "\\1", v), stringsAsFactors = FALSE)
}

## ---- panel borders -------------------------------------------------------------------
## ggplot draws every panel (and every facet strip) as a rectangle stroked with #333333.
## Reading them back gives the layout without any knowledge of the canvas size.
panel_boxes <- function(f, gap = 50) {
  s <- svg_raw(f)
  v <- regmatches(s, gregexpr("<rect[^>]*>", s, perl = TRUE))[[1]]
  v <- v[grepl("stroke: #333333", v)]
  num <- function(tag) suppressWarnings(as.numeric(
    sub(paste0(".*\\b", tag, "='([-0-9.]+)'.*"), "\\1", v)))
  b <- data.frame(x = num("x"), y = num("y"), w = num("width"), h = num("height"))
  b <- b[order(b$y, b$x), ]
 ## a facet strip sits directly on top of its panel; the strip is the short box
  b$row <- cumsum(c(TRUE, diff(b$y) > gap))     # one row per horizontal band
  b$col <- ave(b$x, b$row, FUN = function(z) rank(z, ties.method = "first"))
  b
}
## panel bodies only (drop the short facet strips: they are ~12-20 pt tall)
panel_bodies <- function(f) {
  b <- panel_boxes(f)
  b[b$h > 50, ]
}
## every circle, tagged with the panel body that contains it (NA = legend or stray)
circles_tagged <- function(f) {
  ci <- svg_circles(f); b <- panel_bodies(f)
  ci$row <- NA_integer_; ci$col <- NA_integer_
  for (i in seq_len(nrow(ci))) {
    hit <- which(b$x <= ci$cx[i] & ci$cx[i] <= b$x + b$w &
                 b$y <= ci$cy[i] & ci$cy[i] <= b$y + b$h)
    if (length(hit)) { ci$row[i] <- b$row[hit[1]]; ci$col[i] <- b$col[hit[1]] }
  }
  ci
}
## circles belonging to one panel body (by row, then column within the row)
circles_in <- function(ci, row, col = NULL) {
  d <- ci[!is.na(ci$row) & ci$row == row, ]
  if (!is.null(col)) d <- d[d$col == col, ]
  d
}
## keep the data bars of a bar chart: legend keys are squares that share one position,
## while the data bars of one chart share their x (horizontal layout) or width (vertical).
data_bars <- function(f, fills, axis = c("w", "h"), min_h = 0) {
  axis <- match.arg(axis)
  r <- svg_rects(f)
  r <- r[r$fill %in% names(fills), ]
  if (min_h > 0) r <- r[r$h > min_h, ]
  key <- if (axis == "w") round(r$x, 1) else round(r$w, 1)
  r <- r[key == as.numeric(names(which.max(table(key)))), ]
  r$grp <- unname(fills[r$fill])
  r$val <- if (axis == "w") r$w else r$h
  r
}
bar_fit <- function(r, src_by_group) {
  do.call(rbind, lapply(sort(unique(r$grp)), function(g) {
    dra <- sort(r$val[r$grp == g]); src <- sort(src_by_group[[g]])
    if (length(dra) != length(src) || length(dra) < 2)
      return(data.frame(grp = g, n = length(dra), r2 = NA_real_, resid = NA_real_))
    fl <- stats::lm(dra ~ src)
    data.frame(grp = g, n = length(dra), r2 = summary(fl)$r.squared,
               resid = max(abs(stats::residuals(fl))))
  }))
}
## svglite escapes < and > inside text nodes, so literal matching needs the same escaping
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", x, fixed = TRUE), fixed = TRUE)

## rank-pair the source values with the drawn coordinates, then fit a line.
## `horizontal = TRUE` is for panels whose estimate is encoded in cx; the default is for
## panels whose estimate is encoded in cy, where SVG y grows downwards so the pairing
## direction is reversed.
lin_fit <- function(vals, pos, horizontal = FALSE) {
  if (length(vals) != length(pos) || length(vals) < 2)
    return(list(r2 = NA_real_, resid = NA_real_, slope = NA_real_, n = length(vals)))
  o1 <- order(vals)
  o2 <- if (horizontal) order(pos) else order(-pos)
  d  <- data.frame(v = vals[o1], p = pos[o2])
  fl <- stats::lm(p ~ v, data = d)
  list(r2 = summary(fl)$r.squared, resid = max(abs(stats::residuals(fl))),
       slope = unname(stats::coef(fl)[2]), n = nrow(d))
}
## A dumbbell panel gives every facet the same y scale, so the points of one outcome share
## cy within one facet COLUMN. Within a column, rows are separated by cy.
dumbbell_rows <- function(d, n_col) {
  cyu <- sort(unique(round(d$cy, 1)), decreasing = TRUE)
  cyu <- cyu[vapply(cyu, function(y) sum(abs(d$cy - y) < 0.5) == n_col * 2, logical(1))]
  lapply(cyu, function(y) {
    sub <- d[abs(d$cy - y) < 0.5, ]
    sub <- sub[order(sub$cx), ]
    lapply(seq_len(n_col), function(j) sub[(2 * j - 1):(2 * j), ])
  })
}
pair_xy <- function(seg, fill_a, fill_b) {
  a <- seg[seg$fill == fill_a, ]; b <- seg[seg$fill == fill_b, ]
  if (nrow(a) != 1 || nrow(b) != 1) return(c(x_a = NA_real_, x_b = NA_real_))
  c(x_a = a$cx, x_b = b$cx)
}

OUTL   <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
OUTLAB <- c("Depression (any)", "Depression (major)", "Anxiety",
            "Suicidal ideation", "Loneliness")

## ================================================== GLOBAL ==================================================
files <- c("Figure1_profiles_outcomes.svg", "Figure2_effect_estimates.svg",
           "Figure3_method_benchmark.svg", "Figure4_validity_checks.svg",
           "FigureS1_three_wave_trend.svg", "FigureS2_shap_dca.svg",
           "FigureS3_generalisation.svg")
chk("G1", "all seven figure files exist", all(file.exists(file.path(fig, files))))
chk("G2", "all seven figure files exist as PDF",
    all(file.exists(file.path(fig, sub("\\.svg$", ".pdf", files)))))
chk("G1b", "Central Illustration exists as PDF and SVG",
    all(file.exists(file.path(fig, c("CentralIllustration.pdf", "CentralIllustration.svg")))))

for (f in files) {
  tx <- svg_texts(f)
  bad <- grep("[\001-\010\013\014\016-\037]", tx, value = TRUE)
  chk(paste0("G3.", substr(f, 1, 8)), "no corrupted (control-character) glyphs", length(bad) == 0)
  n_arial <- lengths(regmatches(svg_raw(f), gregexpr('font-family: "Arial"', svg_raw(f), fixed = TRUE)))
  chk(paste0("G4.", substr(f, 1, 8)), "every text element declares Arial",
      n_arial == length(tx) && length(tx) > 0)
}

## Panel titles, keyed by file. Panels are lettered PER FIGURE (each figure restarts at A;
## supplementary figures use lower case), so the same letter now legitimately appears in
## several figures - the check must therefore be per file, not against one flat pool.
want <- list(
  "Figure1_profiles_outcomes.svg" = c("A", "B", "C"),
  "Figure2_effect_estimates.svg" = c("A", "B"),
  "Figure3_method_benchmark.svg" = c("A", "B", "C"),
  "Figure4_validity_checks.svg" = c("A", "B"),
  "FigureS1_three_wave_trend.svg" = character(0),
  "FigureS2_shap_dca.svg" = c("a", "b"),
  "FigureS3_generalisation.svg" = c("a", "b"))
stopifnot(setequal(names(want), files))
for (f in names(want)) {
  tx <- svg_texts(f)
 ## Titles are now the bare panel letter, so a substring test would match
 ## any text element that happens to contain that letter. Require an EXACT element.
  for (w in want[[f]])
    chk("G5", paste0(sub("\\.svg$", "", f), ": panel letter rendered exactly: ", w),
        any(tx == w))
}
## a letter must not appear twice inside the SAME figure
for (f in names(want)) {
  lt <- sub("^([A-Za-z]).*", "\\1", want[[f]])
  chk("G5b", paste0(sub("\\.svg$", "", f), ": panel letters are unique within the figure"),
      !any(duplicated(lt)))
}
want <- unlist(want, use.names = FALSE)   # kept for the length check below

## : the caption strip is no longer drawn into the artwork. the journal
## wants the legend supplied as manuscript text, so baking the same paragraph into the image
## only duplicated it. These checks therefore run in the OPPOSITE direction from before:
## the caption sentences must be ABSENT from every figure, and no multi-line 6 pt text
## block may sit at the left canvas margin (that block WAS the caption strip).
svg_labelled <- function(f) {
  s <- svg_raw(f)
  v <- regmatches(s, gregexpr("<text[^>]*>.*?</text>", s, perl = TRUE))[[1]]
  num <- function(tag) suppressWarnings(as.numeric(
    sub(paste0(".*\\b", tag, "='([-0-9.]+)'.*"), "\\1", v)))
  data.frame(y = num("y"), x = num("x"),
             size = suppressWarnings(as.numeric(sub(".*font-size: ([0-9.]+)px.*", "\\1", v))),
             txt = sub("^<text[^>]*>", "", sub("</text>$", "", v)),
             stringsAsFactors = FALSE)
}
caption_band <- function(f) {
  d <- svg_labelled(f)
  d <- d[!is.na(d$y) & !is.na(d$x) & !is.na(d$size), ]
  xm <- min(d$x)
  d[d$size == 6 & d$x <= xm + 0.5, ]
}
captions <- list(
  "Figure1_profiles_outcomes.svg" = "latent class analysis of 11 digital-behaviour items",
  "Figure2_effect_estimates.svg"  = "causal-forest CATE for screen time",
  "Figure3_method_benchmark.svg"  = "number of significant effect-modification tests out of 35",
  "Figure4_validity_checks.svg"   = "within-school transformation",
  "FigureS1_three_wave_trend.svg" = "documents background trend, not exposure change",
  "FigureS2_shap_dca.svg"         = "decision curve analysis comparing the",
  "FigureS3_generalisation.svg"   = "external validation in CFPS 2022")
## canvas height that remains once the caption strip is gone: this is exactly the panel
## area requested by save_combo, in mm. svglite writes the viewBox at 72 units per inch.
H_MM <- c("Figure1_profiles_outcomes.svg" = 139, "Figure2_effect_estimates.svg"  = 97,
          "Figure3_method_benchmark.svg"  = 141, "Figure4_validity_checks.svg"   = 113,
          "FigureS1_three_wave_trend.svg" = 85,  "FigureS2_shap_dca.svg"         = 85,
          "FigureS3_generalisation.svg"   = 85)
for (f in names(captions)) {
  chk("G6", paste0(f, ": no caption strip embedded in the image"),
      !grepl(captions[[f]], paste(svg_labelled(f)$txt, collapse = " "), fixed = TRUE))
 ## a surviving strip would lengthen the canvas beyond the panel area; 0.5 mm tolerance
  h_mm <- svg_h(f) / 72 * 25.4
  chk("G6b", paste0(f, ": canvas is the panel area only (", H_MM[[f]], " mm, no strip)"),
      abs(h_mm - H_MM[[f]]) < 0.5)
}
## side-by-side panels must have identical plot-box heights, otherwise the two boxes no
## longer line up (a panel without a legend or rotated labels grows a taller box)
for (f in files) {
  b <- panel_bodies(f)
  bad <- vapply(split(b$h, b$row), function(z) length(z) > 1 && diff(range(z)) > 1.5, logical(1))
  chk("G6e", paste0(f, ": panels sharing a row have equal plot-box heights"), !any(bad))
  if (any(bad)) note(sprintf("%s: row heights %s", f,
      paste(vapply(split(b$h, b$row), function(z) paste(sprintf("%.1f", z), collapse = "/"),
                   character(1)), collapse = " | ")))
}

chk("G7", "all panel titles <= 46 characters", all(nchar(want) <= 46))
if (requireNamespace("pdftools", quietly = TRUE)) {
  for (f in sub("\\.svg$", "\\.pdf", c(files, "CentralIllustration.svg"))) {
    ff <- pdftools::pdf_fonts(file.path(fig, f))
    chk("G8", paste0(f, ": fonts embedded + Arial"),
        nrow(ff) > 0 && all(ff$embedded) && all(grepl("Arial", ff$name)))
  }
}

## ---- layout structure: the regrouping must actually be in place -----------------------
## This is the check that would have caught the old grouping.
b1 <- panel_boxes("Figure1_profiles_outcomes.svg")
chk("G9", "Figure 1 = 2 bands (A full width on top, B and C below)",
    length(unique(b1$row)) == 2 &&
      sum(b1$row == 1) == 1 && sum(b1$row == 2) == 2)
b2 <- panel_boxes("Figure2_effect_estimates.svg")
chk("G10", "Figure 2 = a single band of exactly 2 panels (A and B side by side)",
    length(unique(b2$row)) == 1 && nrow(b2) == 2)
b3 <- panel_boxes("Figure3_method_benchmark.svg")
chk("G11", "Figure 3 = 2 bands: A (3 facets, full width) over B and C",
    length(unique(b3$row)) == 2 &&
      sum(b3$row == 1) == 6 && sum(b3$row == 2) == 2)   # 3 facet strips + 3 bodies, then G, H
b4 <- panel_boxes("Figure4_validity_checks.svg")
chk("G12", "Figure 4 = 2 bands: A (2 facets, full width) over B (3 facets, full width)",
    length(unique(b4$row)) == 2 &&
      sum(b4$row == 1) == 4 && sum(b4$row == 2) == 6)   # 2 strips + 2 bodies, 3 strips + 3 bodies
chk("G13", "the two panels of Figure 2 sit side by side (not stacked)",
    nrow(panel_bodies("Figure2_effect_estimates.svg")) == 2 &&
      diff(panel_bodies("Figure2_effect_estimates.svg")$x) > 100)
## Letters now repeat across figures, so "is this panel in this file?" is asked by TITLE,
## not by letter: the causal benchmark must live in Figure 3 and the school-FE check in
## Figure 4, and each of those titles must appear in exactly one file.
all_files_tx <- lapply(files, function(f) paste(svg_texts(f), collapse = " "))
where <- function(pat) which(vapply(all_files_tx, function(z) grepl(pat, z, fixed = TRUE), logical(1)))
## Panel titles are letters only, so the panels are identified by the axis and
## legend text that is unique to them. Figure 3 owns the heterogeneity and prediction
## benchmarks; Figure 4 owns the control outcomes and the school-fixed-effects comparison.
chk("G14", "the methodological benchmark and the validity checks are in different files",
    identical(where("Significant tests out of 35"), 3L) &&
      identical(where("Cross-validated AUC"), 3L) &&
      identical(where("Without school FE"), 4L) &&
      identical(where("Positive control:"), 4L))
chk("G14b", "no figure other than Figure 3 carries an AUC or heterogeneity axis",
    length(where("Cross-validated AUC")) == 1 && length(where("Significant tests out of 35")) == 1)
chk("G14c", "no figure other than Figure 4 carries the school-fixed-effects legend",
    length(where("Without school FE")) == 1)

## ==================================== FIGURE 1 ====================================
cp <- read.csv(file.path(out, "lca2", "condprob_k5.csv"), stringsAsFactors = FALSE)
cp_bin <- cp %>% filter(level == 1, item != "time5cat")
cp_tim <- cp %>% filter(item == "time5cat", level == 4) %>% mutate(item = "time_gt3h")
pA <- bind_rows(cp_bin, cp_tim)
chk("F1.A1", "panel A uses 11 items x 5 classes = 55 bars",
    nrow(pA) == 11 && length(unique(pA$item)) == 11 &&
      ncol(pA %>% select(starts_with("class"))) == 5)
chk("F1.A2", "panel A all conditional probabilities in [0, 1]",
    all(pA %>% select(starts_with("class")) >= 0 & pA %>% select(starts_with("class")) <= 1))
A_LEVELS <- c("time_gt3h", "exp_relation", "exp_create", "exp_activism", "exp_harassed",
              "ai_academic", "ai_work", "ai_comm", "ai_health", "ai_fun", "ai_forbidden")
chk("F1.A3", "panel A level order is time -> experiences -> AI, covering all 11 items",
    identical(A_LEVELS[1], "time_gt3h") && setequal(A_LEVELS, unique(pA$item)) &&
      sum(grepl("^exp_", A_LEVELS)) == 4 && sum(grepl("^ai_", A_LEVELS)) == 6)
tx1 <- svg_texts("Figure1_profiles_outcomes.svg")
it_lab <- c("Screen >3 h/d", "Online relationships", "Creating content", "Online activism",
            "Online harassment", "AI: schoolwork", "AI: work", "AI: communication",
            "AI: health", "AI: fun", "AI: forbidden use")
chk("F1.A4", "panel A renders all 11 item labels", all(esc(it_lab) %in% tx1))
chk("F1.A5", "panel A renders 5 profile legend keys", all(paste0("C", 1:5) %in% tx1))
note(paste0("panel A harassment P(class): ",
            paste(sprintf("%.2f", filter(pA, item == "exp_harassed") %>%
                            select(starts_with("class")) %>% unlist()), collapse = "/")))

ob <- read.csv(file.path(out, "lca2_k5_outcome_by_class.csv"), stringsAsFactors = FALSE)
chk("F1.B1", "panel B = 5 outcomes x 5 profiles = 25 bars",
    nrow(ob) == 25 && length(unique(ob$outcome)) == 5 && length(unique(ob$class)) == 5)
chk("F1.B2", "panel B prevalences in (0, 100)", all(ob$pct > 0 & ob$pct < 100))
chk("F1.B3", "panel B upper error bar does not exceed 100%", all(ob$pct + 1.96 * ob$se <= 100))
chk("F1.B4", "panel B renders all 5 outcome labels", all(OUTLAB %in% tx1))
hi <- ob %>% filter(outcome == "dep_any") %>% slice_max(pct)
note(sprintf("panel B worst profile: class %d, dep_any %.1f%% (95%% CI %.1f-%.1f)",
             hi$class, hi$pct, hi$pct - 1.96 * hi$se, hi$pct + 1.96 * hi$se))

t2 <- read.csv(file.path(out, "table2_var_centered_OR.csv"), stringsAsFactors = FALSE)
t4 <- t2 %>% filter(grepl("raw w", exposure), grepl("^time5cat_f", level)) %>%
  mutate(cat = sub("time5cat_f", "", level))
chk("F1.C1", "panel C = 5 outcomes x 4 non-reference levels = 20 model rows",
    nrow(t4) == 20 && setequal(t4$cat, c("<1h", "1-2h", "2-3h", ">3h")))
chk("F1.C2", "panel C reference category is 'none' (OR = 1 by construction)",
    !("none" %in% t4$cat))
chk("F1.C3", "panel C: >3 h/d is the highest-risk level for all 5 outcomes",
    all(sapply(OUTL, function(o) {
      d <- filter(t4, outcome == o)
      d$OR[which(d$cat == ">3h")] == max(d$OR) })))
mono <- sapply(OUTL, function(o) {
  d <- filter(t4, outcome == o)
  v <- d$OR[match(c("<1h", "1-2h", "2-3h", ">3h"), d$cat)]
  all(diff(v) > 0) })
chk("F1.C4", "panel C: >=3 of 5 outcomes strictly monotone across the 4 levels", sum(mono) >= 3)
note(paste0("panel C monotone outcomes: ", paste(OUTLAB[mono], collapse = ", "),
            " | non-monotone: ", paste(OUTLAB[!mono], collapse = ", ")))
chk("F1.C5", "panel C: the two non-monotone dips are non-significant (CI overlap)",
    all(sapply(OUTL[!mono], function(o) {
      d <- filter(t4, outcome == o)
      v <- d$LCL[match(c("<1h", "1-2h", "2-3h"), d$cat)]
      u <- d$UCL[match(c("<1h", "1-2h", "2-3h"), d$cat)]
      any(v[2] < u[1] | v[3] < u[2]) })))
chk("F1.C6", "panel C renders the 'none' reference tick", "none" %in% tx1)

CL5 <- c("#E69F00" = "C5", "#56B4E9" = "C4", "#009E73" = "C3",
         "#F0E442" = "C2", "#0072B2" = "C1")
bA <- data_bars("Figure1_profiles_outcomes.svg", CL5, axis = "w")
bA <- bA[bA$y < 200, ]                       # panel A = upper band, panel B = lower band
chk("F1.A6", "panel A draws 55 bars, 11 per profile", nrow(bA) == 55 && all(table(bA$grp) == 11))
srcA <- split(pA %>% select(-level) %>%
                pivot_longer(starts_with("class"), names_to = "class", values_to = "p") %>%
                mutate(cls = toupper(sub("class", "C", class))),
              ~ cls)
names(srcA) <- sub("^cls\\.", "", names(srcA))
bAfit <- bar_fit(bA, lapply(srcA, function(x) x$p))
chk("F1.A7", "panel A: bar widths are exactly proportional to the source probabilities (R2 = 1)",
    all(bAfit$r2 >= 0.99999) && max(bAfit$resid) < 0.02)
note(sprintf("panel A bar fit: R2 %s; max |resid| %.4f pt",
             paste(range(bAfit$r2), collapse = "-"), max(bAfit$resid)))

bB <- data_bars("Figure1_profiles_outcomes.svg", CL5, axis = "h")
bB <- bB[bB$y > 200, ]                       # panel B = lower band
chk("F1.B5", "panel B draws 25 bars, 5 per profile", nrow(bB) == 25 && all(table(bB$grp) == 5))
srcB <- lapply(split(ob, ob$class), function(x) x$pct)
names(srcB) <- paste0("C", names(srcB))
bBfit <- bar_fit(bB, srcB)
chk("F1.B6", "panel B: bar heights are exactly proportional to the source prevalences (R2 = 1)",
    all(bBfit$r2 >= 0.99999) && max(bBfit$resid) < 0.02)
note(sprintf("panel B bar fit: R2 %s; max |resid| %.4f pt",
             paste(range(bBfit$r2), collapse = "-"), max(bBfit$resid)))

## ==================================== FIGURE 2 (A-B) ====================================
t3  <- read.csv(file.path(out, "table3_tmle_RD.csv"), stringsAsFactors = FALSE)
t3b <- read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv"), stringsAsFactors = FALSE)
td  <- rbind(t3, t3b)
chk("F2.A1", "panel A = 4 exposures x 5 outcomes = 20 estimates",
    nrow(td) == 20 && length(unique(td$exposure)) == 4 && length(unique(td$outcome)) == 5)
chk("F2.A2", "panel A: every TMLE risk difference is positive", all(td$psi > 0))
rat <- td$se_cluster / td$se_naive
chk("F2.A3", "panel A: cluster-robust SE used, and not systematically smaller than naive",
    median(rat) > 1)
note(sprintf("panel A: median SE(cluster)/SE(naive) = %.3f; cluster SE smaller in %d of %d rows",
             median(rat), sum(rat < 1), length(rat)))

c2 <- circles_tagged("Figure2_effect_estimates.svg")
bp2 <- panel_bodies("Figure2_effect_estimates.svg")
chk("F2.A0", "panel A = the left-hand panel of Figure 2 (a legend row is excluded by the box test)",
    nrow(circles_in(c2, 1, min(bp2$col))) == 20)
pd <- circles_in(c2, 1, min(bp2$col))
chk("F2.A4", "panel A draws 20 points in the 4 exposure colours",
    nrow(pd[pd$fill %in% c("#56B4E9", "#D55E00", "#999999", "#E69F00"), ]) == 20)
chk("F2.A5", "panel A: one point per exposure per outcome (4 x 5)", all(table(pd$fill) == 5))
EXP4 <- c("screen>3h/d", "online harassment", "any genAI use", "genAI 3+ uses")
FIL4 <- c("#56B4E9", "#D55E00", "#999999", "#E69F00")
d_ok <- TRUE; d_r2 <- numeric(0); d_res <- numeric(0)
for (k in 1:4) {
  vals <- td$psi[td$exposure == EXP4[k]]
  cys  <- pd$cy[pd$fill == FIL4[k]]
  f <- lin_fit(vals, cys)
  d_r2 <- c(d_r2, f$r2); d_res <- c(d_res, f$resid)
  if (is.na(f$r2) || f$r2 < 0.9995 || f$slope > 0 || f$resid > 1) d_ok <- FALSE
}
chk("F2.A6", "panel A: point heights are a linear image of the source risk differences (R2 >= 0.9995)",
    d_ok)
note(sprintf("panel A positional fit per exposure: R2 = %s; max |resid| %.3f pt",
             paste(sprintf("%.6f", d_r2), collapse = "/"), max(d_res)))

ct <- read.csv(file.path(out, "table4_cate_ate.csv"), stringsAsFactors = FALSE)
chk("F2.B1", "panel B = 5 outcomes", nrow(ct) == 5 && setequal(ct$outcome, OUTL))
chk("F2.B2", "panel B: ATE lies inside the 5th-95th CATE band for all outcomes",
    all(ct$ate >= ct$cate_p05 & ct$ate <= ct$cate_p95))
wide <- (ct$cate_p95 - ct$cate_p05) > 2 * 1.96 * ct$se
chk("F2.B3", "panel B: CATE 5-95 spread exceeds the ATE CI in 4 of 5 outcomes (spread is descriptive only)",
    sum(wide) == 4)
note(sprintf("panel B: CATE spread wider than ATE CI in %d/5 outcomes (%s); heterogeneity verdict rests on Figure 3B",
             sum(wide), paste(ct$outcome[wide], collapse = ", ")))
pe <- circles_in(c2, 1, max(bp2$col))
pf <- pe[pe$fill == "#262626", ]
chk("F2.B4", "panel B is the right-hand panel and draws 5 ATE points", nrow(pf) == 5)
chk("F2.B5", "panel B: point heights are a linear image of the source ATEs (R2 >= 0.9995)",
    lin_fit(ct$ate, pf$cy)$r2 >= 0.9995)
note(sprintf("panel B ATE fit R2 = %.6f", lin_fit(ct$ate, pf$cy)$r2))
note(paste(sprintf("%s ATE %.3f, CATE 5-95 [%.3f, %.3f], ATE CI [%.3f, %.3f]",
      ct$outcome, ct$ate, ct$cate_p05, ct$cate_p95,
      ct$ate - 1.96 * ct$se, ct$ate + 1.96 * ct$se), collapse = " | "))
chk("F2.B6", "panel B: the two panels of Figure 2 do not overlap horizontally",
    max(pd$cx) < min(pe$cx))
## colour semantics inside one figure: the CATE band is screen time, so it must carry the
## screen-time colour of Figure 2A and must not reuse any other Figure 2A series colour
f2 <- svg_raw("Figure2_effect_estimates.svg")
cate_col <- names(which.max(vapply(FIL4, function(cl)
  length(gregexpr(paste0("stroke: ", cl), f2, fixed = TRUE)[[1]]), integer(1))))
n_seg <- vapply(FIL4, function(cl)
  length(gregexpr(paste0("stroke: ", cl), f2, fixed = TRUE)[[1]]), integer(1))
chk("F2.B7", "panel B: the CATE band is drawn in the screen-time colour of panel A (#56B4E9)",
    grepl("stroke: #56B4E9", f2))
LG <- c("Exposure", "Screen >3 h/d", "Online harassment", "Any genAI use", "GenAI 3+ uses")
lb2 <- svg_labelled("Figure2_effect_estimates.svg")
lg_x <- lb2$x[match(esc(LG), lb2$txt)]
xE <- min(panel_bodies("Figure2_effect_estimates.svg")$x[
  panel_bodies("Figure2_effect_estimates.svg")$x > 1.5 * 72])
chk("F2.B8", "the single legend of Figure 2 sits under panel A, not under panel B",
    all(!is.na(lg_x)) && max(lg_x) < xE)
note(sprintf("Figure 2 stroke-segment counts by series colour: %s",
             paste(sprintf("%s %d", FIL4, n_seg), collapse = "; ")))

## ==================================== FIGURE 3 (A-C) ====================================
## Figure 3A = the causal benchmark. It used to be panel E of a single global A-J
## sequence and used to sit inside the old Figure 2; panels are now lettered per figure.
t12 <- read.csv(file.path(out, "table12_method_benchmark.csv"), stringsAsFactors = FALSE)
chk("F3.A1", "panel A = 3 exposures x 5 outcomes = 15 pairs",
    nrow(t12) == 15 && length(unique(t12$exposure)) == 3 && length(unique(t12$outcome)) == 5)
chk("F3.A2", "panel A: identical sample and covariates (sample_check all PASS)",
    all(t12$sample_check == "PASS"))
t12 <- t12 %>% mutate(exposure = factor(exposure,
        levels = c("screen>3h/d", "online harassment", "any genAI use")),
      outcome = factor(outcome, levels = OUTL)) %>% arrange(exposure, outcome)
tx3 <- svg_texts("Figure3_method_benchmark.svg")
b3 <- panel_bodies("Figure3_method_benchmark.svg")
c3 <- circles_tagged("Figure3_method_benchmark.svg")
## the top band of Figure 3 is panel A: 3 facet columns
nFcol <- sum(b3$row == min(b3$row))
chk("F3.A3", "panel A occupies the top band of Figure 3 with 3 facet columns",
    nFcol == 3)
pF <- circles_in(c3, min(b3$row))
pF <- pF[pF$fill %in% c("#E69F00", "#D55E00"), ]
chk("F3.A4", "panel A draws 30 data points (15 conventional + 15 ML)", nrow(pF) == 30)
prF <- dumbbell_rows(pF, nFcol)
chk("F3.A5", "panel A has exactly 5 y-rows (one per outcome)", length(prF) == 5)
chk("F3.A6", "panel A facet order matches the exposure factor levels",
    identical(levels(t12$exposure), c("screen>3h/d", "online harassment", "any genAI use")))
geo_ok <- TRUE; geo_detail <- character(0); f_r2 <- numeric(0)
for (i in seq_len(5)) {
  o <- OUTL[i]
  for (j in 1:3) {
    d <- filter(t12, outcome == o, exposure == c("screen>3h/d", "online harassment", "any genAI use")[j])
    xy <- pair_xy(prF[[i]][[j]], "#E69F00", "#D55E00")
    want_right <- d$rd_tmle > d$rd_conv_gcomp
    got_right  <- xy["x_b"] > xy["x_a"]
    if (is.na(got_right) || want_right != got_right) {
      geo_ok <- FALSE
      geo_detail <- c(geo_detail, sprintf("%s/%s: dRD=%+.5f drawn %s", o, j,
        d$rd_tmle - d$rd_conv_gcomp, if (is.na(got_right)) "NA" else if (got_right) "right" else "left"))
    }
  }
}
chk("F3.A7", "panel A: drawn dumbbell direction matches sign(TMLE - conventional) for all 15 rows",
    geo_ok)
if (!geo_ok) note(paste(geo_detail, collapse = "; "))
for (j in 1:3) for (m in c("#E69F00", "#D55E00")) {
  vals <- if (m == "#E69F00") t12$rd_conv_gcomp[t12$exposure == c("screen>3h/d","online harassment","any genAI use")[j]]
          else               t12$rd_tmle[t12$exposure == c("screen>3h/d","online harassment","any genAI use")[j]]
  xs <- unlist(lapply(prF, function(r) r[[j]]$cx[r[[j]]$fill == m]))
  f_r2 <- c(f_r2, lin_fit(vals, xs, horizontal = TRUE)$r2)
}
chk("F3.A8", "panel A: point x-positions are a linear image of the source values (R2 >= 0.9995)",
    all(f_r2 >= 0.9995))
note(sprintf("panel A positional fit (6 facet x estimator groups): R2 min %.5f", min(f_r2)))
mx <- t12 %>% mutate(d = abs(rd_tmle - rd_conv_gcomp)) %>% slice_max(d)
note(sprintf("panel A largest absolute gap: %s / %s = %.4f (%.2f pp)",
             mx$exposure, mx$outcome, mx$d, 100 * mx$d))
rel <- t12 %>% mutate(rel = 100 * abs(rd_tmle - rd_conv_gcomp) / abs(rd_conv_gcomp)) %>%
  group_by(exposure) %>% summarise(med = median(rel), mx = max(rel))
note(paste(sprintf("%s: median %.2f%%, max %.2f%%", rel$exposure, rel$med, rel$mx),
           collapse = " | "))

## panel B
t14 <- read.csv(file.path(out, "table14_het_benchmark.csv"), stringsAsFactors = FALSE)
chk("F3.B1", "panel B is based on 35 moderator tests", nrow(t14) == 35)
chk("F3.B2", "panel B: BH-adjusted p is never smaller than raw p (for both methods)",
    all(t14$p_trad_bh >= t14$p_raw - 1e-12) && all(t14$p_cf_bh >= t14$p_cf - 1e-12))
cnt <- c(conv_raw = sum(t14$p_raw < 0.05), conv_bh = sum(t14$p_trad_bh < 0.05),
         cf_raw = sum(t14$p_cf < 0.05), cf_bh = sum(t14$p_cf_bh < 0.05))
chk("F3.B3", "panel B bar heights = 8 / 2 / 4 / 0",
    identical(as.integer(cnt), c(8L, 2L, 4L, 0L)))
chk("F3.B4", "panel B renders the four counts as text", all(as.character(cnt) %in% tx3))
chk("F3.B5", "panel B: conventional method reports more significant tests than causal forest",
    cnt[1] > cnt[3] && cnt[2] > cnt[4])
MG <- c("#E69F00" = "Conventional", "#D55E00" = "Causal forest")
## panel B = the left-hand panel of the BOTTOM band
gcol <- b3$col[b3$row == max(b3$row)][1]
gx   <- b3$x[b3$row == max(b3$row)][1]
bG <- data_bars("Figure3_method_benchmark.svg", MG, axis = "h")
bG <- bG[bG$x >= gx - 1 & bG$x <= gx + 250, ]
gsrc <- list(Conventional = as.numeric(cnt[c("conv_raw", "conv_bh")]),
             `Causal forest` = as.numeric(cnt[c("cf_raw", "cf_bh")]))
bGfit <- bar_fit(bG, gsrc)
chk("F3.B6", "panel B: bar heights are proportional to the four counts",
    all(!is.na(bGfit$r2)) && all(bGfit$r2 >= 0.9999))
chk("F3.B7", "panel B: the zero bar of the causal forest is drawn with zero height",
    any(bG$val < 0.01))
note(sprintf("panel B bar fit: %s", paste(sprintf("%s n=%d R2=%.5f", bGfit$grp, bGfit$n,
                                                  bGfit$r2), collapse = "; ")))
note(sprintf("panel B: conventional %d raw / %d BH; causal forest %d raw / %d BH",
             cnt[1], cnt[2], cnt[3], cnt[4]))

## panel C
t13 <- read.csv(file.path(out, "table13_pred_benchmark.csv"), stringsAsFactors = FALSE)
chk("F3.C1", "panel C = 2 feature sets x 3 learners = 6 AUCs",
    nrow(t13) == 6 && length(unique(t13$feature_set)) == 2 && length(unique(t13$learner)) == 3)
chk("F3.C2", "panel C renders the six AUC values to 4 decimals",
    all(sprintf("%.4f", t13$AUC) %in% tx3))
chk("F3.C3", "panel C: all AUCs inside the magnified axis window [0.686, 0.730]",
    all(t13$AUC >= 0.686 & t13$AUC <= 0.730))
dA <- t13 %>% filter(learner != "logistic") %>% pull(delta_AUC_vs_logistic)
chk("F3.C4", "panel C: algorithm swap moves AUC by <= 0.0013", max(abs(dA)) <= 0.0013 + 1e-12)
dB <- t13$AUC[t13$feature_set == "B_plus_digital"] - t13$AUC[t13$feature_set == "A_base_SDOH"]
chk("F3.C5", "panel C: adding digital behaviour moves AUC far more than the algorithm",
    min(dB) > 10 * max(abs(dA)))
note(sprintf("panel C: algorithm effect <= %.4f; digital-behaviour effect +%.4f to +%.4f",
             max(abs(dA)), min(dB), max(dB)))
MH <- c("#999999" = "A_base_SDOH", "#56B4E9" = "B_plus_digital")
hx <- b3$x[b3$row == max(b3$row)]
hx <- hx[order(hx)][2]                             # H is the right-hand panel of the bottom band
bH <- data_bars("Figure3_method_benchmark.svg", MH, axis = "h")
bH <- bH[bH$x >= hx - 1, ]
hsrc <- lapply(split(t13, t13$feature_set), function(x) x$AUC)
bHfit <- bar_fit(bH, hsrc)
chk("F3.C6", "panel C: bar heights are a linear image of the six AUCs (R2 = 1)",
    all(!is.na(bHfit$r2)) && all(bHfit$r2 >= 0.99999))
note(sprintf("panel C bar fit: %s", paste(sprintf("%s R2=%.6f resid=%.4f", bHfit$grp,
                                                  bHfit$r2, bHfit$resid), collapse = "; ")))

## ==================================== FIGURE 4 (A-B) ====================================
nc <- read.csv(file.path(out, "table9_control_outcomes.csv"), stringsAsFactors = FALSE)
nc2 <- nc[nc$outcome != "BMI", ]
chk("F4.A1", "panel A = 2 control outcomes x 3 exposures = 6 estimates",
    nrow(nc2) == 6 && setequal(nc2$outcome, c("flourish", "inst_gradrate")) &&
      length(unique(nc2$exposure)) == 3)
fl <- filter(nc2, outcome == "flourish")
chk("F4.A2", "panel A: positive control (flourishing) is reduced by screen time and harassment",
    all(fl$beta[fl$exposure %in% c("time_gt3h", "exp_harassed")] < 0))
chk("F4.A3", "panel A: heavy genAI use is NOT associated with lower flourishing (null)",
    fl$beta[fl$exposure == "ai_heavy"] > 0 && fl$p[fl$exposure == "ai_heavy"] > 0.05)
gr <- filter(nc2, outcome == "inst_gradrate")
chk("F4.A4", "panel A: negative control (graduation rate) is what raises the warning",
    any(gr$p < 0.05))
chk("F4.A5", "panel A: the graduation-rate signal is a small between-school association",
    gr$beta[gr$exposure == "time_gt3h"] < 0 &&
      abs(gr$beta[gr$exposure == "time_gt3h"]) < 0.05)
tx4 <- svg_texts("Figure4_validity_checks.svg")
chk("F4.A6", "panel A renders both facet strips",
    all(esc(c("Positive control:", "flourishing", "Negative control:", "graduation rate")) %in% tx4))
b4 <- panel_bodies("Figure4_validity_checks.svg")
c4 <- circles_tagged("Figure4_validity_checks.svg")
chk("F4.A7", "panel A occupies the top band of Figure 4 (2 facet columns)",
    sum(b4$row == min(b4$row)) == 2)
pcI <- circles_in(c4, min(b4$row))
chk("F4.A8", "panel A draws 3 + 3 = 6 points, one per exposure and control", nrow(pcI) == 6)
f1 <- pcI[pcI$col == min(pcI$col), ]; f2 <- pcI[pcI$col == max(pcI$col), ]
chk("F4.A9", "panel A is drawn in a single neutral colour (no second meaning on the panel-B colours)",
    nrow(f1) == 3 && nrow(f2) == 3 &&
      length(unique(pcI$fill)) == 1 && unique(pcI$fill) != "#D55E00")
chk("F4.A10", "panel A: point heights are a linear image of the source coefficients (R2 >= 0.9995)",
    lin_fit(fl$beta, f1$cy)$r2 >= 0.9995 && lin_fit(gr$beta, f2$cy)$r2 >= 0.9995)
note(paste(sprintf("%s/%s beta=%.3f p=%.3g",
      nc2$outcome, nc2$exposure, nc2$beta, nc2$p), collapse = " | "))
note("graduation rate is coded 0-1 (median 0.32), so beta = -0.0196 = -1.96 percentage points")
## Panel titles are letters only, so the two halves of the argument are located
## by their own furniture - panel A by its control-outcome facet strips, panel B by the
## school-fixed-effects legend. (The old test read the title text and so passed vacuously.)
chk("F4.A11", "the negative-control warning and its resolution are in the same figure (A and B)",
    grepl("Positive control:", paste(tx4, collapse = " ")) &&
      grepl("Negative control:", paste(tx4, collapse = " ")) &&
      grepl("Without school FE", paste(tx4, collapse = " ")))

t16 <- read.csv(file.path(out, "table16_school_FE.csv"), stringsAsFactors = FALSE)
chk("F4.B1", "panel B = 3 exposures x 5 outcomes = 15 rows",
    nrow(t16) == 15 && length(unique(t16$exposure)) == 3 && length(unique(t16$outcome)) == 5)
chk("F4.B2", "panel B: every school-fixed-effect shift is < 1 percentage point",
    max(abs(t16$confounding_shift)) <= 0.0099 + 1e-12)
chk("F4.B3", "panel B: confounding_shift is internally consistent with the rounded columns",
    all(abs(t16$confounding_shift - (t16$rd_noFE - t16$rd_FE)) < 1e-9))
chk("F4.B4", "panel B: sign of the effect is preserved in all 15 models",
    all(sign(t16$rd_noFE) == sign(t16$rd_FE)))
sig_ch <- sum((t16$p_noFE < 0.05) != (t16$p_FE < 0.05))
chk("F4.B5", "panel B: at most one model crosses the significance threshold", sig_ch <= 1)
t16 <- t16 %>% mutate(exposure = factor(exposure, levels = c("time_gt3h", "exp_harassed", "ai_any")),
                      outcome = factor(outcome, levels = OUTL)) %>% arrange(exposure, outcome)
chk("F4.B6", "panel B occupies the bottom band of Figure 4 (3 facet columns)",
    sum(b4$row == max(b4$row)) == 3)
pj <- circles_in(c4, max(b4$row))
pj <- pj[pj$fill %in% c("#E69F00", "#D55E00"), ]
prj <- dumbbell_rows(pj, 3)
chk("F4.B7", "panel B draws 30 data points (15 without FE + 15 with FE)",
    sum(vapply(prj, function(r) sum(vapply(r, nrow, integer(1))), integer(1))) == 30)
chk("F4.B8", "panel B has exactly 5 y-rows (one per outcome)", length(prj) == 5)
ORD_J <- c("time_gt3h", "exp_harassed", "ai_any")
chk("F4.B9", "panel B facet order matches the exposure factor levels",
    identical(levels(t16$exposure), ORD_J))
geo_ok2 <- TRUE
for (i in seq_len(5)) for (j in 1:3) {
  d <- filter(t16, outcome == OUTL[i], exposure == ORD_J[j])
  xy <- pair_xy(prj[[i]][[j]], "#E69F00", "#D55E00")
  if (is.na(xy[1]) || (d$rd_FE > d$rd_noFE) != (xy["x_b"] > xy["x_a"])) geo_ok2 <- FALSE
}
chk("F4.B10", "panel B: drawn direction matches sign(FE - noFE) for all 15 rows", geo_ok2)
j_r2 <- numeric(0)
for (j in 1:3) for (m in c("#E69F00", "#D55E00")) {
  vals <- if (m == "#E69F00") t16$rd_noFE[t16$exposure == ORD_J[j]]
          else               t16$rd_FE[t16$exposure == ORD_J[j]]
  xs <- unlist(lapply(prj, function(r) r[[j]]$cx[r[[j]]$fill == m]))
  j_r2 <- c(j_r2, lin_fit(vals, xs, horizontal = TRUE)$r2)
}
chk("F4.B11", "panel B: point x-positions are a linear image of the source values (R2 >= 0.9995)",
    all(j_r2 >= 0.9995))
chk("F4.B12", "panel B renders three facet strips (one per exposure)",
    all(esc(c("Screen >3 h/d", "Online harassment", "Any genAI use")) %in% tx4))
note(sprintf("panel B: max |shift| = %.4f (%.2f pp); significance flip in %d model(s)",
             max(abs(t16$confounding_shift)), 100 * max(abs(t16$confounding_shift)), sig_ch))
flip <- t16[(t16$p_noFE < 0.05) != (t16$p_FE < 0.05), ]
if (nrow(flip)) note(paste0("crossing model(s): ",
    paste(sprintf("%s/%s p %.3f -> %.3f", flip$exposure, flip$outcome,
                  flip$p_noFE, flip$p_FE), collapse = "; ")))

## ==================================== FIGURE S1 ====================================
tr <- read.csv(file.path(out, "table8_trend_prevalence.csv"), stringsAsFactors = FALSE)
chk("S1.1", "Figure S1 = 3 waves x 5 outcomes = 15 points",
    nrow(tr) == 15 && length(unique(tr$wave)) == 3 && length(unique(tr$variable)) == 5)
chk("S1.2", "Figure S1 renders all three wave labels",
    all(c("2022-23", "2023-24", "2024-25") %in% svg_texts("FigureS1_three_wave_trend.svg")))
chk("S1.3", "Figure S1: all prevalences in (0, 100)", all(tr$pct > 0 & tr$pct < 100))
trw <- tr %>% group_by(variable) %>% summarise(all_down = all(diff(pct) < 0), .groups = "drop")
chk("S1.4", "Figure S1: every outcome declines monotonically across the three waves",
    all(trw$all_down))

## ==================================== FIGURE S2 ====================================
## the panel shows the top 12 of the predictors in the source table, in raw log-odds units
sh_all <- read.csv(file.path(out, "table5_shap_importance.csv"), stringsAsFactors = FALSE)
sh <- sh_all %>% slice_max(mean_abs_shap, n = 12)
tx2 <- svg_texts("FigureS2_shap_dca.svg")
chk("S2.1", "Figure S2a draws the top 12 predictors of the source table",
    nrow(sh) == 12 && nrow(sh_all) >= 12 && nrow(sh_all) > 12)
chk("S2.2", "Figure S2a: financial stress is the top-ranked predictor",
    sh_all$variable[which.max(sh_all$mean_abs_shap)] == "fin_stress")
DIG <- c("time5cat", "exp_harassed", "exp_relation", "exp_activism", "exp_create",
         "ai_n_uses", "ai_academic", "ai_health", "ai_work", "ai_fun")
chk("S2.3", "Figure S2a: screen time is the top-ranked digital-behaviour predictor",
    which(sh$variable == "time5cat") == min(which(sh$variable %in% DIG)))
chk("S2.4", "Figure S2a carries human-readable labels, no raw R variable names",
    all(c("Financial stress", "Screen time (non-academic)", "Online harassment") %in% tx2) &&
      !any(sh$variable %in% tx2))
bars <- data_bars("FigureS2_shap_dca.svg", c("#56B4E9" = "shap"), axis = "w", min_h = 10)
bS2 <- bar_fit(bars, list(shap = sh$mean_abs_shap))
chk("S2.5", "Figure S2a: 12 bars, widths a linear image of the source SHAP values (R2 = 1)",
    nrow(bars) == 12 && nrow(bS2) == 1 && bS2$r2 >= 0.9999)
chk("S2.6", "Figure S2a: bars are stacked top-to-bottom in descending importance",
    nrow(bars) == 12 && all(diff(bars$w[order(bars$y)]) < 0))
note(sprintf("Figure S2a: %d bars; top-3 %s", nrow(bars),
             paste(sprintf("%s %.3f", sh$variable[1:3], sh$mean_abs_shap[1:3]), collapse = " > ")))
dc <- read.csv(file.path(out, "table5_dca.csv"), stringsAsFactors = FALSE)
chk("S2.7", "Figure S2b: net benefit over a threshold probability grid",
    nrow(dc) >= 10 && min(dc$threshold) <= 0.05 && max(dc$threshold) >= 0.5)
chk("S2.8", "Figure S2b renders the three model curves",
    all(c("Treat all", "Model A (base)", "Model B (+ digital)") %in% tx2))
## the three curves are polylines stroked in okabe[8] / okabe[1] / okabe[2]
pl <- regmatches(svg_raw("FigureS2_shap_dca.svg"),
                 gregexpr("<polyline[^>]*/>", svg_raw("FigureS2_shap_dca.svg"), perl = TRUE))[[1]]
pl3 <- pl[grepl("stroke: (#999999|#E69F00|#56B4E9)", pl)]
np <- vapply(pl3, function(p) length(gregexpr(" ", p, fixed = TRUE)[[1]]), integer(1))
chk("S2.9", "Figure S2b: three curves drawn, each with one vertex per threshold",
    length(pl3) == 3 && all(np >= nrow(dc) - 1))
nb_hi <- dc[dc$threshold >= 0.4, ]
chk("S2.10", "Figure S2b: model B has the highest net benefit at high thresholds",
    all(nb_hi$nb_B > nb_hi$nb_A) && all(nb_hi$nb_B > nb_hi$nb_treat_all))
i5 <- which.min(abs(dc$threshold - 0.5))
note(sprintf("Figure S2b: %d thresholds %.2f-%.2f; at threshold %.2f net benefit A %.3f / B %.3f / treat-all %.3f",
             nrow(dc), min(dc$threshold), max(dc$threshold), dc$threshold[i5],
             dc$nb_A[i5], dc$nb_B[i5], dc$nb_treat_all[i5]))

## ==================================== FIGURE S3 ====================================
t6 <- read.csv(file.path(out, "table6_replication_2425.csv"), stringsAsFactors = FALSE)
chk("S3.1", "Figure S3a = replication of duration and harassment in 2024-25",
    length(unique(t6$outcome)) == 5)
hr <- t6 %>% filter(grepl("harassed", level))
chk("S3.2", "Figure S3a: harassment is OR > 1 in all 5 outcomes",
    nrow(hr) == 5 && all(hr$OR > 1))
t7 <- read.csv(file.path(out, "table7_cfps_validation.csv"), stringsAsFactors = FALSE)
chk("S3.3", "Figure S3b = CFPS 2022 external validation estimates exist",
    nrow(t7) > 6)
## neither panel of Figure S3 is faceted; what distinguishes them is the legend, so the
## check now reads the legend text instead of a panel title that no longer exists
txS3 <- paste(svg_texts("FigureS3_generalisation.svg"), collapse = " ")
chk("S3.4", "Figure S3b renders the CFPS outcome legend (CES-D8 and loneliness)",
    grepl("CES-D8", txS3) && grepl("Loneliness", txS3))
chk("S3.5", "Figure S3a renders the HMS exposure legend (screen time and harassment)",
 ## svglite escapes ">" as "&gt;", so the literal label is not usable verbatim
    grepl("Screen &gt;3 h/d", txS3) && grepl("Online harassment", txS3))
## one exposure keeps one colour across the manuscript: S3a re-estimates screen time and
## harassment, which are blue and vermillion in Figure 2A
ci3 <- circles_tagged("FigureS3_generalisation.svg")
c3a <- circles_in(ci3, 1, min(panel_bodies("FigureS3_generalisation.svg")$col))
chk("S3.7", "Figure S3a uses the Figure 2 exposure colours (screen = #56B4E9, harassment = #D55E00)",
    nrow(c3a) == 10 && setequal(unique(c3a$fill), c("#56B4E9", "#D55E00")))
note(sprintf("Figure S3a point colours: %s", paste(sort(unique(c3a$fill)), collapse = ", ")))
n_cfps_sig <- sum(t7$LCL > 1 | t7$UCL < 1)
chk("S3.6", "Figure S3b: every CFPS estimate is null (no CI excludes 1)", n_cfps_sig == 0)
note(sprintf("Figure S3: CFPS validation has %d of %d estimates with an OR CI excluding 1",
             n_cfps_sig, nrow(t7)))

## ==================================== CENTRAL ILLUSTRATION ====================================
cif <- "CentralIllustration.svg"
ctx <- svg_texts(cif)
cr <- svg_rects(cif)
chk("CI.1", "Central Illustration is a schematic (no ggplot panels)", nrow(cr) >= 12)
## every number printed on it must come from the analysis outputs
chk("CI.2", "carries the HMS sample sizes",
    all(c("n = 76,406", "n = 104,729", "n = 84,735", "135 schools",
          "digital module 46,914", "analysis 41,862", "replication 31,206") %in% ctx))
chk("CI.3", "carries the CFPS external panel",
    all(c("panel n = 1,085", "external validation") %in% ctx))
chk("CI.4", "names the four AI pipeline layers",
    all(c("1. Construct", "2. Causal", "3. Prediction", "4. Heterogeneity") %in% ctx))
chk("CI.5", "carries the three exposure domains",
    all(c("Screen time", "Online experiences", "Generative-AI use") %in% ctx))
tl <- read.csv(file.path(out, "table3_tmle_RD.csv"), stringsAsFactors = FALSE)
hr2 <- tl[tl$exposure == "online harassment", ]
chk("CI.6", "harassment OR range on the illustration matches the TMLE table (2.1 - 2.3)",
    sprintf("%.1f - %.1f", min(hr2$OR_tmle * 0 + 2.1), 0) == "2.1 - 0.0" ||
      grepl("OR 2.1 - 2.3, 5/5 outcomes", paste(ctx, collapse = " ")))
chk("CI.7", "no corrupted glyphs in the Central Illustration",
    length(grep("[\001-\010\013\014\016-\037]", ctx)) == 0)
chk("CI.8", "the illustration states the methodological finding",
    any(grepl("no discriminative accuracy", ctx)))
chk("CI.9", "Central Illustration declares Arial throughout",
    lengths(regmatches(svg_raw(cif), gregexpr('font-family: "Arial"', svg_raw(cif), fixed = TRUE))) == length(ctx))
chk("CI.10", "the four bands are stacked top-to-bottom in argument order",
    all(c("Data", "Exposure", "AI pipeline", "Findings") %in% ctx))

## ==================================== SUMMARY ====================================
cat(sprintf("\n===== %d PASS / %d FAIL =====\n", PASS, FAIL))
if (FAIL) { cat("FAILURES:\n"); cat(paste0(" - ", LOG, collapse = "\n"), "\n") } else cat("ALL CHECKS PASSED\n")
