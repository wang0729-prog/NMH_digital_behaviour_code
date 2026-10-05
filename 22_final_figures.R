source("00_config.R")
## 22_final_figures.R
## Reorganised figure set for the manuscript.
##
## Problems with the previous set of 13 single-panel figures:
## * too many standalone files for a journal article;
## * Fig5a (AUC of two nested models) is entirely subsumed by the prediction-benchmark
## figure, which shows the same two feature sets crossed with three learners -> deleted;
## * the negative-control figure (raises a warning) and the school-FE figure (resolves it)
## lived in separate files although they are one argument -> now adjacent panels;
## * the methodological-benchmark layer had been bolted on as four extra files -> now
## integrated into the main figures rather than appended.
##
## Result (retargeted at the journal):
## Figure 1 (A-C) digital-behaviour profiles and their outcomes [substance]
## Figure 2 (A-B) effect estimates: TMLE risk differences and CATE spread [substance]
## Figure 3 (A-C) methodological benchmark: causal / heterogeneity / prediction
## Figure 4 (A-B) validity: control outcomes raise a warning, school FE resolve it
## Figure S1 three-wave outcome trend
## Figure S2 (a-b) SHAP importance and decision curve
## Figure S3 (a-b) internal replication and external validation
##
## PANEL LETTERING (changed on request):
## every figure starts at A again - Figure 2 is A/B, Figure 3 is A/B/C, Figure 4 is A/B.
## A single A-J sequence across the whole manuscript is NOT how journals letter panels.
## Main figures use CAPITALS, supplementary figures lower case (current journal practice).
## Cross-figure references are written "Figure 3B", never a bare letter.
## The ggplot objects below are named after the panel they draw (pF2a = Figure 2 panel A),
## so the object name and the printed letter can never drift apart.
##
## Mapping from the previous global A-J scheme:
## old D -> Figure 2A | old E -> Figure 2B | old F -> Figure 3A
## old G -> Figure 3B | old H -> Figure 3C | old I -> Figure 4A | old J -> Figure 4B
##
## Why Figure 3 was split (it used to hold G/H/I/J):
## the journal allows at most 6 display items (figures AND tables combined) and
## explicitly states "we do not use schemes". Four main figures + two main tables fills
## that budget exactly, and the Central Illustration is submitted separately as a
## graphical abstract precisely BECAUSE schemes are not accepted as display items
## (see 27_central_illustration.R).
## The old Figure 3 mixed two unrelated arguments: F/G/H answer "does the ML pipeline
## matter?" (the methodological claim) while I/J answer "are the estimates credible?"
## (validity). One caption cannot serve both. They are now separate figures.
##
## Layout conventions (after a visual pass):
## * panel titles are the PANEL LETTER ONLY (A, B, C for main figures, a, b for
## supplementary figures), 8 pt bold. The descriptive sentence for each panel lives in
## the figure legend supplied with the manuscript, never inside the image. This is the
## author's house style and applies to every figure in this project;
## a single-panel supplementary figure carries no title at all;
## * legends with more than two entries are wrapped to two rows;
## * panels that carry facets get the full 180 mm width.
## Assembled with base grid (no patchwork/cowplot available in this library).

if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(grid) })

dir <- ANALYSIS_DIR
out <- file.path(dir, "outputs"); fig <- file.path(dir, "figures")
dir.create(fig, showWarnings = FALSE)

SZ_PANEL <- 8; SZ_AXIS <- 7; SZ_TICK <- 6
okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
theme_ms <- function(base = 7) {
  theme_bw(base_size = base) +
    theme(text = element_text(family = "Arial", size = SZ_TICK),
          axis.title = element_text(size = SZ_AXIS),
          axis.text = element_text(size = SZ_TICK),
          legend.title = element_text(size = SZ_TICK),
          legend.text = element_text(size = SZ_TICK),
          strip.text = element_text(size = SZ_AXIS, face = "bold"),
          plot.title = element_text(size = SZ_PANEL, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
          legend.background = element_blank(),
          legend.key.size = unit(3.2, "mm"),
          legend.spacing.x = unit(1.6, "mm"),
          plot.margin = margin(3, 4, 3, 3))
}
## legends with 3+ entries are too wide for a half-width panel -> wrap to two rows
two_row <- function(...) guides(..., colour = guide_legend(nrow = 2, byrow = TRUE),
                                fill = guide_legend(nrow = 2, byrow = TRUE))

## ---- grid composer: layout matrix as in base graphics, 1 panel per integer ----
## : the caption strip has been REMOVED from the image. the journal
## requires figure legends to be supplied as text with the manuscript, so baking the same
## paragraph into the artwork duplicated it and left the submitted files carrying a legend
## they should not have. The canvas is now exactly the panel area; panel titles (A, B, C)
## stay inside the panels so a figure is still self-identifying.
save_combo <- function(panels, lay, name, w_mm, h_mm, heights = NULL) {
 ## svglite/cairo mangle non-ASCII glyphs (delta, <= were corrupted on a previous run).
  for (p in panels) {
    ti <- c(p$labels$title, p$labels$x, p$labels$y, p$labels$colour, p$labels$fill)
    for (s in ti) if (!is.null(s) && any(as.integer(charToRaw(as.character(s))) > 127L))
      stop("non-ASCII label in ", name, ": ", s)
  }
 ## panel titles must fit the panel width, otherwise they are silently clipped
  for (p in panels) {
    ti <- p$labels$title
    if (!is.null(ti) && nchar(ti) > 46) stop("panel title too long (", nchar(ti), "): ", ti)
  }
  grobs <- lapply(panels, ggplot2::ggplotGrob)
  nr <- nrow(lay); nc <- ncol(lay)
  if (is.null(heights)) heights <- unit(rep(1, nr), "null")
  stopifnot(length(heights) == nr)

 ## ---- equalise the plot boxes of panels that sit side by side in one layout row ----
 ## Every ggplot panel row is unit(1, "null"), so a panel that carries less fixed-height
 ## furniture underneath (no legend, no rotated axis labels) silently grows a TALLER plot
 ## box than the panel next to it. The two boxes then no longer line up. Padding the
 ## shorter grob with a fixed spacer row below its panel restores the alignment; nothing
 ## is added inside the panel, so no data element moves.
  grid_rows <- function(gs) sort(unique(gs$layout$t[gs$layout$name == "panel"]))
  fixed_h <- function(g, rows) {
    if (!length(rows)) return(0)
    grid::convertHeight(sum(g$heights[rows]), "pt", valueOnly = TRUE)
  }
  for (r in seq_len(nr)) {
    ids <- which(vapply(seq_along(lay), function(i) {
      pos <- which(lay == i, arr.ind = TRUE)
      length(unique(pos[, 1])) == 1 && unique(pos[, 1]) == r
    }, logical(1)))
    if (length(ids) < 2) next
    gs <- grobs[ids]
    pr <- lapply(gs, grid_rows)
    top <- vapply(seq_along(gs), function(i)
      fixed_h(gs[[i]], if (min(pr[[i]]) > 1) seq_len(min(pr[[i]]) - 1) else integer(0)), 0)
    bot <- vapply(seq_along(gs), function(i)
      { n <- length(gs[[i]]$heights); if (max(pr[[i]]) < n)
          fixed_h(gs[[i]], (max(pr[[i]]) + 1):n) else 0 }, 0)
    for (i in seq_along(gs)) {
      if (max(top) - top[i] > 0.05)
        gs[[i]] <- gtable::gtable_add_rows(
          gs[[i]], grid::unit(max(top) - top[i], "pt"), pos = min(pr[[i]]))
      n <- length(gs[[i]]$heights)
      if (max(bot) - bot[i] > 0.05)
        gs[[i]] <- gtable::gtable_add_rows(
          gs[[i]], grid::unit(max(bot) - bot[i], "pt"), pos = n)
    }
    grobs[ids] <- gs
  }

  draw <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc,
                                               widths = unit(rep(1, nc), "null"),
                                               heights = heights)))
    for (i in seq_along(grobs)) {
      pos <- which(lay == i, arr.ind = TRUE)
      r <- min(pos[, 1]); r2 <- max(pos[, 1])
      c1 <- min(pos[, 2]); c2 <- max(pos[, 2])
      pushViewport(viewport(layout.pos.row = r:r2, layout.pos.col = c1:c2))
      grid.draw(grobs[[i]])
      popViewport()
    }
    popViewport()
  }
  grDevices::cairo_pdf(file.path(fig, paste0(name, ".pdf")),
                       width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
  suppressWarnings(draw()); dev.off()
  svglite::svglite(file.path(fig, paste0(name, ".svg")),
                   width = w_mm / 25.4, height = h_mm / 25.4,
                   system_fonts = list(sans = "Arial"))
  suppressWarnings(draw()); dev.off()
  cat("wrote", name, "\n")
}

OUTL <- c("dep_any", "dep_maj", "anx_any", "sui_idea", "lonely")
OUTLAB <- c("Depression (any)", "Depression (major)", "Anxiety", "Suicidal ideation", "Loneliness")
setO <- function(x) factor(x, levels = OUTL, labels = OUTLAB)

## ============================ Figure 1: profiles and outcomes ============================
cp <- read.csv(file.path(out, "lca2", "condprob_k5.csv"))
cp_bin <- cp %>% filter(level == 1, item != "time5cat")
cp_tim <- cp %>% filter(item == "time5cat", level == 4) %>% mutate(item = "time_gt3h")
pF1a <- bind_rows(cp_bin, cp_tim) %>%
  mutate(item = factor(item,
                       levels = c("time_gt3h","exp_relation","exp_create","exp_activism",
                                  "exp_harassed","ai_academic","ai_work","ai_comm",
                                  "ai_health","ai_fun","ai_forbidden"),
                       labels = c("Screen >3 h/d","Online relationships","Creating content",
                                  "Online activism","Online harassment","AI: schoolwork",
                                  "AI: work","AI: communication","AI: health","AI: fun",
                                  "AI: forbidden use"))) %>%
  select(-level) %>% pivot_longer(starts_with("class"), names_to = "class", values_to = "p") %>%
  mutate(class = factor(sub("class", "C", class), levels = c("C5","C4","C3","C2","C1")))
pF1a <- ggplot(pF1a, aes(x = item, y = p, fill = class)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.15) +
  coord_flip() + scale_fill_manual(values = okabe, name = "Profile") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Conditional probability",
       title = "A") + theme_ms()

ob <- read.csv(file.path(out, "lca2_k5_outcome_by_class.csv"))
ob$outcome <- setO(ob$outcome); ob$class <- factor(ob$class, levels = 5:1)
pF1b <- ggplot(ob, aes(x = outcome, y = pct, fill = class)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.15) +
  geom_errorbar(aes(ymin = pct - 1.96 * se, ymax = pct + 1.96 * se),
                width = 0.12, position = position_dodge(0.8), linewidth = 0.2) +
  scale_fill_manual(values = okabe, name = "Profile") +
  labs(x = NULL, y = "Weighted prevalence (%)",
       title = "B") +
  theme_ms() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

t2 <- read.csv(file.path(out, "table2_var_centered_OR.csv"))
t4 <- t2 %>% filter(grepl("raw w", exposure), grepl("^time5cat_f", level)) %>%
  mutate(cat = factor(sub("time5cat_f", "", level), levels = c("none","<1h","1-2h","2-3h",">3h")),
         outcome = setO(outcome))
## the reference category ("none") has OR = 1 by construction and is not in the model output;
## add it explicitly so the dose-response starts from the reference line.
ref <- data.frame(outcome = setO(OUTL), cat = factor("none", levels = levels(t4$cat)),
                  OR = 1, LCL = NA_real_, UCL = NA_real_)
t4 <- bind_rows(t4, ref)
pF1c <- ggplot(t4, aes(x = cat, y = OR, group = outcome, colour = outcome)) +
  geom_line(linewidth = 0.35) + geom_point(size = 1.3) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.1, linewidth = 0.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.25) +
  scale_colour_manual(values = okabe, name = NULL) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3)) +
  labs(x = "Non-academic screen time per day", y = "Adjusted odds ratio (log scale)",
       title = "C") +
  theme_ms() + theme(legend.position = "bottom") + two_row()
save_combo(list(pF1a, pF1b, pF1c), matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
           "Figure1_profiles_outcomes", 180, 139)

## ============ Figure 2: effects, causal benchmark, CATE spread ============
t3 <- read.csv(file.path(out, "table3_tmle_RD.csv"))
t3b <- read.csv(file.path(out, "table3b_tmle_aiheavy_RD.csv"))
td <- rbind(t3, t3b)
td$outcome <- setO(td$outcome)
td$exposure <- factor(td$exposure,
                      levels = c("screen>3h/d","online harassment","any genAI use","genAI 3+ uses"),
                      labels = c("Screen >3 h/d","Online harassment","Any genAI use","GenAI 3+ uses"))
pF2a <- ggplot(td, aes(x = outcome, y = 100 * psi, colour = exposure)) +
  geom_point(position = position_dodge(0.5), size = 1.4) +
  geom_errorbar(aes(ymin = 100 * (psi - 1.96 * se_cluster),
                    ymax = 100 * (psi + 1.96 * se_cluster)),
                width = 0.15, position = position_dodge(0.5), linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_colour_manual(values = okabe[c(2, 6, 8, 1)], name = "Exposure") +
  labs(x = NULL, y = "Risk difference (percentage points)",
       title = "A") +
  theme_ms() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     legend.position = "bottom") + two_row()

t12 <- read.csv(file.path(out, "table12_method_benchmark.csv"))
t12$outcome <- setO(t12$outcome)
t12$exposure <- factor(t12$exposure,
                       levels = c("screen>3h/d","online harassment","any genAI use"),
                       labels = c("Screen >3 h/d","Online harassment","Any genAI use"))
dE <- t12 %>% select(exposure, outcome, rd_conv_gcomp, rd_tmle) %>%
  pivot_longer(c(rd_conv_gcomp, rd_tmle), names_to = "method", values_to = "rd") %>%
  mutate(method = factor(method, levels = c("rd_conv_gcomp", "rd_tmle"),
                         labels = c("Conventional", "Machine learning")))
pF3a <- ggplot(dE, aes(x = 100 * rd, y = outcome, colour = method)) +
  geom_line(aes(group = outcome), colour = "grey55", linewidth = 0.35) +
  geom_point(size = 1.7) +
  facet_wrap(~ exposure, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = okabe[c(1, 6)], name = "Estimator") +
  labs(x = "Marginal risk difference (percentage points)", y = NULL,
       title = "A") +
  theme_ms() + theme(legend.position = "bottom",
                     panel.spacing.x = unit(6, "mm"))

ct <- read.csv(file.path(out, "table4_cate_ate.csv"))
ct$outcome <- setO(ct$outcome)
## the CATE band is the screen-time exposure, so it uses the screen-time colour of panel A
## (okabe[2]); okabe[6] is "Online harassment" in panel A and reusing it here would put two
## different meanings on the same colour inside one figure
pF2b <- ggplot(ct, aes(x = outcome, y = ate)) +
  geom_errorbar(aes(ymin = cate_p05, ymax = cate_p95), width = 0.22,
                linewidth = 0.9, colour = okabe[2]) +
  geom_errorbar(aes(ymin = ate - 1.96 * se, ymax = ate + 1.96 * se),
                width = 0.34, linewidth = 0.25, colour = "grey25") +
  geom_point(size = 1.7, colour = "grey15") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  labs(x = NULL, y = "Risk difference (proportion)",
       title = "B") +
  theme_ms() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
## A and B are both effect estimates for the SAME exposures: A is the average effect,
## B is how much individual effects vary around it. They belong together; the estimator
## comparison belongs with the other methodological benchmarks in Figure 3.
save_combo(list(pF2a, pF2b), matrix(c(1, 2), 1, 2, byrow = TRUE),
           "Figure2_effect_estimates", 180, 97)

## ====== Figure 3: benchmark (heterogeneity, prediction) + validity checks ======
t14 <- read.csv(file.path(out, "table14_het_benchmark.csv"))
dG <- data.frame(
  method = factor(c("Conventional", "Conventional", "Causal forest", "Causal forest"),
                  levels = c("Conventional", "Causal forest")),
  correction = factor(c("Raw", "BH-adjusted", "Raw", "BH-adjusted"),
                      levels = c("Raw", "BH-adjusted")),
  n = c(sum(t14$p_raw < 0.05), sum(t14$p_trad_bh < 0.05),
        sum(t14$p_cf < 0.05), sum(t14$p_cf_bh < 0.05)))
pF3b <- ggplot(dG, aes(x = correction, y = n, fill = method)) +
  geom_col(position = position_dodge(0.75), width = 0.6, colour = "black", linewidth = 0.15) +
  geom_text(aes(label = n, group = method), position = position_dodge(0.75),
            vjust = -0.7, size = 6 / .pt, family = "Arial") +
  geom_hline(yintercept = 0.05 * 35, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = okabe[c(1, 6)], name = "Method") +
  scale_y_continuous(limits = c(0, 9.5)) +
  labs(x = NULL, y = "Significant tests out of 35",
       title = "B") +
  theme_ms() + theme(legend.position = "bottom")

t13 <- read.csv(file.path(out, "table13_pred_benchmark.csv"))
t13$learner <- factor(t13$learner, levels = c("logistic", "lasso", "xgboost"),
                      labels = c("Logistic", "LASSO", "XGBoost"))
t13$feature_set <- factor(t13$feature_set, levels = c("A_base_SDOH", "B_plus_digital"),
                          labels = c("Base: demographics + SDOH", "Base + digital behaviour"))
pF3c <- ggplot(t13, aes(x = learner, y = AUC, fill = feature_set)) +
  geom_col(position = position_dodge(0.75), width = 0.65, colour = "black", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.4f", AUC), group = feature_set),
            position = position_dodge(0.75), vjust = -0.7, size = 5.5 / .pt, family = "Arial") +
  scale_fill_manual(values = okabe[c(8, 2)], name = NULL) +
  scale_y_continuous(limits = c(0.686, 0.730), oob = scales::squish) +
  labs(x = NULL, y = "Cross-validated AUC (magnified)",
       title = "C") +
 ## the two feature-set labels fill the full panel width on one row; stack them instead
  theme_ms() + theme(legend.position = "bottom") + two_row()

nc <- read.csv(file.path(out, "table9_control_outcomes.csv"))
nc <- nc[nc$outcome != "BMI", ]
nc$outcome <- factor(nc$outcome, levels = c("flourish", "inst_gradrate"),
                     labels = c("Positive control:\nflourishing",
                                "Negative control:\ngraduation rate"))
nc$exposure <- factor(nc$exposure, levels = c("time_gt3h", "exp_harassed", "ai_heavy"),
                      labels = c("Screen >3 h/d", "Online harassment", "GenAI 3+ uses"))
## Panel I carries no colour encoding: the facet strips already name the two controls, and
## the only legend in Figure 4 belongs to panel B, where colour means "with / without school
## fixed effects". Colouring panel A would put a second meaning on those same colours.
pF4a <- ggplot(nc, aes(x = exposure, y = beta)) +
  geom_point(size = 1.8, colour = "grey20") +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.12, linewidth = 0.25,
                colour = "grey20") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Adjusted difference",
       title = "A") +
  theme_ms() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

t16 <- read.csv(file.path(out, "table16_school_FE.csv"))
t16$outcome <- setO(t16$outcome)
t16$exposure <- factor(t16$exposure, levels = c("time_gt3h", "exp_harassed", "ai_any"),
                       labels = c("Screen >3 h/d", "Online harassment", "Any genAI use"))
dJ <- t16 %>% select(exposure, outcome, rd_noFE, rd_FE) %>%
  pivot_longer(c(rd_noFE, rd_FE), names_to = "model", values_to = "rd") %>%
  mutate(model = factor(model, levels = c("rd_noFE", "rd_FE"),
                        labels = c("Without school FE", "With school FE")))
pF4b <- ggplot(dJ, aes(x = rd, y = outcome, colour = model)) +
  geom_line(aes(group = outcome), colour = "grey55", linewidth = 0.35) +
  geom_point(size = 1.7) +
  facet_wrap(~ exposure, scales = "free_x", nrow = 1) +
  scale_x_continuous(n.breaks = 2, labels = function(x) sprintf("%.02f", x)) +
  scale_colour_manual(values = okabe[c(1, 6)], name = "Model") +
  labs(x = "Risk difference (proportion)", y = NULL,
       title = "B") +
 ## every facet has its own x range, so the outermost tick labels of neighbouring
 ## facets would otherwise touch; widen the facet gap
  theme_ms() + theme(legend.position = "bottom",
                     panel.spacing.x = unit(7, "mm"))
## ---------------------------------------------------------------------------
## Figure 3 = the methodological claim only. Panels F/G/H are the three layers of the
## benchmark (causal / heterogeneity / prediction), i.e. the answer to "why does the
## machine-learning pipeline matter?". F carries facets, so it takes the full width on
## row 1; G and H are single-panel plots and share row 2.
## This used to be bundled with the validity checks (now Figure 4) in one figure, which forced a
## single caption to serve two unrelated arguments.
## ---------------------------------------------------------------------------
save_combo(list(pF3a, pF3b, pF3c), matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
           "Figure3_method_benchmark", 180, 141,
           heights = unit(c(1.0, 1.05), "null"))

## ---------------------------------------------------------------------------
## Figure 4 = validity, as one self-contained "warning and resolution" pair: panel A
## raises a confounding warning (the negative control), panel B resolves it (school fixed
## effects). Both carry facets, so each takes a full-width row.
## ---------------------------------------------------------------------------
save_combo(list(pF4a, pF4b), matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE),
           "Figure4_validity_checks", 180, 113,
           heights = unit(c(1.0, 1.0), "null"))

## ============================ Supplementary ============================
tr <- read.csv(file.path(out, "table8_trend_prevalence.csv"))
tr$wave <- factor(tr$wave, levels = c("2022-23", "2023-24", "2024-25"))
tr$outcome <- setO(tr$variable)
pS1 <- ggplot(tr, aes(x = wave, y = pct, group = outcome, colour = outcome)) +
  geom_line(linewidth = 0.35) + geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = pct - 1.96 * se, ymax = pct + 1.96 * se),
                width = 0.12, linewidth = 0.2) +
  scale_colour_manual(values = okabe, name = NULL) +
  labs(x = NULL, y = "Weighted prevalence (%)",
       title = NULL) +
  theme_ms() + theme(legend.position = "bottom")
save_combo(list(pS1), matrix(1, 1, 1), "FigureS1_three_wave_trend", 150, 85)

## raw R names (fin_stress, race_catBlack, ai_n_uses ...) must never reach the axis of a
## submitted figure; every predictor is mapped to a label, and an unmapped name is an error
sh_lab <- c(fin_stress = "Financial stress", time5cat = "Screen time (non-academic)",
            food_insec = "Food insecurity", exp_harassed = "Online harassment",
            female = "Female sex", age_num = "Age (years)",
            exp_relation = "Online relationship", exp_activism = "Online activism",
            ai_n_uses = "No. of genAI uses", ai_academic = "GenAI: schoolwork",
            grad_student = "Graduate student", race_catBlack = "Race: Black",
            ai_health = "GenAI: health", undergrad = "Undergraduate",
            ai_work = "GenAI: work", ai_fun = "GenAI: entertainment",
            race_catMultiracial = "Race: multiracial", firstgen = "First-generation",
            exp_create = "Content creation")
sh <- read.csv(file.path(out, "table5_shap_importance.csv")) %>%
  slice_max(mean_abs_shap, n = 12) %>%
  mutate(lab = unname(ifelse(variable %in% names(sh_lab), sh_lab[variable], NA_character_)))
if (any(is.na(sh$lab)))
  stop("unmapped SHAP predictor(s): ", paste(sh$variable[is.na(sh$lab)], collapse = ", "))
sh <- mutate(sh, lab = factor(lab, levels = rev(lab)))
pS2a <- ggplot(sh, aes(x = lab, y = mean_abs_shap)) +
  geom_col(fill = okabe[2], colour = "black", width = 0.7, linewidth = 0.15) +
  coord_flip() + labs(x = NULL, y = "Mean |SHAP| (log-odds)",
                      title = "a") + theme_ms()

dca <- read.csv(file.path(out, "table5_dca.csv"))
dc <- dca %>% pivot_longer(starts_with("nb_"), names_to = "model", values_to = "nb") %>%
  mutate(model = factor(model, levels = c("nb_treat_all", "nb_A", "nb_B"),
                        labels = c("Treat all", "Model A (base)", "Model B (+ digital)")))
pS2b <- ggplot(dc, aes(x = threshold, y = nb, colour = model)) +
  geom_line(linewidth = 0.35) +
  scale_colour_manual(values = okabe[c(8, 1, 2)], name = NULL) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "b") +
  theme_ms() + theme(legend.position = "bottom")
save_combo(list(pS2a, pS2b), matrix(c(1, 2), 1, 2, byrow = TRUE),
           "FigureS2_shap_dca", 180, 85)

## ---- Figure S3: generalisation (internal replication + external validation) ----
## These two results previously had no figure at all even though they answer the
## first external-validity question a reader will ask.
rep6 <- read.csv(file.path(out, "table6_replication_2425.csv"))
repA <- rep6 %>% filter(level == "time5cat_f>3h") %>%
  mutate(exposure = "Screen >3 h/d", outcome = setO(outcome))
repB <- rep6 %>% filter(level == "exp_harassed") %>%
  mutate(exposure = "Online harassment", outcome = setO(outcome))
rp <- bind_rows(repA, repB) %>%
  mutate(exposure = factor(exposure, levels = c("Screen >3 h/d", "Online harassment")))
pS3a <- ggplot(rp, aes(x = outcome, y = OR, colour = exposure)) +
  geom_point(position = position_dodge(0.4), size = 1.5) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.14,
                position = position_dodge(0.4), linewidth = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.25) +
 ## the same two exposures are coloured in Figure 2, panel A: screen time is blue and
 ## harassment vermillion there too, so an exposure keeps one colour across the manuscript
  scale_colour_manual(values = okabe[c(2, 6)], name = NULL) +
  scale_y_log10(breaks = c(1, 1.5, 2, 2.5)) +
  labs(x = NULL, y = "Adjusted odds ratio (log scale)",
       title = "a") +
  theme_ms() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     legend.position = "bottom")

cf <- read.csv(file.path(out, "table7_cfps_validation.csv"))
cf <- cf %>% filter(level != "mobile_h") %>%
  mutate(exposure = factor(level,
                           levels = c("sv_watch", "mobile_g2-4h", "mobile_g4-6h", "mobile_g>6h"),
                           labels = c("Short video", "Mobile 2-4 h/d", "Mobile 4-6 h/d",
                                      "Mobile >6 h/d")),
         outcome = factor(outcome, levels = c("dep8", "dep10", "lonely"),
                          labels = c("Depression (CES-D8 >= 8)",
                                     "Depression (CES-D8 >= 10)", "Loneliness")))
pS3b <- ggplot(cf, aes(x = exposure, y = OR, colour = outcome)) +
  geom_point(position = position_dodge(0.4), size = 1.5) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.14,
                position = position_dodge(0.4), linewidth = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.25) +
  scale_colour_manual(values = okabe[c(1, 6, 2)], name = NULL) +
  scale_y_log10(breaks = c(0.5, 1, 2)) +
  labs(x = NULL, y = "Adjusted odds ratio (log scale)",
       title = "b") +
  theme_ms() + theme(axis.text.x = element_text(angle = 25, hjust = 1),
                     legend.position = "bottom")
save_combo(list(pS3a, pS3b), matrix(c(1, 2), 1, 2, byrow = TRUE),
           "FigureS3_generalisation", 180, 85)

cat("FINAL FIGURES DONE\n")
