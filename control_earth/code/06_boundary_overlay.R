# 06_boundary_overlay.R -- the headline figure (METHODOLOGY.md 7.1, Fig 4).
#
# REBUILT 2026-08-09 for the post-amendment framing (D051-D059). The previous
# version drew the boundary on three planets, a loess recovery curve on the
# graded planet, and the boundary mapped into truth space through a single
# attenuation factor. That mapping is retired: attenuation is layout-dependent
# (0.43 / 0.33 / 0.23 / 0.16), so no scalar converts a reported trend into a
# true one, and the 0.75 cm/yr crossing is a property of the depletion-
# dominated planet rather than a universal floor.
#
# a) The R = 1 detectability boundary on every planet, against the real-Earth
#    value, all computed by the same estimator (D020). Two balanced planets
#    have NO crossing -- almost no aquifer escapes the preprocessing spread --
#    and are labelled as such rather than given an interval.
# b) Model-free recovery: median P(correct sign) per band of true trend
#    magnitude, depletion-dominated planet against the three balanced ones.
#    This is the amendment's primary presentation; no fitted curve.
# c) Sign recovery and attenuation against the degree of sign interleaving.
# d) Groundwater-sign recovery under a matched land model, a mismatched one
#    and no correction at all (D060).
#
#   Rscript code/06_boundary_overlay.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "ggplot2", "cowplot",
                          "scales"))

# Project root (here() resolves to $HOME in this repo; use absolute paths) -----

ce_root <- here::here("control_earth")

# Source all .R files in the "functions" folder -------------------------------

r_functions <- list.files(path = file.path(ce_root, "functions"),
                          pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_functions, source))

# theme_AP lives with the repository functions -----------------------

source(here::here("functions/theme_AP.R"))

# DATA #########################################################################

summ <- fread(file.path(ce_root, "datasets/output/balanced_layouts_summary.csv"))
bins <- fread(file.path(ce_root, "datasets/output/balanced_layouts_bins.csv"))

# Empirical boundary (real data) ----------------------------------------------

b_real <- c(est = 0.347, lo = 0.286, hi = 0.423)

# Panel-a table: real Earth, the three original planets (banked D029/D032) and
# the three balanced layouts (banked by 17_balanced_layouts.R). The pilot and
# mixed-sign values predate the summary table and are carried as constants.

lab <- c(graded_balanced_v1 = "Balanced signs 1",
         graded_balanced_v2 = "Balanced signs 2",
         graded_balanced_v3 = "Balanced signs 3")
bal <- summ[planet %in% names(lab),
            .(planet = lab[planet], est = boundary,
              lo = boundary_lo, hi = boundary_hi)]

overlay_dt <- rbindlist(list(
  data.table(planet = "Real Earth", est = b_real[["est"]],
             lo = b_real[["lo"]], hi = b_real[["hi"]]),
  data.table(planet = "Uniform depletion", est = 0.339, lo = 0.263, hi = 0.408),
  data.table(planet = "Mixed signs", est = 0.412, lo = 0.270, hi = 0.685),
  data.table(planet = "Graded magnitudes", est = 0.299, lo = 0.254, hi = 0.342),
  bal))
overlay_dt[, real := planet == "Real Earth"]

# Explicit top-to-bottom reading order; ggplot puts the FIRST level at the
# bottom, so the vector is reversed.

row_order <- c("Real Earth", "Uniform depletion", "Mixed signs",
               "Graded magnitudes", "Balanced signs 1", "Balanced signs 2",
               "Balanced signs 3")
stopifnot(setequal(row_order, overlay_dt$planet))
overlay_dt[, planet := factor(planet, levels = rev(row_order))]

# No crossing exists on layout 3: only 7 of 73 aquifers fall below R = 1, and
# just 759 of 2,000 bootstrap draws find a crossing at all. Showing that
# interval would imply an estimate that does not exist.

no_cross <- overlay_dt[is.na(est), planet]

# PANEL A: boundary across planets #############################################

# All rows go to the scale (na.rm drops the no-crossing point, not its row), so
# a single factor drives the axis order.

plot_overlay <- ggplot(overlay_dt, aes(x = est, y = planet, color = real)) +
  geom_vline(xintercept = b_real[["est"]], linetype = "dotted",
             linewidth = 0.3, color = "#0072B2") +
  geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 0.4, size = 0.28,
                  na.rm = TRUE) +
  scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey30"),
                     guide = "none") +
  scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.35, 0.7)) +
  labs(x = bquote("cm" ~ yr^-1), y = NULL) +
  theme_AP() +
  theme(plot.margin = margin(2, 4, 2, 2))

if (length(no_cross)) {
  txt <- data.table(planet = factor(no_cross, levels = levels(overlay_dt$planet)),
                    est = 0.02)
  plot_overlay <- plot_overlay +
    geom_text(data = txt, aes(x = est, y = planet), label = "no crossing",
              hjust = 0, size = 2.1, color = "grey30", fontface = "italic",
              inherit.aes = FALSE)
}

# PANEL B: model-free recovery by magnitude band ###############################

# Individual balanced layouts stay visible as thin lines (their per-bin cells
# hold only ~15 aquifers and are correspondingly noisy); their median carries
# the comparison against the depletion-dominated planet. ASCII bin labels: the
# PDF device drops a literal "<=" glyph.

bal_ids <- c("graded_balanced_v1", "graded_balanced_v2", "graded_balanced_v3")
bd <- bins[planet %in% c("graded_floor_v3", bal_ids)]
bd[, bin := factor(bin, levels = c("<=0.25", "0.25-0.5", "0.5-0.75",
                                   "0.75-1", ">1"),
                   labels = c("0-0.25", "0.25-0.5", "0.5-0.75",
                              "0.75-1", ">1"))]
bal_med <- bd[planet %in% bal_ids,
              .(median_p_sign = median(median_p_sign)), by = bin]

plot_bins <- ggplot() +
  geom_hline(yintercept = 0.9, linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.25,
             color = "grey60") +
  geom_line(data = bd[planet %in% bal_ids],
            aes(x = bin, y = median_p_sign, group = planet),
            color = "grey70", linewidth = 0.25) +
  geom_line(data = bal_med, aes(x = bin, y = median_p_sign, group = 1,
                                color = "Balanced signs"), linewidth = 0.5) +
  geom_point(data = bal_med, aes(x = bin, y = median_p_sign,
                                 color = "Balanced signs"), size = 1.2) +
  geom_line(data = bd[planet == "graded_floor_v3"],
            aes(x = bin, y = median_p_sign, group = 1,
                color = "Graded magnitudes"), linewidth = 0.5) +
  geom_point(data = bd[planet == "graded_floor_v3"],
             aes(x = bin, y = median_p_sign, color = "Graded magnitudes"),
             size = 1.2) +
  scale_color_manual(values = c("Graded magnitudes" = "#D55E00",
                                "Balanced signs" = "grey15"), name = NULL) +
  # Full data range: one layout dips to 0.34 in the 0.5-0.75 band (n = 7),
  # which a 0.4 floor would silently clip.
  scale_y_continuous(limits = c(0.3, 1)) +
  labs(x = bquote("True trend magnitude (cm" ~ yr^-1 * ")"),
       y = "Fraction with\ncorrect sign") +
  theme_AP() +
  theme(legend.position = c(0.97, 0.03), legend.justification = c(1, 0),
        axis.text.x = element_text(size = 5.4, angle = 30, hjust = 1),
        plot.margin = margin(2, 4, 2, 7))

# PANEL C: recovery and attenuation against sign interleaving (D055, D059) ####

# The sign-imbalanced planet anchors this panel at zero interleaving: all its
# planted signs are negative, so no aquifer has an opposite-signed neighbour.
# It shares its magnitude field with the balanced planets, so the four points
# differ in sign arrangement alone.

il <- fread(file.path(ce_root, "datasets/output/sign_interleaving.csv"))
pc_ids <- c("graded_floor_v3", bal_ids)
stopifnot(il[planet == "graded_floor_v3", n_within_500km] == 0L)
pc <- merge(summ[planet %in% pc_ids,
                 .(planet, median_p_sign, attenuation)],
            il[, .(planet, n_within_500km)], by = "planet")
pc <- melt(pc, id.vars = c("planet", "n_within_500km"),
           measure.vars = c("median_p_sign", "attenuation"),
           variable.name = "metric", value.name = "value")
pc[, metric := factor(metric, levels = c("median_p_sign", "attenuation"),
                      labels = c("Correct sign", "Recovered / true"))]

plot_interleave <- ggplot(pc, aes(x = n_within_500km, y = value,
                                  color = metric, shape = metric)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.4) +
  scale_color_manual(values = c("Correct sign" = "grey15",
                                "Recovered / true" = "#0072B2"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(breaks = seq(0, 50, by = 25)) +
  scale_y_continuous(limits = c(0, 1)) +
  # ASCII only: the PDF device silently drops the masculine-ordinal glyph
  # (it rendered as "N.."), so spell out "No.".
  labs(x = "Number of aquifers with opposite-signed\nneighbour (< 500 km)",
       y = "Fraction") +
  theme_AP() +
  theme(legend.position = c(0.97, 0.03), legend.justification = c(1, 0),
        plot.margin = margin(2, 4, 2, 2))

# PANEL D: groundwater subtraction on the balanced planets (D060) #############

aux <- fread(file.path(ce_root, "datasets/output/aux_subtraction_balanced.csv"))
pd_dt <- aux[, .(p_ok = mean(median_sign_ok)), by = .(planet, estimator)]
pd_dt[, layout := factor(planet, levels = bal_ids, labels = c("1", "2", "3"))]
pd_dt[, estimator := factor(estimator, levels = c("gldas", "none", "era5"),
                            labels = c("Matched (GLDAS-Noah)",
                                       "No correction",
                                       "Mismatched (ERA5-Land)"))]

plot_aux <- ggplot(pd_dt, aes(x = layout, y = p_ok, color = estimator,
                              shape = estimator, group = estimator)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.4) +
  scale_color_manual(values = c("Matched (GLDAS-Noah)" = "#0072B2",
                                "No correction" = "grey55",
                                "Mismatched (ERA5-Land)" = "#D55E00"),
                     name = NULL) +
  scale_shape_manual(values = c(16, 15, 17), name = NULL) +
  # Headroom for the legend, which would otherwise sit on the descending lines.
  scale_y_continuous(limits = c(0.52, 1.12), breaks = c(0.6, 0.7, 0.8, 0.9)) +
  labs(x = "Balanced planet", y = "Fraction with\ncorrect sign") +
  theme_AP() +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        plot.margin = margin(2, 4, 2, 7))

# ASSEMBLE #####################################################################

# The control Earth. a) The R = 1 detectability boundary recomputed by the same
# estimator on six synthetic planets against the real-Earth value (blue, dotted
# line at 0.35); error bars are 95% confidence intervals. On the third balanced
# planet no crossing exists: only 7 of 73 aquifers fall below R = 1. b) Median
# fraction of the 1,920 recipes recovering an aquifer's planted sign, by band of
# true trend magnitude, on the depletion-dominated planet (orange) and the three
# balanced planets (greys). Dotted 0.90, dashed 0.50 (chance). Weak signals
# collapse towards chance once signs are balanced; above 1 cm/yr the planets
# agree. c) Per-recipe sign-accuracy against true magnitude on the
# depletion-dominated planet: 1,920 thin curves (grey), the 180 fully-corrected
# recipes in orange, and the ensemble median as an estimator in bold black;
# dotted 0.90. No recipe class reaches 0.90 at low magnitude on the planet whose
# sign imbalance most flatters recovery.

# Total width kept under the MS text block (15 cm ~ 5.9 in, A4 + 3 cm margins).
# Panel a keeps its share for the long planet labels; b and c compressed
# relatively (rel_widths sum to 1).

plot_main <- plot_grid(plot_overlay, plot_bins, plot_interleave, plot_aux,
                       ncol = 2L, rel_widths = c(0.5, 0.5),
                       rel_heights = c(0.5, 0.5), labels = "auto")
ggsave(file.path(ce_root, "figures/fig_headline_overlay.pdf"), plot_main,
       width = 5.2, height = 3.35)
cat("wrote figures/fig_headline_overlay.pdf\n")
