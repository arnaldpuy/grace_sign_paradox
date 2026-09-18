# 08b_val_of_val.R -- validation-of-validation (Phase 7, D034).
#
# The field's gold-standard check validates GRACE groundwater trends by sign
# agreement with in-situ wells (MS16: 68% agreement above the detectability
# boundary, 51% below). On a planet with known truth we can ask what that
# check is worth: synthetic wells sample the TRUE groundwater layer at the
# real Jasechko well locations, converted to depth-to-water via a prescribed
# specific yield, with per-well noise calibrated to the real within-aquifer
# trend spread, subsampled to the real study's per-aquifer well counts. The
# well-agreement test is then re-run with the real-side code (the
# well-sign-validation chunk of code_robustness_analysis_clean.Rmd), and
# every verdict is scored against the known truth.
#
#   Rscript code/08b_val_of_val.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "ggplot2", "cowplot",
                          "scales"))

# Paths (here() resolves to $HOME in this repo -- use absolute paths) ----------

ce_root <- here::here("control_earth")
ps_repo <- here::here()
ps_code <- here::here()
vv_dir <- file.path(ce_root, "datasets/output/val_of_val")
rn_dir <- file.path(ce_root, "datasets/output/rolling_null")

# Source the paper's own metric + boundary functions (D020 invariant) ---------

source(file.path(ps_repo, "functions/recipe_subspace_metrics_fun.R"))
source(file.path(ps_repo, "functions/detectability_boundary_fun.R"))
source(file.path(ps_repo, "functions/theme_AP.R"))

# Fixed design (registered, D034) ---------------------------------------------

planets <- c("pilot_uniform_depletion_v2", "graded_floor_v2", "graded_floor_v3")
sy_levels <- c(0.05, 0.15, 0.25)
sigma_well <- fread(file.path(vv_dir, "well_noise_calibration.csv"))[
  n_wells >= 10, median(sd_trend_m_yr)]
n_mc <- 200L
seed_headline <- 71L

# Real per-aquifer well counts (the Jasechko table the real test used) --------

cohort <- fread(file.path(ps_code, "datasets/jasechko_2024/cohort_50K_n73.csv"))

# ONE MONTE-CARLO DRAW ########################################################

# Synthesize the well table for one (planet, Sy, seed): per-well depth-to-
# water trend = -(true gw trend)/(100 Sy) + noise, aggregated per aquifer as
# the median over a subsample of the real study's n_wells, rounded to the
# 0.01 m/yr precision of the real table. Then the real-side verdicts.

draw_wells_fun <- function(wells, sy, seed) {
  set.seed(seed)
  w <- wells[, .(true_gw = true_gw_trend_cm_yr,
                 n_avail = .N,
                 i = seq_len(.N)), by = Study_area]
  w <- merge(w, cohort[, .(Study_area, n_target = n_wells)], by = "Study_area")
  w[, keep := i %in% sample(n_avail[1], min(n_target[1], n_avail[1])),
    by = Study_area]
  w <- w[keep == TRUE]
  w[, dtw_trend := -(true_gw / 100) / sy + rnorm(.N, 0, sigma_well)]
  w[, .(trend_m_yr = round(median(dtw_trend), 2), n_wells = .N),
    by = Study_area]
}

score_draw_fun <- function(well_aq, grace, truth, b_w) {
  v <- merge(grace, well_aq, by = "Study_area")
  v <- merge(v, truth, by = "Study_area")
  v <- v[n_wells >= 10 & trend_m_yr != 0]
  v[, well_storage_sign := -sign(trend_m_yr)]
  v[, agree := sign(median_trend) == well_storage_sign]
  v[, grace_detectable := abs_median >= b_w]
  v[, grace_true := sign(median_trend) == sign(tau_true_cm_yr)]
  v[, well_true := well_storage_sign == sign(tau_true_cm_yr)]
  v[, well_true_gw := well_storage_sign == sign(tau_gw_cm_yr)]
  v
}

summarise_draw_fun <- function(v) {
  data.table(
    n_aq = nrow(v),
    agree_above_pct = 100 * v[grace_detectable == TRUE, mean(agree)],
    agree_below_pct = 100 * v[grace_detectable == FALSE, mean(agree)],
    n_above = v[grace_detectable == TRUE, .N],
    false_valid_pct = 100 * v[agree == TRUE, mean(!grace_true)],
    missed_valid_pct = 100 * v[agree == FALSE, mean(grace_true)],
    well_wrong_pct = 100 * v[, mean(!well_true)],
    well_wrong_gw_pct = 100 * v[, mean(!well_true_gw)])
}

# RUN #########################################################################

res <- list()
res_num <- list()
for (p in planets) {
  wells <- fread(file.path(vv_dir, sprintf("synth_wells_%s.csv", p)))
  truth <- fread(file.path(vv_dir, sprintf("truth_gw_aquifer_%s.csv", p)))
  prt <- fread(file.path(rn_dir, sprintf("window_trends_%s_w_full.csv", p)))
  grace <- recipe_subspace_metrics_fun(prt)
  b_w <- detectability_boundary_fun(grace$abs_median, grace$dominance_R,
                                    type = "R1")$est
  for (sy in sy_levels) {
    mc <- rbindlist(lapply(seq_len(n_mc), function(k) {
      well_aq <- draw_wells_fun(wells, sy, seed = 1000L * k + round(100 * sy))
      summarise_draw_fun(score_draw_fun(well_aq, grace, truth, b_w))
    }))
    q <- function(x) sprintf("%.0f [%.0f, %.0f]", mean(x, na.rm = TRUE),
                             quantile(x, 0.025, na.rm = TRUE),
                             quantile(x, 0.975, na.rm = TRUE))
    res[[paste(p, sy)]] <- data.table(
      planet = p, sy = sy, boundary = round(b_w, 3),
      n_aq = round(mean(mc$n_aq)), n_above = round(mean(mc$n_above)),
      agree_above = q(mc$agree_above_pct),
      agree_below = q(mc$agree_below_pct),
      false_valid = q(mc$false_valid_pct),
      missed_valid = q(mc$missed_valid_pct),
      well_wrong = q(mc$well_wrong_pct),
      well_wrong_gw = q(mc$well_wrong_gw_pct))
    mcl <- melt(mc, measure.vars = setdiff(names(mc), c("n_aq", "n_above")),
                variable.name = "metric")
    res_num[[paste(p, sy)]] <- mcl[, .(
      planet = p, sy = sy,
      mean = mean(value, na.rm = TRUE),
      lo = quantile(value, 0.025, na.rm = TRUE),
      hi = quantile(value, 0.975, na.rm = TRUE)), by = metric]
  }
}
res <- rbindlist(res)
res_num <- rbindlist(res_num)
fwrite(res, file.path(vv_dir, "val_of_val_summary.csv"))
fwrite(res_num, file.path(vv_dir, "val_of_val_mc_numeric.csv"))

cat("=== Validation-of-validation: synthetic wells vs known truth ===\n")
cat(sprintf("sigma_well = %.3f m/yr (real within-aquifer SD, calibrated); %d MC draws\n\n",
            sigma_well, n_mc))
print(res)

# Headline per-aquifer table (v3, Sy = 0.15, fixed seed) ----------------------

p <- "graded_floor_v3"
wells <- fread(file.path(vv_dir, sprintf("synth_wells_%s.csv", p)))
truth <- fread(file.path(vv_dir, sprintf("truth_gw_aquifer_%s.csv", p)))
prt <- fread(file.path(rn_dir, sprintf("window_trends_%s_w_full.csv", p)))
grace <- recipe_subspace_metrics_fun(prt)
b_w <- detectability_boundary_fun(grace$abs_median, grace$dominance_R,
                                  type = "R1")$est
v <- score_draw_fun(draw_wells_fun(wells, 0.15, seed_headline),
                    grace, truth, b_w)
fwrite(v, file.path(vv_dir, "val_of_val_aquifer_graded_floor_v3.csv"))

cat("\n=== Headline draw (graded_floor_v3, Sy = 0.15, seed 71) ===\n")
ab <- v[grace_detectable == TRUE]; bl <- v[grace_detectable == FALSE]
bt_a <- binom.test(sum(ab$agree), nrow(ab), 0.5, alternative = "greater")
bt_b <- binom.test(sum(bl$agree), nrow(bl), 0.5, alternative = "greater")
cat(sprintf("above boundary: %d/%d agree (%.0f%%), binomial p = %.3g\n",
            sum(ab$agree), nrow(ab), 100 * mean(ab$agree), bt_a$p.value))
cat(sprintf("below boundary: %d/%d agree (%.0f%%), binomial p = %.3g\n",
            sum(bl$agree), nrow(bl), 100 * mean(bl$agree), bt_b$p.value))
cat(sprintf("wells validated a wrong-sign GRACE trend: %d of %d agreements\n",
            v[agree == TRUE & !grace_true, .N], v[agree == TRUE, .N]))
cat(sprintf("wells rejected a right-sign GRACE trend:  %d of %d disagreements\n",
            v[agree == FALSE & grace_true, .N], v[agree == FALSE, .N]))
cat(sprintf("well aggregate itself wrong on total-truth sign: %d/%d aquifers\n",
            v[!(well_true), .N], nrow(v)))

# FIGURE (METHODOLOGY Fig 6) ##################################################

# Validation-of-validation. a) Sign agreement between the recipe-median GRACE
# trend and the synthetic wells, above and below each planet's own R = 1
# boundary, as a function of the assumed specific yield (points: MC mean,
# bars: 2.5-97.5% over 200 draws). Horizontal dashed lines: the REAL-Earth
# values (68% above, 51% below the boundary, MS16). b) Error decomposition
# on the flagship graded planet: when GRACE and wells disagree, GRACE is
# usually the one that is right (false rejection); when they agree, GRACE is
# almost never wrong (false validation); the well aggregate itself carries
# the wrong storage sign for a growing share of aquifers as the specific
# yield weakens the head signal.

pd_a <- res_num[metric %in% c("agree_above_pct", "agree_below_pct")]
pd_a[, side := fifelse(metric == "agree_above_pct", "above boundary",
                       "below boundary")]
pd_a[, planet_lab := factor(planet,
  levels = c("pilot_uniform_depletion_v2", "graded_floor_v2",
             "graded_floor_v3"),
  labels = c("Uniform depletion", "Mixed signs", "Graded magnitudes"))]

plot_agree <- ggplot(pd_a, aes(x = factor(sy), y = mean, color = side)) +
  geom_hline(yintercept = 68, linetype = "dashed", linewidth = 0.3,
             color = "#0072B2") +
  geom_hline(yintercept = 51, linetype = "dashed", linewidth = 0.3,
             color = "#D55E00") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.25, linewidth = 0.4,
                  position = position_dodge(width = 0.45)) +
  facet_wrap(~planet_lab) +
  scale_color_manual(values = c("above boundary" = "#0072B2",
                                "below boundary" = "#D55E00"),
                     name = NULL) +
  labs(x = "Assumed specific yield", y = "GRACE-well sign agreement (%)") +
  theme_AP() +
  theme(legend.position = "top")

pd_b <- res_num[planet == "graded_floor_v3" &
                metric %in% c("false_valid_pct", "missed_valid_pct",
                              "well_wrong_pct")]
pd_b[, metric_lab := factor(metric,
  levels = c("missed_valid_pct", "well_wrong_pct", "false_valid_pct"),
  labels = c("GRACE right when\nwells disagree",
             "well aggregate wrong\non true sign",
             "GRACE wrong when\nwells agree"))]

plot_err <- ggplot(pd_b, aes(x = factor(sy), y = mean, color = metric_lab)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.25, linewidth = 0.4,
                  position = position_dodge(width = 0.45)) +
  scale_color_manual(values = c("#6a3d9a", "#e31a1c", "#33a02c"),
                     name = NULL, guide = guide_legend(ncol = 1)) +
  labs(x = "Assumed specific yield", y = "Rate (%)") +
  theme_AP() +
  theme(legend.position = "top", legend.justification = "left")

plot_vv <- plot_grid(plot_agree, plot_err, ncol = 2L,
                     rel_widths = c(0.62, 0.38), labels = "auto")
ggsave(file.path(ce_root, "figures/fig_val_of_val.pdf"), plot_vv,
       width = 5.5, height = 2.9)
cat("\nwrote figures/fig_val_of_val.pdf\n")
