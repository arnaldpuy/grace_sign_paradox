# 07b_rolling_null.R -- the rolling-window null and the mission-window
# paradox decomposition (D033).
#
# The synthetic planets have LINEAR truths (trend + seasonal), so any
# window-to-window sign flip of the recipe-median trend, and any change in
# the sign-paradox share between mission windows, is pure noise + recipe
# arithmetic. This script classifies the synthetic planets with EXACTLY the
# real-side code (the `rolling-windows` chunk of
# code_grace_uncertainty_clean.Rmd; recipe_subspace_metrics_fun for the
# mission windows) and puts the real numbers alongside, recomputed from the
# real artefacts restricted to the same 73 aquifers.
#
#   Rscript code/07b_rolling_null.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "ggplot2", "cowplot",
                          "scales"))

# Paths (here() resolves to $HOME in this repo -- use absolute paths) ----------

ce_root <- here::here("control_earth")
ps_repo <- here::here()
ps_code <- here::here()
rn_dir <- file.path(ce_root, "datasets/output/rolling_null")

# Source the paper's own metric functions (same-code-path invariant, D020) ----

source(file.path(ps_repo, "functions/recipe_subspace_metrics_fun.R"))
source(file.path(ps_repo, "functions/theme_AP.R"))

# ROLLING-WINDOW CLASSIFICATION ###############################################

# Mirror of the real-side per-aquifer classification: win_flip when the sign
# of the recipe-MEDIAN trend differs across the 15 windows; spans_full when
# the full-record ensemble spans zero.

classify_rolling_fun <- function(rw, full_trends) {
  rw_aq <- rw[, .(
    n_win = .N,
    win_flip = uniqueN(sign(med_trend_cmyr[med_trend_cmyr != 0])) > 1L,
    frac_win_spans = mean(spans),
    n_win_consensus = sum(!spans)), by = Study_area]
  rw_aq <- merge(rw_aq,
                 full_trends[!is.na(trend_cm_yr), .(
                   spans_full = (min(trend_cm_yr) < 0) & (max(trend_cm_yr) > 0)),
                   by = Study_area], by = "Study_area")
  rw_aq
}

summarise_rolling_fun <- function(rw_aq, label) {
  data.table(
    dataset = label,
    n_aq = nrow(rw_aq),
    win_flip = rw_aq[win_flip == TRUE, .N],
    both = rw_aq[win_flip & spans_full, .N],
    win_only = rw_aq[win_flip & !spans_full, .N],
    recipe_only = rw_aq[!win_flip & spans_full, .N],
    stable = rw_aq[!win_flip & !spans_full, .N],
    med_frac_win_spans = median(rw_aq$frac_win_spans))
}

# The 73 aquifer names, taken from the synthetic series ------------------------

syn_names <- unique(fread(file.path(rn_dir,
  "rolling_window_signs_graded_floor_v3.csv"))$Study_area)
stopifnot(length(syn_names) == 73L)

# Real Earth -------------------------------------------------------------------

rw_real <- fread(file.path(ps_repo, "datasets/output/tier_a",
                           "rolling_window_signs.csv"))
rw_real73 <- rw_real[Study_area %in% syn_names]
stopifnot(uniqueN(rw_real73$Study_area) == 73L)
real_full <- fread(file.path(ps_code, "datasets/output/tier_a_full",
                             "jasechko_per_recipe_trends.csv"))
real_full73 <- real_full[Study_area %in% syn_names]
res_rolling <- summarise_rolling_fun(
  classify_rolling_fun(rw_real73, real_full73), "real_earth")

# Synthetic planets ------------------------------------------------------------

planets <- c("pilot_uniform_depletion_v2", "graded_floor_v2", "graded_floor_v3")
for (p in planets) {
  rw_p <- fread(file.path(rn_dir, sprintf("rolling_window_signs_%s.csv", p)))
  full_p <- fread(file.path(rn_dir, sprintf("window_trends_%s_w_full.csv", p)))
  stopifnot(uniqueN(rw_p$Study_area) == 73L)
  res_rolling <- rbind(res_rolling,
                       summarise_rolling_fun(classify_rolling_fun(rw_p, full_p), p))
}
fwrite(res_rolling, file.path(rn_dir, "rolling_null_summary.csv"))

cat("=== Rolling-window sign stability: real vs synthetic (73 aquifers x 15 windows) ===\n")
print(res_rolling, digits = 3)

# MISSION-WINDOW PARADOX DECOMPOSITION ########################################

# Sign-paradox share per mission window, real (recomputed on the same 73)
# vs each synthetic planet, all through recipe_subspace_metrics_fun.

real_win_src <- list(
  w_full = file.path(ps_code, "datasets/output/tier_a_full",
                     "jasechko_per_recipe_trends.csv"),
  w_grace = file.path(ps_code, "datasets/output/tier_a_grace_only",
                      "jasechko_per_recipe_trends.csv"),
  w_fo = file.path(ps_code, "datasets/output/tier_a_276",
                   "jasechko_per_recipe_trends.csv"))

window_metrics_fun <- function(prt, label, win) {
  m <- recipe_subspace_metrics_fun(prt[Study_area %in% syn_names])
  stopifnot(nrow(m) == 73L)
  data.table(dataset = label, window = win, n_aq = nrow(m),
             n_paradox = sum(m$spans_zero),
             spans_zero_pct = 100 * mean(m$spans_zero),
             paradox10_pct = 100 * mean(m$paradox_gt10))
}

res_win <- rbindlist(lapply(names(real_win_src), function(w) {
  f <- real_win_src[[w]]
  if (!file.exists(f)) f <- paste0(f, ".gz")
  window_metrics_fun(fread(f), "real_earth", w)
}))
for (p in planets) {
  res_win <- rbind(res_win, rbindlist(lapply(
    c("w_full", "w_grace", "w_fo"), function(w) {
      prt <- fread(file.path(rn_dir, sprintf("window_trends_%s_%s.csv", p, w)))
      window_metrics_fun(prt, p, w)
    })))
}
fwrite(res_win, file.path(rn_dir, "window_paradox_summary.csv"))

cat("\n=== Sign-paradox share by mission window: real vs synthetic ===\n")
print(res_win, digits = 3)

# WINDOW DISPERSION DIAGNOSTIC ################################################

# Separates the two possible drivers of a paradox-share difference on a
# window: ensemble spread (median per-aquifer recipe IQR) vs signal strength
# (median per-aquifer |ensemble median|). The calibration gate (G1-G4)
# matched FULL-RECORD spread only; window-level spread was never a target,
# so this diagnostic is the honesty check on the null.

dispersion_fun <- function(prt, label, win) {
  d <- prt[Study_area %in% syn_names][!is.na(trend_cm_yr),
    .(iqr = IQR(trend_cm_yr), amed = abs(median(trend_cm_yr))),
    by = Study_area]
  data.table(dataset = label, window = win,
             med_iqr = median(d$iqr), med_absmed = median(d$amed))
}

res_disp <- rbindlist(lapply(names(real_win_src), function(w) {
  f <- real_win_src[[w]]
  if (!file.exists(f)) f <- paste0(f, ".gz")
  dispersion_fun(fread(f), "real_earth", w)
}))
for (p in planets) {
  res_disp <- rbind(res_disp, rbindlist(lapply(
    c("w_full", "w_grace", "w_fo"), function(w) {
      prt <- fread(file.path(rn_dir, sprintf("window_trends_%s_%s.csv", p, w)))
      dispersion_fun(prt, p, w)
    })))
}
fwrite(res_disp, file.path(rn_dir, "window_dispersion_summary.csv"))

cat("\n=== Window dispersion: spread (IQR) vs signal (|median|) ===\n")
print(res_disp, digits = 3)
