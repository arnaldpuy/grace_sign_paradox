# 03_calibration.R -- Phase 3: reproduce-reality gate (METHODOLOGY.md 5,
# PREREGISTRATION.md gates G1-G4).
#
# Compares the SYNTHETIC per-recipe trend ensemble (recovered by the
# unchanged pipeline from the Control Earth pilot) against the REAL one on
# four pre-registered signature statistics. Recovery scores must not be
# interpreted unless all four gates pass (gate thresholds frozen in
# PREREGISTRATION.md before the recovery run).
#
# Run AFTER 04_recover.py completes:
#   Rscript code/03_calibration.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table"))

# PATHS ########################################################################

control_root <- here::here("control_earth")
real_csv <- file.path(here::here(),
                      "datasets/output/tier_a_full",
                      "jasechko_per_recipe_trends.csv")
synth_dir <- file.path(control_root, "datasets/output/tier_a_synth")
out_csv <- file.path(control_root, "datasets/output/calibration_report.csv")

# LOAD ENSEMBLES ###############################################################

# Real 73-aquifer x 1,920-recipe ensemble ------------------------------------

real_dt <- fread(real_csv)

# Synthetic ensemble: one CSV per centre (D006) -------------------------------

synth_files <- list.files(synth_dir, pattern = "jasechko_per_recipe_trends",
                          recursive = TRUE, full.names = TRUE)
stopifnot(length(synth_files) == 3)
synth_dt <- rbindlist(lapply(synth_files, fread))
setnames(synth_dt, tolower(names(synth_dt)))
setnames(real_dt, tolower(names(real_dt)))

# SIGNATURE STATISTICS #########################################################

# Shared metric code path: the same function runs on both cubes (METHODOLOGY
# 11 design invariant), so the gate is apples-to-apples by construction.

signature_stats_fun <- function(dt) {

  # G1: per-aquifer recipe IQR --------------------------------------------

  iqr_dt <- dt[is.finite(trend_cm_yr),
               .(iqr = IQR(trend_cm_yr), n_cells = .N),
               by = aquifer_idx]

  # G2: filter-axis variance share (between-filter / total, per aquifer) ---

  g2_dt <- dt[is.finite(trend_cm_yr),
              .(share = {
                total_var <- var(trend_cm_yr)
                grp_mean <- .SD[, mean(trend_cm_yr), by = filter]$V1
                grp_n <- .SD[, .N, by = filter]$N
                between <- sum(grp_n * (grp_mean - mean(trend_cm_yr))^2) /
                  (.N - 1)
                between / total_var
              }), by = aquifer_idx]

  # G3: inter-centre spread (SD across centre medians, per aquifer) --------

  g3_dt <- dt[is.finite(trend_cm_yr),
              .(med = median(trend_cm_yr)), by = .(aquifer_idx, centre)][
                , .(centre_sd = sd(med)), by = aquifer_idx]

  list(iqr = iqr_dt, filter_share = g2_dt, centre_sd = g3_dt)
}

real_sig <- signature_stats_fun(real_dt)
synth_sig <- signature_stats_fun(synth_dt)

# G4: area-scaling exponent of the recipe spread ------------------------------

# n_cells per aquifer proxies area identically on both sides; the exponent is
# the OLS slope of log(IQR) on log(n_cells). A pure-leakage world gives -0.5;
# reality gives about -0.085 on the area axis.

cells_dt <- fread(file.path(here::here(),
                            "datasets/jasechko_2024/cell_masks_0p5deg.csv"))[
                              , .(n_cells = .N), by = aquifer_idx]
setkey(cells_dt, aquifer_idx)

expo_fun <- function(iqr_dt) {
  dt <- cells_dt[iqr_dt[, .(aquifer_idx, iqr)], on = .(aquifer_idx)]
  coef(lm(log(iqr) ~ log(n_cells), data = dt))[2]
}

# GATES ########################################################################

g1_ratio <- median(synth_sig$iqr$iqr) / median(real_sig$iqr$iqr)
g2_real <- mean(real_sig$filter_share$share)
g2_synth <- mean(synth_sig$filter_share$share)
g3_ratio <- median(synth_sig$centre_sd$centre_sd) /
  median(real_sig$centre_sd$centre_sd)
g4_synth <- expo_fun(synth_sig$iqr)
g4_real <- expo_fun(real_sig$iqr)

report_dt <- data.table(
  gate = c("G1_iqr_ratio", "G2_filter_share", "G3_centre_ratio",
           "G4_area_exponent"),
  real = c(1, g2_real, 1, g4_real),
  synthetic = c(g1_ratio, g2_synth, g3_ratio, g4_synth),
  criterion = c("0.5-1.5", sprintf("%.2f +/- 0.15", g2_real), "0.5-1.5",
                "-0.20 to -0.02"),
  pass = c(g1_ratio >= 0.5 & g1_ratio <= 1.5,
           abs(g2_synth - g2_real) <= 0.15,
           g3_ratio >= 0.5 & g3_ratio <= 1.5,
           g4_synth >= -0.20 & g4_synth <= -0.02))
fwrite(report_dt, out_csv)
print(report_dt)
cat("\nREPRODUCE-REALITY GATE:",
    if (all(report_dt$pass)) "PASS" else "FAIL -- retune D015, re-emit, re-run",
    "\n")
