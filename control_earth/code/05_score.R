# 05_score.R -- Phase 5: score the recovered ensemble against known truth
# (METHODOLOGY.md 6.2; predictions P1/P3 in PREREGISTRATION.md).
#
# Only interpretable AFTER 03_calibration.R reports PASS.
#
#   Rscript code/05_score.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table"))

# PATHS ########################################################################

control_root <- here::here("control_earth")

# Planet from the command line; the pilot's outputs predate the suffix convention.
.arg <- commandArgs(trailingOnly = TRUE)
planet <- if (length(.arg)) .arg[1] else "pilot_uniform_depletion_v2"
synth_dir <- file.path(control_root, "datasets/output",
                       if (planet == "pilot_uniform_depletion_v2") "tier_a_synth"
                       else paste0("tier_a_synth_", planet))
truth_csv <- file.path(control_root, "datasets/truth",
                       paste0("truth_aquifer_trends_", planet, ".csv"))
out_csv <- file.path(control_root, "datasets/output",
                     paste0("score_", planet, ".csv"))

# LOAD #########################################################################

synth_files <- list.files(synth_dir, pattern = "jasechko_per_recipe_trends",
                          recursive = TRUE, full.names = TRUE)
stopifnot(length(synth_files) == 3)
recipes_dt <- rbindlist(lapply(synth_files, fread))
setnames(recipes_dt, tolower(names(recipes_dt)))
truth_dt <- fread(truth_csv)
setkey(truth_dt, aquifer_idx)

# SCORE (METHODOLOGY 6.2) ######################################################

# Per aquifer: sign-recovery probability, bias, bracketing ---------------------

score_dt <- recipes_dt[is.finite(trend_cm_yr),
                       .(n_recipes = .N,
                         median_rec = median(trend_cm_yr),
                         q25 = quantile(trend_cm_yr, 0.25),
                         q75 = quantile(trend_cm_yr, 0.75),
                         frac_neg = mean(trend_cm_yr < 0),
                         frac_pos = mean(trend_cm_yr > 0)),
                       by = aquifer_idx]
setkey(score_dt, aquifer_idx)
score_dt <- truth_dt[, .(aquifer_idx, study_area = Study_area,
                         tau_true = tau_true_cm_yr)][score_dt]
score_dt[, p_sign := fifelse(tau_true < 0, frac_neg, frac_pos)]
score_dt[, bias := median_rec - tau_true]
score_dt[, brackets := tau_true >= q25 & tau_true <= q75]
score_dt[, sign_recovered_by_median := sign(median_rec) == sign(tau_true)]
fwrite(score_dt, out_csv)

# SUMMARY vs pre-registered predictions ###############################

cat(sprintf("planet: %s  (n = %d aquifers, %d recipes/aquifer)\n",
            planet, nrow(score_dt), score_dt$n_recipes[1]))
cat(sprintf("P1  median P_sign            : %.3f  (predicted >= 0.95)\n",
            median(score_dt$p_sign)))
cat(sprintf("    aquifers with P_sign>=.95: %d / %d\n",
            score_dt[p_sign >= 0.95, .N], nrow(score_dt)))
cat(sprintf("    median-sign recovered    : %d / %d aquifers\n",
            score_dt[sign_recovered_by_median == TRUE, .N], nrow(score_dt)))
cat(sprintf("P3  bracketing rate (IQR)    : %.2f  (predicted >= 0.8)\n",
            mean(score_dt$brackets)))
cat(sprintf("    median bias              : %+.3f cm/yr\n",
            median(score_dt$bias)))
cat(sprintf("    median |bias|            : %.3f cm/yr\n",
            median(abs(score_dt$bias))))
