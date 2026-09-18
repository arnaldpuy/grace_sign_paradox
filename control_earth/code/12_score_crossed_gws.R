# 12_score_crossed_gws.R -- score the GROUNDWATER SUBTRACTION on the
# crossed-truth planet (D045).
#
# The planet's background hydrology is ERA5-Land, so removing a GLDAS-Noah
# auxiliary is no longer an inverse crime: the subtraction error is a real
# inter-model discrepancy and the recovered groundwater sign can be scored
# against the planted answer.
#
# Three estimators of the same quantity, all on the same 250 months:
#   none    recovered TWS, no subtraction        (what the aquifer looks like
#                                                 before any correction)
#   era5    recovered TWS - ERA5-Land auxiliary  (RIGHT model; residual error
#                                                 is the filter/attenuation
#                                                 mismatch alone)
#   gldas   recovered TWS - GLDAS-Noah auxiliary (WRONG model; what the real
#                                                 analysis actually does)
# scored against tau_gw, the tapered groundwater layer's own aquifer-mean
# trend, from code/forward/12_crossed_aux_trends.py.
#
#   Rscript code/12_score_crossed_gws.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("data.table"))

# PATHS ########################################################################

control_root <- here::here("control_earth")
planet <- "crossed_era5_v1"
synth_dir <- file.path(control_root, "datasets/output",
                       paste0("tier_a_synth_", planet))
truth_csv <- file.path(control_root, "datasets/truth",
                       paste0("truth_gw_aux_trends_", planet, ".csv"))
out_csv <- file.path(control_root, "datasets/output",
                     paste0("score_gws_", planet, ".csv"))

# LOAD #########################################################################

synth_files <- list.files(synth_dir, pattern = "jasechko_per_recipe_trends",
                          recursive = TRUE, full.names = TRUE)
stopifnot(length(synth_files) == 3)
recipes_dt <- rbindlist(lapply(synth_files, fread))
setnames(recipes_dt, tolower(names(recipes_dt)))
truth_dt <- fread(truth_csv)
setnames(truth_dt, tolower(names(truth_dt)))
stopifnot(uniqueN(recipes_dt$aquifer_idx) == 73L,
          recipes_dt[, .N, by = aquifer_idx][, all(N == 1920L)],
          nrow(truth_dt) == 73L)

# SUBTRACT #####################################################################

# The auxiliary is recipe-independent, so removing it at trend level is exact
# (the linear-subtraction framing of the main analysis).

recipes_dt <- merge(recipes_dt,
                    truth_dt[, .(aquifer_idx, tau_gw = tau_gw_cm_yr,
                                 aux_era5 = aux_era5_cm_yr,
                                 aux_gldas = aux_gldas_cm_yr)],
                    by = "aquifer_idx")
recipes_dt[, `:=`(gws_none = trend_cm_yr,
                  gws_era5 = trend_cm_yr - aux_era5,
                  gws_gldas = trend_cm_yr - aux_gldas)]

# SCORE ########################################################################

score_one <- function(dt, col) {
  s <- dt[is.finite(get(col)),
          .(median_rec = median(get(col)),
            q25 = quantile(get(col), 0.25),
            q75 = quantile(get(col), 0.75),
            frac_neg = mean(get(col) < 0),
            frac_pos = mean(get(col) > 0)),
          by = .(aquifer_idx, tau_gw)]
  s[, `:=`(estimator = sub("^gws_", "", col),
           p_sign = fifelse(tau_gw < 0, frac_neg, frac_pos),
           bias = median_rec - tau_gw,
           brackets = tau_gw >= q25 & tau_gw <= q75,
           median_sign_ok = sign(median_rec) == sign(tau_gw))]
  s[]
}

score_dt <- rbindlist(lapply(c("gws_none", "gws_era5", "gws_gldas"),
                             function(x) score_one(recipes_dt, x)))
score_dt <- merge(score_dt, truth_dt[, .(aquifer_idx, study_area,
                                         tau_stipulated_cm_yr)],
                  by = "aquifer_idx")
setcolorder(score_dt, c("estimator", "aquifer_idx", "study_area", "tau_gw"))
setorder(score_dt, estimator, aquifer_idx)
fwrite(score_dt, out_csv)

# SUMMARY ######################################################################

cat(sprintf("planet %s: %d aquifers x %d recipes; groundwater truth spans %+.3f to %+.3f cm/yr\n",
            planet, nrow(truth_dt), 1920L,
            min(truth_dt$tau_gw_cm_yr), max(truth_dt$tau_gw_cm_yr)))
cat(sprintf("auxiliary discrepancy (gldas - era5): median |%.3f|, max |%.3f| cm/yr\n\n",
            median(abs(truth_dt$aux_gldas_cm_yr - truth_dt$aux_era5_cm_yr)),
            max(abs(truth_dt$aux_gldas_cm_yr - truth_dt$aux_era5_cm_yr))))

cat("                    median-sign   median   mean    bracketing\n")
cat("estimator            recovered    p_sign   bias      rate\n")
for (e in c("none", "era5", "gldas")) {
  d <- score_dt[estimator == e]
  cat(sprintf("%-18s   %2d / 73     %.3f   %+.3f     %.2f\n",
              e, d[median_sign_ok == TRUE, .N], median(d$p_sign),
              mean(d$bias), mean(d$brackets)))
}

# Does the WRONG auxiliary flip the verdict? -----------------------------------

w <- dcast(score_dt, aquifer_idx + tau_gw ~ estimator,
           value.var = c("median_rec", "median_sign_ok", "p_sign"))
flip <- w[sign(median_rec_gldas) != sign(median_rec_era5)]
cat(sprintf("\nSign of the recipe-median groundwater trend differs between the\n"))
cat(sprintf("GLDAS-Noah and ERA5-Land auxiliaries for %d of 73 aquifers.\n", nrow(flip)))
cat(sprintf("Of those, the wrong (GLDAS) auxiliary is the one that errs in %d.\n",
            flip[median_sign_ok_era5 == TRUE & median_sign_ok_gldas == FALSE, .N]))
if (nrow(flip)) print(flip[order(abs(tau_gw)),
                           .(aquifer_idx, tau_gw = round(tau_gw, 3),
                             era5 = round(median_rec_era5, 3),
                             gldas = round(median_rec_gldas, 3))])

# Where the residual error comes from -----------------------------------------

cat(sprintf("\nWith the RIGHT model (era5) the subtraction still misses the sign for\n"))
cat(sprintf("%d aquifers -- attenuation, not model error: the recovered trend is\n",
            score_dt[estimator == "era5" & median_sign_ok == FALSE, .N]))
cat(sprintf("filtered but the model auxiliary is not.\n"))
att <- score_dt[estimator == "none"][truth_dt, on = "aquifer_idx"]
cat(sprintf("Recovered/true magnitude ratio, |tau_gw| > 0.5: %.2f (median)\n",
            att[abs(tau_gw) > 0.5, median(median_rec / tau_true_cm_yr)]))
cat(sprintf("\nwrote %s\n", basename(out_csv)))
