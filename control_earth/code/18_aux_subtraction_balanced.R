# 18_aux_subtraction_balanced.R -- does land-model error interact with sign
# geometry? (extension to AMENDMENT_01; D060)
#
# The crossed planet (D045) scored the groundwater subtraction on an
# ALL-DEPLETING key (73 negative / 0 positive), where "recover the sign" is
# arithmetically "return a negative number" and a more-negative estimator wins
# for free. This script repeats the test on the three BALANCED layouts, whose
# groundwater keys are ~37/36, so a directional bias no longer pays.
#
# Setup differs from the crossed planet and the asymmetry matters:
#   crossed planet   background ERA5-Land -> subtracting GLDAS is MISMATCHED
#   balanced planets background GLDAS     -> subtracting ERA5  is MISMATCHED
# so here the GLDAS arm is the inverse crime (zero model error by
# construction, verified to 1.3e-09 cm/yr) and the informative comparison is
#   none   recovered TWS, no subtraction
#   gldas  minus the matched auxiliary  (upper bound: filter mismatch only)
#   era5   minus the MISMATCHED auxiliary (what land-model error costs)
#
# The auxiliary trend fields are properties of the land-model datasets over the
# fixed 250 recovery months, NOT of the planet: aux_gldas is identical to
# 0.00e+00 cm/yr across the crossed and balanced planets, so the ERA5-Land
# field banked for the crossed planet is reused here rather than recomputed.
#
#   Rscript code/18_aux_subtraction_balanced.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("data.table"))

control_root <- here::here("control_earth")
setwd(control_root)

# The ERA5-Land auxiliary, from the crossed planet's banked table -------------

era5 <- fread("datasets/truth/truth_gw_aux_trends_crossed_era5_v1.csv")
setnames(era5, tolower(names(era5)))
era5 <- era5[, .(aquifer_idx, aux_era5 = aux_era5_cm_yr)]

planets <- c("graded_balanced_v1", "graded_balanced_v2", "graded_balanced_v3")

score_planet_fun <- function(planet) {
  aux <- fread(file.path("datasets/truth",
                         paste0("truth_gw_aux_trends_", planet, ".csv")))
  setnames(aux, tolower(names(aux)))
  # On a GLDAS-background planet both banked columns are GLDAS; keep one.
  aux <- aux[, .(aquifer_idx, tau_gw = tau_gw_cm_yr,
                 aux_gldas = aux_gldas_cm_yr)]
  aux <- merge(aux, era5, by = "aquifer_idx")
  stopifnot(nrow(aux) == 73L)

  files <- list.files(file.path("datasets/output",
                                paste0("tier_a_synth_", planet)),
                      pattern = "jasechko_per_recipe_trends",
                      recursive = TRUE, full.names = TRUE)
  stopifnot(length(files) == 3)
  rec <- rbindlist(lapply(files, fread))
  setnames(rec, tolower(names(rec)))
  rec <- merge(rec, aux, by = "aquifer_idx")
  rec[, `:=`(gws_none = trend_cm_yr,
             gws_gldas = trend_cm_yr - aux_gldas,
             gws_era5 = trend_cm_yr - aux_era5)]

  one <- function(col) {
    s <- rec[is.finite(get(col)),
             .(median_rec = median(get(col)),
               frac_neg = mean(get(col) < 0),
               frac_pos = mean(get(col) > 0)),
             by = .(aquifer_idx, tau_gw)]
    s[, `:=`(planet = planet, estimator = sub("^gws_", "", col),
             p_sign = fifelse(tau_gw < 0, frac_neg, frac_pos),
             bias = median_rec - tau_gw,
             median_sign_ok = sign(median_rec) == sign(tau_gw))]
    s[]
  }
  rbindlist(lapply(c("gws_none", "gws_gldas", "gws_era5"), one))
}

res <- rbindlist(lapply(planets, score_planet_fun))
fwrite(res, "datasets/output/aux_subtraction_balanced.csv")

# REPORT #######################################################################

cat("=== Groundwater key balance (contrast with the crossed planet's 73/0) ===\n")
k <- res[estimator == "none", .(neg = sum(tau_gw < 0), pos = sum(tau_gw > 0),
                                always_depleting = round(mean(tau_gw < 0), 3)),
         by = planet]
print(k)

cat("\n=== Median-sign recovery out of 73, by estimator ===\n")
tab <- dcast(res[, .(ok = sum(median_sign_ok)), by = .(planet, estimator)],
             planet ~ estimator, value.var = "ok")
setcolorder(tab, c("planet", "none", "gldas", "era5"))
print(tab)
cat("\n(gldas = MATCHED here, an inverse crime; era5 = MISMATCHED land model)\n")

cat("\n=== What the mismatched auxiliary costs (era5 - none, era5 - gldas) ===\n")
tab[, `:=`(cost_vs_none = era5 - none, cost_vs_matched = era5 - gldas)]
print(tab[, .(planet, none, gldas, era5, cost_vs_none, cost_vs_matched)])

cat("\n=== Median bias by estimator (directional-bias check) ===\n")
print(res[, .(median_bias = round(median(bias), 3)),
          by = .(planet, estimator)][order(planet, estimator)])

cat("\n=== Does the cost concentrate at low |tau_gw|? ===\n")
w <- dcast(res[estimator %in% c("gldas", "era5")],
           planet + aquifer_idx + tau_gw ~ estimator, value.var = "median_sign_ok")
w[, band := cut(abs(tau_gw), c(0, 0.25, 0.5, 0.75, 1, Inf),
                labels = c("<=0.25", "0.25-0.5", "0.5-0.75", "0.75-1", ">1"))]
print(w[, .(n = .N, matched_ok = sum(gldas), mismatched_ok = sum(era5),
            lost = sum(gldas) - sum(era5)), by = band][order(band)])

cat("\nwrote datasets/output/aux_subtraction_balanced.csv\n")
