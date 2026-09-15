recipe_subspace_metrics_fun <- function(per_recipe_dt) {

  # PER-AQUIFER METRICS OVER A RECIPE SUBSPACE #################################

  # Reduce a per-(aquifer, recipe) trend table -- already filtered to the
  # recipe subspace of interest -- to one row per aquifer, carrying the

  stopifnot(data.table::is.data.table(per_recipe_dt))
  stopifnot(all(c("Study_area", "trend_cm_yr") %in% names(per_recipe_dt)))

  out <- per_recipe_dt[!is.na(trend_cm_yr),
    .(n_recipes    = .N,
      median_trend = stats::median(trend_cm_yr),
      iqr_trend    = stats::IQR(trend_cm_yr),
      min_trend    = min(trend_cm_yr),
      max_trend    = max(trend_cm_yr),
      n_pos        = sum(trend_cm_yr > 0),
      n_neg        = sum(trend_cm_yr < 0)),
    by = Study_area]

  out[, abs_median := abs(median_trend)]
  out[, dominance_R := iqr_trend / abs_median]
  out[, sign_agree := pmax(n_pos, n_neg) / n_recipes]
  out[, sign_flips := pmin(n_pos, n_neg) / n_recipes]
  out[, spans_zero := (min_trend < 0) & (max_trend > 0)]
  out[, paradox_gt05 := sign_flips >= 0.05]
  out[, paradox_gt10 := sign_flips >= 0.10]
  out[]
}
