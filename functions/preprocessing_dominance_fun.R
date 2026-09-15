preprocessing_dominance_fun <- function(recipe_trends) {

  # PREPROCESSING DOMINANCE R = IQR / |median| #################################

  # Input: array (recipe, lat, lon) of per-recipe per-cell trends.
  # Returns list with per-cell median, IQR, R, and the global p_above_one.

  med <- apply(recipe_trends, c(2L, 3L), median, na.rm = TRUE)
  q25 <- apply(recipe_trends, c(2L, 3L),
               function(x) quantile(x, 0.25, na.rm = TRUE))
  q75 <- apply(recipe_trends, c(2L, 3L),
               function(x) quantile(x, 0.75, na.rm = TRUE))
  iqr <- q75 - q25
  med_safe <- med
  med_safe[abs(med_safe) < 1e-8] <- NA_real_
  R <- iqr / abs(med_safe)

  finite <- is.finite(R)
  p_above <- mean(R[finite] > 1, na.rm = TRUE)

  list(median = med, iqr = iqr, dominance_R = R, p_above_one = p_above)
}
