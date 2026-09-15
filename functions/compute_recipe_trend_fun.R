compute_recipe_trend_fun <- function(recipes, time_slice = NULL) {

  # PER-RECIPE PER-CELL LINEAR TRENDS ###########################################

  # arr (time, recipe, lat, lon) -> slopes (recipe, lat, lon), cm/yr. NaN-aware:
  # pre-2015 BB01 months are all-NaN, so slopes fit each cell's finite months.

  arr <- recipes$arr
  times <- recipes$times
  if (!is.null(time_slice)) {
    idx <- which(times >= time_slice[1] & times <= time_slice[2])
    arr <- arr[idx, , , , drop = FALSE]
    times <- times[idx]
  }

  t_years <- as.numeric(times - times[1]) / 365.25

  dim_arr <- dim(arr)
  flat <- matrix(arr, nrow = dim_arr[1], ncol = prod(dim_arr[-1]))
  col_means <- colMeans(flat, na.rm = TRUE)
  flat_dem <- sweep(flat, 2L, col_means, FUN = "-")

  mask <- is.finite(flat_dem)
  # Centre time per recipe-cell's own months (not the global mean): avoids a
  # toward-zero bias for cells with missing months. Mirrors the Python _ols_slopes.
  x_masked <- ifelse(mask, t_years, 0)
  n_ok <- colSums(mask)
  x_bar <- colSums(x_masked) / pmax(n_ok, 1L)
  x_c <- ifelse(mask, t_years - rep(x_bar, each = dim_arr[1]), 0)
  y_filled <- ifelse(mask, flat_dem, 0)
  num <- colSums(x_c * y_filled)
  den <- colSums(x_c * x_c)
  slopes <- ifelse(den > 0, num / den, NA_real_)
  array(slopes, dim = dim_arr[-1])
}
