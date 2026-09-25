fit_recipes_to_wells_fun <- function(recipes, wris_data,
                                       baseline_start = "2015-04-01",
                                       baseline_end = "2020-12-31",
                                       min_months = 6L,
                                       deseasonalize = TRUE) {

  # SCORE EACH (RECIPE × WELL) PAIR BY CORRELATION AND MAE ######################

  # Returns a list with:

  arr <- recipes$arr
  n_t <- dim(arr)[1L]
  n_rec <- dim(arr)[2L]
  lats <- recipes$lat
  lons <- recipes$lon
  times <- recipes$times

  if (isTRUE(deseasonalize)) {
    # Subtract per-calendar-month climatology over the full GRACE record
    # at each (recipe, cell) before joining with the well.
    months <- as.integer(format(times, "%m"))
    arr_des <- arr
    for (m in unique(months)) {
      idx <- which(months == m)
      if (length(idx) <= 1L) next
      arr_des[idx, , , ] <- sweep(
        arr_des[idx, , , , drop = FALSE],
        c(2L, 3L, 4L),
        apply(arr[idx, , , , drop = FALSE], c(2L, 3L, 4L), mean,
              na.rm = TRUE),
        FUN = "-",
      )
    }
    arr <- arr_des
  }

  stations <- wris_data$stations
  records <- wris_data$records
  n_wells <- nrow(stations)
  cors <- matrix(NA_real_, nrow = n_wells, ncol = n_rec)
  maes <- matrix(NA_real_, nrow = n_wells, ncol = n_rec)
  n_used <- integer(n_wells)

  for (k in seq_len(n_wells)) {
    sc <- stations$station_code[k]
    well <- station_storage_anomaly_fun(
      records, station_id = sc,
      baseline_start = baseline_start,
      baseline_end = baseline_end,
      deseasonalize = deseasonalize)
    if (nrow(well) < min_months) next

    i <- which.min(abs(lats - stations$lat[k]))
    j <- which.min(abs(lons - stations$lon[k]))
    cell <- arr[, , i, j]
    cell_dt <- data.table::data.table(date = times, cell)
    merged <- well[cell_dt, on = "date", nomatch = NULL]
    if (nrow(merged) < min_months) next

    n_used[k] <- nrow(merged)
    storage <- merged$storage_anom_cm
    recipe_cols <- setdiff(names(merged), c("date", "storage_anom_cm"))
    sd_storage <- sd(storage, na.rm = TRUE)
    if (!is.finite(sd_storage) || sd_storage < 1e-8) next

    for (r in seq_len(n_rec)) {
      rcol <- merged[[recipe_cols[r]]]
      sd_r <- sd(rcol, na.rm = TRUE)
      if (!is.finite(sd_r) || sd_r < 1e-8) next
      cors[k, r] <- stats::cor(rcol, storage, use = "pairwise.complete.obs")
      maes[k, r] <- mean(abs(rcol - storage), na.rm = TRUE)
    }
  }

  ok <- n_used >= min_months
  fits <- data.table::data.table(
    station_code = stations$station_code[ok],
    lat = stations$lat[ok],
    lon = stations$lon[ok],
    state = stations$state[ok],
    reference_sy = stations$reference_sy[ok],
    n_months_used = n_used[ok]
  )
  list(
    correlations = cors[ok, , drop = FALSE],
    mae_cm = maes[ok, , drop = FALSE],
    fits = fits,
    recipe_keys = recipes$recipe_keys
  )
}
