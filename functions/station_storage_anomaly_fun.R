station_storage_anomaly_fun <- function(records, station_id,
                                          baseline_start = "2015-04-01",
                                          baseline_end = "2020-12-31",
                                          deseasonalize = TRUE) {

  # WELL STORAGE ANOMALY (cm), OPTIONALLY DE-SEASONALISED ######################

  # Returns a 2-column data.table (date, storage_anom_cm) for the
  # requested station; empty if no observations exist or baseline empty.

  sub <- records[station_code == station_id]
  if (nrow(sub) == 0L) {
    return(data.table::data.table(date = as.Date(character()),
                                   storage_anom_cm = numeric()))
  }
  sy <- sub$reference_sy[1L]

  monthly <- sub[, .(depth_m = mean(depth_to_water_m, na.rm = TRUE)),
                  by = .(date)]
  data.table::setorder(monthly, date)

  baseline <- monthly[date >= as.Date(baseline_start)
                      & date <= as.Date(baseline_end)]
  if (nrow(baseline) == 0L) {
    return(data.table::data.table(date = as.Date(character()),
                                   storage_anom_cm = numeric()))
  }
  baseline_mean <- mean(baseline$depth_m, na.rm = TRUE)
  monthly[, storage_anom_cm := -(depth_m - baseline_mean) * sy * 100]

  if (isTRUE(deseasonalize)) {
    monthly[, month := data.table::month(date)]
    clim <- monthly[, .(clim_m = mean(storage_anom_cm, na.rm = TRUE)),
                     by = month]
    monthly <- clim[monthly, on = "month"]
    monthly[, storage_anom_cm := storage_anom_cm - clim_m]
    monthly[, c("month", "clim_m") := NULL]
  }
  monthly[, .(date, storage_anom_cm)]
}
