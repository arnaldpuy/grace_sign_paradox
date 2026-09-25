load_wris_wells_fun <- function(wris_csv,
                                  aoi_lat_min = 24.0, aoi_lat_max = 32.0,
                                  aoi_lon_min = 73.0, aoi_lon_max = 88.0) {

  # LOAD AND RESHAPE THE CGWB FIGSHARE WELL CSV #################################

  # Input: path to CGWB_India_filtered_GWLs_ref_sy_2000_2022.csv
  # Returns list with:

  raw <- data.table::fread(wris_csv)
  setnames(raw, c("Station Code", "Latitude", "Longitude", "Reference_Sy"),
           c("station_code", "lat", "lon", "reference_sy"),
           skip_absent = TRUE)
  in_aoi <- raw[, lat >= aoi_lat_min & lat <= aoi_lat_max
                  & lon >= aoi_lon_min & lon <= aoi_lon_max]
  sub <- raw[in_aoi]

  meta_cols <- c("station_code", "Station Name", "Station Type",
                  "Agency Name", "State", "District", "Tehsil", "Block",
                  "Village", "lat", "lon", "Type of Well", "Aquifer Type",
                  "Well Depth", "reference_sy")
  meta_cols <- intersect(meta_cols, names(sub))
  stations <- sub[, ..meta_cols]
  data.table::setnames(stations, c("State", "District"),
                        c("state", "district"), skip_absent = TRUE)

  obs_cols <- grep("^(Jan|May|Aug|Nov)-\\d{2}$", names(sub), value = TRUE)
  long <- data.table::melt(
    sub,
    id.vars = c("station_code", "lat", "lon", "reference_sy"),
    measure.vars = obs_cols,
    variable.name = "obs_col",
    value.name = "depth_to_water_m",
    variable.factor = FALSE,
  )
  long <- long[!is.na(depth_to_water_m)]

  month_lookup <- c(Jan = 1L, May = 5L, Aug = 8L, Nov = 11L)
  parts <- data.table::tstrsplit(long$obs_col, "-", fixed = TRUE)
  long[, month := month_lookup[parts[[1L]]]]
  long[, year := 2000L + as.integer(parts[[2L]])]
  long[, date := as.Date(sprintf("%04d-%02d-15", year, month))]
  long[, obs_col := NULL]
  long[, c("year", "month") := NULL]

  list(stations = stations, records = long)
}
