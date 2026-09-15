load_grace_recipes_fun <- function(variant_dir, var_name = "ewh_anom",
                                     suffix = "") {

  # LOAD GRACE RECIPE CUBES #####################################################

  # Reads the per-centre NetCDFs into arr (time, recipe, lat, lon) + times/lat/lon +
  # recipe_keys (one row per recipe, seven axis columns; order matches arr).
  # suffix = "" for the 89-month POC cubes, "_appendix" for the 272-month cubes.
  # Recipe index mirrors Python row-major flatten (geocenter varies fastest).

  centres <- c("CSR", "JPL", "GFZ")
  paths <- file.path(variant_dir, paste0("grace_ewh_variants_",
                                         tolower(centres), suffix, ".nc"))
  stopifnot(all(file.exists(paths)))

  arrays <- vector("list", length(centres))
  meta <- list()
  for (k in seq_along(centres)) {
    nc <- ncdf4::nc_open(paths[k])
    a <- ncdf4::ncvar_get(nc, var_name)

    # ncdf4 returns the NetCDF axis order reversed (lon first, time last).
    if (k == 1) {
      # Parse the CF time epoch from the units attribute. xarray writes
      # something like "days since 2018-06-15 00:00:00"; pull the date.
      time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
      epoch_match <- regmatches(
        time_units,
        regexpr("\\d{4}-\\d{2}-\\d{2}", time_units))
      time_origin <- if (length(epoch_match)) as.Date(epoch_match)
                       else as.Date("1970-01-01")
      meta$time <- as.Date(ncdf4::ncvar_get(nc, "time"),
                           origin = time_origin)
      meta$lat <- as.numeric(ncdf4::ncvar_get(nc, "lat"))
      meta$lon <- as.numeric(ncdf4::ncvar_get(nc, "lon"))
      meta$truncation <- as.integer(ncdf4::ncvar_get(nc, "truncation"))
      meta$filter <- as.character(ncdf4::ncvar_get(nc, "filter"))
      meta$gia_model <- as.character(ncdf4::ncvar_get(nc, "gia_model"))
      meta$c20_treatment <- as.character(
        ncdf4::ncvar_get(nc, "c20_treatment"))
      meta$c30_treatment <- as.character(
        ncdf4::ncvar_get(nc, "c30_treatment"))
      meta$geocenter <- as.character(ncdf4::ncvar_get(nc, "geocenter"))
    }
    arrays[[k]] <- a
    ncdf4::nc_close(nc)
  }

  # ncdf4 returns dims reversed; after abind over centre, permute so geocenter is the
  # fastest non-time axis, matching Python's product() order exactly.
  stacked <- abind::abind(arrays, along = length(dim(arrays[[1]])) + 1L)
  stacked <- aperm(
    stacked,
    perm = c(9, 3, 4, 5, 6, 7, 8, 10, 2, 1),
  )

  # Canonical month-15 time axis.
  ym <- format(meta$time, "%Y-%m")
  times <- as.Date(paste0(ym, "-15"))

  n_t <- length(times)
  n_c <- length(centres)
  n_tr <- length(meta$truncation)
  n_f <- length(meta$filter)
  n_g <- length(meta$gia_model)
  n_c20 <- length(meta$c20_treatment)
  n_c30 <- length(meta$c30_treatment)
  n_gc <- length(meta$geocenter)
  n_lat <- length(meta$lat)
  n_lon <- length(meta$lon)
  # After the c(9,3,4,5,6,7,8,10,2,1) aperm the dim order is
  # (time, geo, c30, c20, gia, filter, trunc, centre, lat, lon).
  stopifnot(dim(stacked)[1] == n_t,
            dim(stacked)[2] == n_gc,
            dim(stacked)[3] == n_c30,
            dim(stacked)[4] == n_c20,
            dim(stacked)[5] == n_g,
            dim(stacked)[6] == n_f,
            dim(stacked)[7] == n_tr,
            dim(stacked)[8] == n_c,
            dim(stacked)[9] == n_lat,
            dim(stacked)[10] == n_lon)

  # expand.grid varies the FIRST arg fastest, Python the LAST: pass axes reversed,
  # then reverse the column order.
  recipe_keys <- expand.grid(
    geocenter = meta$geocenter,
    c30_treatment = meta$c30_treatment,
    c20_treatment = meta$c20_treatment,
    gia_model = meta$gia_model,
    filter = meta$filter,
    truncation = meta$truncation,
    centre = centres,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )
  data.table::setDT(recipe_keys)
  data.table::setcolorder(recipe_keys, c("centre", "truncation", "filter",
                                          "gia_model", "c20_treatment",
                                          "c30_treatment", "geocenter"))

  # Collapse the seven recipe axes into one: (time, recipe, lat, lon).
  flat <- array(
    stacked,
    dim = c(n_t,
            n_c * n_tr * n_f * n_g * n_c20 * n_c30 * n_gc,
            n_lat, n_lon),
  )
  stopifnot(dim(flat)[2] == nrow(recipe_keys))

  list(
    arr = flat,
    times = times,
    lat = meta$lat,
    lon = meta$lon,
    recipe_keys = recipe_keys,
    centres = centres
  )
}
