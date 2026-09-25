# 19_sign_interleaving.R -- sign-interleaving diagnostic of the planted layouts
# (D055): for every planet, the number of aquifers whose nearest aquifer of the
# OPPOSITE planted sign lies within 500 km (the smoothing scale) and the median
# of those nearest distances. Spec-only: it reads the planted trends of each
# truth specification and the cohort cell masks, no recovered cube.
#
# Distance: haversine (R = 6371 km) between aquifer centroids, each centroid
# the mean latitude and longitude of the aquifer's 0.5-degree cells. The median
# is taken over the aquifers that have an opposite-sign neighbour at all (NA on
# a planet with a single sign).
#
#   Rscript control_earth/code/19_sign_interleaving.R
#
# Output: control_earth/datasets/output/sign_interleaving.csv

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "jsonlite"))

# PATHS ########################################################################

ce_root <- here::here("control_earth")
masks_csv <- here::here("datasets", "jasechko_2024", "cell_masks_0p5deg.csv")
out_csv <- file.path(ce_root, "datasets", "output", "sign_interleaving.csv")
planets <- c("graded_floor_v3", "graded_balanced_v1", "graded_balanced_v2",
             "graded_balanced_v3", "graded_balanced_v1_reversed",
             "graded_balanced_v3_noise2")
radius_km <- 6371
scale_km <- 500

# CENTROIDS ####################################################################

centroids <- fread(masks_csv)[, .(lon = mean(lon), lat = mean(lat)),
                              by = aquifer_idx]
setkey(centroids, aquifer_idx)
stopifnot(nrow(centroids) == 73L)

haversine_fun <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * radius_km * asin(sqrt(a))
}

# DIAGNOSTIC PER PLANET ########################################################

interleaving_fun <- function(planet) {
  spec <- fromJSON(file.path(ce_root, "datasets", "truth",
                             sprintf("truth_spec_%s.json", planet)))
  tr <- spec$aquifer_trends_cm_yr
  tr <- tr[names(tr) != "__default__"]
  dt <- data.table(aquifer_idx = as.integer(names(tr)), tau = unlist(tr))
  dt <- centroids[dt, on = .(aquifer_idx)]
  stopifnot(nrow(dt) == 73L, !anyNA(dt$lon))
  n <- nrow(dt)
  d <- outer(seq_len(n), seq_len(n), function(i, j)
    haversine_fun(dt$lon[i], dt$lat[i], dt$lon[j], dt$lat[j]))
  opposite <- outer(sign(dt$tau), sign(dt$tau), function(a, b) a * b < 0)
  d[!opposite] <- Inf
  nearest <- apply(d, 1, min)
  has_opposite <- is.finite(nearest)
  data.table(planet = planet,
             n_within_500km = sum(has_opposite & nearest < scale_km),
             median_km = if (any(has_opposite)) median(nearest[has_opposite])
                         else NA_real_)
}

result <- rbindlist(lapply(planets, interleaving_fun))
fwrite(result, out_csv)
print(result)
