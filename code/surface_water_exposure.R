## ---------------------------------------------------------------------------
## Surface-water exposure screening.
##
## The groundwater residual GWS = TWS - (soil moisture + snow + canopy) does
## not remove surface water. This script bounds, per Tier-1 aquifer, the
## largest OLS trend surface water could have imprinted over the 2002-2025
## record: the summed HydroLAKES volume of every lake and reservoir whose
## centroid falls inside the aquifer footprint, converted to aquifer-mean
## equivalent water height and multiplied by the extremal OLS factor of the
## real GRACE month sampling. For a storage history bounded in [0, C], the
## largest attainable OLS slope is C times the sum of the positive OLS time
## weights (attained by a step at the record's temporal mean, ~1.5 C/T for
## uniform sampling, not C/T). This is a stock-based EXPOSURE bound on the
## registered water bodies, not an estimate of the surface-water trend.
##
## Inputs (external downloads, not shipped in this repository):
##   datasets/jasechko_2024/aquifers_shp/jasechko_et_al_2024_aquifers.shp
##   datasets/jasechko_2024/cohort_50K_n73.csv
##   datasets/hydrolakes/HydroLAKES_points_v10_shp/HydroLAKES_points_v10.shp
##     (HydroLAKES v1.0, Messager et al. 2016, doi:10.1038/ncomms13603;
##      Vol_total in million m^3; Lake_type 1 = natural, 2 = reservoir,
##      3 = regulated lake; reservoir volumes incorporate GRanD)
##
## Output (banked in this repository):
##   datasets/output/tier_a/surface_water_exposure.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
})
sf::sf_use_s2(FALSE)

# Extremal OLS factor from the real GRACE month midpoints (250 months common
# to the three centres; the sampling the recipe trends are fitted on).

mons <- data.table::fread(
  "datasets/output/tier_a/grace_month_midpoints.csv")
t_i <- mons$mid_decimal_year
w_i <- (t_i - mean(t_i)) / sum((t_i - mean(t_i))^2)
k_ols <- sum(w_i[w_i > 0])
stopifnot(nrow(mons) == 250L, abs(k_ols - 0.0618) < 0.001)

# 1. Cohort polygons ----------------------------------------------------------

coh <- data.table::fread("datasets/jasechko_2024/cohort_50K_n73.csv")
stopifnot(nrow(coh) == 73L)

shp <- sf::st_read("datasets/jasechko_2024/aquifers_shp/jasechko_et_al_2024_aquifers.shp",
                   quiet = TRUE)
shp$Study_area <- ifelse(is.na(shp$Broader) | shp$Broader %in% c("-", ""),
                         shp$Aquifer,
                         paste0(shp$Aquifer, " (", shp$Broader, ")"))
shp <- shp[shp$Study_area %in% coh$Study_area, ]
shp <- sf::st_make_valid(shp)
stopifnot(length(unique(shp$Study_area)) == 73L)

# 2. HydroLAKES points-in-polygon ---------------------------------------------

hl <- sf::st_read("datasets/hydrolakes/HydroLAKES_points_v10_shp/HydroLAKES_points_v10.shp",
                  quiet = TRUE)
hl <- hl[, c("Hylak_id", "Lake_type", "Lake_area", "Vol_total")]
hl <- sf::st_transform(hl, sf::st_crs(shp))
jn <- sf::st_join(hl, shp["Study_area"], left = FALSE)
lakes <- data.table::as.data.table(sf::st_drop_geometry(jn))

# 3. Per-aquifer exposure bound -----------------------------------------------

exp_dt <- lakes[, .(n_waterbodies = .N,
                    n_reservoirs = sum(Lake_type %in% c(2L, 3L)),
                    vol_total_km3 = sum(Vol_total) / 1e3,
                    vol_reservoir_km3 = sum(Vol_total[Lake_type %in% c(2L, 3L)]) / 1e3),
                by = Study_area]
exp_dt <- merge(coh[, .(Study_area, area_km2)], exp_dt,
                by = "Study_area", all.x = TRUE)
for (j in c("n_waterbodies", "n_reservoirs", "vol_total_km3", "vol_reservoir_km3"))
  data.table::set(exp_dt, which(is.na(exp_dt[[j]])), j, 0)
exp_dt[, stock_ewh_cm := vol_total_km3 / area_km2 * 1e5]
exp_dt[, bound_cmyr := stock_ewh_cm * k_ols]
data.table::setorder(exp_dt, -bound_cmyr)

# 4. Bank ---------------------------------------------------------------------

data.table::fwrite(exp_dt, "datasets/output/tier_a/surface_water_exposure.csv")

cat(sprintf("73 aquifers; %d with >= 1 water body; total volume %.0f km^3\n",
            exp_dt[n_waterbodies > 0, .N], sum(exp_dt$vol_total_km3)))
cat(sprintf("bound >= 0.35 cm/yr: %d aquifers; >= 0.10: %d; median bound %.3f\n",
            exp_dt[bound_cmyr >= 0.35, .N], exp_dt[bound_cmyr >= 0.10, .N],
            stats::median(exp_dt$bound_cmyr)))
print(exp_dt[1:10, .(Study_area = substr(Study_area, 1, 45), n_waterbodies,
                     vol_total_km3 = round(vol_total_km3, 1),
                     bound_cmyr = round(bound_cmyr, 3))])
