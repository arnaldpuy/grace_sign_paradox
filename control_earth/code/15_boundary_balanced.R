# 15_boundary_balanced.R -- data-space boundary of the sign-balanced planet
# (AMENDMENT_01_sign_balance.md section 4; D051).
#
# The claim that synthetic planets reproduce the real-Earth detectability
# boundary (~0.35 cm/yr) was established on three planets before sign balance
# was imposed. This script tests it on graded_balanced_v1 with the UNCHANGED
# empirical protocol (D020): per-aquifer dominance ratio R = IQR/|median| over
# the stacked three-centre 1,920-recipe ensemble, loess R = 1 crossing,
# aquifer bootstrap. As a protocol check it first re-computes graded_floor_v3,
# which must reproduce the banked 0.299 [0.254, 0.342].
#
#   Rscript code/15_boundary_balanced.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("data.table"))

ce_root <- here::here("control_earth")
ps_repo <- here::here()
setwd(ce_root)

# The paper's own metric + boundary functions (D020 invariant) -----------------

source(file.path(ps_repo, "functions/recipe_subspace_metrics_fun.R"))
source(file.path(ps_repo, "functions/detectability_boundary_fun.R"))

# BOUNDARY PER PLANET ##########################################################

boundary_planet_fun <- function(planet) {
  files <- list.files(file.path("datasets/output",
                                paste0("tier_a_synth_", planet)),
                      pattern = "jasechko_per_recipe_trends",
                      recursive = TRUE, full.names = TRUE)
  stopifnot(length(files) == 3)
  prt <- rbindlist(lapply(files, fread))
  m <- recipe_subspace_metrics_fun(prt)
  b <- detectability_boundary_fun(m$abs_median, m$dominance_R, type = "R1")
  data.table(planet = planet, n_aquifers = b$n, n_recipes = prt[, .N] / b$n,
             est = b$est, lo = b$lo, hi = b$hi)
}

res <- rbindlist(lapply(c("graded_floor_v3", "graded_balanced_v1"),
                        boundary_planet_fun))

# Protocol check: v3 must reproduce its banked boundary ------------------------

stopifnot(abs(res[planet == "graded_floor_v3", est] - 0.299) < 0.005)

cat("=== Data-space boundary (R = 1), unchanged empirical protocol ===\n")
cat("real Earth: 0.347 [0.286, 0.423] (MS16, fixed reference)\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-20s %.3f [%.3f, %.3f] cm/yr\n",
              res$planet[i], res$est[i], res$lo[i], res$hi[i]))
}

fwrite(res, "datasets/output/boundary_balanced.csv")
cat("wrote datasets/output/boundary_balanced.csv\n")
