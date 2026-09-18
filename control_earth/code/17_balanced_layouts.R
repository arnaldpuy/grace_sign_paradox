# 17_balanced_layouts.R -- the amendment's fixed analyses, per balanced planet
# (AMENDMENT_01_sign_balance.md section 4; D051). Idempotent: analyses every
# balanced planet whose score CSV exists, plus graded_floor_v3 as the
# depletion-dominated negative control, and pools s90 across the balanced
# LAYOUTS (1-3; the reversed planet is assessed on its own, amendment 5).
#
#   Rscript code/17_balanced_layouts.R      # rerun after each planet lands

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("data.table"))

control_root <- here::here("control_earth")
ps_repo <- here::here()
setwd(control_root)
source(file.path(ps_repo, "functions/recipe_subspace_metrics_fun.R"))
source(file.path(ps_repo, "functions/detectability_boundary_fun.R"))

# Same estimator as 09_parametric_floor.R / 14_floor_sign_balance.R -----------

fit_cross <- function(d, q) {
  d <- d[is.finite(p_sign) & abs(tau_true) > 0]
  if (nrow(d) < 10L) return(NA_real_)
  cf <- coef(suppressWarnings(
    glm(p_sign ~ log(abs(tau_true)), family = quasibinomial(), data = d)))
  if (!is.finite(cf[2L]) || cf[2L] <= 0) return(NA_real_)
  exp(unname((qlogis(q) - cf[1L]) / cf[2L]))
}
boot_cross <- function(d, q, n = 2000L, seed = 71L) {
  set.seed(seed)
  bs <- replicate(n, fit_cross(d[sample.int(nrow(d), replace = TRUE)], q))
  bs <- bs[is.finite(bs)]
  if (length(bs) < 2L) return(c(NA_real_, NA_real_))
  stats::quantile(bs, c(0.025, 0.975))
}
s90_row <- function(d, label) {
  e <- fit_cross(d, 0.90); ci <- boot_cross(d, 0.90)
  data.table(subset = label, n = nrow(d), s90 = e, lo = ci[[1L]], hi = ci[[2L]])
}

# PER-PLANET ANALYSES ##########################################################

planets <- c("graded_floor_v3", "graded_balanced_v1", "graded_balanced_v2",
             "graded_balanced_v3", "graded_balanced_v1_reversed",
             "graded_balanced_v3_noise2")
planets <- planets[file.exists(
  file.path("datasets/output", paste0("score_", planets, ".csv")))]
cat("planets with banked scores:", paste(planets, collapse = ", "), "\n\n")

bin_edges <- c(0, 0.25, 0.5, 0.75, 1, Inf)
bin_labels <- c("<=0.25", "0.25-0.5", "0.5-0.75", "0.75-1", ">1")

summary_list <- list(); bins_list <- list(); floors_list <- list()
for (p in planets) {
  sc <- fread(file.path("datasets/output", paste0("score_", p, ".csv")))

  # Fitted crossings (s95 recorded but NOT estimable within the planted range)
  fl <- rbind(s90_row(sc, p),
              if (sc[tau_true < 0, .N] >= 10) s90_row(sc[tau_true < 0], paste0(p, " depleting")),
              if (sc[tau_true > 0, .N] >= 10) s90_row(sc[tau_true > 0], paste0(p, " recovering")))
  floors_list[[p]] <- fl

  # Model-free binned table (primary presentation)
  b <- sc[, .(median_p_sign = median(p_sign), n = .N),
          by = .(bin = cut(abs(tau_true), bin_edges, labels = bin_labels))]
  setorder(b, bin)
  bins_list[[p]] <- data.table(planet = p, b)

  # Class-neutral metrics + attenuation + data-space boundary
  acc_neg <- sc[tau_true < 0, mean(p_sign)]
  acc_pos <- sc[tau_true > 0, mean(p_sign)]
  files <- list.files(file.path("datasets/output", paste0("tier_a_synth_", p)),
                      pattern = "jasechko_per_recipe_trends",
                      recursive = TRUE, full.names = TRUE)
  bnd <- if (length(files) == 3) {
    m <- recipe_subspace_metrics_fun(rbindlist(lapply(files, fread)))
    detectability_boundary_fun(m$abs_median, m$dominance_R, type = "R1")
  } else list(est = NA_real_, lo = NA_real_, hi = NA_real_)
  summary_list[[p]] <- data.table(
    planet = p, n_neg = sc[tau_true < 0, .N], n_pos = sc[tau_true > 0, .N],
    median_p_sign = median(sc$p_sign),
    n_psign_95 = sc[p_sign >= 0.95, .N],
    acc = sc[, mean(p_sign)], acc_neg = acc_neg, acc_pos = acc_pos,
    balanced_acc = (acc_neg + acc_pos) / 2,
    attenuation = sc[abs(tau_true) > 0.5, median(median_rec / tau_true)],
    s90 = fl$s90[1], s90_lo = fl$lo[1], s90_hi = fl$hi[1],
    boundary = bnd$est, boundary_lo = bnd$lo, boundary_hi = bnd$hi)
}

# POOLED s90 across the balanced layouts (amendment 4; reversed excluded) ------

layout_ids <- intersect(planets, c("graded_balanced_v1", "graded_balanced_v2",
                                   "graded_balanced_v3"))
pooled <- NULL
if (length(layout_ids) >= 2) {
  stacked <- rbindlist(lapply(layout_ids, function(p)
    fread(file.path("datasets/output", paste0("score_", p, ".csv")))[
      , .(tau_true, p_sign)]))
  pooled <- s90_row(stacked, sprintf("POOLED (%d layouts)", length(layout_ids)))
  floors_list[["pooled"]] <- pooled
}

# REPORT + BANK ################################################################

res_summary <- rbindlist(summary_list)
res_bins <- dcast(rbindlist(bins_list), bin ~ planet, value.var = "median_p_sign")
res_floors <- rbindlist(floors_list)

cat("=== Binned median P_sign by |tau_true| (primary presentation) ===\n")
print(res_bins, digits = 3)
cat("\n=== Fitted 90% crossings (secondary summary) ===\n")
res_floors[, sprintf("%-42s n=%3d  s90 = %.2f [%.2f, %.2f]",
                     subset, n, s90, lo, hi)] |> paste(collapse = "\n") |> cat()
cat("\n\n=== Per-planet summary ===\n")
print(res_summary[, .(planet, median_p_sign, n_psign_95, balanced_acc,
                      attenuation, boundary)], digits = 3)

fwrite(res_summary, "datasets/output/balanced_layouts_summary.csv")
fwrite(rbindlist(bins_list), "datasets/output/balanced_layouts_bins.csv")
fwrite(res_floors, "datasets/output/balanced_layouts_floors.csv")
cat("\nwrote balanced_layouts_{summary,bins,floors}.csv\n")
