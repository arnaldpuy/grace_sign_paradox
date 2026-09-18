# 10_recipe_skill.R -- per-recipe sign skill on known truth (D036; leverage N1+N2).
#
# The strongest remaining objection to the control-Earth result is that the
# ensemble statistics dilute skilled recipes: "a good analyst with the right
# recipe escapes the floor". Because the graded planet's truth is known, this
# is directly testable: for EACH of the 1,920 recipes, compute its individual
# sign-accuracy as a function of the planted magnitude, and ask whether any
# recipe (or the 180 fully-corrected subset, or the ensemble median as an
# estimator) reliably recovers signs below the floor. Selection effects are
# guarded by split-half replication: recipes selected as "best" on a random
# half of the aquifers are re-scored on the held-out half.
#
#   Rscript code/10_recipe_skill.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "ggplot2", "cowplot",
                          "scales"))

# Paths (here() resolves to $HOME in this repo -- use absolute paths) ----------

ce_root <- here::here("control_earth")
source(here::here("functions/theme_AP.R"))
out_dir <- file.path(ce_root, "datasets/output/recipe_skill")
dir.create(out_dir, showWarnings = FALSE)

# DATA #########################################################################

trends <- rbindlist(lapply(c("csr", "jpl", "gfz"), function(cc) {
  fread(file.path(ce_root, "datasets/output/tier_a_synth_graded_floor_v3", cc,
                  "jasechko_per_recipe_trends.csv"))
}))
sc <- fread(file.path(ce_root, "datasets/output/score_graded_floor_v3.csv"))
trends <- merge(trends, sc[, .(aquifer_idx, tau_true)], by = "aquifer_idx")
trends[, recipe_id := paste(centre, recipe_idx_per_centre, sep = "_")]
trends[, correct := sign(trend_cm_yr) == sign(tau_true)]
trends[, abs_tau := abs(tau_true)]
trends[, fully_corrected := c20_treatment == "replaced" &
         c30_treatment == "replaced" & geocenter != "none" &
         gia_model != "none"]
stopifnot(uniqueN(trends$recipe_id) == 1920L,
          trends[fully_corrected == TRUE, uniqueN(recipe_id)] == 180L)

floor_est <- 0.75

# N1: PER-RECIPE SKILL #########################################################

# Sub-floor / supra-floor accuracy per recipe ---------------------------------

skill <- trends[, .(
  acc_all = mean(correct),
  acc_sub = mean(correct[abs_tau < floor_est]),
  acc_supra = mean(correct[abs_tau >= floor_est]),
  n_sub = sum(abs_tau < floor_est),
  fully_corrected = fully_corrected[1],
  filter = filter[1], gia_model = gia_model[1], centre = centre[1],
  truncation = truncation[1], geocenter = geocenter[1]), by = recipe_id]

# Per-recipe logistic floor (s90 crossing, mirroring 09_parametric_floor) -----

fit_one_fun <- function(d) {
  fit <- suppressWarnings(glm(correct ~ log(abs_tau), family = binomial(),
                              data = d))
  cf <- coef(fit)
  if (!is.finite(cf[2]) || cf[2] <= 0) return(NA_real_)
  exp(unname((qlogis(0.9) - cf[1]) / cf[2]))
}
skill_floor <- trends[, .(s90 = fit_one_fun(.SD)), by = recipe_id]
skill <- merge(skill, skill_floor, by = "recipe_id")

# Split-half replication guard (seed 71) ---------------------------------------

set.seed(71)
aq <- unique(trends$aquifer_idx)
half_a <- sample(aq, floor(length(aq) / 2))
rep_check <- merge(
  trends[abs_tau < floor_est & aquifer_idx %in% half_a,
         .(acc_select = mean(correct)), by = recipe_id],
  trends[abs_tau < floor_est & !(aquifer_idx %in% half_a),
         .(acc_holdout = mean(correct)), by = recipe_id], by = "recipe_id")
setorder(rep_check, -acc_select)
top20 <- rep_check[1:20]

# N2: THE ENSEMBLE MEDIAN AS AN ESTIMATOR ######################################

med_est <- sc[, .(aquifer_idx, abs_tau = abs(tau_true),
                  correct = sign(median_rec) == sign(tau_true))]
med_sub <- med_est[abs_tau < floor_est, mean(correct)]
med_supra <- med_est[abs_tau >= floor_est, mean(correct)]
med_s90 <- fit_one_fun(med_est)
med_pct <- mean(skill$acc_sub <= med_sub)

# OUTPUT #######################################################################

fwrite(skill, file.path(out_dir, "per_recipe_skill_graded_floor_v3.csv"))
fwrite(rep_check, file.path(out_dir, "split_half_replication.csv"))

cat("=== N1: per-recipe sign skill on known truth (graded planet, 1,920 recipes) ===\n")
cat(sprintf("sub-floor aquifers (|tau_true| < %.2f): n = %d\n",
            floor_est, skill$n_sub[1]))
cat(sprintf("sub-floor accuracy: median %.3f | best %.3f | worst %.3f\n",
            median(skill$acc_sub), max(skill$acc_sub), min(skill$acc_sub)))
cat(sprintf("recipes with sub-floor accuracy >= 0.90: %d of 1,920\n",
            skill[acc_sub >= 0.90, .N]))
cat(sprintf("recipes with sub-floor accuracy >= 0.95: %d of 1,920\n",
            skill[acc_sub >= 0.95, .N]))
cat(sprintf("fully-corrected 180: median sub-floor accuracy %.3f | best %.3f\n",
            skill[fully_corrected == TRUE, median(acc_sub)],
            skill[fully_corrected == TRUE, max(acc_sub)]))
cat(sprintf("per-recipe s90 floor: median %.2f | min %.2f | max %.2f cm/yr (%d of 1,920 defined)\n",
            median(skill$s90, na.rm = TRUE), min(skill$s90, na.rm = TRUE),
            max(skill$s90, na.rm = TRUE), sum(is.finite(skill$s90))))
best <- skill[which.max(acc_sub)]
cat(sprintf("best recipe: %s trunc %s filter %s GIA %s geoc %s | sub-floor %.3f, s90 %.2f\n",
            best$centre, best$truncation, best$filter, best$gia_model,
            best$geocenter, best$acc_sub, best$s90))

cat("\n=== split-half replication (top-20 recipes selected on half A) ===\n")
cat(sprintf("selection-half accuracy of top-20: mean %.3f\n", mean(top20$acc_select)))
cat(sprintf("held-out-half accuracy of top-20 : mean %.3f\n", mean(top20$acc_holdout)))
cat(sprintf("ensemble-wide held-out mean      : %.3f\n", mean(rep_check$acc_holdout)))

cat("\n=== N2: the ensemble median as an estimator ===\n")
cat(sprintf("median-estimator sub-floor accuracy: %.3f (beats %.0f%% of single recipes)\n",
            med_sub, 100 * med_pct))
cat(sprintf("median-estimator supra-floor accuracy: %.3f | s90 = %.2f cm/yr\n",
            med_supra, med_s90))

# FIGURE ######################################################################

# Per-recipe skill on known truth. Thin lines: per-recipe logistic fits of
# sign-accuracy against true magnitude (grey: all 1,920; orange: the 180
# fully-corrected recipes). Bold black: the ensemble median as an estimator.
# Dotted line: 0.90. Vertical dashed line: the ensemble recovery floor (0.75
# cm/yr). No recipe class crosses 0.90 materially below the floor.

grid_tau <- data.table(ltau = seq(log(0.04), log(3.5), length.out = 120))
grid_tau[, tau := exp(ltau)]
coefs <- trends[, {
  fit <- suppressWarnings(glm(correct ~ log(abs_tau), family = binomial()))
  .(a = coef(fit)[1], b = coef(fit)[2],
    fully_corrected = fully_corrected[1])
}, by = recipe_id]
curves <- coefs[, .(tau = grid_tau$tau,
                    p = plogis(a + b * grid_tau$ltau),
                    fully_corrected = fully_corrected), by = recipe_id]
med_fit <- suppressWarnings(glm(correct ~ log(abs_tau), family = binomial(),
                                data = med_est))
med_curve <- data.table(tau = grid_tau$tau,
                        p = plogis(coef(med_fit)[1] +
                                     coef(med_fit)[2] * grid_tau$ltau))

# Export the fitted curves so 06_boundary_overlay.R can draw Fig 4c cheaply ----

fwrite(coefs, file.path(out_dir, "per_recipe_logistic_coefs.csv"))
fwrite(data.table(a = unname(coef(med_fit)[1]), b = unname(coef(med_fit)[2])),
       file.path(out_dir, "median_estimator_coefs.csv"))

plot_skill <- ggplot() +
  geom_line(data = curves[fully_corrected == FALSE],
            aes(x = tau, y = p, group = recipe_id),
            color = "grey75", linewidth = 0.08, alpha = 0.25) +
  geom_line(data = curves[fully_corrected == TRUE],
            aes(x = tau, y = p, group = recipe_id),
            color = "#D55E00", linewidth = 0.12, alpha = 0.4) +
  geom_line(data = med_curve, aes(x = tau, y = p),
            color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0.9, linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = floor_est, linetype = "dashed", linewidth = 0.35,
             color = "#0072B2") +
  scale_x_log10(breaks = c(0.05, 0.1, 0.35, 0.75, 3),
                labels = c("0.05", "0.1", "0.35", "0.75", "3")) +
  labs(x = bquote("True trend magnitude (cm" ~ yr^-1 * ")"),
       y = "P(correct sign)") +
  theme_AP()

ggsave(file.path(ce_root, "figures/fig_recipe_skill.pdf"), plot_skill,
       width = 3.4, height = 2.4)
cat("\nwrote figures/fig_recipe_skill.pdf\n")
