# 09_parametric_floor.R -- parametric fit of the sign-recovery floor (D035).
#
# D032 found the loess 0.95-crossing ill-conditioned on graded_floor_v3: the
# recovery curve plateaus at ~0.95 across 0.6-1.5 cm/yr, so the crossing
# point estimate (2.38) fell outside its own bootstrap CI [0.52, 2.26]. This
# script replaces the single crossing with a standard parametric description:
# a logistic regression of the per-aquifer correct-sign fraction on
# log|tau_true| (quasibinomial GLM, one observation per aquifer),
#
#   logit(P_sign) = a + b log|tau|,   s*_q = exp((logit(q) - a) / b).
#
# A three-parameter logistic with a free lower asymptote was tried first and
# is UNIDENTIFIABLE on these data (no left tail below |tau| = 0.05: the
# asymptote ran to -24 and the bootstrap CI spanned [0.7, 67]) -- recorded
# here so it is not re-tried. Because the curve plateaus near 0.95, the 0.95
# crossing is intrinsically borderline; the 0.90 crossing is reported
# alongside as the well-determined summary. CIs from a 2,000-draw aquifer
# bootstrap (seed 71, project convention).
#
#   Rscript code/09_parametric_floor.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("here", "data.table", "ggplot2", "cowplot",
                          "scales"))

# Paths (here() resolves to $HOME in this repo -- use absolute paths) ----------

ce_root <- here::here("control_earth")
source(here::here("functions/theme_AP.R"))

sc <- fread(file.path(ce_root, "datasets/output/score_graded_floor_v3.csv"))
sc[, ltau := log(abs(tau_true))]

# FIT #########################################################################

fit_floor_fun <- function(d) {
  fit <- suppressWarnings(
    glm(p_sign ~ ltau, family = quasibinomial(), data = d))
  cf <- coef(fit)
  if (!is.finite(cf[2]) || cf[2] <= 0) return(NULL)
  cross <- function(q) exp(unname((qlogis(q) - cf[1]) / cf[2]))
  c(a = unname(cf[1]), b = unname(cf[2]),
    s50 = cross(0.50), s90 = cross(0.90), s95 = cross(0.95))
}

full_fit <- fit_floor_fun(sc)
cat("=== Parametric floor fit (graded_floor_v3, 73 aquifers) ===\n")
cat("logit(P_sign) = a + b log|tau|, quasibinomial GLM\n")
print(round(full_fit, 4))

# Aquifer bootstrap ------------------------------------------------------------

set.seed(71)
n_boot <- 2000L
boot <- rbindlist(lapply(seq_len(n_boot), function(k) {
  f <- fit_floor_fun(sc[sample(.N, .N, replace = TRUE)])
  if (is.null(f)) return(NULL)
  as.data.table(as.list(f))
}))
cat(sprintf("\nbootstrap: %d of %d draws with a defined upward crossing\n",
            nrow(boot), n_boot))
ci <- boot[, .(s50_lo = quantile(s50, 0.025), s50_hi = quantile(s50, 0.975),
               s90_lo = quantile(s90, 0.025), s90_hi = quantile(s90, 0.975),
               s95_lo = quantile(s95, 0.025), s95_hi = quantile(s95, 0.975))]
cat(sprintf("s*_0.90 = %.2f [%.2f, %.2f] cm/yr\n",
            full_fit[["s90"]], ci$s90_lo, ci$s90_hi))
cat(sprintf("s*_0.95 = %.2f [%.2f, %.2f] cm/yr\n",
            full_fit[["s95"]], ci$s95_lo, ci$s95_hi))

res <- data.table(planet = "graded_floor_v3", t(as.list(full_fit)),
                  s50_lo = ci$s50_lo, s50_hi = ci$s50_hi,
                  s90_lo = ci$s90_lo, s90_hi = ci$s90_hi,
                  s95_lo = ci$s95_lo, s95_hi = ci$s95_hi,
                  n_boot_ok = nrow(boot))
fwrite(res, file.path(ce_root, "datasets/output/parametric_floor_v3.csv"))

# FIGURE ######################################################################

# The parametric alternative to the loess crossing. Points: the 73 aquifers
# of graded_floor_v3. Line: the fitted logistic; ribbon: the 2.5-97.5%
# envelope of the bootstrap fits. Orange band: the fitted 0.90 crossing with
# its bootstrap CI (the well-determined floor summary); the 0.95 crossing is
# reported in text but not drawn -- the plateau makes it borderline by
# construction. Dotted lines: P = 0.90 and 0.95.

grid_tau <- data.table(ltau = seq(log(0.04), log(3.5), length.out = 200))
grid_tau[, tau := exp(ltau)]
grid_tau[, fit := plogis(full_fit[["a"]] + full_fit[["b"]] * ltau)]
env <- rbindlist(lapply(seq_len(nrow(boot)), function(k) {
  data.table(ltau = grid_tau$ltau,
             p = plogis(boot$a[k] + boot$b[k] * grid_tau$ltau))
}))
env_q <- env[, .(lo = quantile(p, 0.025), hi = quantile(p, 0.975)), by = ltau]
grid_tau <- merge(grid_tau, env_q, by = "ltau")

plot_floor <- ggplot(sc, aes(x = abs(tau_true), y = p_sign)) +
  annotate("rect", xmin = ci$s90_lo, xmax = ci$s90_hi,
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#D55E00") +
  geom_vline(xintercept = full_fit[["s90"]], linetype = "dashed",
             linewidth = 0.35, color = "#D55E00") +
  geom_hline(yintercept = 0.90, linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = 0.95, linetype = "dotted", linewidth = 0.3,
             color = "grey55") +
  geom_ribbon(data = grid_tau, aes(x = tau, ymin = lo, ymax = hi),
              inherit.aes = FALSE, alpha = 0.25, fill = "grey60") +
  geom_line(data = grid_tau, aes(x = tau, y = fit), inherit.aes = FALSE,
            linewidth = 0.5, color = "grey20") +
  geom_point(size = 1.1, alpha = 0.7) +
  annotate("label", x = full_fit[["s90"]], y = 0.15,
           label = "0.90 crossing", size = 2.2, color = "#D55E00",
           fill = "white", alpha = 0.9) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.35, 1, 3),
                labels = c("0.05", "0.1", "0.35", "1", "3")) +
  labs(x = bquote("True trend magnitude |" * tau * "| (cm" ~ yr^-1 * ")"),
       y = "P(correct sign) across 1,920 recipes") +
  theme_AP()

ggsave(file.path(ce_root, "figures/fig_parametric_floor.pdf"), plot_floor,
       width = 4.2, height = 2.8)
cat("\nwrote figures/fig_parametric_floor.pdf\n")
