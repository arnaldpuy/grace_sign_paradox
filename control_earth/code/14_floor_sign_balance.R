# 14_floor_sign_balance.R -- does the recovery floor survive sign balance? (D048)
#
# The 0.75 cm/yr floor of D035 is fitted on graded_floor_v3, where 69 of 73
# aquifers carry a negative total trend: an always-depleting rule scores 95%
# there, so coherent same-sign leakage could have flattered recovery at low
# magnitudes. graded_balanced_v1 carries the SAME per-aquifer |tau| with the
# signs re-randomised (43 negative / 30 positive in total-field terms), so a
# difference in the fitted floor is attributable to sign balance alone.
#
# Also refits separately on the depleting and the recovering aquifers, which
# only a balanced planet can support and which an external review requested.
#
#   Rscript code/14_floor_sign_balance.R

# Load packages ---------------------------------------------------------------

if (!requireNamespace("sensobol", quietly = TRUE)) install.packages("sensobol")
sensobol::load_packages(c("data.table"))

control_root <- here::here("control_earth")
setwd(control_root)

# FIT ##########################################################################

# Same estimator as 09_parametric_floor.R: quasibinomial logistic in log|tau|.
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
report <- function(d, lab, q = 0.90) {
  e <- fit_cross(d, q); ci <- boot_cross(d, q)
  cat(sprintf("%-34s n=%2d  s%.0f = %.2f [%.2f, %.2f] cm/yr\n",
              lab, nrow(d), 100 * q, e, ci[1L], ci[2L]))
  data.table(subset = lab, n = nrow(d), q = q, est = e,
             lo = ci[[1L]], hi = ci[[2L]])
}

bal <- fread("datasets/output/score_graded_balanced_v1.csv")
v3  <- fread("datasets/output/score_graded_floor_v3.csv")

cat("=== Recovery floor under sign balance (D048) ===\n")
cat(sprintf("balanced planet: %d negative / %d positive total trends; an\n",
            bal[tau_true < 0, .N], bal[tau_true > 0, .N]))
cat(sprintf("always-depleting rule scores %.0f%% here against %.0f%% on graded_floor_v3\n\n",
            100 * bal[, mean(tau_true < 0)], 100 * v3[, mean(tau_true < 0)]))

res <- rbindlist(list(
  report(v3,  "graded_floor_v3 (69 neg / 4 pos)"),
  report(bal, "graded_balanced_v1 (all)"),
  report(bal[tau_true < 0], "  balanced, depleting only"),
  report(bal[tau_true > 0], "  balanced, recovering only")))

# 95% crossings, for the same reason the triplet is reported in the SM ---------

cat("\n")
res <- rbindlist(list(res,
  report(v3,  "graded_floor_v3", q = 0.95),
  report(bal, "graded_balanced_v1", q = 0.95)))

fwrite(res, "datasets/output/floor_sign_balance.csv")

# VERDICT ######################################################################

# CLASS-NEUTRAL PERFORMANCE (AMENDMENT_01 section 4) ##########################

# Ordinary accuracy (mean per-aquifer P_sign) is misleading under sign
# imbalance: an always-depleting rule scores 95% on v3. Balanced accuracy
# averages the per-direction means, so that rule scores 0.5 by construction.

class_neutral <- rbindlist(lapply(
  list(list(d = v3, planet = "graded_floor_v3"),
       list(d = bal, planet = "graded_balanced_v1")),
  function(p) {
    acc_neg <- p$d[tau_true < 0, mean(p_sign)]
    acc_pos <- p$d[tau_true > 0, mean(p_sign)]
    data.table(planet = p$planet,
               n_neg = p$d[tau_true < 0, .N], n_pos = p$d[tau_true > 0, .N],
               acc = p$d[, mean(p_sign)],
               acc_neg = acc_neg, acc_pos = acc_pos,
               balanced_acc = (acc_neg + acc_pos) / 2,
               always_depleting = p$d[, mean(tau_true < 0)])
  }))

cat("\n=== Class-neutral performance ===\n")
print(class_neutral[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(class_neutral, "datasets/output/class_neutral_recovery.csv")

# VERDICT continued ------------------------------------------------------------

v3_lo <- res[subset == "graded_floor_v3 (69 neg / 4 pos)", lo]
v3_hi <- res[subset == "graded_floor_v3 (69 neg / 4 pos)", hi]
bal_e <- res[subset == "graded_balanced_v1 (all)" & q == 0.90, est]
cat(sprintf("\nbalanced s90 (%.2f) inside the v3 interval [%.2f, %.2f]: %s\n",
            bal_e, v3_lo, v3_hi,
            if (bal_e >= v3_lo && bal_e <= v3_hi) "YES" else "NO"))
cat("wrote datasets/output/floor_sign_balance.csv\n")
