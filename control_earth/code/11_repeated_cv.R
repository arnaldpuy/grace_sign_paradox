# Repeated split-half cross-validation of recipe selection, replacing the
# single split of D036, plus joint uncertainty propagation for
# boundary / attenuation.

suppressMessages(library(data.table))
setwd(here::here("control_earth"))
floor_est <- 0.75

sc <- fread("datasets/output/score_graded_floor_v3.csv")
tr <- rbindlist(lapply(c("csr", "gfz", "jpl"), function(c)
  fread(sprintf("datasets/output/tier_a_synth_graded_floor_v3/%s/jasechko_per_recipe_trends.csv", c))))
tr[, recipe_id := paste0(centre, "_", recipe_idx_per_centre)]
tr <- merge(tr, sc[, .(aquifer_idx, tau_true)], by = "aquifer_idx")
tr[, correct := sign(trend_cm_yr) == sign(tau_true)]
sub <- tr[abs(tau_true) < floor_est]
aqs <- sort(unique(sub$aquifer_idx))
cat(sprintf("sub-floor aquifers: %d | recipes: %d\n", length(aqs),
            uniqueN(sub$recipe_id)))

# --- repeated stratified split-half -----------------------------------------
set.seed(71)
res <- rbindlist(lapply(seq_len(200L), function(k) {
  half <- sample(aqs, floor(length(aqs) / 2))
  a <- sub[aquifer_idx %in% half, .(acc = mean(correct)), by = recipe_id]
  b <- sub[!aquifer_idx %in% half, .(acc = mean(correct)), by = recipe_id]
  setorder(a, -acc)
  data.table(k = k,
             top20_holdout = b[recipe_id %in% a$recipe_id[1:20], mean(acc)],
             top50_holdout = b[recipe_id %in% a$recipe_id[1:50], mean(acc)],
             baseline = b[, mean(acc)],
             best_holdout = b[, max(acc)])
}))
res[, `:=`(adv20 = top20_holdout - baseline, adv50 = top50_holdout - baseline)]
q <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, .025), quantile(x, .975))
cat(sprintf("held-out advantage of the top 20 : %s\n", q(res$adv20)))
cat(sprintf("held-out advantage of the top 50 : %s\n", q(res$adv50)))
cat(sprintf("held-out accuracy, top 20        : %s\n", q(res$top20_holdout)))
cat(sprintf("held-out accuracy, all recipes   : %s\n", q(res$baseline)))
cat(sprintf("splits where the advantage <= 0  : %.0f%%\n", 100 * mean(res$adv20 <= 0)))

# --- boundary / attenuation with both uncertainties --------------------------
set.seed(71)
b_real <- c(est = 0.347, lo = 0.286, hi = 0.423)
# boundary draws: lognormal matched to the published bootstrap interval
lb <- log(b_real); sdb <- (lb[["hi"]] - lb[["lo"]]) / (2 * 1.96)
bd <- exp(rnorm(20000, lb[["est"]], sdb))
# attenuation draws: aquifer bootstrap of median(recovered/true), |tau|>0.5
supra <- sc[abs(tau_true) > 0.5]
at <- replicate(20000, {
  i <- sample.int(nrow(supra), replace = TRUE)
  stats::median(supra$median_rec[i] / supra$tau_true[i])
})
ratio <- bd / at
cat(sprintf("\nattenuation: %.3f [%.3f, %.3f]\n", median(at),
            quantile(at, .025), quantile(at, .975)))
cat(sprintf("boundary / attenuation: %.2f [%.2f, %.2f] cm/yr\n",
            median(ratio), quantile(ratio, .025), quantile(ratio, .975)))
cat(sprintf("overlap with the 0.90 floor CI [0.47, 1.50]: %s\n",
            ifelse(quantile(ratio, .025) < 1.50 & quantile(ratio, .975) > 0.47,
                   "yes", "no")))

# --- bank the manuscript quantities ------------------------------------------

# One row per published statistic, so the release notebook can pin them.
smry <- data.table(
  quantity = c("top20_holdout", "baseline_holdout", "adv20",
               "attenuation", "boundary_over_attenuation"),
  median = c(median(res$top20_holdout), median(res$baseline), median(res$adv20),
             median(at), median(ratio)),
  lo = c(quantile(res$top20_holdout, .025), quantile(res$baseline, .025),
         quantile(res$adv20, .025), quantile(at, .025), quantile(ratio, .025)),
  hi = c(quantile(res$top20_holdout, .975), quantile(res$baseline, .975),
         quantile(res$adv20, .975), quantile(at, .975), quantile(ratio, .975)))
fwrite(smry, "datasets/output/recipe_skill/repeated_cv_summary.csv")
fwrite(res, "datasets/output/recipe_skill/repeated_cv_splits.csv")
cat("banked repeated_cv_summary.csv + repeated_cv_splits.csv\n")
