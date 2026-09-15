cluster_recipe_families_fun <- function(well_fit, k = 4L) {

  # WARD HIERARCHICAL CLUSTERING OF RECIPES BY PER-WELL CORRELATION ############

  # Each recipe is a point in n_wells-dim space (its per-well r vector).

  mat <- well_fit$correlations  # (n_wells, n_recipes)
  # Transpose: rows = recipes (the points to cluster).
  x <- t(mat)
  col_means <- colMeans(x, na.rm = TRUE)
  nan_frac <- mean(is.na(x))
  # Mean-impute per column; if a column was entirely NA, col_means is
  # NaN — fall back to 0 so hclust doesn't break. The Python port has
  col_means[!is.finite(col_means)] <- 0
  x[is.na(x)] <- col_means[col(x)[is.na(x)]]
  x[!is.finite(x)] <- 0

  d <- stats::dist(x)
  hc <- stats::hclust(d, method = "ward.D2")
  raw_labels <- as.integer(stats::cutree(hc, k = k))

  mean_score_per_recipe <- colMeans(mat, na.rm = TRUE)

  # Canonical relabel (family 0 = best mean correlation) to keep hclust IDs in
  # parity with scipy.fcluster in the Python pipeline.
  fam_score <- vapply(sort(unique(raw_labels)), function(g) {
    mean(mean_score_per_recipe[raw_labels == g], na.rm = TRUE)
  }, numeric(1))
  rank_order <- order(fam_score, decreasing = TRUE)
  relabel_map <- integer(length(rank_order))
  relabel_map[rank_order] <- seq_along(rank_order) - 1L
  labels <- relabel_map[raw_labels]

  clusters <- data.table::copy(well_fit$recipe_keys)
  clusters[, family := labels]
  clusters[, mean_score := mean_score_per_recipe]
  list(clusters = clusters, hc = hc, nan_fraction_imputed = nan_frac)
}
