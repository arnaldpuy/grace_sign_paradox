# recipe_skill — per-recipe sign skill on known truth

Computed on the sign-imbalanced graded planet's banked per-recipe trend
tables. Tests whether a single well-chosen recipe recovers weak-signal signs
better than the ensemble.

Files:

- `per_recipe_skill_graded_floor_v3.csv` — one row per recipe (1,920):
  sign-accuracy below and above 0.75 cm/yr in the planted trend (the weak
  stratum holds 49 of the 73 aquifers), per-recipe logistic s90 (unstable, do
  not headline), axis metadata, fully-corrected flag.
- `split_half_replication.csv` — per-recipe weak-stratum accuracy on the
  selection half vs the held-out half (seed 71).
- `per_recipe_logistic_coefs.csv`, `median_estimator_coefs.csv` — logistic
  fits behind the recovery curves.

Headline: no recipe reaches 0.95 weak-stratum accuracy (best 0.939, median
0.816); the 180 fully-corrected recipes are indistinguishable from the full
ensemble; recipes selected on half the aquifers keep a modest held-out
advantage over 200 stratified splits (0.854 [0.754, 0.928] vs 0.799
[0.743, 0.864] for the ensemble at large; `code/11_repeated_cv.R` in the
control_earth repository) but never approach the 0.95 level at which a sign
is dependable — the single split banked in `split_half_replication.csv`
understates the transfer. The ensemble median is the best available sign
estimator (weak stratum 0.898, beating 94% of single recipes).
