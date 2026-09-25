# Robustness outputs

Banked inputs for the robustness section of `code/code_main_analysis.Rmd`.

One supersession note. The `Published (attested)` rows of
`nested_recipe_space_summary.csv` (72 recipes, 27/73 sign-paradoxical
aquifers, 0.09 cm/yr boundary) predate the human adjudication of the
bibliometric corpus and are retained only as a historical record. The
literature-attested recipe space used in the manuscript is rebuilt from the
adjudicated corpus by the notebook itself (126 recipes, 33/73, 0.14 cm/yr;
chunk `gws-strata`), which exports the exact recipe set to
`literature_attested_recipes.csv` in this folder and asserts its size,
counts and boundary at every knit. All other rows of
`nested_recipe_space_summary.csv` remain current.

Generators. `leakage_mask_sensitivity.csv` is written by
`code/grace_pipeline/tier_a_leakage.py` and `leakage_ring_per_recipe.csv.gz` by
`tier_a_leakage_ring.py`, both from the recipe cubes;
`common_months_per_recipe_trends.csv.gz` by `tier_a_common_months.py` from the
per-aquifer monthly series; `trend_estimator_global.csv.gz` by
`code/estimator_global.py`. The remaining tables were produced by analysis
scripts of the private working repository from the per-recipe tables shipped
here; `PROVENANCE.md` at the repository root names each one.
