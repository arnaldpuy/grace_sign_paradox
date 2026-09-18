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
