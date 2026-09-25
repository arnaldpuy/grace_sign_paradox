# Provenance manifest

One row per data file in the repository: the script or notebook chunk that writes it,
the reproduction track it needs (README, *How to run*: 1 = the notebook and the shipped
tables; 2 = the per-aquifer monthly series, rebuilt from the recipe cubes; 3 = the
recipe cubes or another raw input; 4 = the control-Earth rebuild; input = a curated or
third-party input with no generator here), and the notebook chunks that read it. Every
number in the paper is asserted from these files by the `assert-numbers` chunk of
`code/code_main_analysis.Rmd`; the Supplementary Materials' *Code* boxes map each
supplementary item to its chunks. Regenerate this file with the command in the
README (*Provenance*).


## `control_earth/datasets/output/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `aux_subtraction_balanced.csv` | `control_earth/code/18_aux_subtraction_balanced.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `balanced_layouts_bins.csv` | `control_earth/code/17_balanced_layouts.R` | 1 | `ce-load` |
| `balanced_layouts_floors.csv` | `control_earth/code/17_balanced_layouts.R` | 1 | `ce-load` |
| `balanced_layouts_summary.csv` | `control_earth/code/17_balanced_layouts.R` | 1 | `ce-load` |
| `boundary_balanced.csv` | `control_earth/code/15_boundary_balanced.R` | 4 (synthetic per-recipe tables) | - |
| `calibration_report.csv` | `control_earth/code/03_calibration.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `calibration_report_2026-07-14_real276.csv` | preserved first gate report (real ensemble on 276 aquifers), superseded by calibration_report.csv | input | - |
| `class_neutral_recovery.csv` | `control_earth/code/14_floor_sign_balance.R` | 1 | - |
| `floor_sign_balance.csv` | `control_earth/code/14_floor_sign_balance.R` | 1 | - |
| `parametric_floor_v3.csv` | `control_earth/code/09_parametric_floor.R` | 1 | `ce-load` |

## `control_earth/datasets/output/recipe_skill/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `median_estimator_coefs.csv` | `control_earth/code/10_recipe_skill.R` | 4 (synthetic per-recipe tables) | - |
| `per_recipe_logistic_coefs.csv` | `control_earth/code/10_recipe_skill.R` | 4 (synthetic per-recipe tables) | - |
| `per_recipe_skill_graded_floor_v3.csv` | `control_earth/code/10_recipe_skill.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `repeated_cv_summary.csv` | `control_earth/code/11_repeated_cv.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `split_half_replication.csv` | `control_earth/code/10_recipe_skill.R` | 4 (synthetic per-recipe tables) | `ce-load` |

## `control_earth/datasets/output/rolling_null/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `rolling_attribution_planets.csv` | `control_earth/code/07c_rolling_attribution.py` | 4 (truth fields) | `ce-rolling-attribution` |
| `rolling_attribution_replicates.csv` | `control_earth/code/07c_rolling_attribution.py` | 4 (truth fields) | `ce-rolling-attribution` |
| `rolling_null_summary.csv` | `control_earth/code/07b_rolling_null.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `rolling_truth_window_trends.csv` | `control_earth/code/07c_rolling_attribution.py` | 4 (truth fields) | - |
| `rolling_window_signs_graded_balanced_v3.csv` | `control_earth/code/07_rolling_windows_null.py` | 4 (synthetic cubes) | - |
| `rolling_window_signs_graded_balanced_v3_noise2.csv` | `control_earth/code/07_rolling_windows_null.py` | 4 (synthetic cubes) | - |
| `rolling_window_signs_graded_floor_v2.csv` | `control_earth/code/07_rolling_windows_null.py` | 4 (synthetic cubes) | - |
| `rolling_window_signs_graded_floor_v3.csv` | `control_earth/code/07_rolling_windows_null.py` | 4 (synthetic cubes) | - |
| `rolling_window_signs_pilot_uniform_depletion_v2.csv` | `control_earth/code/07_rolling_windows_null.py` | 4 (synthetic cubes) | - |
| `window_dispersion_summary.csv` | `control_earth/code/07b_rolling_null.R` | 4 (synthetic per-recipe tables) | - |
| `window_paradox_summary.csv` | `control_earth/code/07b_rolling_null.R` | 4 (synthetic per-recipe tables) | - |

## `control_earth/datasets/output/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `score_graded_balanced_v1.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | - |
| `score_graded_balanced_v1_reversed.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | - |
| `score_graded_balanced_v2.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | - |
| `score_graded_balanced_v3.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | - |
| `score_graded_balanced_v3_noise2.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | - |
| `score_graded_floor_v2.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | `ce-load` |
| `score_graded_floor_v3.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | `ce-load` |
| `score_gws_crossed_era5_v1.csv` | `control_earth/code/12_score_crossed_gws.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `score_pilot_uniform_depletion_v2.csv` | `control_earth/code/05_score.R <planet>` | 4 (synthetic per-recipe tables) | `ce-load` |
| `sign_interleaving.csv` | `control_earth/code/19_sign_interleaving.R` | 1 | `ce-load` |
| `truth_gw_aux_trends_crossed_era5_v1.csv` | `control_earth/code/forward/12_crossed_aux_trends.py (per planet)` |  | `ce-load` |

## `control_earth/datasets/output/val_of_val/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `truth_gw_aquifer_graded_floor_v2.csv` | `control_earth/code/08a_synth_wells.py` | 4 (synthetic cubes; wells) | - |
| `truth_gw_aquifer_graded_floor_v3.csv` | `control_earth/code/08a_synth_wells.py` | 4 (synthetic cubes; wells) | - |
| `truth_gw_aquifer_pilot_uniform_depletion_v2.csv` | `control_earth/code/08a_synth_wells.py` | 4 (synthetic cubes; wells) | - |
| `val_of_val_aquifer_graded_floor_v3.csv` | `control_earth/code/08b_val_of_val.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `val_of_val_mc_numeric.csv` | `control_earth/code/08b_val_of_val.R` | 4 (synthetic per-recipe tables) | `ce-load` |
| `val_of_val_summary.csv` | `control_earth/code/08b_val_of_val.R` | 4 (synthetic per-recipe tables) | - |
| `well_noise_calibration.csv` | `control_earth/code/08a_synth_wells.py` | 4 (synthetic cubes; wells) | - |

## `control_earth/datasets/truth/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `params_iter2.json` | planet specification, answer key or manifest written by control_earth/code/forward/ (track 4) | input | - |
| `params_iter3.json` | planet specification, answer key or manifest written by control_earth/code/forward/ (track 4) | input | - |
| `truth_aquifer_trends_crossed_era5_v1.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_balanced_v1.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_balanced_v1_reversed.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_balanced_v2.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_balanced_v3.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_balanced_v3_noise2.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_floor_v1.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_floor_v2.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_graded_floor_v3.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_pilot_uniform_depletion.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_aquifer_trends_pilot_uniform_depletion_v2.csv` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_gw_aux_trends_crossed_era5_v1.csv` | `control_earth/code/forward/12_crossed_aux_trends.py (per planet)` |  | `ce-load` |
| `truth_gw_aux_trends_graded_balanced_v1.csv` | `control_earth/code/forward/12_crossed_aux_trends.py (per planet)` |  | - |
| `truth_gw_aux_trends_graded_balanced_v2.csv` | `control_earth/code/forward/12_crossed_aux_trends.py (per planet)` |  | - |
| `truth_gw_aux_trends_graded_balanced_v3.csv` | `control_earth/code/forward/12_crossed_aux_trends.py (per planet)` |  | - |
| `truth_manifest_crossed_era5_v1.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_balanced_v1.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_balanced_v1_reversed.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_balanced_v2.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_balanced_v3.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_balanced_v3_noise2.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_floor_v1.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_floor_v2.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_graded_floor_v3.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_pilot_uniform_depletion.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_manifest_pilot_uniform_depletion_v2.json` | `control_earth/code/forward/01_truth_builder.py (per planet)` |  | - |
| `truth_spec_crossed_era5_v1.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_balanced_v1.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_balanced_v1_reversed.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_balanced_v2.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_balanced_v3.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_balanced_v3_noise2.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_floor_v1.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_floor_v2.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_graded_floor_v3.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_pilot.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |
| `truth_spec_pilot_v2.json` | `control_earth/code/forward/13_make_balanced_spec.py (per planet)`; `control_earth/code/forward/16_make_balanced_specs_batch.py (per planet)` |  | - |

## `datasets/bibliometric/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `bibliometric_master_all.csv` | notebook chunk `gws-strata` | 1 | `robustness-lit-weights`, `gws-strata`, `corpus-recipe-prevalence`, `literature-stats`, `policy-synthesis-table` |
| `claim_signs.csv` | LLM claim-sign extraction plus 24 hand-coded papers (Supplementary Materials) | input | `literature-stats`, `worked-example-murray` |
| `corpus_provenance.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | - |
| `literature_candidate_basins.txt` | superseded basin list of the earlier region-level matching | input | - |
| `literature_footprint_basins.txt` | the 37 cohort polygons of the coded footprints; input of tier_a_window_trends.py --windows | input | - |
| `method_adjudication.csv` | curated: stated processing choices outside the recipe space (unsupported) and family-level restrictions, one row per record and axis | input | `literature-stats`, `assert-numbers` |

## `datasets/bibliometric/option_b/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `d_regime_aquifers.csv` | notebook chunk `aquifer-regimes` | 1 | `aquifer-regimes` |
| `d_regime_reclassification.csv` | notebook chunk `aquifer-regimes` | 1 | `aquifer-regimes` |

## `datasets/bibliometric/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `paper_footprints.csv` | hand-coded study footprints, one row per paper, with the coding note on coverage and claim type | input | `literature-stats` |
| `paper_study_windows.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `literature-stats`, `fig-consequences` |
| `policy_citations_fullcorpus.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-uptake-stats` |
| `policy_context_final_classification.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-context-stats`, `fig-consequences` |
| `policy_context_retrieval.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-context-stats` |
| `policy_documents_fullcorpus.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-uptake-stats` |
| `policy_fullcorpus_by_region.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-uptake-stats`, `fig-consequences`, `policy-synthesis-table` |
| `policy_triangulation_direct.csv` | curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods) | input | `policy-triangulation-direct` |
| `record_adjudication.csv` | curated: duplicate article records dropped and claims excluded from trend risk, one reason per row | input | `literature-stats`, `assert-numbers` |

## `datasets/jasechko_2024/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `cell_masks_0p5deg.csv` | derived from Jasechko et al. (2024), doi:10.5281/zenodo.10003697 | input | - |
| `cell_masks_0p5deg_276.csv` | derived from Jasechko et al. (2024), doi:10.5281/zenodo.10003697 | input | `literature-stats` |
| `cell_masks_1p0deg.csv` | derived from Jasechko et al. (2024), doi:10.5281/zenodo.10003697 | input | - |
| `cohort_50K_n73.csv` | `code/grace_pipeline/polygon_masks.py` |  | `prelim` |
| `cohort_climate_zones.csv` | `code/grace_pipeline/climate_zones.py` | 3 (climate raster) | `prelim` |

## `datasets/output/robustness/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `aquifer_size_bins.csv` | size-bin definition (hand-set input) | input | `full_first_figure2` |
| `common_months_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_common_months.py` | 2 (monthly series) | `common-months` |
| `leakage_mask_sensitivity.csv` | `code/grace_pipeline/tier_a_leakage.py` | 3 (recipe cubes) | `leakage-median` |
| `leakage_ring_per_recipe.csv.gz` | `code/grace_pipeline/tier_a_leakage_ring.py` | 3 (recipe cubes) | `leakage-median` |
| `leakage_summary.csv` | `code_robustness_analysis.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository; summarises `leakage_mask_sensitivity.csv` | 1 | `full_first_figure2`, `fig-tier2-size-check` |
| `literature_attested_recipes.csv` | notebook chunk `gws-strata` | 1 | `gws-strata` |
| `mascon_vs_recipe_ensemble.csv` | `code_robustness_analysis.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `fig-disagreement`, `assert-numbers` |
| `nested_recipe_space_summary.csv` | `code_robustness_analysis.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository (the attested-space rows are superseded by chunk `gws-strata`, see the folder README) | 1 | `full_first_figure2`, `gws-strata` |
| `trend_estimator_anchorcell.csv` | `code_robustness_analysis.R` (private working repository, not distributed): estimator comparison on the anchor-cell series of the IGP recipe cube; the cohort-wide comparison is `trend_estimator_global.csv.gz` | 3 (recipe cubes) | `tab-paradox-robustness` |
| `trend_estimator_global.csv.gz` | `code/estimator_global.py` | 2 (monthly series) | `tab-paradox-robustness` |
| `trend_window_sensitivity.csv` | `code_robustness_analysis.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `tab-paradox-robustness` |
| `well_grace_sign_validation.csv` | `code_robustness_analysis.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository (`functions/fit_recipes_to_wells_fun.R`) | 1 | `full_first_figure2`, `well-binomial` |

## `datasets/output/section1/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `fig1_flip_timeseries.csv` | `reduce_section1.R` (private working repository, not distributed): per-recipe time series of the anchor cell, read from the IGP recipe cube | 3 (recipe cubes) | `load-section1` |
| `fig1_recipe_families.csv` | `reduce_section1.R` (private working repository, not distributed): Ward clustering of the anchor-cell recipe series (`functions/cluster_recipe_families_fun.R`) | 3 (recipe cubes) | `load-section1` |

## `datasets/output/tier_a/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `fig_gradient_data.csv` | `tier_a_gradient_figure.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `load-tier-a`, `sa-boundary-ci` |
| `grace_month_midpoints.csv` | `code/grace_pipeline/tier_a_window_trends.py --emit-monthly` (the month midpoints of the monthly series) | 2 (monthly series) | - |
| `gws_aux_trends.csv` | `code/grace_pipeline/gws_aux_aquifer.py` | 2 (monthly series; GLDAS, ERA5-Land) | `robustness-gws` |
| `gws_aux_trends_era5land.csv` | `code/grace_pipeline/gws_aux_aquifer.py --source era5land` | 2 (monthly series; ERA5-Land) | `robustness-gws` |
| `gws_aux_window_trends_gldas.csv` | `code/grace_pipeline/gws_aux_aquifer.py --windows` | 2 (monthly series; GLDAS) | `claims-gws` |
| `jasechko_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_per_aquifer.py` (one run per cohort and window, see `functions/resolve_tier_a_window.R`) | 3 (recipe cubes) | `load-tier-a`, `sa-boundary-ci`, `fig-tier2-sign-flip`, `worked-example-murray` |
| `per_recipe_trend_se.csv.gz` | `code/grace_pipeline/tier_a_trend_se.py` | 2 (monthly series) | `noise-boundary` |
| `recipe_stripe_index.csv` | `code/grace_pipeline/tier_a_stripe_index.py` | 2 (monthly series) | `robustness-stripes` |
| `rolling_window_signs.csv` | `code/grace_pipeline/tier_a_rolling_windows.py` | 2 (monthly series) | `rolling-windows` |
| `sensitivity_1deg_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_1deg_gain.py` | 3 (recipe cubes; CLM4 gain grid) | `robustness-1deg-gain` |
| `sensitivity_gain_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_1deg_gain.py` | 3 (recipe cubes; CLM4 gain grid) | `robustness-1deg-gain` |
| `surface_water_exposure.csv` | `code/surface_water_exposure.R` | 1 | `gws-exposure` |

## `datasets/output/tier_a_full/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `gradient_data_276.csv` | `tier_a_two_tier_stats.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `boundary-ci`, `sa-boundary-ci`, `full_first_figure2`, `load-tier-a-276` |
| `jasechko_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_per_aquifer.py` (one run per cohort and window, see `functions/resolve_tier_a_window.R`) | 3 (recipe cubes) | `load-tier-a`, `sa-boundary-ci`, `fig-tier2-sign-flip`, `worked-example-murray` |
| `sign_paradox_276.csv` | `tier_a_two_tier_stats.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `sa-boundary-ci`, `full_first_figure2`, `load-tier-a-276` |
| `size_check_per_aquifer.csv` | `tier_a_size_check.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `merge_sign_flip`, `load-tier-a-276` |
| `variance_decomp_276.csv` | `tier_a_two_tier_stats.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository (`functions/aquifer_variance_decomp_fun.R`) | 1 | `load-tier-a-276`, `assert-numbers` |

## `datasets/output/tier_a_full_73/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `aquifer_regimes.csv` | notebook chunk `aquifer-regimes` | 1 | `aquifer-regimes` |
| `fig_gradient_data.csv` | `tier_a_gradient_figure.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository | 1 | `load-tier-a`, `sa-boundary-ci` |
| `global_top_bottom_recipes.csv` | `tier_a_well_skill.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository (`functions/fit_recipes_to_wells_fun.R`) | 1 | `load-tier-a` |
| `jasechko_aquifer_metrics.csv` | `code/grace_pipeline/tier_a_per_aquifer.py` | 3 (recipe cubes) | `load-tier-a` |
| `jasechko_per_recipe_trends.csv.gz` | `code/grace_pipeline/tier_a_per_aquifer.py` (one run per cohort and window, see `functions/resolve_tier_a_window.R`) | 3 (recipe cubes) | `load-tier-a`, `sa-boundary-ci`, `fig-tier2-sign-flip`, `worked-example-murray` |
| `mascon_aquifer_metrics.csv` | `code/grace_pipeline/tier_a_mascon.py` | 3 (mascon products; wells) | `load-tier-a` |
| `mascon_per_aquifer_trends.csv` | `code/grace_pipeline/tier_a_mascon.py` | 3 (mascon products; wells) | `load-tier-a` |
| `mascon_well_level_skill.csv` | `code/grace_pipeline/tier_a_mascon.py` | 3 (mascon products; wells) | `load-tier-a` |
| `per_recipe_global_skill.csv` | `tier_a_well_skill.R`: analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository (`functions/fit_recipes_to_wells_fun.R`) | 1 | `load-tier-a` |
| `well_level_global_best_worst.csv` | `code/grace_pipeline/tier_a_well_level.py` | 3 (wells) | `load-tier-a` |
| `well_level_per_recipe_skill.csv` | `code/grace_pipeline/tier_a_well_level.py` | 3 (wells) | `load-tier-a` |
| `wells_in_cohort.csv` | `code/grace_pipeline/well_aquifer_join.py` | 3 (Jasechko et al. 2024 wells) | `load-tier-a` |

## `datasets/output/tier_a_windows/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `jasechko_per_recipe_trends__2002-04_2002-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2004-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2008-01.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2008-06.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2008-08.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2009-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2009-08.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2010-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2012-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2013-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2014-03.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2015-04.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2016-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2016-08.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2016-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2017-06.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2019-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2022-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-04_2023-06.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-07_2011-06.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-07_2017-06.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2007-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2007-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2008-01.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2008-04.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2008-10.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2014-09.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-08_2016-05.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2002-09_2004-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2005-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2006-05.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2006-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2007-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2009-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2010-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2012-09.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2012-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2014-05.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2014-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2015-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2016-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2019-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2021-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-01_2021-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-02_2012-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-05_2009-04.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-07_2005-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2003-08_2010-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2004-01_2015-01.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2005-01_2013-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2005-01_2015-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2005-06_2012-10.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2009-01_2016-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2013-01_2015-12.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |
| `jasechko_per_recipe_trends__2016-01_2016-07.csv.gz` | `code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons) | 2 (monthly series) | `literature-stats` |

## `figures/`

| File | Written by | Track | Read by notebook chunk |
|---|---|---|---|
| `preprocessing_choice_space.pdf` | static TikZ source and its PDF | input | `fig-paper-1-flip` |
| `preprocessing_choice_space.tex` | static TikZ source and its PDF | input | `fig-paper-1-flip` |
