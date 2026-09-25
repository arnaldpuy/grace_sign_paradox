<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX) — minted at release -->

# One satellite record yields opposite groundwater trends

[Arnald Puy](https://www.arnaldpuy.com), Samuel Flinders, Seth N. Linga, Carmen Aguiló-Rivera, Ana-Maria Crețu

R and Python code accompanying the paper *One satellite record yields opposite
groundwater trends*, whose abstract is the following:

## Abstract

*The sign of satellite-gravity groundwater trends (depletion or recovery) anchors global water security. Whether this sign is set by the aquifers or by the processing that produces it remains unknown. Here we evaluate 1,920 pre-processing recipes across the full GRACE record (2002–2025) for 73 of the world's largest aquifer systems and show that these choices can change whether three in four aquifers are reported as gaining or losing water. A detectability boundary near 0.35 cm yr⁻¹ equivalent water height marks the regime where the sign of a GRACE trend is defined by pre-processing choices rather than by the aquifer. Eight control Earths with known trends reproduce this boundary wherever it can be estimated and show that there is no universal threshold guaranteeing correct recovery of the storage trend. Such underdetermination matters for policy: although most published estimates are consistent with the trend of the admissible recipe ensemble, one third of aquifers lack GRACE studies and admit both depletion and recovery narratives under admissible processing choices. Our results distinguish robust trends from those conditioned by pre-processing choices and provide a framework for carrying this uncertainty into groundwater assessments and policy making.*

## Contents

The repository combines the reproduction notebook for every figure, table and
in-text number of the paper and its Supplementary Materials with the Python
pipeline that built the underlying 1,920-recipe trend ensembles.

- [Background](#background)
- [Repository structure](#repository-structure)
- [Dependencies](#dependencies)
- [How to run](#how-to-run)
- [Functions](#functions)
- [Key analytical concepts](#key-analytical-concepts)
- [Outputs](#outputs)
- [Provenance](#provenance)
- [Not distributed](#not-distributed)
- [Citation](#citation)
- [License](#license)

## Background

Converting GRACE/GRACE-FO Level-2 spherical-harmonic coefficients into a
water-storage trend requires a chain of defensible choices: processing centre,
truncation degree, C20/C30 treatment, geocentre restoration, GIA model and
destriping/smoothing filter. The full factorial gives 1,920 recipes. The paper
propagates all of them into per-aquifer trends for the world's largest
aquifers, identifies the detectability boundary below which the sign of the
trend is set by the recipe rather than the aquifer, and confirms the boundary
on synthetic planets with stipulated truths ("control Earth") observed
through a forward model of the GRACE measurement system. Three of these planets
balance depletion against recovery at fixed trend magnitudes, and show that the
true trend needed to recover a direction depends on how depleting and
recovering aquifers are arranged in space.

## Repository structure

```
grace_sign_paradox/
├── README.md, LICENSE, requirements.txt
├── PROVENANCE.md                    # every data file: generating script, reproduction track, reading chunk
├── code/
│   ├── code_main_analysis.Rmd       # THE notebook: every figure, table and number
│   ├── code_main_analysis.pdf       # its knitted certificate (assertions passed)
│   ├── grace_pipeline/              # Python: Level-2 archives -> 1,920-recipe cubes -> tables
│   ├── estimator_global.py, surface_water_exposure.R, make_provenance.py, check_imports.py
├── functions/                       # 14 R helpers sourced by the notebook
├── datasets/
│   ├── jasechko_2024/               # aquifer cohorts, climate zones, cell masks
│   ├── bibliometric/                # curated literature/policy corpus tables
│   └── output/                      # the pre-reduced analysis tables the notebook reads
├── control_earth/                   # the synthetic-truth experiment (see control_earth/README.md)
│   ├── code/forward/                # truth builder, forward model, recovery driver, planet specs, launch scripts
│   ├── code/                        # calibration gate, scoring and every analysis behind the summary tables
│   ├── functions/                   # R helper sourced by the control-Earth scripts
│   ├── datasets/truth/              # per planet: specification, SHA-256 manifest and answer key
│   └── datasets/output/             # control-Earth summary tables (gate, crossings, scores)
└── figures/                         # the choice-space diagram (TikZ source + PDF), the notebook's only static figure input
```

> Note: the three per-centre recipe cubes (about 110 GiB each), the
> per-aquifer monthly series reduced from them, the synthetic Level-2 records
> and the control-Earth planet fields are not distributed (see
> [Not distributed](#not-distributed)): the code here rebuilds them (tracks
> 2–4 below) and they are available from the corresponding author. Every
> number in the paper is asserted from the tables in `datasets/` and
> `control_earth/datasets/`, and `PROVENANCE.md` names the script behind each
> table. Every planet is rebuildable from its specification
> in `control_earth/datasets/truth/` with the code in `control_earth/code/`. The frozen design, pre-registration and decision register of
> the control-Earth experiment are deposited at
> [doi:10.5281/zenodo.21454588](https://doi.org/10.5281/zenodo.21454588), and the
> amendment specifying the balanced planets, their analyses and their stopping
> rules — frozen before those planets were run — at
> [doi:10.5281/zenodo.21819495](https://doi.org/10.5281/zenodo.21819495).

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | all table operations |
| `ggplot2`, `cowplot`, `scales`, `ggrepel`, `viridisLite`, `RColorBrewer` | figures |
| `here` | repository-rooted paths |
| `maps` | base maps |
| `sensobol` | variance decomposition and package loading |
| `mgcv` | loess/GAM fits |
| `knitr`, `rmarkdown` | knitting |
| `R.utils` | transparent reads of the gzipped tables |

Minimum R version: 4.5. The notebook installs `sensobol` if absent and loads
the rest through `sensobol::load_packages()`. The Python pipeline (tracks 2–3
below) needs Python ≥ 3.12 and the packages pinned in `requirements.txt`
(`pip install -r requirements.txt`), plus a free NASA Earthdata account.
`python code/check_imports.py` verifies the installation.

### Versions used for the published results

Later versions are expected to work; these produced the numbers in the paper
(from the recorded session information of the knits and the pipeline
environment).

| | Version |
|---|---|
| R | 4.5.2 (aarch64-apple-darwin20) |
| `data.table` | 1.18.2.1 |
| `ggplot2` | 4.0.3.9000 |
| `sensobol` | 1.1.9 |
| `mgcv` | 1.9-3 |
| `scales` | 1.4.0 · `cowplot` 1.2.0.9000 · `ggrepel` 0.9.6 |
| Python | 3.12.4 |
| `numpy` | 2.4.4 · `scipy` 1.17.1 |
| `xarray` | 2026.4.0 · `netCDF4` 1.7.4 |
| `gravity-toolkit` | 1.2.5 |
| `geopandas` | 1.1.3 · `shapely` 2.1.2 · `h5py` 3.16.0 |

### Random seeds

Every stochastic step is seeded, so all reported numbers are exactly
reproducible:

| Step | Seed |
|---|---|
| Detectability-boundary and recovery-crossing bootstraps (2,000 resamples) | 71 |
| Split-half recipe-skill replication | 71 |
| Synthetic planet noise — uniform depletion | 20260714 |
| Synthetic planet noise — mixed sign | 20260716 |
| Synthetic planet noise — graded magnitude | 20260717 |
| Synthetic planet noise — balanced signs 1 | 20260805 |
| Synthetic planet noise — balanced signs 2 | 20270805 |
| Synthetic planet noise — balanced signs 3 | 20280805 |
| Synthetic planet noise — balanced signs 1, sign-reversed | 20260805 (unchanged, so the noise field is identical) |
| Synthetic planet noise — balanced signs 3, noise replicate | 20290805 |

The planet seeds are also fixed in the truth specifications deposited with the
pre-registration, so a planet can be rebuilt from its specification alone.

## How to run

Four tracks, by increasing depth:

1. **Fast (under an hour).** Set the working directory to the repository root
   and knit the notebook:

   ```r
   rmarkdown::render("code/code_main_analysis.Rmd")
   ```

   This regenerates every computed figure, table and statistic of the paper
   and its supplement from `datasets/`, and ends with a `Numbers in the text`
   section that asserts every published value with `stopifnot()`: **the knit
   fails if any number does not reproduce.** The shipped
   `code_main_analysis.pdf` is the passing certificate.

2. **Intermediate (hours, ~330 GB disk).** From the three per-centre recipe
   cubes (`grace_ewh_global_<centre>_appendix.nc` under
   `$GRACE_DATA_ROOT/processed/grace/global/`, built by track 3 or obtained
   from the authors), rebuild `datasets/output/`. Run every command from
   `code/` (`cd code`). `python -m grace_pipeline.tier_a_per_aquifer` reduces
   a cube to the per-aquifer, per-recipe trend tables;
   `python -m grace_pipeline.tier_a_window_trends --emit-monthly` (one pass
   per cube, about 1.5 h) writes the per-aquifer monthly series
   `datasets/output/tier_a_windows/monthly_<centre>.npz`, from which
   `tier_a_window_trends --windows` (the per-study window tables),
   `tier_a_common_months`, `tier_a_trend_se`, `tier_a_rolling_windows`,
   `tier_a_stripe_index`, `gws_aux_aquifer` and `code/estimator_global.py`
   rebuild their tables in seconds to minutes; `tier_a_leakage`,
   `tier_a_leakage_ring` and `tier_a_1deg_gain` read the cubes directly
   (about 35 min each). The window, standard-error, gain, common-month,
   leakage and ring scripts validate their output against the shipped table
   it replaces and print the maximum difference.

3. **Full (days).** Rebuild the cubes from the public Level-2 archives:
   `code/grace_pipeline/download_grace_shm.py` fetches the CSR/JPL/GFZ
   releases and correction series; `code/grace_pipeline/pipeline_global.py`
   synthesises the 1,920-recipe cubes. The DDK decorrelation kernels are not
   vendored: `python -m grace_pipeline.download_ddk` fetches the five levels
   the pipeline applies from a pinned commit of
   [strawpants/GRACE-filter](https://github.com/strawpants/GRACE-filter)
   (MIT) into `$GRACE_DATA_ROOT/raw/grace/auxiliary/ddk/` and verifies each
   file against the MD5 checksum of the copy used for the published results
   (`--check` verifies without downloading). The other auxiliary inputs
   (GLDAS, ERA5-Land, mascons, SLR C20/C30, degree-1, GIA models) are fetched
   by the `download_*.py` scripts or documented, with their sources and file
   names, in `config.py` and `io_aux.py`.

4. **Control Earth (about 14.5 hours per planet).** Rebuild a synthetic
   planet from its specification and score it: the launch scripts in
   `control_earth/code/forward/` run truth builder, forward model and
   recovery in sequence, `control_earth/code/05_score.R` scores the result
   against the answer key, and `control_earth/code/17_balanced_layouts.R`
   regenerates the summary tables. Synthetic Level-2 files and cubes go to
   `GRACE_CE_CUBE_DIR` (default `control_earth/datasets/synthetic_l2/`).
   Details in `control_earth/README.md`.

The pipeline reads and writes a GRACE data tree (`raw/`, `processed/`)
located by the environment variable `GRACE_DATA_ROOT`; it defaults to `data/`
in this repository. All other paths are relative to the repository.

## Functions

```r
r_functions <- list.files(path = here("functions"),
                          pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_functions, source))
```

| Function | Purpose |
|---|---|
| `load_grace_recipes_fun` | loads a per-recipe NetCDF cube into a long `data.table` |
| `compute_recipe_trend_fun` | OLS trend per recipe, optionally on a time slice |
| `preprocessing_dominance_fun` | per-aquifer dominance R = IQR/\|median\| across recipes |
| `detectability_boundary_fun` | the R = 1 (and 95%-sign-agreement) crossing of the loess fit, with bootstrap CI |
| `recipe_subspace_metrics_fun` | headline statistics for an arbitrary recipe subset |
| `aquifer_variance_decomp_fun` | factorial variance decomposition of the trends into the seven recipe axes |
| `cluster_recipe_families_fun` | Ward clustering of recipes into families |
| `fit_recipes_to_wells_fun` | recipe-vs-well skill metrics (W1 per-well correlation, W2 trend rank) |
| `load_wris_wells_fun`, `station_storage_anomaly_fun` | well-record ingestion and per-station storage anomalies |
| `literature_recipe_match_fun` | maps each paper's stated methods to its compatible recipe set (unstated axes match all levels) |
| `window_match_fun` | assigns each paper its exact study window and slices the matching per-window ensembles |
| `resolve_tier_a_window` | single source of truth for the analysis window |
| `theme_AP` | the ggplot2 theme used throughout |

The control-Earth scripts additionally source
`control_earth/functions/recoverability_floor_fun.R`, which fits the
sign-recovery crossings on a scored planet.

## Key analytical concepts

### Recipe

One complete path from the Level-2 spherical-harmonic coefficients to a
water-storage trend: processing centre × truncation degree × C20 and C30
treatment × geocentre restoration × GIA model × destriping-and-smoothing
filter. The full factorial gives 1,920 recipes.

### Pre-processing dominance R

Per aquifer, the inter-quartile range of the 1,920 recipe trends divided by
the absolute recipe-median trend. R < 1: the signal dominates the recipe
spread. R > 1: pre-processing choices dominate the signal.

### Sign paradox

An aquifer whose recipe ensemble contains both positive and negative trends:
whether it is gaining or losing water depends on the recipe.

### Detectability boundary

The recipe-median trend magnitude at which the loess fit of R against signal
strength crosses R = 1 (0.35 cm yr⁻¹ EWH, with bootstrap CI), below which the
sign of a GRACE trend is recipe-dependent. Cross-checked by the magnitude at
which 95% of recipes agree on the sign (0.46 cm yr⁻¹), by ~16,000 in-situ
wells and by synthetic planets with known truths, on which it is reproduced
wherever it is estimable.

## Outputs

### The per-aquifer recipe-trend tables

Users who only want the trend ensembles can read them directly, without
knitting anything. One row per (centre, aquifer, recipe) with the seven
recipe axes spelled out:

| Table | Contents |
|---|---|
| `datasets/output/tier_a_full_73/jasechko_per_recipe_trends.csv.gz` | the 73-aquifer Tier-1 cohort, full 2002–2025 record (140,160 rows) |
| `datasets/output/tier_a_full/jasechko_per_recipe_trends.csv.gz` | the 276-aquifer two-tier cohort (529,920 rows) |
| `datasets/output/tier_a_windows/jasechko_per_recipe_trends__*.csv.gz` | the same ensembles sliced to each published study's window (55 windows), restricted to the 37 cohort polygons of the coded study footprints |
| `datasets/output/robustness/common_months_per_recipe_trends.csv.gz` | the 73-aquifer ensemble refitted on the 241 months valid in every recipe of every centre |
| `datasets/output/robustness/leakage_ring_per_recipe.csv.gz` | per (centre, aquifer, recipe): the trend over the aquifer's own cells and over its one-cell exterior ring |
| `datasets/output/tier_a/per_recipe_trend_se.csv.gz` | per-recipe trend standard errors (plain and AR(1)) on the 73-aquifer cohort |

Columns of the trend tables: `centre`, `aquifer_idx`, `Study_area`,
`recipe_idx_per_centre`, `trend_cm_yr`, `truncation`, `filter`, `gia_model`,
`c20_treatment`, `c30_treatment`, `geocenter`. The literature tables the
window analysis joins to are `datasets/bibliometric/paper_study_windows.csv`
(each paper's study window), `claim_signs.csv` (the direction each paper
claims) and `paper_footprints.csv` (the cohort polygons of each paper's own
study area, with the coding note; `literature_footprint_basins.txt` lists the
37 polygons). Two curated tables precede the sign-flip risk and carry one
written reason per row: `record_adjudication.csv` (duplicate article records
dropped, claims that are not window trends excluded from the risk) and
`method_adjudication.csv` (stated processing choices the recipe space cannot
represent, which are compared against the full ensemble instead, and
family-level restrictions).

```r
library(data.table)
trends <- fread("datasets/output/tier_a_full_73/jasechko_per_recipe_trends.csv.gz")
trends[, .(spans_zero = min(trend_cm_yr) < 0 & max(trend_cm_yr) > 0),
       by = Study_area][spans_zero == TRUE]
```

### Control-Earth specifications and answer keys

| File | Contents |
|---|---|
| `control_earth/datasets/truth/truth_spec_<planet>.json` | the planted per-aquifer trends, seasonal cycle, noise seed and taper of each planet |
| `control_earth/datasets/truth/truth_manifest_<planet>.json` | SHA-256 chain of the specification and the fields built from it |
| `control_earth/datasets/truth/truth_aquifer_trends_<planet>.csv` | the answer key: true total-storage trend per aquifer that every score is computed against |
| `control_earth/datasets/truth/truth_gw_aux_trends_<planet>.csv` | planted groundwater-layer trend and the two auxiliary trends, for the subtraction test |
| `control_earth/datasets/truth/params_iter3.json` | the accepted forward-operator noise parameters (iteration 3 of the calibration) |

### Control-Earth summary tables

| Table | Contents |
|---|---|
| `control_earth/datasets/output/calibration_report.csv` | the reproduce-reality gate values G1–G4 on the identical 73-aquifer support (all pass); `calibration_report_2026-07-14_real276.csv` preserves the first report, whose real side was the 276-aquifer ensemble |
| `control_earth/datasets/output/rolling_null/rolling_attribution_planets.csv`, `rolling_attribution_replicates.csv`, `rolling_truth_window_trends.csv` | every rolling window of five planets scored against its truth window (`07c_rolling_attribution.py`); the noise replicate pair splits wrong-sign windows shared by both members from those in one |
| `control_earth/datasets/output/parametric_floor_v3.csv` | the 50/90/95% recovery crossings on the sign-imbalanced graded planet (69 depleting, 4 recovering), with bootstrap CIs |
| `control_earth/datasets/output/score_graded_floor_v3.csv` | per-aquifer sign recovery, recovered magnitudes and bracketing on the sign-imbalanced graded planet |
| `control_earth/datasets/output/score_graded_balanced_*.csv` | the same, on the three balanced planets, their sign-reversed counterpart and a noise replicate |
| `control_earth/datasets/output/balanced_layouts_{summary,bins,floors}.csv` | per planet: recovery by magnitude band (the primary presentation), fitted crossings, attenuation, class-neutral accuracy and the data-space boundary |
| `control_earth/datasets/output/boundary_balanced.csv` | the data-space boundary per balanced planet, with the graded-planet protocol check |
| `control_earth/datasets/output/sign_interleaving.csv` | aquifers with an opposite-signed neighbour within 500 km, per planet (`19_sign_interleaving.R`, from the specifications alone) |
| `control_earth/datasets/output/aux_subtraction_balanced.csv` | the groundwater subtraction on the balanced planets under a matched auxiliary, a mismatched one and none |
| `control_earth/datasets/output/score_gws_crossed_era5_v1.csv` | the groundwater subtraction scored on the crossed (ERA5-Land-background) planet, whose planted groundwater layer is depleting on all 73 aquifers |
| `control_earth/datasets/output/truth_gw_aux_trends_crossed_era5_v1.csv` | the crossed planet's answer key: planted groundwater-layer trend and the two auxiliary trends per aquifer |

## Provenance

`PROVENANCE.md` lists every data file in the repository with the script or
notebook chunk that writes it, the reproduction track it needs and the
notebook chunks that read it. Regenerate it from the repository root with

```bash
python code/make_provenance.py
```

The mapping from each published number to the tables runs through the
notebook: the `Numbers in the text` section asserts every value with
`stopifnot()`, and the *Code* boxes of the Supplementary Materials name the
chunk and script behind each supplementary item.

## Not distributed

| Item | Size | How to obtain |
|---|---|---|
| Recipe cubes `grace_ewh_global_<centre>_appendix.nc` (three centres) | about 110 GiB each | track 3, or from the corresponding author |
| Per-aquifer monthly series `datasets/output/tier_a_windows/monthly_<centre>.npz` | about 170 MB each | `tier_a_window_trends.py --emit-monthly` on the cubes (track 2) |
| DDK kernels | 9.8 MB each, five files | `python -m grace_pipeline.download_ddk` |
| Synthetic Level-2 records, recovered cubes, truth fields and per-recipe tables of the control Earths | about 113 GB per centre and planet | track 4 (about 14.5 h per planet) |
| Raw GRACE/GRACE-FO Level-2 archives, mascon products, GLDAS, ERA5-Land, SLR and GIA series | tens of GB | the `download_*.py` scripts and the sources named in `config.py` and `io_aux.py` |
| Well records of Jasechko et al. (2024) beyond the shipped cohort join | | their deposit (doi:10.5281/zenodo.10003697) |
| Policy-document PDFs and the bibliometric retrieval and extraction scripts | | publisher copyright; `policy_context_retrieval.csv` carries every source URL, the extraction procedure is declared in the Supplementary Materials |

## Citation

If you use this material, please cite the paper and, for the control-Earth
design record:

> Puy, A. (2026). *Control Earth: frozen design, pre-registration and decision
> register of a synthetic-truth test of the GRACE groundwater processing
> chain*. Zenodo. doi: [10.5281/zenodo.21454588](https://doi.org/10.5281/zenodo.21454588).

and, for the amendment that specifies the sign-balanced planets:

> Puy, A. (2026). *Control Earth: analysis amendment 01, sign-balanced planets*.
> Zenodo. doi: [10.5281/zenodo.21819495](https://doi.org/10.5281/zenodo.21819495).

<!-- Repo citation: add the Zenodo release DOI when minted. -->

## Licenses and acknowledgements

This project is released under the **MIT License** (see `LICENSE`).

- `datasets/jasechko_2024/` derives from Jasechko et al. (2024), *Rapid
  groundwater decline and some cases of recovery in aquifers globally*
  (data: CC BY 4.0, https://doi.org/10.5281/zenodo.10003697). Their public
  deposit covers the ~59% of well records they had permission to repost.
- The policy-document corpus was discovered through Altmetric
  (https://www.altmetric.com); the shipped tables keep only factual, curated
  fields. Policy PDFs are not redistributed (publisher copyright);
  `datasets/bibliometric/policy_context_retrieval.csv` carries the source URL
  of every document.
- GRACE/GRACE-FO Level-2 products: CSR, JPL and GFZ via NASA PO.DAAC and GFZ
  ISDC. Mascon comparators: CSR (UT Austin), JPL (PO.DAAC,
  doi:10.5067/TEMSC-3JC634) and GSFC.
- The bibliometric tables carry provenance columns from the LLM-assisted
  corpus extraction (model, verbatim evidence quote, rationale). The
  procedure — extraction, an independent second LLM reader and human
  adjudication of every disputed record — is declared in the paper's
  Supplementary Materials (bibliometric methods).
