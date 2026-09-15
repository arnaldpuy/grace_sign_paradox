<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX) — minted at release -->

# One satellite record yields opposite groundwater trends

[Arnald Puy](https://www.arnaldpuy.com), Samuel Flinders, Seth N. Linga, Carmen Aguiló-Rivera, Ana-Maria Crețu

R and Python code accompanying the paper *One satellite record yields opposite
groundwater trends*, whose abstract is the following:

## Abstract

*The sign of satellite-gravity groundwater trends (depletion or recovery)
anchors global water security. Whether this sign is set by the aquifers or by
the processing that produces it remains unknown. Here we evaluate 1,920
pre-processing recipes across the full GRACE record (2002–2025) for 73 of the
world's largest aquifer systems, and show that these choices can change whether
three in four aquifers are reported as gaining or losing water. A detectability
boundary near 0.35 cm yr⁻¹ equivalent water height marks the regime where the
sign of a GRACE trend is defined by pre-processing choices rather than by the
aquifer. Eight control Earths with known trends reproduce this boundary
wherever it can be estimated and show that there is no universal threshold
guaranteeing correct recovery of the groundwater trend. Such underdetermination
matters for policy: on contested aquifers, published claims concentrate on
depletion relative to the admissible recipe ensemble and these verdicts enter
policy reports as single estimates. Our results define where satellite
gravimetry can adjudicate groundwater depletion, where it cannot and how
uncertainty lost in pre-processing becomes certainty in the scientific record.*

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
├── README.md, LICENSE
├── code/
│   ├── code_main_analysis.Rmd       # THE notebook: every figure, table and number
│   ├── code_main_analysis.pdf       # its knitted certificate (assertions passed)
│   └── grace_pipeline/              # Python: Level-2 archives -> 1,920-recipe cubes -> tables
├── functions/                       # 14 R helpers sourced by the notebook
├── datasets/
│   ├── jasechko_2024/               # aquifer cohorts, climate zones, cell masks
│   ├── bibliometric/                # curated literature/policy corpus tables
│   └── output/                      # the pre-reduced analysis tables the notebook reads
├── control_earth/
│   └── datasets/output/             # control-Earth summary tables (gate, crossings, scores)
└── figures/                         # the choice-space diagram (TikZ source + PDF), the notebook's only static figure input
```

> Note: the three ~110 GiB per-centre recipe cubes, the synthetic Level-2
> record and the control-Earth planet fields are archived on the records
> listed in the paper's Data availability statement — they are not in this
> repository. The frozen design, pre-registration and decision register of
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

Three tracks, by increasing depth:

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

2. **Intermediate (hours, ~330 GB disk).** Download the three per-centre
   recipe cubes from the records in the Data availability statement, then
   re-run the per-aquifer reductions (`code/grace_pipeline/tier_a_per_aquifer.py`
   and companions) to rebuild `datasets/output/` from the cubes.

3. **Full (days).** Rebuild the cubes from the public Level-2 archives:
   `code/grace_pipeline/download_grace_shm.py` fetches the CSR/JPL/GFZ
   releases and correction series; `code/grace_pipeline/pipeline_global.py`
   synthesises the 1,920-recipe cubes. DDK coefficients are vendored from
   [strawpants/GRACE-filter](https://github.com/strawpants/GRACE-filter)
   (MIT); auxiliary sources (GLDAS, ERA5-Land, mascons, SLR C20/C30,
   degree-1, GIA models) are documented in the pipeline scripts.

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
| `datasets/output/tier_a_windows/jasechko_per_recipe_trends__*.csv` | the same ensembles sliced to each published study's window |

Columns: `centre`, `aquifer_idx`, `Study_area`, `recipe_idx_per_centre`,
`trend_cm_yr`, `truncation`, `filter`, `gia_model`, `c20_treatment`,
`c30_treatment`, `geocenter`.

```r
library(data.table)
trends <- fread("datasets/output/tier_a_full_73/jasechko_per_recipe_trends.csv.gz")
trends[, .(spans_zero = min(trend_cm_yr) < 0 & max(trend_cm_yr) > 0),
       by = Study_area][spans_zero == TRUE]
```

### Control-Earth summary tables

| Table | Contents |
|---|---|
| `control_earth/datasets/output/calibration_report.csv` | the reproduce-reality gate values G1–G4 (all pass) |
| `control_earth/datasets/output/parametric_floor_v3.csv` | the 50/90/95% recovery crossings on the sign-imbalanced graded planet (69 depleting, 4 recovering), with bootstrap CIs |
| `control_earth/datasets/output/score_graded_floor_v3.csv` | per-aquifer sign recovery, recovered magnitudes and bracketing on the sign-imbalanced graded planet |
| `control_earth/datasets/output/score_graded_balanced_*.csv` | the same, on the three balanced planets, their sign-reversed counterpart and a noise replicate |
| `control_earth/datasets/output/balanced_layouts_{summary,bins,floors}.csv` | per planet: recovery by magnitude band (the primary presentation), fitted crossings, attenuation, class-neutral accuracy and the data-space boundary |
| `control_earth/datasets/output/boundary_balanced.csv` | the data-space boundary per balanced planet, with the graded-planet protocol check |
| `control_earth/datasets/output/sign_interleaving.csv` | aquifers with an opposite-signed neighbour within 500 km, per planet |
| `control_earth/datasets/output/aux_subtraction_balanced.csv` | the groundwater subtraction on the balanced planets under a matched auxiliary, a mismatched one and none |
| `control_earth/datasets/output/score_gws_crossed_era5_v1.csv` | the groundwater subtraction scored on the crossed (ERA5-Land-background) planet, whose planted groundwater layer is depleting on all 73 aquifers |
| `control_earth/datasets/output/truth_gw_aux_trends_crossed_era5_v1.csv` | the crossed planet's answer key: planted groundwater-layer trend and the two auxiliary trends per aquifer |

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
