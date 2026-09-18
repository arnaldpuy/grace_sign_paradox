# Control Earth experiment

Code, planet specifications and answer keys of the synthetic-truth experiment
(Methods, "The control Earth experiment with known truth"). The banked outputs
that the notebook reads and asserts are in `datasets/output/`.

## Layout

| Path | Contents |
|---|---|
| `code/forward/` | Truth builder, forward operator (water to Stokes coefficients, Wahr et al. 1998, via `gravity-toolkit`), recovery driver, planet-specification generators, launch scripts |
| `code/` | Calibration gate, scoring, and every analysis behind the tables in `datasets/output/` |
| `functions/` | R helpers sourced by the analysis scripts |
| `datasets/truth/` | One `truth_spec_*.json` (the planted trends and seed), `truth_manifest_*.json` (SHA-256 chain) and `truth_aquifer_trends_*.csv` (answer key) per planet; forward-operator parameters `params_iter*.json` |
| `datasets/output/` | Scores and analysis tables, one README per subdirectory |

## Paths and environment

Scripts locate the repository from their own position; nothing needs editing.
The recovery driver writes synthetic Level-2 files and cubes (about 113 GB per
centre) to `GRACE_CE_CUBE_DIR`, default `control_earth/datasets/synthetic_l2/`.
Launch scripts use the interpreter in `PYTHON`, default `python3`. The Python
dependencies are those of the main pipeline (`requirements.txt`); the recovery
imports `code/grace_pipeline` unchanged and redirects only its input paths.

## Not distributed

Truth NetCDFs, synthetic Level-2 files and recovered cubes are regenerated
from the specifications (about 14.5 hours per planet). Two inputs of the
well-validation test (`08a_synth_wells.py`), the per-well table and the annual
depth-to-water file, are derived from Jasechko et al. (2024) and are not
redistributed; the script's paths mark them.

## Decision references

Comments of the form `Dnnn` cite entries of the decision register deposited
with the pre-registration (doi:10.5281/zenodo.21454588, amended at
doi:10.5281/zenodo.21819495).
