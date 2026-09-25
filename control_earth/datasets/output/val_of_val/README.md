# val_of_val — validation-of-validation with synthetic wells

Monte-Carlo well synthesis at the real Jasechko well locations, scored
against known truth: what can well validation actually see?

Files:

- `truth_gw_aquifer_<planet>.csv` — per-aquifer mask-averaged
  groundwater-layer trend (the well-accessible truth) alongside the
  total-storage answer key tau_true.
- `well_noise_calibration.csv` — per-aquifer SD of real per-well
  depth-to-water OLS trends (2002-2025, >= 10 annual values); the cohort
  median (0.223 m/yr) is the per-well noise used in the synthesis.
- `val_of_val_summary.csv` — formatted MC summary (200 draws) per planet x
  specific yield: above/below-boundary agreement, false-validation rate,
  false-rejection rate, well-aggregate error rates.
- `val_of_val_mc_numeric.csv` — the same in numeric long format
  (mean, 2.5%, 97.5%).
- `val_of_val_aquifer_graded_floor_v3.csv` — per-aquifer verdicts for the
  headline draw (graded planet, Sy = 0.15, seed 71).

Headline: the real above/below-boundary agreement pattern (68%/51%)
reproduces on known truth (69%/54% headline draw); the check's failure mode
is false REJECTION — when GRACE and wells disagree, GRACE is right ~72% of
the time at realistic parameters, because the well aggregate itself carries
the wrong storage sign for ~1 in 5 aquifers — while false validation is rare
(~1%). Figure: `figures/fig_val_of_val.pdf`.
