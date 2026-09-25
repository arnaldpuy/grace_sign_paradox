# rolling_null — rolling-window null on the synthetic planets

Truth on every planet is LINEAR (trend + seasonal), so window-to-window sign
changes here are pure noise + recipe arithmetic — the null distribution for
the real-Earth rolling-window and mission-window results.

Files:

- `rolling_window_signs_<planet>.csv` — one row per aquifer x 10-yr window
  (starts 2002-2016, Apr-Mar, MIN_MONTHS = 84, mirror of
  `tier_a_rolling_windows.py`): n_recipes, frac_pos, med_trend_cmyr, spans.
- `rolling_null_summary.csv` — per-aquifer classification counts (win_flip;
  both / win_only / recipe_only / stable; median share of windows spanning
  zero), real Earth (recomputed: 68/73; 55/13/1/4; 0.60) vs the pilot,
  mixed-sign and graded planets.
- `window_paradox_summary.csv` — sign-paradox share per mission window
  (real: 77/74/26%).
- `window_dispersion_summary.csv` — per-window median ensemble IQR and median
  |ensemble median|: separates spread effects from signal effects.

Headline: the linear-truth null produces 62-72 of 73 window flips (real: 68)
— window instability is not by itself evidence of non-monotonic storage; the
real 77% -> 26% paradox drop on the GRACE-FO window is a signal effect (real
w_fo median |trend| 1.47 cm/yr vs 0.31-0.45 under linear truth). Caveat: the
synthetic 10-yr/w_grace window spread is about 1.5x the real one (the
calibration gate constrains full-record spread only).
