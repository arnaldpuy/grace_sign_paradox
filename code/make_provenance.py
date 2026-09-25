"""Build PROVENANCE.md: every data file in the repository -> the script or
notebook chunk that writes it, the reproduction track it needs, and the notebook
chunks that read it. Writers are found by scanning the scripts for the file name
next to a write call; the hand tables below cover the names built at run time.

Usage: python code/make_provenance.py   (from the repository root)
"""
import re, subprocess, pathlib, json, collections
ROOT = pathlib.Path(__file__).resolve().parents[1]
def _git(*args):
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True)
    return r.stdout.decode() if r.returncode == 0 else None
ls = _git("ls-files", "-z")
if ls is not None:
    ot = _git("ls-files", "-z", "--others", "--exclude-standard") or ""
    dl = set((_git("ls-files", "--deleted") or "").split("\n"))
    files = sorted(set(f for f in ls.split("\0") + ot.split("\0") if f and f not in dl))
else:  # no git checkout (e.g. an unpacked archive): walk the tree
    files = sorted(str(q.relative_to(ROOT)) for q in ROOT.rglob("*") if q.is_file()
                   and not any(part.startswith(".") or part.endswith(("_cache", "_files")) or part == "__pycache__" for part in q.relative_to(ROOT).parts))
data = [f for f in files if (f.startswith("datasets/") or f.startswith("control_earth/datasets/") or f.startswith("figures/")) and not f.endswith("README.md")]
scripts = [f for f in files if re.search(r"\.(py|R|sh)$", f) and f != "code/check_imports.py"]
rmd = (ROOT / "code/code_main_analysis.Rmd").read_text()
# chunks
chunks = []
for m in re.finditer(r"^```\{r ([^,}]+)[^\n]*\}\n(.*?)^```", rmd, re.S | re.M):
    chunks.append((m.group(1).strip(), m.group(2)))
planets = sorted({p.stem.replace("truth_spec_", "") for p in (ROOT / "control_earth/datasets/truth").glob("truth_spec_*.json")}, key=len, reverse=True)
planets += ["pilot_uniform_depletion", "crossed_era5"]
VERBS = re.compile(r"fwrite\(|write\.csv\(|to_csv\(|savez|np\.save|ggsave\(|save_plot\(|saveRDS\(|writeLines\(|file\.copy\(|to_netcdf|\.save\(|json\.dump|open\(.*[\"']w")
OUTASSIGN = re.compile(r"^\s*\w*(out|OUT)\w*\s*(=|<-)")
OUTHDR = re.compile(r"^\s*(Outputs?|Artefacts?|Writes?)\b.*:")
src_cache = {s: (ROOT / s).read_text(errors="ignore").split("\n") for s in scripts}
def writers_of(bn):
    pats = [bn]
    for p in planets:
        if p in bn:
            pats.append(bn.split(p)[0]); break
    hits = []
    for s, lines in src_cache.items():
        for pat in pats:
            if len(pat) < 6: continue
            named = [pat in l for l in lines]
            if not any(named): continue
            verb = [bool(VERBS.search(l)) for l in lines]
            near = [verb[i] or (i > 0 and verb[i-1]) or (i+1 < len(verb) and verb[i+1]) for i in range(len(lines))]
            hdr = [False]*len(lines)
            for i, l in enumerate(lines):
                if OUTHDR.match(l):
                    for j in range(i, min(len(lines), i+6)): hdr[j] = True
            if any(named[i] and (near[i] or OUTASSIGN.match(lines[i]) or hdr[i]) for i in range(len(lines))):
                hits.append(s + ("" if pat == bn else " (per planet)")); break
    return hits
def chunk_writers(bn):
    return [c for c, body in chunks if bn in body and re.search(r"fwrite\(|ggsave\(|save_plot\(", body) and re.search(re.escape(bn) + r"[\"')]", body)]
def chunk_readers(bn, rel):
    d = str(pathlib.Path(rel).parent)
    out = []
    for c, body in chunks:
        if bn in body: out.append(c)
        elif d.endswith("tier_a_windows") and "tier_a_windows" in body: out.append(c)
        elif d.endswith("section1") and '"section1"' in body: out.append(c)
    return out
TRACK = {  # script -> what it needs beyond the repository
 "code/grace_pipeline/tier_a_per_aquifer.py": "3 (recipe cubes)", "code/grace_pipeline/tier_a_window_trends.py": "3 (recipe cubes; monthly series)",
 "code/grace_pipeline/tier_a_leakage.py": "3 (recipe cubes)", "code/grace_pipeline/tier_a_leakage_ring.py": "3 (recipe cubes)",
 "code/grace_pipeline/tier_a_1deg_gain.py": "3 (recipe cubes; CLM4 gain grid)", "code/grace_pipeline/tier_a_stripe_index.py": "2 (monthly series)",
 "code/grace_pipeline/tier_a_trend_se.py": "2 (monthly series)", "code/grace_pipeline/tier_a_rolling_windows.py": "2 (monthly series)",
 "code/grace_pipeline/tier_a_common_months.py": "2 (monthly series)", "code/grace_pipeline/gws_aux_aquifer.py": "2 (monthly series; GLDAS, ERA5-Land)",
 "code/grace_pipeline/tier_a_mascon.py": "3 (mascon products; wells)", "code/grace_pipeline/tier_a_well_level.py": "3 (wells)",
 "code/grace_pipeline/well_aquifer_join.py": "3 (Jasechko et al. 2024 wells)", "code/grace_pipeline/climate_zones.py": "3 (climate raster)",
 "code/estimator_global.py": "2 (monthly series)", "code/surface_water_exposure.R": "1",
 "control_earth/code/03_calibration.R": "4 (synthetic per-recipe tables)", "control_earth/code/05_score.R": "4 (synthetic per-recipe tables)",
 "control_earth/code/07_rolling_windows_null.py": "4 (synthetic cubes)", "control_earth/code/07b_rolling_null.R": "4 (synthetic per-recipe tables)",
 "control_earth/code/07c_rolling_attribution.py": "4 (truth fields)", "control_earth/code/08a_synth_wells.py": "4 (synthetic cubes; wells)",
 "control_earth/code/08b_val_of_val.R": "4 (synthetic per-recipe tables)", "control_earth/code/09_parametric_floor.R": "1",
 "control_earth/code/10_recipe_skill.R": "4 (synthetic per-recipe tables)", "control_earth/code/11_repeated_cv.R": "4 (synthetic per-recipe tables)",
 "control_earth/code/12_score_crossed_gws.R": "4 (synthetic per-recipe tables)", "control_earth/code/14_floor_sign_balance.R": "1",
 "control_earth/code/15_boundary_balanced.R": "4 (synthetic per-recipe tables)", "control_earth/code/17_balanced_layouts.R": "1",
 "control_earth/code/18_aux_subtraction_balanced.R": "4 (synthetic per-recipe tables)", "control_earth/code/19_sign_interleaving.R": "1",
 "control_earth/code/06_boundary_overlay.R": "1",
}
INPUT_NOTE = {
 "datasets/bibliometric": "curated corpus table (retrieval and extraction pipeline not distributed; Supplementary Materials, bibliometric methods)",
 "datasets/jasechko_2024": "derived from Jasechko et al. (2024), doi:10.5281/zenodo.10003697",
 "control_earth/datasets/truth": "planet specification, answer key or manifest written by control_earth/code/forward/ (track 4)",
 "figures": "static TikZ source and its PDF",
}
HAND = {
 "datasets/bibliometric/claim_signs.csv": "LLM claim-sign extraction plus 24 hand-coded papers (Supplementary Materials)",
 "datasets/bibliometric/paper_footprints.csv": "hand-coded study footprints, one row per paper, with the coding note on coverage and claim type",
 "datasets/bibliometric/record_adjudication.csv": "curated: duplicate article records dropped and claims excluded from trend risk, one reason per row",
 "datasets/bibliometric/method_adjudication.csv": "curated: stated processing choices outside the recipe space (unsupported) and family-level restrictions, one row per record and axis",
 "datasets/bibliometric/literature_footprint_basins.txt": "the 37 cohort polygons of the coded footprints; input of tier_a_window_trends.py --windows",
 "datasets/bibliometric/literature_candidate_basins.txt": "superseded basin list of the earlier region-level matching",
 "control_earth/datasets/output/calibration_report_2026-07-14_real276.csv": "preserved first gate report (real ensemble on 276 aquifers), superseded by calibration_report.csv",
 "datasets/output/robustness/aquifer_size_bins.csv": "size-bin definition (hand-set input)",
}
PRIV = "analysis script of the private working repository (not distributed) run on the per-recipe tables of this repository"
OVERRIDE = {  # file (or basename prefix) -> ([writers], track)
 "jasechko_per_recipe_trends.csv.gz": (["`code/grace_pipeline/tier_a_per_aquifer.py` (one run per cohort and window, see `functions/resolve_tier_a_window.R`)"], "3 (recipe cubes)"),
 "jasechko_aquifer_metrics.csv": (["`code/grace_pipeline/tier_a_per_aquifer.py`"], "3 (recipe cubes)"),
 "jasechko_per_recipe_trends__": (["`code/grace_pipeline/tier_a_window_trends.py --windows` (one table per distinct study window, restricted to the 37 footprint polygons)"], "2 (monthly series)"),
 "gws_aux_trends_era5land.csv": (["`code/grace_pipeline/gws_aux_aquifer.py --source era5land`"], "2 (monthly series; ERA5-Land)"),
 "gws_aux_window_trends_gldas.csv": (["`code/grace_pipeline/gws_aux_aquifer.py --windows`"], "2 (monthly series; GLDAS)"),
 "mascon_aquifer_metrics.csv": (["`code/grace_pipeline/tier_a_mascon.py`"], "3 (mascon products; wells)"),
 "score_gws_crossed_era5_v1.csv": (["`control_earth/code/12_score_crossed_gws.R`"], "4 (synthetic per-recipe tables)"),
 "score_": (["`control_earth/code/05_score.R <planet>`"], "4 (synthetic per-recipe tables)"),
 "rolling_window_signs_": (["`control_earth/code/07_rolling_windows_null.py`"], "4 (synthetic cubes)"),
 "truth_gw_aquifer_": (["`control_earth/code/08a_synth_wells.py`"], "4 (synthetic cubes; wells)"),
 "grace_month_midpoints.csv": (["`code/grace_pipeline/tier_a_window_trends.py --emit-monthly` (the month midpoints of the monthly series)"], "2 (monthly series)"),
 "fig1_flip_timeseries.csv": (["`reduce_section1.R` (private working repository, not distributed): per-recipe time series of the anchor cell, read from the IGP recipe cube"], "3 (recipe cubes)"),
 "fig1_recipe_families.csv": (["`reduce_section1.R` (private working repository, not distributed): Ward clustering of the anchor-cell recipe series (`functions/cluster_recipe_families_fun.R`)"], "3 (recipe cubes)"),
 "fig_gradient_data.csv": (["`tier_a_gradient_figure.R`: " + PRIV], "1"),
 "gradient_data_276.csv": (["`tier_a_two_tier_stats.R`: " + PRIV], "1"),
 "sign_paradox_276.csv": (["`tier_a_two_tier_stats.R`: " + PRIV], "1"),
 "variance_decomp_276.csv": (["`tier_a_two_tier_stats.R`: " + PRIV + " (`functions/aquifer_variance_decomp_fun.R`)"], "1"),
 "size_check_per_aquifer.csv": (["`tier_a_size_check.R`: " + PRIV], "1"),
 "global_top_bottom_recipes.csv": (["`tier_a_well_skill.R`: " + PRIV + " (`functions/fit_recipes_to_wells_fun.R`)"], "1"),
 "per_recipe_global_skill.csv": (["`tier_a_well_skill.R`: " + PRIV + " (`functions/fit_recipes_to_wells_fun.R`)"], "1"),
 "leakage_summary.csv": (["`code_robustness_analysis.R`: " + PRIV + "; summarises `leakage_mask_sensitivity.csv`"], "1"),
 "mascon_vs_recipe_ensemble.csv": (["`code_robustness_analysis.R`: " + PRIV], "1"),
 "nested_recipe_space_summary.csv": (["`code_robustness_analysis.R`: " + PRIV + " (the attested-space rows are superseded by chunk `gws-strata`, see the folder README)"], "1"),
 "trend_estimator_anchorcell.csv": (["`code_robustness_analysis.R` (private working repository, not distributed): estimator comparison on the anchor-cell series of the IGP recipe cube; the cohort-wide comparison is `trend_estimator_global.csv.gz`"], "3 (recipe cubes)"),
 "trend_window_sensitivity.csv": (["`code_robustness_analysis.R`: " + PRIV], "1"),
 "well_grace_sign_validation.csv": (["`code_robustness_analysis.R`: " + PRIV + " (`functions/fit_recipes_to_wells_fun.R`)"], "1"),
}
def override_for(bn):
    if bn in OVERRIDE: return OVERRIDE[bn]
    for k, v in OVERRIDE.items():
        if k.endswith("_") and bn.startswith(k): return v
    return None
rows = collections.OrderedDict()
for f in data:
    bn = pathlib.Path(f).name
    ov = override_for(bn)
    w = [] if ov else writers_of(bn); cw = chunk_writers(bn)
    gen = list(ov[0]) if ov else []; track = ov[1] if ov else ""
    for s in w:
        base = s.split(" (")[0]
        gen.append(f"`{s}`"); track = TRACK.get(base, track or "")
    if cw: gen.append("notebook chunk " + ", ".join(f"`{c}`" for c in cw)); track = track or "1"
    if f in HAND: gen.append(HAND[f]); track = track or "input"
    if not gen:
        for k, v in INPUT_NOTE.items():
            if f.startswith(k): gen.append(v); track = "input"
    if not gen: gen.append("**no generator in the release**"); track = "?"
    rd = chunk_readers(bn, f)
    rows[f] = (gen, track, rd)
out = ["# Provenance manifest", "",
 "One row per data file in the repository: the script or notebook chunk that writes it,",
 "the reproduction track it needs (README, *How to run*: 1 = the notebook and the shipped",
 "tables; 2 = the per-aquifer monthly series, rebuilt from the recipe cubes; 3 = the",
 "recipe cubes or another raw input; 4 = the control-Earth rebuild; input = a curated or",
 "third-party input with no generator here), and the notebook chunks that read it. Every",
 "number in the paper is asserted from these files by the `assert-numbers` chunk of",
 "`code/code_main_analysis.Rmd`; the Supplementary Materials' *Code* boxes map each",
 "supplementary item to its chunks. Regenerate this file with the command in the",
 "README (*Provenance*).", ""]
cur = None
for f, (gen, track, rd) in rows.items():
    d = str(pathlib.Path(f).parent)
    if d != cur:
        cur = d; out += ["", f"## `{d}/`", "", "| File | Written by | Track | Read by notebook chunk |", "|---|---|---|---|"]
    out.append(f"| `{pathlib.Path(f).name}` | {'; '.join(gen)} | {track} | {', '.join(f'`{c}`' for c in rd) if rd else '-' } |")
(ROOT / "PROVENANCE.md").write_text("\n".join(out) + "\n")
print(len(rows), "files;", sum(1 for g,_,_ in rows.values() if "no generator" in g[0]), "without generator")
for f, (gen, track, rd) in rows.items():
    if "no generator" in gen[0] or track == "?": print("  ??", f)
