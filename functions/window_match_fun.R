##  window_match_fun.R
##  ---------------------------------------------------------------------------
## Score every paper on the recipe ensemble sliced to ITS OWN study window
## (mission-era scoring mis-scored window-specific claims, e.g. Leblanc 2009).
## Windows: paper_study_windows.csv; tables: jasechko_per_recipe_trends__<win>.csv;
## the window key is "<win_start>_<win_end>".
##  ---------------------------------------------------------------------------

# null-coalescing helper (base R has none)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Assign each cohort paper its EXACT study window.
#'
#' @param bib_in         data.table with a `paper_id` column, restricted to the
#'                        cohort (in_cohort == TRUE & paper_category == "own-SHM").
#' @param study_windows  data.table read from paper_study_windows.csv; must carry
#'                        `paper_id`, `win_start`, `win_end` (YYYY-MM).
#' @return  bib_in plus `win_start`, `win_end`, and `data_window`
#'           (= "<win_start>_<win_end>"). A cohort paper with no window row
#'           (should not occur) defaults to the full record and is warned about.
classify_paper_windows <- function(bib_in, study_windows) {
  stopifnot(data.table::is.data.table(bib_in),
            "paper_id" %in% names(bib_in),
            data.table::is.data.table(study_windows),
            all(c("paper_id", "win_start", "win_end") %in% names(study_windows)))

  sw <- unique(study_windows[, .(paper_id, win_start, win_end)])
  sw[, data_window := sprintf("%s_%s", win_start, win_end)]

  out <- merge(bib_in, sw, by = "paper_id", all.x = TRUE)
  missing <- out[is.na(data_window), unique(paper_id)]
  if (length(missing)) {
    warning("classify_paper_windows: no study window for paper_id(s) ",
            paste(missing, collapse = ", "),
            " -- defaulting to full record 2002-04..2025-12")
    out[is.na(data_window),
        `:=`(win_start = "2002-04", win_end = "2025-12",
             data_window = "2002-04_2025-12")]
  }
  out[]
}


#' Load per-window per-recipe trend tables produced by
#' tier_a_window_trends.py --windows.
#'
#' @param windows_dir  directory of jasechko_per_recipe_trends__<start>_<end>.csv.
#' @param recipe_grid  data.table mapping (centre + 7 axes) -> recipe_idx
#'                      (the Rmd's recipe_grid). Joined so the returned tables
#'                      carry the same recipe_idx the paper matcher uses.
#' @param windows      optional character vector of "<start>_<end>" keys to load
#'                      (default: every CSV in windows_dir). Pass the set of
#'                      data_window values actually present in the cohort to
#'                      avoid loading unused windows.
#' @param basins       optional character vector of Study_area names; rows are
#'                      filtered to these on load (the literature analysis only
#'                      touches the candidate basins, so this keeps the in-memory
#'                      tables small).
#' @return  named list (key = "<start>_<end>") of data.tables keyed on
#'           (recipe_idx, Study_area) with column `trend = trend_cm_yr`.
load_per_recipe_trends_by_window <- function(windows_dir, recipe_grid,
                                             windows = NULL, basins = NULL) {
  stopifnot(data.table::is.data.table(recipe_grid))
  files <- list.files(
    windows_dir,
    pattern = "^jasechko_per_recipe_trends__.*\\.csv(\\.gz)?$", full.names = TRUE)
  res <- list()
  for (f in files) {
    key <- sub("\\.csv(\\.gz)?$", "",
               sub("^jasechko_per_recipe_trends__", "", basename(f)))
    if (!is.null(windows) && !(key %in% windows)) next
    pr <- if (grepl("\\.gz$", f)) {
      data.table::fread(cmd = paste("gunzip -c", shQuote(f)))
    } else data.table::fread(f)
    if (!is.null(basins)) pr <- pr[Study_area %chin% basins]
    keyed <- merge(pr, recipe_grid,
                   by = c("centre", "truncation", "filter", "gia_model",
                          "c20_treatment", "c30_treatment", "geocenter"))
    res[[key]] <- keyed[, .(recipe_idx, Study_area, trend = trend_cm_yr)]
  }
  missing <- setdiff(windows %||% character(0), names(res))
  if (length(missing)) {
    warning("load_per_recipe_trends_by_window: no trend table for window(s) ",
            paste(missing, collapse = ", "),
            " -- run tier_a_window_trends.py --windows")
  }
  res
}


#' Window-aware sign-flip risk per (paper, basin).
#'
#' For each (paper, basin) pair, join the paper's compatible recipes to the
#' per-recipe trend table for THAT paper's exact window and compute the
#' depleting/recharging fractions. Two statistics follow. `sign_flip_risk` is
#' the fraction of compatible recipes whose sign is OPPOSITE to the paper's
#' claim: frac_pos for a depletion claim, frac_neg for a recovery claim, NA
#' when the paper makes no directional claim (mixed, stable, no-sign).
#' `minority_share` is the smaller of the two fractions; it ignores the claim
#' and is bounded by 0.5. The caller derives the pair's sign-paradox status
#' from frac_pos & frac_neg (spans zero iff both are > 0) on the paper's own
#' window.
#'
#' @param papers          matches data.table with paper_id, recipe_idx,
#'                          data_window, n_matched and claim_sign
#'                          ("depleting", "recovering" or other); without a
#'                          claim_sign column every sign_flip_risk is NA.
#' @param basins          character vector of Study_area values to test.
#' @param trends_by_win   list from load_per_recipe_trends_by_window().
#' @return  data.table: paper_id, candidate_basin, data_window, claim_sign,
#'           frac_pos, frac_neg, sign_flip_risk, minority_share, spans_zero,
#'           n_matched.
compute_sign_flip_window_aware <- function(papers, basins, trends_by_win) {
  if (!"claim_sign" %in% names(papers)) {
    papers <- data.table::copy(papers)[, claim_sign := NA_character_]
  }
  data.table::rbindlist(lapply(basins, function(b) {
    data.table::rbindlist(lapply(names(trends_by_win), function(win) {
      pw <- papers[data_window == win]
      if (nrow(pw) == 0L) return(NULL)
      bt <- trends_by_win[[win]][Study_area == b, .(recipe_idx, trend)]
      if (nrow(bt) == 0L) return(NULL)
      pm <- merge(pw, bt, by = "recipe_idx", allow.cartesian = TRUE)
      r <- pm[, .(frac_pos = mean(trend > 0, na.rm = TRUE),
                  frac_neg = mean(trend < 0, na.rm = TRUE),
                  n_matched = unique(n_matched),
                  claim_sign = unique(claim_sign)[1L]),
              by = paper_id]
      r[, `:=`(candidate_basin = b,
               data_window     = win,
               sign_flip_risk  = data.table::fcase(
                 !is.na(claim_sign) & claim_sign == "depleting",  frac_pos,
                 !is.na(claim_sign) & claim_sign == "recovering", frac_neg,
                 default = NA_real_),
               minority_share  = pmin(frac_pos, frac_neg),
               spans_zero      = (frac_pos > 0 & frac_neg > 0))]
      r
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}


#' Append one row per footprint to every per-window trend table.
#'
#' A paper's footprint is the set of cohort polygons that cover its own study
#' area (paper_footprints.csv). Its per-recipe trend is the cell-count-weighted
#' mean of the member polygons' trends, which equals the trend of the
#' cell-mean series over the union because every polygon shares the recipe's
#' time design. Polygons with a missing trend are dropped from the weights.
#'
#' @param trends_by_win  list from load_per_recipe_trends_by_window().
#' @param fp_members     data.table (fp_key, Study_area), one row per member.
#' @param cell_w         data.table (Study_area, n_cells).
#' @return  the same list with footprint rows appended (Study_area = fp_key).
add_footprint_rows_fun <- function(trends_by_win, fp_members, cell_w) {
  fm <- merge(fp_members, cell_w, by = "Study_area")
  stopifnot(nrow(fm) == nrow(fp_members))
  lapply(trends_by_win, function(tw) {
    x <- merge(tw, fm, by = "Study_area", allow.cartesian = TRUE)
    fp <- x[is.finite(trend), .(Study_area = fp_key[1L],
                                trend = sum(n_cells * trend) / sum(n_cells)),
            by = .(fp_key, recipe_idx)][, fp_key := NULL]
    data.table::rbindlist(list(tw, fp), use.names = TRUE)
  })
}


#' Append one row per footprint to the per-window auxiliary trend table.
#'
#' Same weights as add_footprint_rows_fun(), per (data_window, centre,
#' truncation), so that trend(GWS) = trend(TWS) - trend(AUX) holds on the
#' footprint exactly as on a single polygon.
add_footprint_aux_fun <- function(aux_win, fp_members, cell_w) {
  fm <- merge(fp_members, cell_w, by = "Study_area")
  x <- merge(aux_win, fm, by = "Study_area", allow.cartesian = TRUE)
  fp <- x[is.finite(aux_trend_cmyr), .(
    Study_area = fp_key[1L],
    aux_trend_cmyr = sum(n_cells * aux_trend_cmyr) / sum(n_cells),
    n_months = min(n_months)),
    by = .(fp_key, data_window, centre, truncation)][, fp_key := NULL]
  data.table::rbindlist(list(aux_win, fp), use.names = TRUE, fill = TRUE)
}


#' Window-aware sign-flip risk of every paper on ITS OWN footprint.
#'
#' Evaluates each region's papers against the region's footprint keys with
#' compute_sign_flip_window_aware() and keeps, for every paper, only the row
#' of its own footprint. Papers without a footprint (fp_key NA) are dropped.
#'
#' @param m              matches data.table with region_label and fp_key.
#' @param fp_by_region   list, region -> data.table with column fp_key.
#' @param tbw            per-window trend tables carrying the footprint rows.
own_footprint_flip_fun <- function(m, fp_by_region, tbw) {
  data.table::rbindlist(lapply(names(fp_by_region), function(reg) {
    pir <- m[region_label == reg & !is.na(fp_key)]
    if (nrow(pir) == 0L) return(NULL)
    sf <- compute_sign_flip_window_aware(pir, unique(fp_by_region[[reg]]$fp_key),
                                         tbw)
    sf[, region_label := reg]
    sf <- merge(sf, unique(pir[, .(paper_id, fp_key)]), by = "paper_id")
    sf <- sf[candidate_basin == fp_key]
    sf[, fp_key := NULL][]
  }), use.names = TRUE, fill = TRUE)
}
