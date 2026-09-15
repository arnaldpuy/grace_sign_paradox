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
    pattern = "^jasechko_per_recipe_trends__.*\\.csv$", full.names = TRUE)
  res <- list()
  for (f in files) {
    key <- sub("\\.csv$", "",
               sub("^jasechko_per_recipe_trends__", "", basename(f)))
    if (!is.null(windows) && !(key %in% windows)) next
    pr <- data.table::fread(f)
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
#' per-recipe trend table for THAT paper's exact window, then compute the
#' depleting/recharging fractions and the minority-sign (flip) risk. The
#' caller derives the pair's sign-paradox status from frac_pos & frac_neg
#' (spans zero iff both are > 0) on the paper's own window.
#'
#' @param papers          matches data.table with paper_id, recipe_idx,
#'                          data_window, n_matched.
#' @param basins          character vector of Study_area values to test.
#' @param trends_by_win   list from load_per_recipe_trends_by_window().
#' @return  data.table: paper_id, candidate_basin, data_window, frac_pos,
#'           frac_neg, sign_flip_risk, spans_zero, n_matched.
compute_sign_flip_window_aware <- function(papers, basins, trends_by_win) {
  data.table::rbindlist(lapply(basins, function(b) {
    data.table::rbindlist(lapply(names(trends_by_win), function(win) {
      pw <- papers[data_window == win]
      if (nrow(pw) == 0L) return(NULL)
      bt <- trends_by_win[[win]][Study_area == b, .(recipe_idx, trend)]
      if (nrow(bt) == 0L) return(NULL)
      pm <- merge(pw, bt, by = "recipe_idx", allow.cartesian = TRUE)
      r <- pm[, .(frac_pos = mean(trend > 0, na.rm = TRUE),
                  frac_neg = mean(trend < 0, na.rm = TRUE),
                  n_matched = unique(n_matched)),
              by = paper_id]
      r[, `:=`(candidate_basin = b,
               data_window     = win,
               sign_flip_risk  = pmin(frac_pos, frac_neg),
               spans_zero      = (frac_pos > 0 & frac_neg > 0))]
      r
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}
