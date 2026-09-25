## resolve_tier_a_window.R -----------------------------------------------------
## Guard: standalone scripts default to the W-FULL headline window and print a
## banner; a non-headline window warns (or errors with require_headline = TRUE),
## so legacy POC numbers cannot be mistaken for manuscript values.
## Override with TIER_A_DIR, e.g. TIER_A_DIR=datasets/output/tier_a_276.
## ---------------------------------------------------------------------------

resolve_tier_a_window <- function(require_headline = FALSE) {
  default_dir <- here::here("datasets", "output", "tier_a_full")   # W-FULL headline
  out_dir     <- Sys.getenv("TIER_A_DIR", unset = default_dir)
  base        <- basename(out_dir)

  windows <- c(
    tier_a_full       = "W-FULL   | GRACE+GRACE-FO 2002-04..2025-12 (272 mo) | HEADLINE (manuscript)",
    tier_a_full_73    = "W-FULL   | GRACE+GRACE-FO 2002-04..2025-12 (272 mo) | HEADLINE, Tier-1 (73 aquifers)",
    tier_a_276        = "W-FO POC | GRACE-FO only 2018-06..2025-12 (89 mo)  | LEGACY -- NOT the manuscript window",
    tier_a            = "W-FO POC | GRACE-FO only 2018-06..2025-12 (89 mo)  | LEGACY -- NOT the manuscript window",
    tier_a_grace_only = "W-GRACE  | GRACE only 2002-04..2017-06            | paper-window comparator, not headline")
  win <- if (base %in% names(windows)) windows[[base]] else
           paste0("UNKNOWN window for directory '", base, "' -- verify before citing")
  is_headline <- base %in% c("tier_a_full", "tier_a_full_73")

  bar <- strrep("=", 86)
  cat(bar, "\n[Tier-A window]  ", win,
      "\n[Tier-A window]  out_dir = ", out_dir, "\n", bar, "\n", sep = "")

  if (!is_headline) {
    msg <- paste0(
      "NON-HEADLINE Tier-A window ('", base, "'). The manuscript reports the ",
      "W-FULL (tier_a_full) numbers; outputs from this window are NOT comparable ",
      "to the headline and must not be cited as manuscript values. Set ",
      "TIER_A_DIR=", default_dir, " for the headline window.")
    if (isTRUE(require_headline)) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE, immediate. = TRUE)
  }

  if (!file.exists(file.path(out_dir, "jasechko_per_recipe_trends.csv")) &&
      !file.exists(file.path(out_dir, "jasechko_per_recipe_trends.csv.gz")))
    stop("resolve_tier_a_window: no jasechko_per_recipe_trends.csv(.gz) in '",
         out_dir, "'.", call. = FALSE)

  list(out_dir = out_dir, base = base, window = win, is_headline = is_headline)
}
