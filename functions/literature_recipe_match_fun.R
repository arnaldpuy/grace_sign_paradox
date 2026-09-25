literature_recipe_match_fun <- function(papers_dt, recipe_grid,
                                        overrides = NULL) {

  # PAPER-METHOD -> RECIPE-SET MATCHER #########################################

  # Match each paper to compatible recipes. An axis constrains the match only
  # when the paper states a level that exists in the recipe space. Four cases
  # are distinguished (method_adjudication.csv carries the adjudicated ones):
  #   unstated      -> the axis does not constrain;
  #   in-space      -> the axis is restricted to the stated level(s);
  #   family stated -> the axis is restricted to the family's in-space levels
  #                    ("restrict" rows of the override table);
  #   out of space  -> the paper's method cannot be placed in this space
  #                    ("unsupported" rows); the paper is compared against the
  #                    FULL ensemble and must be reported separately, never as
  #                    a method-compatible result.
  # Unspecified axes match all levels;
  # multi-valued axes (adjudicated composites such as "multiple: 60 (CSR), 90
  # (GFZ)" or "P4M6 destriping + G300") match EVERY stated level that exists in
  # the recipe space -- levels outside the space do not constrain.
  # Stated C20/C30 -> "replaced"; a named degree-1 method -> "GravIS" (same lineage).
  # ------
  # papers_dt : data.table with columns
  #     paper_id, centre, filter, gia_model, truncation, c20, c30, geocenter
  # ------
  # data.table with one row per (paper_id, recipe_idx) match. Empty rows
  # are dropped (papers whose specified method does not match any recipe).

  stopifnot(data.table::is.data.table(papers_dt))
  stopifnot(data.table::is.data.table(recipe_grid))
  needed_paper_cols <- c("paper_id", "centre", "filter", "gia_model",
                          "truncation", "c20", "c30", "geocenter")
  needed_grid_cols <- c("recipe_idx", "centre", "truncation", "filter",
                         "gia_model", "c20_treatment", "c30_treatment",
                         "geocenter")
  stopifnot(all(needed_paper_cols %in% names(papers_dt)))
  stopifnot(all(needed_grid_cols %in% names(recipe_grid)))

  # Sets of valid axis levels in the recipe space.
  valid_centre <- unique(recipe_grid$centre)
  valid_filter <- unique(recipe_grid$filter)
  valid_gia    <- unique(recipe_grid$gia_model)
  valid_trunc  <- unique(recipe_grid$truncation)

  # --- level-set extractors ---------------------------------------------------
  # Each returns the character vector of recipe-space levels stated in the
  # free-text value (possibly several), or NULL when none is recognisable
  # (NULL -> the axis does not constrain the match).
  set_centre <- function(x) {
    x0 <- toupper(as.character(x))
    hits <- valid_centre[vapply(valid_centre, function(l) grepl(l, x0, fixed = TRUE),
                                logical(1))]
    if (length(hits)) hits else NULL
  }
  set_trunc <- function(x) {
    nums <- unique(as.integer(unlist(regmatches(x, gregexpr("[0-9]+", as.character(x))))))
    hits <- intersect(nums, valid_trunc)
    if (length(hits)) hits else NULL
  }
  set_filter <- function(x) {
    x0 <- tolower(as.character(x))
    hits <- character(0)
    for (k in c("2", "3", "5", "7", "8"))
      if (grepl(paste0("ddk[ -]?", k), x0)) hits <- c(hits, paste0("DDK", k))
    for (n in c("300", "400", "500")) {
      has_n <- grepl(n, x0)
      if (has_n && grepl("swenson", x0)) hits <- c(hits, paste0("Swenson_G", n))
      else if (has_n && grepl(paste0("g", n, "|gauss"), x0)) hits <- c(hits, paste0("G", n))
    }
    hits <- intersect(unique(hits), valid_filter)
    if (length(hits)) hits else NULL
  }
  set_gia <- function(x) {
    x0 <- tolower(as.character(x))
    hits <- character(0)
    if (grepl("caron", x0))                    hits <- c(hits, "Caron_2018")
    if (grepl("ice.?6g", x0))                  hits <- c(hits, "Peltier_ICE6G_D")
    if (grepl("ice.?5g|geruo|vm2", x0))        hits <- c(hits, "Geruo_ICE5G_VM2")
    if (grepl("\\bnone\\b", x0))               hits <- c(hits, "none")
    hits <- intersect(unique(hits), valid_gia)
    if (length(hits)) hits else NULL
  }
  # C20 / C30: a stated treatment is always an SLR replacement (TN-14, TN-11,
  # Loomis 2019, ...); an explicit "original" is the un-replaced field. A value
  # asserting both (multi-solution composites) does not constrain.
  set_c2030 <- function(x) {
    x0 <- tolower(trimws(as.character(x)))
    if (x0 == "" || grepl("unclear|not stated|unspecified", x0)) return(NULL)
    # The negation is tested first: "not replaced" is the un-replaced field,
    # and it also contains the string "replac".
    if (grepl("not replac|un-?replac|without replac|no replacement|original", x0))
      return("original")
    if (grepl("replac", x0)) return("replaced")
    NULL
  }
  # Geocentre: any named degree-1 product / method means a correction was
  # applied -> the "GravIS" level; an explicit "none"/"not applied" -> "none".
  # Composites asserting both do not constrain.
  set_geo <- function(x) {
    x0 <- tolower(trimws(as.character(x)))
    if (x0 == "" || grepl("unclear|not stated|unspecified", x0)) return(NULL)
    # The negation is tested first: "not applied" also contains "applied".
    if (grepl("\\bnone\\b|not applied|not restored|no geo|omitted|ignor|n/a", x0))
      return("none")
    if (grepl("swenson|sun|tn-?13|gravis|landerer|replac|restor|geo-?cent|degree-?1|applied|slr", x0))
      return("GravIS")
    NULL
  }

  # Adjudicated overrides for this paper, by axis.
  ov_for <- function(pid, ax) {
    if (is.null(overrides)) return(NULL)
    r <- overrides[paper_id == pid & axis == ax]
    if (!nrow(r)) return(NULL)
    r
  }

  do_one_paper <- function(p) {
    grid <- recipe_grid
    # A paper with any out-of-space choice is compared against the full
    # ensemble; its stated axes are not used to narrow the set.
    if (!is.null(overrides) &&
        nrow(overrides[paper_id == p$paper_id & action == "unsupported"])) {
      return(data.table::data.table(paper_id = p$paper_id,
                                     recipe_idx = grid$recipe_idx,
                                     n_matched = nrow(grid)))
    }
    apply_axis <- function(grid, axis, col, parsed) {
      r <- ov_for(p$paper_id, axis)
      if (!is.null(r) && r$action[1L] == "restrict") {
        lv <- trimws(strsplit(r$levels[1L], ";", fixed = TRUE)[[1]])
        return(grid[get(col) %in% lv])
      }
      if (!is.null(parsed)) grid[get(col) %in% parsed] else grid
    }
    grid <- apply_axis(grid, "centre", "centre", set_centre(p$centre))
    grid <- apply_axis(grid, "filter", "filter", set_filter(p$filter))
    grid <- apply_axis(grid, "gia_model", "gia_model", set_gia(p$gia_model))
    grid <- apply_axis(grid, "truncation", "truncation", set_trunc(p$truncation))
    grid <- apply_axis(grid, "c20", "c20_treatment", set_c2030(p$c20))
    grid <- apply_axis(grid, "c30", "c30_treatment", set_c2030(p$c30))
    grid <- apply_axis(grid, "geocenter", "geocenter", set_geo(p$geocenter))
    if (nrow(grid) == 0L) return(NULL)
    data.table::data.table(paper_id  = p$paper_id,
                             recipe_idx = grid$recipe_idx,
                             n_matched  = nrow(grid))
  }

  out <- data.table::rbindlist(
    lapply(seq_len(nrow(papers_dt)),
            function(i) do_one_paper(papers_dt[i])),
    use.names = TRUE, fill = TRUE)
  out[]
}
