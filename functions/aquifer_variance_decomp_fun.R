aquifer_variance_decomp_fun <- function(trends_dt,
                                          axes = c("centre", "truncation",
                                                    "filter", "gia_model",
                                                    "c20_treatment",
                                                    "c30_treatment",
                                                    "geocenter")) {

  # PER-AQUIFER FACTORIAL VARIANCE DECOMPOSITION ##############################

  # Per aquifer: decompose across-recipe trend variance into 7 main effects,
  # 21 pairwise interactions and a residual, each as a fraction of total SS.
  # Returns Study_area | term | term_type | frac | ax1 | ax2.

  stopifnot(data.table::is.data.table(trends_dt))
  stopifnot(all(c("Study_area", "trend_cm_yr", axes)
                  %in% names(trends_dt)))

  pair_grid <- data.table::as.data.table(t(utils::combn(axes, 2L)))
  data.table::setnames(pair_grid, c("ax1", "ax2"))

  do_one <- function(aq_dt) {
    y <- aq_dt$trend_cm_yr
    grand <- mean(y, na.rm = TRUE)
    total_ss <- sum((y - grand)^2, na.rm = TRUE)
    if (!is.finite(total_ss) || total_ss == 0) return(NULL)

    main <- vapply(axes, function(ax) {
      m <- aq_dt[, .(mm = mean(trend_cm_yr, na.rm = TRUE), nn = .N),
                  by = c(ax)]
      sum(m$nn * (m$mm - grand)^2)
    }, numeric(1L))

    inter <- vapply(seq_len(nrow(pair_grid)), function(k) {
      a1 <- pair_grid$ax1[k]; a2 <- pair_grid$ax2[k]
      m <- aq_dt[, .(mm = mean(trend_cm_yr, na.rm = TRUE), nn = .N),
                  by = c(a1, a2)]
      full_ss <- sum(m$nn * (m$mm - grand)^2)
      full_ss - main[a1] - main[a2]
    }, numeric(1L))

    main_dt <- data.table::data.table(
      term = axes, term_type = "main",
      frac = main / total_ss,
      ax1 = axes, ax2 = NA_character_)
    inter_dt <- data.table::data.table(
      term = paste(pair_grid$ax1, pair_grid$ax2, sep = " x "),
      term_type = "interaction",
      frac = inter / total_ss,
      ax1 = pair_grid$ax1, ax2 = pair_grid$ax2)
    residual <- max(0, 1 - sum(main_dt$frac) - sum(inter_dt$frac))
    res_dt <- data.table::data.table(
      term = "higher-order", term_type = "residual",
      frac = residual,
      ax1 = NA_character_, ax2 = NA_character_)

    rbind(main_dt, inter_dt, res_dt)
  }

  out <- trends_dt[, do_one(.SD),
                     by = Study_area,
                     .SDcols = c("trend_cm_yr", axes)]
  out[]
}
