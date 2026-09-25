detectability_boundary_fun <- function(x, y, type = c("R1", "SA095"),
                                        n_boot = 2000L, seed = 71L,
                                        span = 0.75, grid_n = 400L) {

  # BOOTSTRAPPED DETECTABILITY-BOUNDARY CROSSING ##############################

  # Boundary estimators: "R1" = loess of log10(R) crossing 0 downward; "SA095" =
  # loess of sign-agreement crossing 0.95 upward; aquifer bootstrap gives the 95% CI.
  # Reproduces r1_crossing_loess / sa_crossing_loess in the main notebook.

  type <- match.arg(type)

  crossing <- function(xx, yy) {
    if (type == "R1") {
      ok <- is.finite(xx) & is.finite(yy) & xx > 0 & yy > 0
      if (sum(ok) < 10L) return(NA_real_)
      lx <- log10(xx[ok]); ly <- log10(yy[ok]); thr <- 0
    } else {
      ok <- is.finite(xx) & is.finite(yy) & xx > 0 & yy >= 0.5 & yy <= 1
      if (sum(ok) < 10L) return(NA_real_)
      lx <- log10(xx[ok]); ly <- yy[ok]; thr <- 0.95
    }
    fit <- tryCatch(
      stats::loess(ly ~ lx, span = span,
                   control = stats::loess.control(surface = "direct")),
      error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    gx <- seq(min(lx), max(lx), length.out = grid_n)
    gy <- suppressWarnings(stats::predict(fit, newdata = data.frame(lx = gx)))
    s <- if (type == "R1") {
      which(gy[-1L] < thr & gy[-length(gy)] >= thr)
    } else {
      which(gy[-1L] >= thr & gy[-length(gy)] < thr)
    }
    if (length(s) == 0L) return(NA_real_)
    i <- s[1L]
    10^(gx[i] + (thr - gy[i]) / (gy[i + 1L] - gy[i]) * (gx[i + 1L] - gx[i]))
  }

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  obs <- crossing(x, y)
  set.seed(seed)
  bs <- replicate(n_boot, {
    idx <- sample.int(length(x), replace = TRUE)
    crossing(x[idx], y[idx])
  })
  bs <- bs[is.finite(bs)]
  ci <- if (length(bs) > 1L) stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE)
        else c(NA_real_, NA_real_)
  list(est = obs, lo = ci[[1L]], hi = ci[[2L]],
       n = length(x), n_boot_ok = length(bs))
}
