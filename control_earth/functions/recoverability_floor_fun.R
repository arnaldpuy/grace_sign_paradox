recoverability_floor_fun <- function(tau_true, p_sign, n_boot = 2000L,
                                      seed = 71L, span = 0.75,
                                      grid_n = 400L) {

  # RECOVERABILITY FLOOR s* (METHODOLOGY 6.2; PREREGISTRATION P2) #############

  # The |tau_true| at which the loess of sign-recovery probability P_sign on
  # log10(|tau_true|) crosses 0.95 upward, with a nonparametric bootstrap
  # over aquifers. Deliberately NOT a re-implementation: it sources and calls
  # the paper's own detectability_boundary_fun() (SA095
  # branch, identical loess span, grid, bootstrap size and seed), so the
  # truth-based synthetic floor and the data-based empirical boundary are
  # estimated by the SAME code path (METHODOLOGY 11 design invariant; D020).
  # The only difference is the meaning of the x-axis: |tau_true| (known truth
  # magnitude) here vs |recipe-median trend| (data) in the paper.

  repo_fun <- file.path(here::here(),
                        "functions/detectability_boundary_fun.R")
  source(repo_fun, local = TRUE)
  detectability_boundary_fun(x = abs(tau_true), y = p_sign, type = "SA095",
                             n_boot = n_boot, seed = seed, span = span,
                             grid_n = grid_n)
}
