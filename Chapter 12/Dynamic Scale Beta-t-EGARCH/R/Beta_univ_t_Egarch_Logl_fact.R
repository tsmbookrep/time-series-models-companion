# Beta-t-EGARCH likelihood together with filtered diagnostic quantities.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
source(file.path(script_dir, "Beta_univ_t_Egarch_Logl.R"))

Beta_univ_t_Egarch_Logl_fact <- function(par, y, inf, I) {
  y <- as.numeric(y)
  p <- unpack_parameters(par, I)
  filtered <- beta_t_egarch_filter(par, y, inf, I)
  fit <- y - filtered$res  # Legacy MATLAB output, kept for compatibility.
  log_density_constant <- lgamma((p$df + 1) / 2) - lgamma(p$df / 2) - 0.5 * log(pi * p$df)
  logl <- length(y) * log_density_constant - sum(filtered$lam) - (p$df + 1) / 2 * sum(log1p(filtered$res^2 / p$df))
  list(Logl = -logl, lam = filtered$lam, res = filtered$res, u = filtered$score, Beta = filtered$beta, fit = fit)
}
