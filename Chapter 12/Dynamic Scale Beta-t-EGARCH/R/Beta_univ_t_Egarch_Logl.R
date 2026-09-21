# Negative log-likelihood for the Beta-t-EGARCH dynamic-scale model.

unpack_parameters <- function(par, I) {
  par <- as.numeric(par)
  kappa_2 <- 0
  has_leverage <- FALSE
  if (I == 1) {
    if (!(length(par) %in% c(4, 5))) stop("The integrated specification requires 4 or 5 parameters.")
    mu <- par[1]; omega <- par[2]; kappa <- par[3]; log_df <- par[4]; phi <- NA_real_
    if (length(par) == 5) { kappa_2 <- par[5]; has_leverage <- TRUE }
  } else {
    if (!(length(par) %in% c(5, 6))) stop("The stationary specification requires 5 or 6 parameters.")
    mu <- par[1]; omega <- par[2]; phi <- par[3]; kappa <- par[4]; log_df <- par[5]
    if (length(par) == 6) { kappa_2 <- par[6]; has_leverage <- TRUE }
  }
  list(mu = mu, omega = omega, phi = phi, kappa = kappa, df = exp(log_df), kappa_2 = kappa_2, has_leverage = has_leverage)
}

beta_t_egarch_filter <- function(par, y, inf, I) {
  y <- as.numeric(y)
  p <- unpack_parameters(par, I)
  T <- length(y)
  lam <- res <- score <- beta <- numeric(T)
  information <- if (inf == 1) 2 * p$df / (p$df + 3) else 1
  next_lam <- p$omega
  for (i in seq_len(T)) {
    lam[i] <- next_lam
    res[i] <- (y[i] - p$mu) * exp(-lam[i])
    beta[i] <- (res[i]^2 / p$df) / (1 + res[i]^2 / p$df)
    score[i] <- ((p$df + 1) * beta[i] - 1) / information
    leverage <- if (p$has_leverage) p$kappa_2 * sign(p$mu - y[i]) * (score[i] + 1) else 0
    next_lam <- if (I == 1) lam[i] + p$kappa * score[i] + leverage else p$omega * (1 - p$phi) + p$phi * lam[i] + p$kappa * score[i] + leverage
  }
  list(lam = lam, res = res, score = score, beta = beta)
}

Beta_univ_t_Egarch_Logl <- function(par, y, inf, I) {
  p <- unpack_parameters(par, I)
  filtered <- beta_t_egarch_filter(par, y, inf, I)
  log_density_constant <- lgamma((p$df + 1) / 2) - lgamma(p$df / 2) - 0.5 * log(pi * p$df)
  logl <- length(y) * log_density_constant - sum(filtered$lam) - (p$df + 1) / 2 * sum(log1p(filtered$res^2 / p$df))
  -logl
}
