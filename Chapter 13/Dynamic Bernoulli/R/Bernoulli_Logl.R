# Negative log-likelihood and filter for the dynamic Bernoulli model.

Bernoulli_Logl <- function(par, y, tsmc = 0) {
  # This is the R counterpart to Bernoulli_Logl.m.
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1))) stop("The Bernoulli response must contain only zeros and ones.")

  omega <- par[1] # Unconstrained mean/intercept parameter.
  if (tsmc == 1) {
    kappa <- par[2] # Dynamic parameter in the simple Markov specification.
  } else {
    phi <- par[2]   # Dynamic AR parameter.
    kappa <- par[3] # Dynamic conditional-score parameter.
  }

  T <- length(y)
  mu <- numeric(T)
  res <- numeric(T)

  # Start the recursion at the sample mean, as in the MATLAB implementation.
  mu_next <- mean(y)
  for (i in seq_len(T)) {
    mu[i] <- mu_next
    res[i] <- y[i] - mu_next # For Bernoulli data, the score is the residual.
    if (tsmc == 1) {
      mu_next <- omega + kappa * y[i]
    } else {
      mu_next <- omega + phi * mu_next + kappa * res[i]
    }
  }

  # Infeasible probabilities receive a large objective value during optimization.
  if (any(mu <= 0) || any(mu >= 1) || any(!is.finite(mu))) {
    return(list(Logl = 1e12, mu = mu, res = res))
  }

  Logl <- -sum(y * log(mu) + (1 - y) * log(1 - mu))
  list(Logl = Logl, mu = mu, res = res)
}
