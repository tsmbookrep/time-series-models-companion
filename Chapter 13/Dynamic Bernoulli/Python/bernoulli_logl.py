"""Negative log-likelihood and filter for the dynamic Bernoulli model."""

import numpy as np


def bernoulli_logl(par, y, tsmc=0):
    """Return the negative log-likelihood, dynamic probabilities, and score residuals.

    This is the Python counterpart to ``Bernoulli_Logl.m``.  With ``tsmc=0``
    it uses the score-driven AR(1) probability recursion from the MATLAB code.
    With ``tsmc=1`` it uses the simpler two-parameter Markov specification.
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    if not np.all(np.isin(y, (0.0, 1.0))):
        raise ValueError("The Bernoulli response must contain only zeros and ones.")

    omega = par[0]  # Unconstrained mean/intercept parameter.
    if tsmc:
        kappa = par[1]  # Dynamic parameter in the simple Markov specification.
    else:
        phi = par[1]    # Dynamic AR parameter.
        kappa = par[2]  # Dynamic conditional-score parameter.

    nobs = y.size
    mu = np.empty(nobs)
    res = np.empty(nobs)

    # Start the recursion at the sample mean, as in the MATLAB implementation.
    mu_next = np.mean(y)
    for index in range(nobs):
        mu[index] = mu_next
        res[index] = y[index] - mu_next  # For Bernoulli data, the score is the residual.
        if tsmc:
            mu_next = omega + kappa * y[index]
        else:
            mu_next = omega + phi * mu_next + kappa * res[index]

    # Invalid probabilities can occur at infeasible trial values during optimization.
    if np.any(mu <= 0.0) or np.any(mu >= 1.0) or not np.all(np.isfinite(mu)):
        return np.inf, mu, res

    nll = -np.sum(y * np.log(mu) + (1.0 - y) * np.log1p(-mu))
    return nll, mu, res
