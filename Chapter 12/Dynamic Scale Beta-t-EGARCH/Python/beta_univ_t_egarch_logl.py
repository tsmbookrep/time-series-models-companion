"""Negative log-likelihood for the Beta-t-EGARCH dynamic-scale model."""

import numpy as np
from scipy.special import gammaln


def _unpack_parameters(par, I):
    """Translate the MATLAB parameter layouts into named Python values."""
    par = np.asarray(par, dtype=float).reshape(-1)
    kappa_2 = 0.0
    has_leverage = False

    if I == 1:
        if par.size not in (4, 5):
            raise ValueError("The integrated specification requires 4 or 5 parameters.")
        mu, omega, kappa, log_df = par[:4]
        phi = None
        if par.size == 5:
            kappa_2 = par[4]
            has_leverage = True
    else:
        if par.size not in (5, 6):
            raise ValueError("The stationary specification requires 5 or 6 parameters.")
        mu, omega, phi, kappa, log_df = par[:5]
        if par.size == 6:
            kappa_2 = par[5]
            has_leverage = True

    return mu, omega, phi, kappa, np.exp(log_df), kappa_2, has_leverage


def beta_t_egarch_filter(par, y, inf, I):
    """Filter log scale, standardized residuals, scores, and the Beta transform."""
    y = np.asarray(y, dtype=float).reshape(-1)
    mu, omega, phi, kappa, df, kappa_2, has_leverage = _unpack_parameters(par, I)
    n_obs = y.size
    lam = np.zeros(n_obs)
    res = np.zeros(n_obs)
    score = np.zeros(n_obs)
    beta = np.zeros(n_obs)
    information = 2 * df / (df + 3) if inf == 1 else 1.0
    next_lam = omega

    for index in range(n_obs):
        lam[index] = next_lam
        res[index] = (y[index] - mu) * np.exp(-lam[index])
        beta[index] = (res[index] ** 2 / df) / (1 + res[index] ** 2 / df)
        score[index] = ((df + 1) * beta[index] - 1) / information
        leverage = kappa_2 * np.sign(mu - y[index]) * (score[index] + 1) if has_leverage else 0.0
        if I == 1:
            next_lam = lam[index] + kappa * score[index] + leverage
        else:
            next_lam = omega * (1 - phi) + phi * lam[index] + kappa * score[index] + leverage

    return lam, res, score, beta


def beta_univ_t_egarch_logl(par, y, inf, I):
    """Return the negative Student's t log-likelihood, as in the MATLAB function."""
    _, _, _, _, df, _, _ = _unpack_parameters(par, I)
    lam, res, _, _ = beta_t_egarch_filter(par, y, inf, I)
    log_density_constant = gammaln((df + 1) / 2) - gammaln(df / 2) - 0.5 * np.log(np.pi * df)
    log_likelihood = (
        res.size * log_density_constant
        - np.sum(lam)
        - (df + 1) / 2 * np.sum(np.log1p(res ** 2 / df))
    )
    return -log_likelihood
