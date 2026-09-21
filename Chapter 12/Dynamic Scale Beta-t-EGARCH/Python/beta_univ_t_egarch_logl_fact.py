"""Beta-t-EGARCH likelihood together with filtered diagnostic quantities."""

import numpy as np
from scipy.special import gammaln

try:
    from .beta_univ_t_egarch_logl import _unpack_parameters, beta_t_egarch_filter
except ImportError:
    from beta_univ_t_egarch_logl import _unpack_parameters, beta_t_egarch_filter


def beta_univ_t_egarch_logl_fact(par, y, inf, I):
    """Return negative log-likelihood, scale path, residuals, scores, Beta, and fit."""
    y = np.asarray(y, dtype=float).reshape(-1)
    _, _, _, _, df, _, _ = _unpack_parameters(par, I)
    lam, res, score, beta = beta_t_egarch_filter(par, y, inf, I)
    fit = y - res  # Legacy MATLAB output, kept for compatibility.
    log_density_constant = gammaln((df + 1) / 2) - gammaln(df / 2) - 0.5 * np.log(np.pi * df)
    log_likelihood = (
        y.size * log_density_constant
        - np.sum(lam)
        - (df + 1) / 2 * np.sum(np.log1p(res ** 2 / df))
    )
    return -log_likelihood, lam, res, score, beta, fit
