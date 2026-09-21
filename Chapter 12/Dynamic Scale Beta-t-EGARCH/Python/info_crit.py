"""Information criteria used by the DCSt model."""

import numpy as np


def info_crit(logl, par, n):
    """Return BIC and AIC using the likelihood convention supplied by the caller."""
    n_parameters = np.asarray(par).size
    bic = -2 * logl + n_parameters * np.log(n)
    aic = 2 * n_parameters - 2 * logl
    return bic, aic
