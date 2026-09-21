"""Estimate the dynamic-scale Beta-t-EGARCH model on an equity price series."""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm, t as student_t, probplot
from statsmodels.graphics.tsaplots import plot_acf

from beta_univ_t_egarch_logl import beta_univ_t_egarch_logl
from beta_univ_t_egarch_logl_fact import beta_univ_t_egarch_logl_fact
from info_crit import info_crit


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "Data"
RESULTS_DIR = PROJECT_ROOT / "Results"


def load_log_returns(filename="NSX Close---.xlsx"):
    """Load one price series and convert consecutive prices into log returns."""
    prices = pd.read_excel(DATA_DIR / filename, header=None).to_numpy(dtype=float).reshape(-1)
    return np.log(prices[1:] / prices[:-1])


def initial_parameters(y, integrated, leverage, constrained, kappa_start):
    """Construct the original MATLAB starting values and parameter bounds."""
    if integrated:
        size = 5 if leverage else 4
        start = np.full(size, kappa_start)
        start[0], start[1], start[3] = np.mean(y), np.log(np.std(y, ddof=1)), np.log(8)
        lower, upper = np.full(size, -np.inf), np.full(size, np.inf)
        lower[2] = 0
    else:
        size = 6 if leverage else 5
        start = np.full(size, kappa_start)
        start[0], start[1], start[2], start[4] = np.mean(y), np.log(np.std(y, ddof=1)), 0.9, np.log(8)
        lower, upper = np.full(size, -np.inf), np.full(size, np.inf)
        lower[2], lower[3], upper[2] = 0, 0, 1
    return start, list(zip(lower, upper)) if constrained else None


def plot_exploratory_data(y):
    """Replicate the distribution, ACF, and Q-Q diagnostics from MATLAB."""
    grid = np.linspace(np.min(y), np.max(y), 1000)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].hist(y, bins=100, density=True, alpha=0.65, label="Empirical density")
    axes[0, 0].plot(grid, norm.pdf(grid, np.mean(y), np.std(y, ddof=1)), "g-", label="Gaussian")
    axes[0, 0].set_title("Unconditional distribution of the data"); axes[0, 0].legend()
    plot_acf(y, lags=100, ax=axes[0, 1]); axes[0, 1].set_title("ACF data")
    plot_acf(y ** 2, lags=100, ax=axes[1, 0]); axes[1, 0].set_title("ACF squared data")
    probplot(y, dist="norm", plot=axes[1, 1]); axes[1, 1].set_title("Q-Q plot against Gaussian distribution")
    fig.tight_layout()


def run_model(filename="NSX Close---.xlsx", maxiter=67_000, make_plots=True):
    """Run the dynamic-scale model corresponding to DCStDynamicScaleModeling.m."""
    y = load_log_returns(filename)
    integrated_scale, information_scale, leverage_scale = 0, 1, 1
    constrained_scale, kappa_start = 1, 0.01
    if make_plots:
        plot_exploratory_data(y)

    start, bounds = initial_parameters(y, integrated_scale, leverage_scale, constrained_scale, kappa_start)
    objective = lambda par: beta_univ_t_egarch_logl(par, y, information_scale, integrated_scale)
    result = minimize(objective, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": maxiter})
    if not result.success:
        print(f"Optimizer message: {result.message}")

    RESULTS_DIR.mkdir(exist_ok=True)
    pd.DataFrame({"estimate": result.x}).to_excel(RESULTS_DIR / "estPar_Python.xlsx", index=False)
    logl, lam, res, score, beta, fit = beta_univ_t_egarch_logl_fact(result.x, y, information_scale, integrated_scale)
    bic, aic = info_crit(logl, result.x, y.size)
    mu = result.x[0]
    df = np.exp(result.x[3] if integrated_scale else result.x[4])

    if make_plots:
        fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        axes[0].plot(y - mu, label="Y - mu"); axes[0].plot(np.exp(lam), "r", label="Fitted conditional scale")
        axes[0].set_title("Scale fit: returns"); axes[0].legend()
        axes[1].plot(np.abs(y - mu), label="|Y - mu|"); axes[1].plot(np.exp(lam), "r", label="Fitted conditional scale")
        axes[1].set_title("Scale fit: absolute returns"); axes[1].legend(); fig.tight_layout()

        grid = np.linspace(np.min(res), np.max(res), 500)
        plt.figure(figsize=(10, 5)); plt.hist(res, bins=100, density=True, alpha=0.65, label="Empirical residual density")
        plt.plot(grid, student_t.pdf(grid, df), "r-", label="Fitted t")
        plt.plot(grid, norm.pdf(grid), "g-", label="Standard Gaussian")
        plt.title("Empirical fitted residuals versus fitted t distribution"); plt.legend(); plt.tight_layout()

        fig, axes = plt.subplots(3, 1, figsize=(10, 9))
        plot_acf(res, lags=100, ax=axes[0]); axes[0].set_title("ACF residuals")
        plot_acf(res ** 2, lags=100, ax=axes[1]); axes[1].set_title("ACF squared residuals")
        plot_acf(score, lags=100, ax=axes[2]); axes[2].set_title("ACF fitted scores"); fig.tight_layout()

        pit = student_t.cdf(np.sort(res), df)
        plt.figure(figsize=(10, 4)); plt.plot(pit, label="Ordered PIT")
        plt.plot(np.linspace(0, 1, pit.size), "k--", label="Uniform benchmark")
        plt.title("Ordered PIT DCS model versus uniform"); plt.xlabel("Residual rank"); plt.ylabel("PIT"); plt.legend(); plt.tight_layout()

    return {"result": result, "logl": logl, "bic": bic, "aic": aic, "lambda": lam, "residuals": res, "scores": score, "beta": beta, "fit": fit}


if __name__ == "__main__":
    run_model()
    plt.show()
