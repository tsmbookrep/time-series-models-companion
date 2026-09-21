"""Estimate the score-driven dynamic Bernoulli model for Boat Race outcomes."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint, minimize

from bernoulli_logl import bernoulli_logl


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "Data"
RESULTS_DIR = PROJECT_ROOT / "Results"


def load_boat_race_outcomes(filename="BoatRace46.xlsx"):
    """Load the annual Boat Race binary outcomes and fill inherited missing values."""
    outcomes = pd.to_numeric(
        pd.read_excel(DATA_DIR / filename, header=None).iloc[:, 0], errors="coerce"
    ).ffill().dropna()
    y = outcomes.to_numpy(dtype=float)
    if not np.all(np.isin(y, (0.0, 1.0))):
        raise ValueError("Boat Race outcomes must be coded as 0 or 1.")
    return y


def model_setup(tsmc=0):
    """Reproduce the MATLAB starting values, bounds, and probability constraints."""
    if tsmc:
        # Simple Markov chain: mu_(t+1) = omega + kappa * y_t.
        start = np.array([0.1, 0.1])
        constraint_matrix = np.array([[1.0, 0.0], [1.0, 1.0]])
        bounds = [(-np.inf, np.inf), (-1.0, 1.0)]
    else:
        # Score-driven AR(1): mu_(t+1) = omega + phi*mu_t + kappa*(y_t-mu_t).
        start = np.array([0.0, 0.99, 0.1])
        constraint_matrix = np.array(
            [[1.0, 0.0, 0.0], [1.0, 0.0, 1.0],
             [1.0, 1.0, 0.0], [1.0, 1.0, -1.0]]
        )
        bounds = [(-np.inf, np.inf), (-1.0, 1.0), (-np.inf, np.inf)]

    # The MATLAB C matrix imposes 0 <= C*par <= 1.
    constraints = [LinearConstraint(constraint_matrix, 0.0, 1.0)]
    return start, bounds, constraints


def run_model(filename="BoatRace46.xlsx", tsmc=0, maxiter=67_000, make_plot=True):
    """Estimate the selected Bernoulli specification and optionally plot fitted probabilities."""
    y = load_boat_race_outcomes(filename)
    start, bounds, constraints = model_setup(tsmc)
    objective = lambda par: bernoulli_logl(par, y, tsmc)[0]

    result = minimize(
        objective, start, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": maxiter, "ftol": 1e-8, "disp": True},
    )
    if not result.success:
        print(f"Optimizer message: {result.message}")

    nll, mu, res = bernoulli_logl(result.x, y, tsmc)
    RESULTS_DIR.mkdir(exist_ok=True)
    pd.DataFrame({"estimate": result.x}).to_excel(
        RESULTS_DIR / "estPar_Python.xlsx", index=False
    )

    if make_plot:
        from matplotlib import pyplot as plt

        years = np.arange(1946, 1946 + y.size)
        plt.figure(figsize=(10, 5))
        plt.plot(years, y, ".", markersize=14, label="Winner (0 Oxf, 1 Cam)")
        plt.plot(years, mu, "-", linewidth=1.2, color="#d95319", label="Dyn. Prob. (AR1)")
        plt.xticks(np.arange(1950, years[-1] + 1, 10))
        plt.xlim(years[0], years[-1])
        plt.xlabel("Year"); plt.ylabel("Result and probability")
        plt.legend(loc="upper right"); plt.tight_layout()

    return {"result": result, "negative_log_likelihood": nll, "mu": mu, "residuals": res}


if __name__ == "__main__":
    run_model()
    from matplotlib import pyplot as plt
    plt.show()
