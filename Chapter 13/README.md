# Dynamic Bernoulli Model

## Overview

This chapter models the annual Boat Race outcome as a Bernoulli variable. Its probability evolves through a score-driven recursion: the previous probability and the Bernoulli score residual update the probability for the following year. The default specification is an AR(1) dynamic Bernoulli model; a simpler Markov alternative is retained as an optional setting.

## Files

- `Data/BoatRace46.xlsx`: annual Boat Race outcomes, coded 0 or 1.
- `Matlab/BernoulliSetUp.m`: data preparation, constrained maximum-likelihood estimation, and plot.
- `Matlab/Bernoulli_Logl.m`: Bernoulli likelihood and dynamic-probability filter.
- `R/DynamicBernoulliModel.R` and `R/Bernoulli_Logl.R`: R translations of the MATLAB workflow and helper function.
- `Python/dynamic_bernoulli_model.py` and `Python/bernoulli_logl.py`: Python translations.
- `Python Notebooks/DynamicBernoulliModel.ipynb`: self-contained Python notebook with explained, separate execution stages.
- `Results/`: output parameter estimates.

## Run

From this chapter folder, run MATLAB with `matlab -batch "run('Matlab/BernoulliSetUp.m')"`; R with `Rscript R/DynamicBernoulliModel.R`; or Python with `pip install -r Python/requirements.txt` followed by `python Python/dynamic_bernoulli_model.py`.

For R, install the required packages once with `install.packages(c("readxl", "nloptr"))`. For the notebook, install the Python requirements and open `Python Notebooks/DynamicBernoulliModel.ipynb` in Jupyter.
