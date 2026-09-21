# Dynamic Scale Beta-t-EGARCH Model

## Overview

This chapter estimates a univariate dynamic conditional scale model with a Student-t conditional distribution. The conditional log scale evolves through a score-driven EGARCH-type recursion, allowing the model to react to new information and, in the active specification, to negative-return leverage effects. Parameters are estimated by maximum likelihood and the workflow reports fitted scales, standardized residuals, score diagnostics, information criteria, and probability-integral-transform diagnostics.

## Files in this chapter

```
Data/                 Equity closing-price series used as model inputs
Matlab/               MATLAB implementation
R/                    R implementation
Python/               Python implementation
Python Notebooks/     Interactive, self-contained Python implementation
Results/              Estimation outputs created when the code runs
```

The language folders contain the following files:

- `DCStDynamicScaleModeling`: the main estimation and diagnostic workflow.
- `Beta_univ_t_Egarch_Logl`: the negative log-likelihood function.
- `Beta_univ_t_Egarch_Logl_fact`: the filter that recovers fitted model quantities.
- `Info_Crit`: the BIC and AIC calculation.

The helper files are called by the corresponding main script and are not intended to be run by themselves. `Python Notebooks/DCStDynamicScaleModeling.ipynb` is self-contained and presents the same workflow in small, ordered cells.

## Tested software versions

| Language | Version |
| --- | --- |
| MATLAB | R2023b |
| Python | 3.9.13 |
| R | 4.5.0 |

## Run the code

Open a terminal in this chapter folder (`Beta-t-EGARCH/`). All paths are relative to this folder.

### MATLAB

```bash
matlab -batch "run('Matlab/DCStDynamicScaleModeling.m')"
```

Alternatively, open and run `Matlab/DCStDynamicScaleModeling.m` in MATLAB. Estimates are written to `Results/estPar.xlsx`.

### R

Install the only non-base dependency once:

```r
install.packages("readxl")
```

Then run:

```bash
Rscript R/DCStDynamicScaleModeling.R
```

R estimates are written to `Results/estPar_R.csv`.

### Python

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
pip install -r Python/requirements.txt
python Python/dcst_dynamic_scale_modeling.py
```

Python estimates are written to `Results/estPar_Python.xlsx`.

### Python notebook

After installing the Python requirements, install and launch Jupyter if needed:

```bash
pip install jupyterlab
jupyter lab
```

Open `Python Notebooks/DCStDynamicScaleModeling.ipynb` and run the cells in order.

## Data and outputs

The default input is `Data/NSX Close---.xlsx`. The other price series in `Data/` can be selected by changing the active input filename in the main script or notebook. Input files contain one chronological column of closing prices without a header.

Full maximum-likelihood estimation may take time. For a quick check, reduce the maximum-iteration setting in the relevant main script or notebook. Output files are saved in `Results/`.
