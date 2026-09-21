# Replication Code and Data

## Overview

This repository contains the code and data accompanying the chapters of this project. Each chapter is self-contained: its data, MATLAB, R, Python, Python-notebook, and result folders live together. A reader can download the repository, select a chapter, and run its code without editing computer-specific paths.

## Repository structure

The top-level repository is organized by chapter. Every chapter follows this pattern:

```
<repository-root>/
  Chapter-01/
    Data/                 Input data for Chapter 01 only
    Matlab/               MATLAB main script and helper functions
    R/                    R main script and helper functions
    Python/               Python main script and helper modules
    Python Notebooks/     Self-contained Python Jupyter notebooks
    Results/              Outputs generated for Chapter 01
    README.md             Chapter-specific overview and instructions
  Chapter-02/
    ...                   Same self-contained structure
  README.md               This general guide
```

Each chapter README describes its model, its files, its dependencies, and the tested software versions.

## General workflow

1. Clone or download the repository.
2. Open a terminal in the folder of the chapter you want to run.
3. Keep that chapter's `Data/` folder beside its language folders.
4. Read the chapter README and follow its language-specific instructions.

```bash
git clone <repository-url>
cd <repository-folder>/Chapter-01
```

The source files use paths relative to the chapter structure. A complete chapter folder can therefore be moved, copied, or cloned to another computer without modifying paths.

## Running a chapter

The exact main filenames and package versions are specified in each chapter README. In general, run the main script from the chapter folder, not from inside its data folder.

### MATLAB

Open the chapter's main `.m` file in MATLAB and run it, or use the command stated in the chapter README. MATLAB helper functions remain in the chapter's `Matlab/` folder, while inputs and outputs are read from `Data/` and written to `Results/`.

### R

Install the R packages listed by the chapter, then run its main script with `Rscript`. For example:

```bash
Rscript R/<main-script>.R
```

### Python scripts

Create and activate a virtual environment if desired, install the chapter's requirements, and run its main script:

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
pip install -r Python/requirements.txt
python Python/<main-script>.py
```

### Python Jupyter notebooks

Install the chapter's Python requirements and a Jupyter interface, then open the requested notebook from that chapter's `Python Notebooks/` folder:

```bash
pip install jupyterlab
jupyter lab
```

Run the notebook cells from top to bottom unless its chapter README says otherwise.

## Reproducibility notes

Do not separate a chapter's `Data/` folder from its language folders, and do not move a main script away from the helper files in its language folder. Within that structure, all paths are portable and no source file refers to a user name, desktop location, or local repository path.

## Authors

**Andrew Harvey

## Textbook or reference

<!-- To be completed. -->
