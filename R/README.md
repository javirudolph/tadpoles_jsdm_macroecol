# R

This folder contains all R scripts for data processing, analysis, and simulations. Scripts here are standalone — they read inputs, do their job, and write outputs. Narrative, interpretation, and exploratory work live in `notebooks/`.

Simulation scripts may be run locally or submitted to HPC (e.g., UF HiPerGator). Their outputs go into `outputs/`.

---

## Organization

Number scripts to indicate the order they should be run:
```
R/
├── 01_clean_data.R
├── 02_prepare_covariates.R
├── 03_simulation.R
├── 04_summarize_results.R
└── 05_figures.R
```

---

## Script index

- `01_clean_data.R` — 
- `02_prepare_covariates.R` — 
- `03_simulation.R` — 
- `04_summarize_results.R` — 
- `05_figures.R` —

## Good habits

- Use `here::here()` for all file paths — never hardcode absolute paths
- Always `set.seed()` before anything stochastic
- Write outputs explicitly to `outputs/` or `figures/` — never to the script directory
- Include `sessionInfo()` at the end of scripts that produce key outputs