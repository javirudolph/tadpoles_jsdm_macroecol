# Outputs

This folder contains the outputs of simulation runs from scripts in `R/`. These files are not tracked by Git — they are reproducible by running the scripts with the documented parameters and seeds.

---

## Organization

As your project grows, organize outputs into subfolders by simulation type, date, or run version. A consistent structure will save you hours of confusion later. Some options that work well:
```
outputs/
├── 2026-01-15_pilot_run/
├── 2026-02-01_full_sweep_R0/
├── 2026-02-20_sensitivity_dispersal/
└── final/
```

Or by simulation type if you have clearly distinct experiments:
```
outputs/
├── 2026-01-15_pilot_run/
├── 2026-02-01_full_sweep/
└── 2026-02-20_sensitivity/nsitivity/
```

Pick one convention and stick to it. When in doubt, date-based folders are hard to argue with.

---

## The research loop

Simulation outputs don't exist in isolation — they are part of an iterative process that connects notebooks, scripts, outputs, and figures:

1. **Notebook** (`notebooks/`) — you explore the data, design the simulation, choose parameters, and document your reasoning
2. **Script** (`R/`) — you write and run the simulation script, which may go to HPC for large jobs
3. **Output** (`outputs/`) — results land here, organized by run
4. **Notebook** (`notebooks/`) — you load the outputs, visualize, and interpret
5. **Script** (`R/`) — you write a figure script that reads outputs and saves polished figures
6. **Figures** (`figures/`) — final figures land here, ready for the manuscript

This README should reflect where you are in that loop for each run. Link back to the notebook that motivated each simulation so the reasoning is never lost.

---

## Run metadata — document every run

For every simulation run, copy and fill out this template:
```
## [Run name or folder]

- **Date:**              YYYY-MM-DD
- **Script:**            R/simulations/sim_main.R
- **Notebook:**          notebooks/02_simulation_design.qmd
- **Parameters:**        R0 = 2.5, dispersal = fat-tailed, n_reps = 100
- **Seed:**              set.seed(42) or array job seeds 1:100
- **HPC job:**           Job ID if submitted to HiPerGator, or "local"
- **Output files:**      Brief description of what was saved and in what format
- **Notes:**             Anything worth remembering — convergence issues, unexpected results, follow-up questions
- **Status:**            Complete / In progress / Archived
```

---

## Runs

<!-- Add your run entries here -->