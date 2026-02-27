# Notebooks

This is where your digital lab notebooks live.

These are Quarto (`.qmd`) documents that tell the story of your research process — your exploratory data analysis, model tests, simulation runs, parameter explorations, and any dead ends worth remembering. Think of them as a research diary that happens to run code.

Unlike the scripts in `R/`, notebooks are not meant to be lean and silent. They are meant to be read. Write in them like you are explaining your work to a collaborator — or to yourself six months from now. Include your interpretations, your questions, your surprises, and your decisions. The goal is that anyone who reads through these notebooks can follow the arc of your thinking from raw data to final results.

---

## What belongs here

- Exploratory data analysis and first looks at new datasets
- Simulation design and pilot runs
- Model fitting attempts, including the ones that didn't work
- Sanity checks and diagnostic plots
- Parameter sweeps and sensitivity explorations
- Notes on decisions made and why

## What does not belong here

- Final, clean production scripts — those go in `R/`
- The manuscript — that lives in `manuscript/`

---

## Suggested naming convention

Number your notebooks so they have a natural reading order:
```
00_setup.qmd
01_data_exploration.qmd
02_simulation_design.qmd
03_pilot_runs.qmd
04_results_summary.qmd
```

---

## A few good habits

- Start every notebook with a short paragraph describing what you are trying to figure out that day
- Always `set.seed()` before anything stochastic
- Always end with `sessionInfo()`
- Commit your notebooks to Git regularly, even when they are messy — that is the point