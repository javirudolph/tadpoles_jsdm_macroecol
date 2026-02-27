# Processed Data

This folder contains data that has been cleaned, transformed, or derived from the raw data in `data/raw/`. Every file here should be fully reproducible by running the relevant scripts in `R/`.

Do not edit these files by hand. If something needs to change, update the script that produced it and re-run it.

---

## How to document your files

For every file you add to this folder, copy and fill out this template below:
```
## [Short descriptive name]

- **File:**            filename.csv
- **Produced by:**     R/scripts/01_clean_data.R
- **Input(s):**        data/raw/original_file.csv
- **Date created:**    YYYY-MM-DD
- **Contents:**        Brief description of what this file contains and what transformations were applied
- **Notes:**           Any important decisions made during processing worth remembering
```

---

## Files

<!-- Add your entries here as you create processed files -->