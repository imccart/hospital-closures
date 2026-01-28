# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Academic research project analyzing the effects of Critical Access Hospital (CAH) designation on hospital closures, mergers, and operational outcomes. The paper title is "Survival versus Scale: The Effects of Critical Access Hospital Designation."

**Research team**: Ian McCarthy (Emory), Mayra Pineda-Torres (Georgia Tech), Katie Leinenbach (Emory PhD)

## Running the Analysis

This is a pure R project with no formal build system. Execute scripts in R/RStudio:

1. **Data build**: Run `data-code/build-data.R` first (requires external data symlinks)
2. **Analysis**: Run `analysis/_run-analysis.r` to load/clean data, then individual analysis scripts (0-5)

**Required R packages** (installed via `pacman::p_load()`):
- Data: tidyverse, lubridate, stringr, modelsummary, janitor, here, fedmatch, zipcodeR
- Analysis: fixest, did, did2s, synthdid, survival, sf, glmnet
- Visualization: ggplot2, scales, dotwhisker, plotly

**API keys**: Census API key required in `data-code/api-keys.R` (not in version control). Request from https://api.census.gov/data/key_signup.html

## Code Architecture

### Data Pipeline (`data-code/`)
- `build-data.R`: Main script merging external data sources (AHA, HCRIS, CAH, Form 990, COPA)
- External data are symlinked into `data/input/` from separate GitHub repos
- Outputs cleaned datasets to `data/output/`

### Analysis Pipeline (`analysis/`)
- `_run-analysis.r`: Entry point - loads data, performs cleaning/merging, creates `est.dat`
- `functions.R`: Stacking functions for difference-in-differences with staggered treatment timing
- Numbered scripts (0-5) run sequentially after `_run-analysis.r`

### Stacked Difference-in-Differences Design

The core econometric approach uses stacking to handle heterogeneous treatment timing:

```
stack_hosp(pre.period, post.period, state.period)  # Hospital-level stacking
stack_state(pre.period, post.period, state.period) # State-level stacking
stack_hosp_balance(...)                             # Balanced panel for closure outcomes
```

Key design elements:
- Treatment: `eff_year` (year hospital received CAH designation)
- State-level treatment availability: `state_treat_year`
- Control groups: "never treated" (no CAH in state) + "not yet treated" (CAH available later)
- Weights by cohort proportions: `weight_any`, `weight_notyet`
- Event time: `stacked_event_time` relative to treatment

### Key Variables
- `eff_year`: Year of individual hospital CAH designation
- `state_treat_year`: Year CAH program became available in state (0 = never)
- `treat_state`: Binary for states that ever had CAH program
- `closed`, `merged`, `closed_merged`: Hospital event outcomes
- Financial: `margin`, `current_ratio`, `net_fixed_assets`, `capex`
- Operational: `beds`, `beds_ob`, `ftes_rn`, `ip_days_per_bed`

## External Data Dependencies

Data in `data/input/` are symlinked from separate repositories. **Do not modify these files.** If a problem is found or suspected with input data, alert the user rather than attempting to write to these files.

Source repositories:
- `imccart/aha-data` - American Hospital Association surveys
- `imccart/cah` - Critical Access Hospital designations (Flex Monitoring Team)
- `imccart/HCRIS` - Healthcare Cost Report Information System
- `imccart/form-990s` - IRS Form 990 nonprofit financials

## Output

- `results/`: PNG figures and LaTeX tables for the manuscript
- `paper.tex`: Main LaTeX manuscript (synced with Overleaf)

## Important Notes

- The manuscript syncs with Overleaf via git. Check for Overleaf merge commits before editing `paper.tex`.
- `est.dat` is the main analysis dataset created by `_run-analysis.r` and used by all analysis scripts.
- The `stack_hosp_balance` function handles "zombie" hospital IDs (IDs reused after closure/merger).

## Workflow for Experimentation

When experimenting with new data management tasks, analysis, or regression specifications, create a temporary script in a scratch folder (e.g., `scratch/analysis_test.R`) instead of modifying the original code files. Only modify original files after explicit user approval.
