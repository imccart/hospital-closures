# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Academic research project analyzing the effects of Critical Access Hospital (CAH) designation on hospital closures, mergers, and operational outcomes. The paper title is "Survival versus Scale: The Effects of Critical Access Hospital Designation."

**Research team**: Ian McCarthy (Emory), Mayra Pineda-Torres (Georgia Tech), Katie Leinenbach (Emory PhD)

## Running the Analysis

This is a pure R project with no formal build system. Execute scripts in R/RStudio:

1. **Data build**: Run `data-code/build-data.R` first (requires external data symlinks)
2. **Analysis**: Run `analysis/_run-analysis.r` to load/clean data, then individual analysis scripts (0-4)

**Required R packages** (installed via `pacman::p_load()`):
- Data: tidyverse, purrr, lubridate, stringr, modelsummary, janitor, here, fedmatch, zipcodeR
- Analysis: fixest, did, did2s, synthdid, survival, sf, glmnet
- Visualization: ggplot2, scales, dotwhisker, plotly

**API keys**: Census API key required in `data-code/api-keys.R` (not in version control). Request from https://api.census.gov/data/key_signup.html

## Code Architecture

### Data Pipeline (`data-code/`)
- `build-data.R`: Main script merging external data sources (AHA, HCRIS, CAH, Form 990, COPA)
- External data are symlinked into `data/input/` from separate GitHub repos
- Outputs cleaned datasets to `data/output/`
- Computes nearest neighbor distances for each hospital-year using vectorized haversine formula

### Analysis Pipeline (`analysis/`)
- `_run-analysis.r`: Entry point - loads data, performs cleaning/merging, creates `est.dat`
- `functions.R`: Stacking functions for difference-in-differences with staggered treatment timing
- Numbered scripts (0-4) run sequentially after `_run-analysis.r`
- `app-dd-diagnostics.R`: Appendix diagnostics (pre-trends, sample comparison, synthetic weight concentration)

### Stacked Difference-in-Differences Design

The core econometric approach uses stacking to handle heterogeneous treatment timing:

```
stack_hosp(pre.period, post.period, state.period)  # Hospital-level stacking
stack_state(pre.period, post.period, state.period) # State-level stacking
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
- `results/diagnostics/`: Pre-trend plots and SDID/CS diagnostic tables (from `app-dd-diagnostics.R`)
- `data/output/`: Main analysis datasets (`aha_final.csv`, `aha_neighbors.csv`, crosswalk files)
- `data/output/fuzzy/`: QA diagnostics from fuzzy matching scripts
- `paper.tex`: Main LaTeX manuscript (synced with Overleaf)

## Important Notes

- The manuscript syncs with Overleaf via git. Check for Overleaf merge commits before editing `paper.tex`.
- `est.dat` is the main analysis dataset created by `_run-analysis.r` and used by all analysis scripts.

## Workflow for Experimentation

When experimenting with new data management tasks, analysis, or regression specifications, create a temporary script in a scratch folder (e.g., `scratch/analysis_test.R`) instead of modifying the original code files. Only modify original files after explicit user approval.

## Current Status

### Unified Analysis Pipeline (completed)
The analysis pipeline in `_run-analysis.r` now uses a unified outcome loop with `results.table`:

- **`outcome_map`** defines all outcomes with their script paths, labels, file stubs, and cohort years
- **Outcome loop** iterates through outcomes, sources the appropriate DD script, and collects ATT estimates
- **`results.table`** collects SDID and CS estimates with 95% CIs for each outcome
- **LaTeX output** writes to `results/att_overall.tex` with formatted table rows

### DD Scripts Architecture

| Script | Level | Outcomes | Notes |
|--------|-------|----------|-------|
| `2-hospital-dd.R` | Hospital | margin, current_ratio, net_fixed, capex, beds, OB beds, FTE RNs, IP days, system | Continuous outcomes using `stack.hosp` |
| `3-changes-state-dd.R` | State | closures, mergers | Count outcomes using `stack.state`, normalized to rates per 100 hospitals |
| `4-changes-hospital-hazard.R` | Hospital | time-to-closure/merger | Survival/hazard models |
| `dd-other.R` | Mixed | Various | Legacy/exploratory code (TWFE, Sun & Abraham, early implementations) |

### Stacking Functions (`functions.R`)

- **`stack_hosp()`**: Hospital-level stacking for continuous outcomes
- **`stack_state()`**: State-level stacking, collapses to state-year counts

### Archived Code (`archived/`)

#### Hospital-Level DD for Closures/Mergers (archived 2026-01-29)

**Files**: `archived/3-changes-hospital-dd.R`, `archived/functions.R` (contains `stack_hosp_balance()`)

**Why abandoned**: Hospital-level difference-in-differences for closure/merger outcomes is not econometrically sound due to a fundamental identification problem:

1. **Treatment is conditional on survival**: A hospital must remain open to receive CAH designation. By construction, every treated hospital survived until its treatment year—we cannot observe hospitals that "would have received CAH but closed first."

2. **Survivorship bias in the treated group**: The comparison becomes "hospitals that survived long enough to get CAH" vs. "everyone else (including those who closed before they could be treated)." This is not a valid counterfactual.

3. **Selection on risk**: Hospitals likely select into CAH based on factors correlated with closure risk. The direction is ambiguous:
   - At-risk hospitals might seek CAH as a lifeline (negative selection)
   - Only financially stable hospitals can navigate the application (positive selection)

4. **Balanced panel doesn't solve it**: The `stack_hosp_balance()` function created balanced panels by filling forward after closure, but this addresses panel attrition, not the selection-into-treatment problem.

**Credible alternatives**:
- **State-level analysis** (`3-changes-state-dd.R`): Treatment is "state enacted CAH program"—a policy shock plausibly exogenous to individual hospital outcomes
- **Hazard models** (`4-changes-hospital-hazard.R`): Explicitly model time-to-event with CAH as time-varying covariate, though still requires careful assumption defense
- **Intent-to-treat on eligibility**: Define treatment as "eligible for CAH" based on pre-determined characteristics rather than actual adoption (not currently implemented)

### Active Outcomes in Current Analysis
- **Hospital continuous**: margin, current_ratio, net_fixed, capex, BDTOT, OBBD, FTERN, IPDTOT, system (cohorts 1999-2001)
- **State counts**: closures, mergers (cohorts 1999-2001)

### Fuzzy Matching Pipeline (`data-code/fuzzy*.R`)

Links external data sources (HCRIS, CAH, Form 990) to AHA hospital IDs via fuzzy name matching.

| Script | Links | Output | Expected Match Rate |
|--------|-------|--------|---------------------|
| `fuzzy-match.R` | Wrapper | Loads data, cleans state/zip, creates `aha.collapsed` | - |
| `fuzzyhcris.R` | HCRIS MCRNUM → AHA ID | `unique_hcris.csv` | ~73% (90% via direct AHA link) |
| `fuzzycah.R` | CAH supplement → AHA ID | `unique_cah.csv` | ~73% of CAH records |
| `fuzzy990.R` | Form 990 EIN → AHA ID | `unique_990.csv` | ~37% (many for-profits don't file) |

#### Core Design Decisions

**1. Time-invariant matching**
- AHA ID and MCRNUM are persistent identifiers, so matching is done at hospital level (not hospital-year)
- `aha.collapsed` takes modal name/city/zip across all years for each ID
- Rationale: A match found in year 1 should apply to all years; year-based blocking would miss valid matches

**2. Geographic data cleaning** (`fuzzy-match.R`)
- Invalid state codes (e.g., "0", "ml", "n") → NA
- Fill missing states: (1) zip-to-state crosswalk, then (2) modal state from other years of same hospital
- Fill missing zips: modal zip from other years of same hospital
- Invalid zips ("00000", non-5-digit) → NA
- Rationale: Recover geographic data where possible; ~135 hospitals remain with NA state after cleaning

**3. Ensemble matching strategies**
Each script uses multiple blocking/matching approaches, combined before selection:
- Stringdist (Jaro-Winkler) with various blocking (state, state+city, state+zip)
- Word-level Jaccard with blocking
- State-only fallback with stricter 0.85 threshold (vs 0.75-0.80 for primary strategies)
- HCRIS also matches against all historical name variants
- Rationale: Different strategies catch different cases; ensemble improves recall without sacrificing precision

**4. Threshold selection**
- HCRIS/CAH primary strategies: 0.75 minimum score
- Form 990 primary strategies: 0.80 minimum score (stricter due to diverse nonprofit population)
- State-only fallback: 0.85 minimum score (stricter because looser blocking increases false positive risk)
- Rationale: Form 990 data contains many non-hospital nonprofits, requiring tighter filtering

**5. Deterministic tie-breaking**
- Sort order: `score (desc) → city_exact (desc) → zip_exact (desc) → alphabetical ID`
- Score is primary; geographic exactness only breaks ties
- Rationale: Avoids arbitrary selection when multiple candidates tie

**6. 1:1 enforcement** (HCRIS only)
- After selecting best match per AHA ID, also select best match per MCRNUM
- Rationale: Prevents many-to-one mappings; each MCRNUM should map to one AHA ID

**7. Many-to-many handling** (Form 990)
- An EIN can have multiple zips over time (address changes)
- Joins use `relationship = "many-to-many"` to avoid warnings
- Downstream deduplication (best match per ID) prevents row multiplication
- Rationale: The goal is EIN→AHA ID; which zip produced the match is irrelevant

**8. Manual match priority** (Form 990)
- Prior manual matches (from Hannah's data) take priority over fuzzy matches
- Fuzzy matches only fill gaps for unmatched AHA IDs
- Rationale: Human-verified matches are more reliable

#### QA Outputs

Each script produces diagnostics in `data/output/fuzzy/` for manual review:
- `qa_*_diagnostics.csv`: Top 3 candidates per ID with score gaps, ambiguous flag (gap < 0.05)
- `qa_*_summary.csv`: Score distributions (mean, median, p10, p90), match counts by source
- `qa_*_by_state.csv`: Match rates by state (useful for identifying systematic issues)

#### Helper Functions (`data-code/functions.R`)

- `modal_value(x)`: Returns most frequent non-NA value (handles all-NA groups safely)
- `haversine_vec()`: Vectorized great-circle distance calculation (returns miles by default to match `zipcodeR::zip_distance`)
- `preprocess_hospital_name()`: Standardizes names (removes system prefixes like "HCA", expands abbreviations like "St." → "Saint", normalizes "Medical Center" → "Hospital")
- `match_jaccard_words()`: Word-level Jaccard matching with blocking variables

#### Form 990 Keyword Filtering

Only Form 990 records matching hospital keywords are considered (reduces noise from non-hospital nonprofits):
- **Include if contains**: hospital, medical, health
- **Exclude patterns**: Simple string patterns in a vector, combined with `paste(collapse = "|")`:
  ```r
  exclude_patterns <- c(
    "foundation", "fdn", "fndtn", "fdtn",
    "trust", " tr$", " tr ",
    "auxiliary", "guild", "volunteer",
    "self.ins", "selfins",
    "employee benefit", "emp bene", "bene tr", "benefit plan",
    "insurance plan", "insurance guaranty", "insurance guarnaty",
    "development fund", "endowment fund", "dues fund", "health fund",
    "retirement"
  )
  ```
- **Threshold**: 0.80 minimum match score

**Pattern notes**:
- `self.ins` uses `.` as regex wildcard to match "self-ins", "self ins", etc.
- ` tr$` matches " tr" at end of string (common abbreviation for "trust")
- ` tr ` matches " tr " mid-string

**Rationale for exclusions**: These are separate legal entities from hospitals with fundamentally different financials:
- *Foundations*: Donations, investment income, grants (not patient revenue)
- *Self-insurance trusts*: Employee benefit reserves and claims payments
- *Development/endowment funds*: Fundraising and investment management

Matching a hospital to any of these entities would produce invalid financial metrics. ~1,282 non-hospital entities excluded from matching pool (2026-01-31).
