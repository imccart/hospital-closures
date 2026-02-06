# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Academic research project analyzing the effects of Critical Access Hospital (CAH) designation on hospital closures, mergers, and operational outcomes. The paper title is "Survival versus Scale: The Effects of Critical Access Hospital Designation."

**Research team**: Ian McCarthy (Emory), Mayra Pineda-Torres (Georgia Tech), Katie Leinenbach (Emory PhD)

## Running the Analysis

This is a pure R project with no formal build system. Execute scripts in R/RStudio:

1. **Fuzzy matching**: Run `data-code/fuzzy-match.R` first (links Form 990, HCRIS, CAH to AHA IDs)
2. **Data build**: Run `data-code/build-data.R` (requires external data symlinks and fuzzy match outputs)
3. **Analysis**: Run `analysis/_run-analysis.r` to load/clean data, then individual analysis scripts (2-4)

**Required R packages** (installed via `pacman::p_load()`):
- Data: tidyverse, purrr, lubridate, stringr, modelsummary, janitor, here, fedmatch, zipcodeR
- Analysis: fixest, did, did2s, synthdid, survival, sf, glmnet
- Visualization: ggplot2, scales, dotwhisker, plotly

**API keys**: Census API key required in `data-code/api-keys.R` (not in version control). Request from https://api.census.gov/data/key_signup.html

## Code Architecture

### Data Pipeline (`data-code/`)
- `fuzzy-match.R`: Wrapper that runs fuzzy matching scripts (fuzzy990.R, fuzzyhcris.R, fuzzycah.R)
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
- `results/diagnostics/`: Pre-trend plots, CSVs, and LaTeX tabular innards (from `app-dd-diagnostics.R`)
- `data/output/`: Main analysis datasets (`aha_final.csv`, `aha_neighbors.csv`, crosswalk files)
- `data/output/fuzzy/`: QA diagnostics from fuzzy matching scripts
- `paper.tex`: Main LaTeX manuscript (synced with Overleaf)
- `appendix.tex`: Standalone supplemental appendix (financial measure construction, robustness, SDID diagnostics)

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
| `dd-other.R` | Mixed | Various | Deleted — legacy/exploratory code (TWFE, Sun & Abraham, early implementations), superseded by stacked SDID + CS approach |

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

### Data Cleaning Notes

**IPDTOT (inpatient days per bed)**: Winsorized at 365 days per bed maximum. Raw IPDTOT values above 365 indicate data errors (impossible utilization) or bed count misreporting, particularly for small hospitals. Four hospitals had extreme values (430-843 days/bed) that were causing unstable SDID estimates with large cohort heterogeneity.

### Diagnostics Script (`app-dd-diagnostics.R`)

Covers both hospital-level (9 outcomes via `stack.hosp`) and state-level (closures, mergers via `stack.state`) diagnostics. Outputs CSVs, pre-trend PNGs, and three `.tex` tabular-innards files to `results/diagnostics/`. These `.tex` files are `\input`'d by `appendix.tex`.

**Critical**: When building data for `panel.matrices()` in synthdid, use `post_treat` (time-varying 0/1) not `treated` (ever-treated indicator). The `panel.matrices` function infers treatment timing from the data and requires a time-varying treatment indicator. Using the ever-treated indicator causes "Treatment adoption is not simultaneous" errors.

**CS control counts**: Must match the `control_group="notyettreated"` specification in the DD scripts. For cohort `c`, CS controls = never-treated (`treat_group == 0`) + not-yet-treated (`treat_group > c`). Hospital-level `treat_group` construction must use `state.cut` to match `2-hospital-dd.R`.

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

### Financial Outcome Missing Data (2026-02-04)

**Problem**: SDID requires balanced panels, but financial outcomes (margin, current_ratio, net_fixed, capex) suffer massive attrition in early cohorts. For the 1999 cohort on margin, only 3 of 84 treated hospitals survive balancing with 5 pre-periods. AHA outcomes (beds, FTEs) lose only ~27% because AHA survey coverage is near-complete.

**Root cause — ownership composition, not matching quality**: Financial data comes from HCRIS (primary, from ~1998) and Form 990 (fallback, for nonprofits only). Diagnostic analysis of the 84 treated hospitals in the 1999 cohort (after `bed.cut <= 50`) revealed:
- **47 (56%) are government-owned** (county, city, state). Government hospitals do not file Form 990 — this is structural and unfixable.
- **36 are nonprofit** (should file 990), but only 25 have any 990 match, and only 4 have pre-1999 990 data. Most 990 data appears in 2009+ (IRS e-filing mandate).
- **1 is for-profit** (no 990 filing).

For these hospitals, the only pre-treatment financial data source is HCRIS, which ramps up slowly: 7 hospitals in 1998, 27 in 1999, near-complete (~77) from 2000. The 3 hospitals surviving the balanced panel at `pre_period=5` are the rare nonprofits with Form 990 matches reaching back to 1994.

**Pre-period sensitivity** (`analysis/app-margin-preperiod.R`): Tests financial outcomes across pre-period lengths 2-5. Key finding: margin ATT is +0.035 at pre=5 (main spec) and +0.041 at pre=4, but collapses to ~0 at pre=3 and pre=2. Net fixed assets is stable (~-0.08) across all pre-periods. Current ratio and capex are too imprecise to be informative at any pre-period length. The margin instability reflects sensitivity to sample composition — which hospitals survive the balanced panel requirement — rather than precision differences.

**Decided approach**: Ridge regression imputation for margin in 1994-1996 using HCRIS PPS data (see "Margin Imputation" below). Pre-period sensitivity analysis in appendix. Other financial outcomes (current_ratio, net_fixed, capex) remain unimputed — report with caveats about small treated samples.

**Prior approaches explored and abandoned**:
- *AHA-based imputation* (archived in `archived/0-impute-financials.R`): Elastic net using AHA predictors + hospital-specific residual correction. Abandoned — AHA operational variables are poor predictors of financial outcomes (cross-sectional R² 0.01-0.39).
- *Manual 990 collection*: Feasible for ~32 nonprofits missing pre-1998 data, but superseded by HCRIS PPS imputation which covers all ownership types.
- *Improving fuzzy matching*: Not a matching quality issue — coverage ceiling is structural (ownership composition).
- *`gsynth`/`fect`*: Handles unbalanced panels natively under factor model identification. Potential future appendix analysis.

**Data additions retained**: EXPTOT and PAYTOT added to `_aha-data.R` (all four data blocks). Available for future use.

### HCRIS PPS Data Extension (2026-02-06)

**New data**: The HCRIS data repo (`imccart/HCRIS`) now includes PPS minimum dataset (1985-1999) from NBER, merged with HCRIS v1996 (1998-2011) and v2010 (2010-2020). Symlink updated.

**Variable availability by year** (key for financial outcomes):
- **1994-1995**: `tot_charges` and `tot_operating_exp` near-complete (~6,200 records), BUT `tot_operating_exp` = Medicare-specific operating costs (not total) and `tot_charges` = cost allocation charges (not revenue statement). `net_pat_rev` only available for ~120-150 records. No balance sheet.
- **1996**: Split — 3,326 records (~57%) have `net_pat_rev` from Form 2552-96 crosswalk; rest are PPS-only.
- **1997+**: Near-complete for all financial variables including `net_pat_rev`, `tot_operating_exp` (total), balance sheet. ~6,000 records/year.

**Key data facts**:
- `tot_operating_exp` jumps ~16% at 1996/1997 boundary (Medicare-specific to total definition change)
- Charge-to-revenue ratio (`tot_charges / net_pat_rev`) is stable within hospitals: correlation 0.85 between 1997 and 1998, mean ~1.72
- Medicare-specific operating costs are NOT separately available in HCRIS for 1997+ (only PPS minimum dataset has `opertots`)

### Margin Imputation (`analysis/supp-impute-margin.R`)

Extends margin coverage to 1994-1996 using HCRIS PPS variables. Sourced after `stack_hosp` in `_run-analysis.r` (line 244). Creates `margin_imp` column in `stack.hosp`: observed margin where available, ridge-imputed where missing.

**Note**: `margin_imp` is currently not wired into the outcome loop — `2-hospital-dd.R` still uses `margin`. Integration pending.

**Calibration approaches explored** (details in `scratch/diag-margin-decomposed.R`):
1. *Decomposed* (discharge-ratio scaling + charge-to-rev projection): RMSE 0.168 on 1998 holdout. Pre-1996 estimates implausibly high (mean +0.22 to +0.28 vs observed -0.03 to -0.05), indicating PPS "charges" and revenue-statement charges are not comparable.
2. *Simple additive* (hospital-specific `margin_true - margin_charges` from 1997): RMSE 0.115, correlation 0.76. Better than decomposed but cannot be applied to 1994-1995 where `tot_operating_exp` has a different definition.
3. *Ridge regression* (chosen): Predicts `margin_true` from `margin_charges`, `mcare_share`, and `charge_to_rev`. Trained on 1997, validated on 1998.

**Production model** (in `supp-impute-margin.R`):
- Ridge regression trained on 1997 within `stack.hosp` (non-CAH, bed-cut sample)
- For 1994-1996: constructs `est_total_opex = tot_operating_exp / mcare_share`, then `margin_charges = (tot_charges - est_total_opex) / tot_charges`, and uses hospital-specific `charge_to_rev_avg` from 1997-1998
- Training-support bounds enforced on all predictors
- Imputed margins recentered by year to match 1997-1998 observed mean (-0.019)

**Validation** (1998 holdout, bed-cut non-CAH sample): RMSE 0.074, MAE 0.048, correlation 0.87, bias 0.009 (n=484)

**Coverage**: Fills ~30% (1994), ~28% (1995), ~16% (1996) of missing margins in the bed-cut sample. Coverage limited by hospitals lacking charge-to-rev ratios in 1997-1998 or having predictor values outside training support.

**Exploratory diagnostics**: `scratch/diag-hcris-coverage.R` (variable availability), `scratch/diag-margin-decomposed.R` (full calibration exploration)

### Next Steps (as of 2026-02-06)

**Finalize margin imputation** — two parts:

1. **Review imputation methodology**: Revisit `supp-impute-margin.R` decisions before committing to production use. Key questions to resolve:
   - Is ridge the right model, or should we consider alternatives (e.g., OLS performed slightly better on 1998 holdout)?
   - Is recentering to the 1997-1998 mean the right normalization, or does it impose too strong an assumption?
   - Are the training-support bounds too restrictive (they limit coverage to ~30% of missing observations)?
   - How sensitive are results to the choice of training year (1997 only vs pooling 1997-1998)?
   - Should the imputation apply only to margin, or also extend to other financial outcomes?

2. **Wire `margin_imp` into the analysis pipeline**: Once the methodology is finalized, integrate imputed margins into the SDID and CS estimators. Options:
   - Add `margin_imp` as a separate outcome in `outcome_map` (run both imputed and non-imputed as robustness)
   - Replace `margin` with `margin_imp` in `stack.hosp` (simpler but less transparent)
   - The CS estimator runs on `est.dat` (not `stack.hosp`), so it would need its own imputation path or a pre-stacking imputation step
