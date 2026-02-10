# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Academic research project analyzing the effects of Critical Access Hospital (CAH) designation on hospital closures, mergers, and operational outcomes. The paper title is "Survival versus Scale: The Effects of Critical Access Hospital Designation."

**Research team**: Ian McCarthy (Emory), Mayra Pineda-Torres (Georgia Tech), Katie Leinenbach (Emory PhD)

## Running the Analysis

This is a pure R project with no formal build system. Execute scripts in R/RStudio:

1. **Fuzzy matching**: Run `data-code/fuzzy-match.R` first (links Form 990, HCRIS, CAH to AHA IDs)
2. **Data build**: Run `data-code/build-data.R` (requires external data symlinks and fuzzy match outputs)
3. **Estimation data**: Run `analysis/_build-estimation-data.r` to load/clean/merge data, construct variables, and write `est.dat` and `state.dat` CSVs
4. **Analysis**: Run `analysis/_run-analysis.r` to read estimation data, stack, and run all estimation scripts

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
- `_build-estimation-data.r`: Data build — loads raw data, merges, constructs variables, writes `est.dat` and `state.dat` CSVs, runs summary stats. Run separately when data inputs change.
- `_run-analysis.r`: Orchestrator — reads saved CSVs, sets globals, builds stacked datasets (`stack.hosp`, `stack.state`, `stack.elig`), then sequentially sources DD scripts and diagnostics. Does not rebuild data and does not contain outcome loops.
- `functions.R`: Stacking functions for difference-in-differences with staggered treatment timing (`stack_hosp`, `stack_hosp_elig`, `stack_state`)
- DD scripts (`2-hospital-dd.R`, `3-hospital-dd-alt.R`, `4-changes-state-dd.R`) are self-contained: each owns its outcome map, loop, results collection, and tex/figure output
- `app-dd-diagnostics.R`: Appendix diagnostics (pre-trends, sample comparison, synthetic weight concentration)

### Stacked Difference-in-Differences Design

The core econometric approach uses stacking to handle heterogeneous treatment timing:

```
stack_hosp(pre.period, post.period, state.period)    # Hospital-level, state-timing controls
stack_hosp_elig(pre.period, post.period, cohort.years) # Hospital-level, no state restriction
stack_state(pre.period, post.period, state.period)    # State-level stacking
```

Key design elements:
- Treatment: `eff_year` (year hospital received CAH designation)
- State-level treatment availability: `state_treat_year`
- `stack_hosp`: Controls from never-treated states + not-yet-treated states (requires `state.period`)
- `stack_hosp_elig`: Controls are all non-CAH or not-yet-CAH hospitals regardless of state (used by `3-hospital-dd-alt.R`)
- Weights by cohort proportions: `weight_any`, `weight_notyet`
- Event time: `stacked_event_time` relative to treatment

### Key Variables
- `eff_year`: Year of individual hospital CAH designation
- `state_treat_year`: Year CAH program became available in state (0 = never)
- `treat_state`: Binary for states that ever had CAH program
- `closed`, `merged`, `closed_merged`: Hospital event outcomes
- Financial: `margin`, `current_ratio`, `net_fixed_assets`, `capex`, `net_pat_rev`, `tot_operating_exp`
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
- `est.dat` and `state.dat` are the main analysis datasets, built by `_build-estimation-data.r` and written to `data/output/`. `_run-analysis.r` reads these CSVs — it does not rebuild the data. Re-run `_build-estimation-data.r` when underlying data inputs change.

## Workflow for Experimentation

When experimenting with new data management tasks, analysis, or regression specifications, create a temporary script in a scratch folder (e.g., `scratch/analysis_test.R`) instead of modifying the original code files. Only modify original files after explicit user approval.

## Current Status

### Analysis Pipeline Structure (completed)

The analysis is split into two entry points:

- **`_build-estimation-data.r`**: Loads raw data, merges, constructs all variables (financial, capacity, event indicators), applies winsorization/interpolation, computes IPW weights, writes `est.dat` and `state.dat` CSVs, and runs summary statistics (`1-sum-stats.R`). Run once when data inputs change.
- **`_run-analysis.r`**: Reads saved CSVs, sets global parameters (`bed.cut`, `post`, `state.cut`, `financial.pre`), builds three stacked datasets (`stack.hosp`, `stack.state`, `stack.elig`), then sequentially sources:
  1. `2-hospital-dd.R` → `hosp.results.table`, `hosp.cohort.results`, `results/att_cohort.tex`
  2. `3-hospital-dd-alt.R` → `elig.results`, `elig.cohort.results`, `results/att_elig_overall.tex`, `results/att_elig_cohort.tex`
  3. `4-changes-state-dd.R` → `state.results.table`
  4. Assembles combined `results/att_overall.tex` (with `$N_{tr}$` column) from `hosp.results.table` + `state.results.table`
  5. Forest plot (`results/forest-sdid.png`): SDID ATTs comparing state-timing vs eligibility-restricted control groups
  6. Diagnostics (`app-dd-diagnostics.R`) and sensitivity (`app-financial-preperiod.R`)

### DD Scripts Architecture

Each DD script is self-contained: owns its outcome map, loop, results collection, and tex/figure output. The orchestrator (`_run-analysis.r`) provides globals and stacked datasets; each script produces a results tibble left in the shared environment.

**Globals contract** — `_run-analysis.r` provides: `bed.cut`, `post`, `state.cut`, `financial.pre`, `est.dat`, `state.dat`, `stack.hosp`, `stack.state`, `stack.elig`

| Script | Level | Outcomes | Produces | Notes |
|--------|-------|----------|----------|-------|
| `2-hospital-dd.R` | Hospital | margin, current_ratio, net_fixed, capex, net_pat_rev, tot_operating_exp, beds, OB beds, FTE RNs, IP days, system | `hosp.results.table`, `hosp.cohort.results`, `att_cohort.tex` | Uses `stack.hosp`, cohorts 1999-2001 |
| `3-hospital-dd-alt.R` | Hospital | All hospital-level | `elig.results`, `elig.cohort.results`, `att_elig_overall.tex`, `att_elig_cohort.tex` | Eligibility-restricted, uses `stack.elig`, cohorts 1999-2005 |
| `4-changes-state-dd.R` | State | closures, mergers | `state.results.table` | Uses `stack.state`, rates per 100 hospitals, cohorts 1999-2001 |
| `5-changes-hospital-hazard.R` | Hospital | time-to-closure/merger | — | Survival/hazard models |

### Stacking Functions (`functions.R`)

- **`stack_hosp()`**: Hospital-level stacking with state-timing controls (main identification)
- **`stack_hosp_elig()`**: Hospital-level stacking without state restriction (alternative identification, cohorts 1999-2005)
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
- **State-level analysis** (`4-changes-state-dd.R`): Treatment is "state enacted CAH program"—a policy shock plausibly exogenous to individual hospital outcomes
- **Hazard models** (`5-changes-hospital-hazard.R`): Explicitly model time-to-event with CAH as time-varying covariate, though still requires careful assumption defense
- **Intent-to-treat on eligibility**: Define treatment as "eligible for CAH" based on pre-determined characteristics rather than actual adoption (not currently implemented)

### Active Outcomes in Current Analysis
- **Hospital continuous**: margin, current_ratio, net_fixed, capex, net_pat_rev, tot_operating_exp, BDTOT, OBBD, FTERN, IPDTOT, system (cohorts 1999-2001)
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

**Pre-period sensitivity** (`analysis/app-financial-preperiod.R`): Tests all 6 financial outcomes (margin, current_ratio, net_fixed, capex, net_pat_rev, tot_operating_exp) across pre-period lengths 2-5. Key findings: margin is consistently null (~0) at all pre-period lengths. Net fixed assets is stable (~-0.08) across all pre-periods. Net patient revenue per bed (~-200 to -219) and operating expenses per bed (~-177 to -208) are remarkably stable — point estimates barely move with pre-period length. Current ratio attenuates from -0.71 (pre=5) to -0.42 (pre=2), losing significance at pre=2. Capex is uninformative at any length. Treated obs: 34 (pre=5) → 43 (pre=2); controls: 139 → 319.

**Decided approach**: Use observed margin data only; HCRIS PPS extension (1985-1999) now provides `net_pat_rev` for ~57% of 1996 records, improving pre-treatment coverage without imputation. Financial outcomes use `pre_period=financial.pre` (set in `_run-analysis.r`) to allow shorter pre-periods than AHA outcomes. Pre-period sensitivity analysis in appendix. Other financial outcomes (current_ratio, net_fixed, capex) reported with caveats about small treated samples.

**Prior approaches explored and abandoned**:
- *HCRIS PPS margin imputation* (`supp-impute-margin.R`, archived 2026-02-09): Two-stage hybrid OLS — Stage 1 translated PPS-measured `margin_charges` to HCRIS-equivalent (1998 training), Stage 2 predicted true margin from HCRIS predictors (1997 training). Abandoned due to poor out-of-sample fit: full-pipeline 1999 holdout RMSE 0.214, correlation 0.24. Stage 1 was the bottleneck — PPS cost-allocation variables weakly predict revenue-statement margins. The fundamental problem is that PPS minimum dataset variables (Medicare-specific costs, cost-allocation charges) are structurally different from HCRIS revenue-statement variables, and no statistical model could reliably bridge this gap. Script remains in repo but is not sourced by either entry point.
- *AHA-based imputation* (archived in `archived/0-impute-financials.R`): Elastic net using AHA predictors + hospital-specific residual correction. Abandoned — AHA operational variables are poor predictors of financial outcomes (cross-sectional R² 0.01-0.39).
- *Manual 990 collection*: Feasible for ~32 nonprofits missing pre-1998 data, but superseded by HCRIS PPS data extension.
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

### Margin Imputation (archived 2026-02-09)

Attempted to extend margin coverage to 1994-1996 using HCRIS PPS data via `analysis/supp-impute-margin.R`. Abandoned due to poor out-of-sample fit — the fundamental problem is that PPS minimum dataset variables (Medicare-specific costs, cost-allocation charges) are structurally different from HCRIS revenue-statement variables needed to compute operating margin.

**Best model achieved**: Two-stage hybrid OLS (Stage 1: PPS→HCRIS measurement bridge; Stage 2: HCRIS predictors→margin). Full-pipeline 1999 holdout: RMSE 0.214, correlation 0.24. Six alternative approaches were also tested (single-stage OLS variants, decomposed calibration, additive correction, ridge, recentering) — none performed adequately. Details preserved in git history (commit `7da3ab8`).

**Current approach**: Use observed margin data only. The HCRIS PPS extension provides `net_pat_rev` for ~57% of 1996 records, improving coverage without imputation. Financial outcomes use shorter pre-periods (`financial.pre` variable) with sensitivity analysis in appendix.

**Exploratory diagnostics retained**: `scratch/diag-hcris-coverage.R`, `scratch/diag-margin-decomposed.R`

### Revenue/Expense Decomposition (2026-02-09)

**New outcomes**: `net_pat_rev` and `tot_operating_exp` added as per-bed, CPI-deflated outcomes (thousands of 2010 dollars per bed, divisor 1e3). Constructed in `_build-estimation-data.r`: winsorized → scaled by beds_base and cpi_deflator → interpolated. Both use `financial.pre` pre-period and run through `2-hospital-dd.R` (outcome-agnostic).

**Motivation**: Decomposing margin into its components reveals whether CAH designation stabilizes hospitals via cost expansion (cost-based reimbursement relaxing cost discipline) or some other mechanism.

**Key results** (from `att_overall.tex`):
- SDID: net_pat_rev ATT = -203 [-583, 177]; tot_operating_exp ATT = -208 [-623, 207]. Both imprecise, CIs include zero.
- CS: net_pat_rev ATT = -522 [-719, -325]; tot_operating_exp ATT = -540 [-731, -348]. Both significant.
- Revenue and expenses decline by nearly identical amounts → explains the null margin result (proportional contraction, not no effect).
- Interpretation: CAH is associated with financial *contraction* (consistent with downsizing to meet 25-bed cap), not cost expansion. Undermines the narrative that cost-based reimbursement inflates spending.

**Pre-trend concerns**: `tot_operating_exp` has significant differential pre-trends in ALL three cohorts (p < 0.002 in each). `net_pat_rev` flags in 1999 (p = 0.003) and marginally in 2001 (p = 0.06); 2000 is clean. However, pre-period sensitivity analysis shows SDID point estimates are stable across pre-period lengths 2-5, suggesting the estimates are not driven by fitting pre-trends.

**Paper/appendix status (2026-02-10)**: paper.tex reflects current results — null margin narrative, revenue/expense decomposition in Section 4.2, 4-panel financial figure (margin, current ratio, revenue, expenses), and current numbers for all outcomes. Organizational changes section reframed around risk reduction rather than margin improvement. appendix.tex includes measure definitions for all six financial outcomes (Section A), expanded pre-trend figures and discussion for revenue/expense, and a new Section D with pre-period sensitivity table. Decomposition framed as descriptive/suggestive given pre-trend concerns. SDID is primary estimator; CS results moving to appendix. Results tables now include $N_{tr}$.

### Forward Plan (2026-02-10)

See `scratch/refreport_202602.md` for detailed progress log and plan.

**Phase 1 — Paper writing** (no new estimation):
- Section 4.X: eligibility-restricted results as robustness subsection, forest plot figure
- Section 5 (Discussion): limitations, policy implications, future directions
- Appendix B: `state.cut` sensitivity (coauthor contributing)
- Expositional: typos, date, $p_c(\delta)$ notation, IPW documentation, $\rho$ justification, CS plots to appendix

**Phase 2 — Additional estimation**:
1. Heterogeneity analysis: ownership (govt vs nonprofit), initial bed size (≤25 vs 26-50), geographic isolation. Existing cohort-specific results are the foundation.
2. Permutation/placebo inference for SDID (addresses small-sample SE concern)
3. gsynth/fect for financial outcomes (handles unbalanced panels + differential trends)
4. ITT specification (state-level treatment for all eligible hospitals)
5. IV feasibility (state availability as instrument)
