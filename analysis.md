# Analysis

This page describes our analysis of the data.

## Overview
The [analysis](/analysis/) folder contains the code used to analyze the data. The code is organized into three files:

- **_run-analysis.R**: This file loads the data from the data folder and cleans and merges the data. It takes as input the 'aha_data' and 'cah_data' created by the data code.
- **0-ipw-weights.R**: This file calculates inverse probability weights (IPW) for the data. It takes as input the cleaned and merged data created by _run-analysis.R. These weights are not used in the main analysis, but are included for completeness.
- **1-sum-stats.R**: This file generates summary statistics and tables for the analysis. It takes as input the cleaned and merged data created by `_run-analysis`.
- **2-hospital-dd.R**: This file performs the main difference-in-differences analysis for hospital-level outcomes. It takes as input the cleaned and merged data created by `_run-analysis`. The file also calls the helper functions defined in the `functions.R` script, in particular for generating a stacked data set by treated cohort.
- **3-changes-state.dd.R**: This file performs the difference-in-differences analysis for state-level outcomes. It takes as input the cleaned and merged data created by `_run-analysis` The file again calls the helper functions defined in the `functions.R` script, but this time for generating the stacked data at the state level.
- **4-changes-hospital-hazard.R**: This file estimates discrete-time hazard rates for closure at the hospital level. It takes as input the cleaned and merged data created by `_run-analysis`.