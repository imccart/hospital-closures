# Analysis

This page describes our analysis of the data.

## Overview
The [analysis](/analysis/) folder contains the code used to analyze the data. The code is organized into three files:

- **_run-analysis.R**: This file loads the data from the data folder and cleans and merges the data. It then calls the **sum-stats.R** code and the **dd-estimates.R** code. It takes as input the 'aha_data' and 'cah_data' created by the data code.