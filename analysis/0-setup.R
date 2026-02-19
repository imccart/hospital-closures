# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  2/18/2026
## Date Edited:   2/18/2026
## Description:   Loads packages at fixed versions via groundhog for
##                reproducibility. Sourced by _build-estimation-data.r
##                and _run-analysis.r.

# Groundhog setup ---------------------------------------------------------

options(repos = c(CRAN = "https://cran.r-project.org/"))
if (!require("groundhog", quietly = TRUE)) {
  install.packages("groundhog")
  library(groundhog)
}

# R 4.5.x CRAN window: 2025-04-11 to ~2026-04
ghog_date <- "2025-10-15"

# Packages ----------------------------------------------------------------

pkgs <- c(
  "tidyverse", "haven", "readxl", "janitor", "here", "zoo",
  "fedmatch", "zipcodeR",
  "fixest", "did", "did2s", "BMisc", "fect",
  "glmnet", "nnet", "mlogit", "survival",
  "scales", "plotly", "panelView", "dotwhisker", "patchwork",
  "sf", "igraph",
  "modelsummary", "kableExtra", "broom"
)
suppressWarnings(groundhog.library(pkgs, ghog_date))
if (!require("synthdid", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("synth-inference/synthdid")
}
library(synthdid)
