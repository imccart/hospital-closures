# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  2/18/2026
## Date Edited:   2/19/2026
## Description:   Activates renv and loads packages. Sourced by
##                _build-estimation-data.r and _run-analysis.r.

# renv activation ---------------------------------------------------------
source("renv/activate.R")

# Packages ----------------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(janitor)
library(here)
library(zoo)
library(fedmatch)
library(zipcodeR)
library(fixest)
library(did)
library(did2s)
library(BMisc)
library(fect)
library(glmnet)
library(nnet)
library(mlogit)
library(survival)
library(scales)
library(plotly)
library(panelView)
library(dotwhisker)
library(patchwork)
library(sf)
library(modelsummary)
library(kableExtra)
library(broom)
library(synthdid)
