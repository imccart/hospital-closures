# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   7/17/2023
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR)

# Read-in data ------------------------------------------------------------
aha.data <- read_csv('data/output/aha_final.csv')
aha.neighbors <- read_csv('data/output/aha_neighbors.csv')


# Source analysis code files -----------------------------------------------

source('analysis/sum-stats.R')
source('analysis/matching.R')