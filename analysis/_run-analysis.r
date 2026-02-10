# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   2/24/2026
## Description:   Run Analysis Files
## Note:          Run _build-estimation-data.r first to generate the CSVs


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here, kableExtra,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s, dotwhisker, mlogit, readxl,
               haven, sf, igraph, plotly, synthdid, BMisc, nnet, glmnet, zoo, purrr, grid, rlang, survival)

source('analysis/functions.R')


# Read estimation data -----------------------------------------------------
est.dat   <- read_csv('data/output/estimation_data.csv')
state.dat <- read_csv('data/output/state_estimation_data.csv')


# Global parameters --------------------------------------------------------
bed.cut   <- 50
post      <- 5
state.cut <- 0
financial.pre <- 5


# Stacked datasets ---------------------------------------------------------
stack.hosp  <- stack_hosp(pre.period=5, post.period=post, state.period=state.cut)
stack.state <- stack_state(pre.period=5, post.period=post, state.period=state.cut)


# Main analysis: state-timing identification --------------------------------

## Unified outcome map
outcome_map <- list(

  # Hospital continuous outcomes (cohorts 1999:2001)
  margin        = list(script="analysis/2-hospital-dd.R", label="Operating margin",             stub="margin",       cohorts=1999:2001, pre_period=financial.pre),
  current_ratio = list(script="analysis/2-hospital-dd.R", label="Current ratio",                stub="currentratio", cohorts=1999:2001, pre_period=financial.pre),
  net_fixed     = list(script="analysis/2-hospital-dd.R", label="Net fixed assets",             stub="netfixed",     cohorts=1999:2001, pre_period=financial.pre),
  capex         = list(script="analysis/2-hospital-dd.R", label="Capital expenditures per bed", stub="capex",        cohorts=1999:2001, pre_period=financial.pre),
  net_pat_rev       = list(script="analysis/2-hospital-dd.R", label="Net patient revenue per bed",  stub="netpatrev",  cohorts=1999:2001, pre_period=financial.pre),
  tot_operating_exp = list(script="analysis/2-hospital-dd.R", label="Operating expenses per bed",   stub="totopexp",   cohorts=1999:2001, pre_period=financial.pre),
  BDTOT         = list(script="analysis/2-hospital-dd.R", label="Total beds",                   stub="beds",         cohorts=1999:2001),
  OBBD          = list(script="analysis/2-hospital-dd.R", label="OB beds",                      stub="beds_ob",      cohorts=1999:2001),
  FTERN         = list(script="analysis/2-hospital-dd.R", label="FTE RNs",                      stub="ftern",        cohorts=1999:2001),
  ip_per_bed    = list(script="analysis/2-hospital-dd.R", label="Inpatient days per bed",       stub="ipdays",       cohorts=1999:2001),
  system        = list(script="analysis/2-hospital-dd.R", label="System membership",            stub="system",       cohorts=1999:2001),

  # State-level count outcomes (cohorts 1999:2001)
  closures = list(script="analysis/4-changes-state-dd.R", label="Closures", stub="closure-rate", cohorts=1999:2001),
  mergers  = list(script="analysis/4-changes-state-dd.R", label="Mergers",  stub="merger-rate",  cohorts=1999:2001)
)

## Loop over outcomes, collect into results table
results.table <- tibble(
  outcome = character(), sdid_att = numeric(), sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  cs_att = numeric(), cs_ci_low = numeric(), cs_ci_high = numeric()
)

for (oname in names(outcome_map)) {
  print(paste0("Running analysis for outcome: ", oname))
  o <- outcome_map[[oname]]
  outcome_var   <- oname
  outcome_sym   <- sym(outcome_var)
  outcome_label <- o$label
  file_stub     <- o$stub
  cohorts       <- o$cohorts
  pre_period    <- if (!is.null(o$pre_period)) o$pre_period else 5

  source(o$script)

  # CS results — apply scale_cs for state-level outcomes
  cs_att_val <- csa.att$overall.att
  cs_se_val  <- csa.att$overall.se
  if (o$script == "analysis/4-changes-state-dd.R") {
    cs_att_val <- cs_att_val * scale_cs
    cs_se_val  <- cs_se_val * scale_cs
  }

  results.table <- bind_rows(results.table, tibble(
    outcome      = outcome_label,
    sdid_att     = as.numeric(att_w),
    sdid_ci_low  = as.numeric(ci_low),
    sdid_ci_high = as.numeric(ci_high),
    cs_att       = cs_att_val,
    cs_ci_low    = cs_att_val - 1.96 * cs_se_val,
    cs_ci_high   = cs_att_val + 1.96 * cs_se_val
  ))
}

## LaTeX output (tabular innards only)
int <- function(lo, hi) sprintf("[%.2f, %.2f]", lo, hi)

tex.lines <- results.table %>%
  mutate(line = sprintf("%s & %.3f & %s & %.3f & %s \\\\",
    outcome, sdid_att, int(sdid_ci_low, sdid_ci_high),
    cs_att, int(cs_ci_low, cs_ci_high))) %>%
  pull(line)

writeLines(c(
  "Outcome & SDID ATT & SDID 95\\% CI & CS ATT & CS 95\\% CI \\\\",
  tex.lines
), "results/att_overall.tex")


# Alternative identification (eligibility-restricted, no state requirement) -
source('analysis/3-hospital-dd-alt.R')


# Diagnostics and sensitivity ----------------------------------------------
source('analysis/app-dd-diagnostics.R')
source('analysis/app-financial-preperiod.R')
