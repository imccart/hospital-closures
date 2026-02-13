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
stack.elig  <- stack_hosp_elig(pre.period=5, post.period=post, cohort.years=1999:2005)


# Main analysis ------------------------------------------------------------
source('analysis/2-hospital-dd.R')        # -> hosp.results.table, hosp.cohort.results
source('analysis/3-hospital-dd-alt.R')     # -> elig.results, elig.cohort.results
source('analysis/4-changes-state-dd.R')    # -> state.results.table


# Combined results table (hospital + state) --------------------------------
results.table <- bind_rows(hosp.results.table, state.results.table)

int <- function(lo, hi) sprintf("[%.2f, %.2f]", lo, hi)

## Main table (SDID only) for paper — complete tabular block
## (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
tex.lines.sdid <- results.table %>%
  mutate(line = sprintf("%s & %.3f & %s & %d \\\\",
    outcome, sdid_att, int(sdid_ci_low, sdid_ci_high), sdid_ntr)) %>%
  pull(line)

writeLines(c(
  "\\begin{tabular}[t]{lclr}",
  "Outcome & SDID ATT & SDID 95\\% CI & $N_{tr}$ \\\\",
  tex.lines.sdid,
  "\\end{tabular}"
), "results/att_overall.tex")

## CS table for appendix — complete tabular block
tex.lines.cs <- results.table %>%
  mutate(line = sprintf("%s & %.3f & %s \\\\",
    outcome, cs_att, int(cs_ci_low, cs_ci_high))) %>%
  pull(line)

writeLines(c(
  "\\begin{tabular}{lcl}",
  "Outcome & CS ATT & CS 95\\% CI \\\\",
  tex.lines.cs,
  "\\end{tabular}"
), "results/att_cs_overall.tex")


# Forest plot: SDID results across control group constructions ---------------

forest_labels <- c(
  "Operating margin" = "Operating\nmargin",
  "Current ratio" = "Current ratio",
  "Net fixed assets" = "Net fixed assets",
  "Capital expenditures per bed" = "Capital expenditures\nper bed",
  "Net patient revenue per bed" = "Net patient revenue\nper bed",
  "Operating expenses per bed" = "Operating expenses\nper bed",
  "Total beds" = "Total beds",
  "OB beds" = "OB beds",
  "FTE RNs" = "FTE RNs",
  "Inpatient days per bed" = "Inpatient days\nper bed",
  "System membership" = "System\nmembership",
  "Closures" = "Closures",
  "Mergers" = "Mergers"
)

forest_dat <- bind_rows(
  results.table %>%
    transmute(outcome, controls = "State-timing",
              att = sdid_att, ci_low = sdid_ci_low, ci_high = sdid_ci_high),
  elig.results %>%
    transmute(outcome, controls = "Eligibility-restricted",
              att = sdid_att, ci_low = sdid_ci_low, ci_high = sdid_ci_high)
) %>%
  mutate(outcome_label = forest_labels[outcome],
         outcome_label = factor(outcome_label,
                                levels = forest_labels[unique(results.table$outcome)]),
         controls = factor(controls, levels = c("Eligibility-restricted", "State-timing")))

n_out <- length(unique(forest_dat$outcome_label))
p_forest <- ggplot(forest_dat, aes(x = att, y = controls, shape = controls)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "gray40") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  size = 0.4, linewidth = 0.4) +
  scale_shape_manual(values = c("State-timing" = 16, "Eligibility-restricted" = 1)) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x", strip.position = "left",
             dir = "v") +
  labs(x = "ATT (95% CI)", y = NULL, shape = "Controls") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                     lineheight = 0.9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )
ggsave("results/forest-sdid.png", p_forest,
       width = 7, height = 0.9 * n_out + 1.2, dpi = 300)


# Heterogeneity ------------------------------------------------------------
source('analysis/6-heterogeneity.R')         # -> het.results


# Diagnostics and sensitivity ----------------------------------------------
source('analysis/app-dd-diagnostics.R')
source('analysis/app-financial-preperiod.R')
