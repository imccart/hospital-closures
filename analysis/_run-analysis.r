# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   2/24/2026
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here, kableExtra,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s, dotwhisker, mlogit, readxl,
               haven, sf, igraph, plotly, synthdid, BMisc, nnet, glmnet, zoo, purrr, grid, rlang, survival)

source('analysis/functions.R')

## Global parameters for estimation and stacking
bed.cut   <- 50
post      <- 5
state.cut <- 0
pre_period <- 5

## Winsorization function
winsor_by_year <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

# Read-in data ------------------------------------------------------------
aha.data <- read_csv('data/output/aha_final.csv')
aha.neighbors <- read_csv('data/output/aha_neighbors.csv')
cah.dates <- read_csv('data/input/cah-states.csv')

state.xwalk <- tibble(
  state = state.name,
  MSTATE = state.abb
)

copa.dat <- read_csv('data/input/copa-states.csv', col_names=c("state", "empty", "copa_status", "lat_lon", "details")) %>%
  separate(details, sep=",", into=c("state2","statute","status","other")) %>%
  mutate(
    copa_start = str_extract(status, "(?<=Enacted )\\d{4}") %>% as.integer(),
    copa_end = str_extract(status, "(?<=Repealed )\\d{4}"),
    copa_end = replace_na(copa_end, "2025") %>% as.integer()
  ) %>%
  mutate(
    copa_start = case_when(
      state == "Maine" ~ 2005,
      state == "New York" ~ 2014,
      TRUE ~ copa_start
    )
  ) %>%
  select(state, copa_status, copa_start, copa_end) %>%
  filter(!is.na(state)) %>%
  left_join(state.xwalk, by="state")

write_csv(copa.dat,'data/output/copa_data.csv')

cpi.data <- read_xlsx("data/input/CPI_1913_2019.xlsx", skip = 11)
cpi.data <- pivot_longer(cpi.data,
                         cols=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                         names_to="month",
                         values_to="index") %>%
  mutate(year=Year) %>%
  group_by(year) %>%
  summarize(index=mean(index, na.rm=TRUE))

cpi.2010 <- cpi.data %>%
  filter(year==2010) %>%
  select(index) %>% as.numeric()

cpi.final <- cpi.data %>%
  mutate(cpi_2010 = cpi.2010,
         cpi_deflator = index/cpi_2010) %>%
  select(year, cpi_deflator)


# Merge and finalize hospital-level data -----------------------------------

data.merge <- aha.data %>%
    left_join(aha.neighbors, by=c("ID", "year")) %>%
    left_join(cah.dates %>% select(cah_date_law=cah_date, cah_year_law=cah_year, MSTATE="abb"), by="MSTATE") %>%
    filter(! MSTATE %in% c("AK","HI","PR","VI","GU","MP","AS", "N","0", "AS", "DC", "DE", "MH", "ML","MD"),
           !is.na(MSTATE), MSTATE!="NA") %>%
    filter(COMMTY=="Y", hosp_type=="General") %>%
    group_by(ID, year) %>%
    mutate(hosp_count=n()) %>%
    filter(hosp_count==1) %>% ungroup () %>%
    select(-hosp_count)

## construct new variables
final.dat <- data.merge %>%
  mutate(closed=case_when(
                !is.na(change_type) & change_type=="Closure" ~ 1,
                !is.na(change_type) & change_type!="Closure" ~ 0,
                is.na(change_type) ~ 0),
         merged=case_when(
                !is.na(change_type) & change_type=="Merger" ~ 1,
                !is.na(change_type) & change_type!="Merger" ~ 0,
                is.na(change_type) ~ 0),
         closed_merged=ifelse(closed==1 | merged==1, 1, 0),
         ) %>%
    group_by(ID) %>% mutate(ever_cah=max(cah, na.rm=TRUE)) %>% ungroup() %>%
    group_by(MSTATE) %>% mutate(state_first_obs=min(eff_year, na.rm=TRUE),
                                state_first_law=min(cah_year_law, na.rm=TRUE)) %>% ungroup() %>%
    mutate(state_first_obs=ifelse(is.infinite(state_first_obs), 0, state_first_obs),
           state_first_law=ifelse(is.infinite(state_first_law), 0, state_first_law),
           state_treat_year=state_first_obs)

## mean hospital beds pre-treatment
hosp.beds <- final.dat %>%
  filter(year < eff_year | is.na(eff_year), year>1995) %>%
  group_by(ID) %>%
  summarize(beds_base=mean(BDTOT, na.rm=TRUE))

est.dat <- final.dat %>%
  left_join(hosp.beds, by="ID") %>%
  left_join(cpi.final , by="year") %>%
  mutate(margin_990=case_when(
            !is.na(margin_1) ~ margin_1,
            is.na(margin_1) & !is.na(margin_2) ~ margin_2,
            is.na(margin_1) & is.na(margin_2) & !is.na(margin_3) ~ margin_3),
         net_fixed_990=case_when(
            !is.na(net_fixed_1) ~ net_fixed_1,
            is.na(net_fixed_1) & !is.na(net_fixed_2) ~ net_fixed_2,
            is.na(net_fixed_1) & is.na(net_fixed_2) & !is.na(net_fixed_3) ~ net_fixed_3),
         depreciation_990=case_when(
            !is.na(depreciation_1) ~ depreciation_1,
            is.na(depreciation_1) & !is.na(depreciation_2) ~ depreciation_2,
            is.na(depreciation_1) & is.na(depreciation_2) & !is.na(depreciation_3) ~ depreciation_3),
         current_ratio_990=case_when(
            !is.na(current_ratio_1) ~ current_ratio_1,
            is.na(current_ratio_1) & !is.na(current_ratio_2) ~ current_ratio_2,
            is.na(current_ratio_1) & is.na(current_ratio_2) & !is.na(current_ratio_3) ~ current_ratio_3)
          ) %>%
  mutate(
    ip_per_bed=IPDTOT/beds_base,
    gross_fixed_hcris = fixed_assets + accum_dep,
    current_ratio_hcris=current_assets/current_liabilities,
    margin_hcris=(net_pat_rev - tot_operating_exp)/net_pat_rev,
    margin=ifelse(!is.na(margin_hcris), margin_hcris, margin_990),
    net_fixed=ifelse(!is.na(fixed_assets), fixed_assets, net_fixed_990),
    current_ratio=ifelse(!is.na(current_ratio_hcris), current_ratio_hcris, current_ratio_990),
    state_event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    hosp_event_time=case_when(
      !is.na(eff_year) ~ year - eff_year,
      TRUE ~ -1),
    aha_id=as.numeric(ID),
    treat_state=ifelse(state_treat_year>0, 1, 0),
    treat_state_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0)
  ) %>%
  group_by(ID) %>%
  arrange(year, .by_group=TRUE) %>%
  mutate(capex_990=net_fixed_990 - lag(net_fixed_990) + depreciation_990,
         capex_hcris=gross_fixed_hcris - lag(gross_fixed_hcris)) %>%
  ungroup() %>%
  mutate(capex=ifelse(!is.na(capex_hcris), capex_hcris, capex_990)) %>%
  group_by(year, ever_cah) %>%
  mutate(
    margin        = winsor_by_year(margin),
    net_fixed     = winsor_by_year(net_fixed),
    current_ratio = winsor_by_year(current_ratio),
    capex         = winsor_by_year(capex),
    BDTOT         = winsor_by_year(BDTOT),
    OBBD          = winsor_by_year(OBBD),
    FTERN         = winsor_by_year(FTERN),
    IPDTOT        = winsor_by_year(IPDTOT),
    ip_per_bed    = if_else(ip_per_bed > (365*.85), 365*.85, ip_per_bed)
  ) %>%
  ungroup() %>%
  mutate(
    net_fixed    = net_fixed / 1e6 / beds_base / cpi_deflator,
    capex        = capex / 1e4 / beds_base / cpi_deflator
  ) %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  fill(state, .direction = "downup") %>%
  mutate(margin=na.approx(margin, x=year, na.rm=FALSE),
         net_fixed=na.approx(net_fixed, x=year, na.rm=FALSE),
         current_ratio=na.approx(current_ratio, x=year, na.rm=FALSE),
         capex=na.approx(capex, x=year, na.rm=FALSE),
         BDTOT = na.approx(BDTOT, x=year, na.rm=FALSE),
         OBBD = na.approx(OBBD, x=year, na.rm=FALSE),
         FTERN = na.approx(FTERN, x=year, na.rm=FALSE),
         IPDTOT = na.approx(IPDTOT, x=year, na.rm=FALSE),
         ip_per_bed = na.approx(ip_per_bed, x=year, na.rm=FALSE)) %>%
  ungroup()


## add ipw weights to estimation data
source('analysis/0-ipw-weights.R')
est.dat <- est.dat %>%
  left_join(id.weights %>% select(ID, ipw, ever_rural), by="ID")

## check missing values by year
missing.dat <- est.dat %>%
  group_by(year) %>%
  summarize(
    n=n(),
    missing_margin=sum(is.na(margin)),
    missing_net_fixed=sum(is.na(net_fixed)),
    missing_current_ratio=sum(is.na(current_ratio)),
    missing_capex=sum(is.na(capex))
  )

# Aggregate to state level --------------------------------------------------

state.dat <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year),
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), changes=sum(closed_merged),
    state_treat_year=first(state_treat_year),
    mean_beds=mean(BDTOT, na.rm=TRUE), mean_distance=mean(distance, na.rm=TRUE)) %>%
  group_by(MSTATE) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(hospitals_lag = lag(hospitals)) %>%
  ungroup() %>%
  mutate(
    event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    treat=ifelse(state_treat_year>0, 1, 0),
    treat_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0),
    state=as.numeric(factor(MSTATE)),
    rate_closed=closures/hospitals_lag,
    rate_merged=mergers/hospitals_lag,
    rate_changes=changes/hospitals_lag)

# source('analysis/supp-extract-manual.R')

# Source analysis code files -----------------------------------------------
write_csv(est.dat, 'data/output/estimation_data.csv')
write_csv(state.dat, 'data/output/state_estimation_data.csv')

source('analysis/1-sum-stats.R')

stack.hosp  <- stack_hosp(pre.period=5, post.period=post, state.period=state.cut)
stack.state <- stack_state(pre.period=5, post.period=post, state.period=state.cut)

## impute missing margin values for very early years
source('analysis/supp-impute-margin.R')

## Unified outcome map
outcome_map <- list(

  # Hospital continuous outcomes (cohorts 1999:2001)
  margin        = list(script="analysis/2-hospital-dd.R", label="Operating margin",             stub="margin",       cohorts=1999:2001),
  current_ratio = list(script="analysis/2-hospital-dd.R", label="Current ratio",                stub="currentratio", cohorts=1999:2001),
  net_fixed     = list(script="analysis/2-hospital-dd.R", label="Net fixed assets",             stub="netfixed",     cohorts=1999:2001),
  capex         = list(script="analysis/2-hospital-dd.R", label="Capital expenditures per bed", stub="capex",        cohorts=1999:2001),
  BDTOT         = list(script="analysis/2-hospital-dd.R", label="Total beds",                   stub="beds",         cohorts=1999:2001),
  OBBD          = list(script="analysis/2-hospital-dd.R", label="OB beds",                      stub="beds_ob",      cohorts=1999:2001),
  FTERN         = list(script="analysis/2-hospital-dd.R", label="FTE RNs",                      stub="ftern",        cohorts=1999:2001),
  ip_per_bed    = list(script="analysis/2-hospital-dd.R", label="Inpatient days per bed",       stub="ipdays",       cohorts=1999:2001),
  system        = list(script="analysis/2-hospital-dd.R", label="System membership",            stub="system",       cohorts=1999:2001),

  # State-level count outcomes (cohorts 1999:2001)
  closures = list(script="analysis/3-changes-state-dd.R", label="Closures", stub="closure-rate", cohorts=1999:2001),
  mergers  = list(script="analysis/3-changes-state-dd.R", label="Mergers",  stub="merger-rate",  cohorts=1999:2001)
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
  pre_period    <- if (!is.null(o$pre_period)) o$pre_period else 5  # default 5, margin uses 3

  source(o$script)

  # CS results — apply scale_cs for state-level outcomes
  cs_att_val <- csa.att$overall.att
  cs_se_val  <- csa.att$overall.se
  if (o$script == "analysis/4-changes-state-dd.r") {
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

## Diagnostics ------------------------------------------------------------
source('analysis/app-dd-diagnostics.R')