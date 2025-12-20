# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   12/16/2025
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here, kableExtra,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s, dotwhisker, mlogit, readxl,
               haven, sf, igraph, plotly, synthdid, BMisc, nnet, glmnet, zoo, purrr, grid, rlang, survival)

source('analysis/functions.R')

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

# construct new variables
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

hosp.beds <- final.dat %>%
  filter(year < eff_year | is.na(eff_year)) %>%
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
            is.na(net_fixed_1) & is.na(net_fixed_2) & !is.na(net_fixed_3) ~ net_fixed_3)
          ) %>%
  mutate(
    margin_hcris=(net_pat_rev - tot_operating_exp)/net_pat_rev,
    margin=ifelse(!is.na(margin_990), margin_990, margin_hcris),
    net_fixed=ifelse(!is.na(net_fixed_990), net_fixed_990, fixed_assets),
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
  group_by(year) %>% 
  mutate(m_top=quantile(margin, probs=0.95, na.rm=TRUE),
         m_bottom=quantile(margin, probs=0.05, na.rm=TRUE),
         fa_top=quantile(net_fixed, probs=0.95, na.rm=TRUE),
         fa_bottom=quantile(net_fixed, probs=0.05, na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate(margin=ifelse(margin>m_top, m_top, margin),
         margin=ifelse(margin<m_bottom, m_bottom, margin),
         net_fixed=ifelse(net_fixed>fa_top, fa_top, net_fixed),
         net_fixed=ifelse(net_fixed<fa_bottom, fa_bottom, net_fixed),
         net_fixed=net_fixed/1000000,
         net_fixed=net_fixed/beds_base,
         net_fixed=net_fixed/cpi_deflator,
         margin=margin/cpi_deflator) %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  fill(state, .direction = "downup") %>%    
  mutate(margin=na.approx(margin, x=year, na.rm=FALSE),
         net_fixed=na.approx(net_fixed, x=year, na.rm=FALSE)) %>%
  ungroup()


## add ipw weights to estimation data
source('analysis/0-ipw-weights.R')
est.dat <- est.dat %>%
  left_join(id.weights %>% select(ID, ipw, ever_rural), by="ID")


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


# Source analysis code files -----------------------------------------------

source('analysis/1-sum-stats.R')
source('analysis/2-margins.R')
source('analysis/3-closures-mergers.R')
