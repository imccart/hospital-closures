# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   8/12/2025
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s,
               haven, sf, igraph, plotly, synthdid, BMisc, nnet, glmnet, zoo)

# Read-in data ------------------------------------------------------------
aha.data <- read_csv('data/output/aha_final.csv')
aha.neighbors <- read_csv('data/output/aha_neighbors.csv')
cah.dates <- read_csv('data/input/cah-states.csv')

non.missing.counts <- aha.data %>%
  group_by(year) %>%
  summarise(across(everything(), ~ sum(!is.na(.) & !(. %in% "")))) %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Count") %>%
  pivot_wider(names_from = year, values_from = Count)

check.distance <- aha.neighbors %>%
  group_by(year) %>%
  summarize(mean_dist=mean(distance, na.rm=TRUE),
            count_dist=sum(!is.na(distance)),
            count_all=n())

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


# Merge and finalize hospital-level data -----------------------------------

data.merge <- aha.data %>% 
    left_join(aha.neighbors, by=c("ID", "year")) %>%
    left_join(cah.dates %>% select(cah_date_law=cah_date, cah_year_law=cah_year, MSTATE="abb"), by="MSTATE") %>%
    filter(! MSTATE %in% c("AK","HI","PR","VI","GU","MP","AS", "N","0", "AS", "DC", "DE", "MH", "ML"), 
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

est.dat <- final.dat %>%
  mutate(margin_990=case_when(
            !is.na(margin_1) ~ margin_1,
            is.na(margin_1) & !is.na(margin_2) ~ margin_2,
            is.na(margin_1) & is.na(margin_2) & !is.na(margin_3) ~ margin_3)
          ) %>%
  mutate(
    margin_hcris=(net_pat_rev - tot_operating_exp)/net_pat_rev,
    margin=ifelse(!is.na(margin_990), margin_990, margin_hcris),
    state_event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    hosp_event_time=case_when(
      !is.na(eff_year) ~ year - eff_year,
      TRUE ~ -1),
    aha_id=as.numeric(ID),
    treat_state=ifelse(state_treat_year>0, 1, 0),
    treat_state_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0),
    compare_hosp=case_when(
      year<1999 & BDTOT<12 & distance>30 ~ 1,
      year>=1999 & BDTOT<25 & distance>30 ~ 1,
      TRUE ~ 0)) %>%
    group_by(year) %>% 
    mutate(m_top=quantile(margin, probs=0.95, na.rm=TRUE),
           m_bottom=quantile(margin, probs=0.05, na.rm=TRUE)) %>%
    ungroup() %>% 
    mutate(margin=ifelse(margin>m_top, m_top, margin),
           margin=ifelse(margin<m_bottom, m_bottom, margin)) %>%
    arrange(ID, year) %>%
    group_by(ID) %>%
    mutate(margin=na.approx(margin, x=year, na.rm=FALSE)) %>%
    ungroup()

## add ipw weights to estimation data
source('analysis/0-ipw-weights.R')
est.dat <- est.dat %>%
  left_join(id.weights %>% select(ID, ipw), by="ID")


## Aggregate to state level --------------------------------------------------

## state.dat1 is the raw count of hospitals, closures, and mergers
state.dat1 <- est.dat %>% 
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    state_treat_year=first(state_treat_year)) %>%
  mutate(
    event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    treat=ifelse(state_treat_year>0, 1, 0),
    treat_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

## state.dat2 is the same, but with ipw weights applied
state.dat2 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=sum(ipw), cah_treat=min(state_treat_year), 
    closures=sum(closed*ipw), sum_cah=sum(cah*ipw, na.rm=TRUE),
    mergers=sum(merged*ipw), any_changes=sum(closed_merged*ipw),
    state_treat_year=first(state_treat_year)) %>%
  mutate(
    event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    treat=ifelse(state_treat_year>0, 1, 0),
    treat_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

## state.dat3 is raw counts but with distance and size restrictions
state.dat3 <- est.dat %>% filter(BDTOT<30, distance>10) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    state_year_treat=first(state_treat_year)) %>%
  mutate(
    event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    treat=ifelse(state_treat_year>0, 1, 0),
    treat_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))


# Source analysis code files -----------------------------------------------

source('analysis/1-sum-stats.R')
source('analysis/2-margins.R')
source('analysis/3-closures-mergers.R')