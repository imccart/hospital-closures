# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   2/5/2025
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s,
               haven, sf, igraph, plotly, synthdid, BMisc)

# Read-in data ------------------------------------------------------------
aha.data <- read_csv('data/output/aha_final.csv')
aha.neighbors <- read_csv('data/output/aha_neighbors.csv')
cah.dates <- read_csv('data/input/cah-states.csv')

non.missing.counts <- aha.data %>%
  group_by(year) %>%
  summarise(across(everything(), ~ sum(!is.na(.) & !(. %in% "")))) %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Count") %>%
  pivot_wider(names_from = year, values_from = Count)

    
# Merge and finalize hospital-level data -----------------------------------

data.merge <- aha.data %>% 
    left_join(aha.neighbors, by=c("ID", "year")) %>%
    left_join(cah.dates %>% select(cah_date_law=cah_date, cah_year_law=cah_year, MSTATE="abb"), by="MSTATE") %>%
    filter(! MSTATE %in% c("AK","HI","PR","VI","GU","MP","AS", "N","0", "AS", "DC", "DE", "MH", "ML"), !is.na(MSTATE)) %>%
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
    group_by(MSTATE) %>% mutate(first_year_obs=min(eff_year, na.rm=TRUE),
                                first_year_law=min(cah_year_law, na.rm=TRUE)) %>% ungroup() %>%
    mutate(first_year_obs=ifelse(is.infinite(first_year_obs), 0, first_year_obs),
           first_year_law=ifelse(is.infinite(first_year_law),0,first_year_law),
           first_year_treat=first_year_obs) %>%
    group_by(ID) %>% mutate(max_treat=max(first_year_treat, na.rm=TRUE), min_treat=min(first_year_treat, na.rm=TRUE)) %>%
    filter(max_treat==min_treat) %>% ungroup() %>% select(-c(max_treat, min_treat))


est.dat <- final.dat %>%
  mutate(
    event_time=case_when(
      first_year_treat>0 ~ year - first_year_treat,
      first_year_treat==0 ~ -1),
    aha_id=as.numeric(ID),
    treat_state=ifelse(first_year_treat>0, 1, 0),
    treat_time=ifelse(year>=first_year_treat & first_year_treat>0, 1, 0),
    compare_hosp=case_when(
      year<1999 & BDTOT<12 & distance>30 ~ 1,
      year>=1999 & BDTOT<25 & distance>30 ~ 1,
      TRUE ~ 0))

## add ipw weights to estimation data
source('analysis/0-ipw-weights.R')
est.dat <- est.dat %>%
  left_join(id.weights %>% select(ID, ipw), by="ID")


## Aggregate to state level --------------------------------------------------

state.dat1 <- est.dat %>% 
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_treat), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    first_year_treat=first(first_year_treat)) %>%
  mutate(
    event_time=case_when(
      first_year_treat>0 ~ year - first_year_treat,
      first_year_treat==0 ~ -1),
    treat_state=ifelse(first_year_treat>0, 1, 0),
    treat_time=ifelse(year>=first_year_treat & first_year_treat>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat2 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=sum(ipw), cah_treat=min(first_year_treat), 
    closures=sum(closed*ipw), sum_cah=sum(cah*ipw, na.rm=TRUE),
    mergers=sum(merged*ipw), any_changes=sum(closed_merged*ipw),
    first_year_treat=first(first_year_treat)) %>%
  mutate(
    event_time=case_when(
      first_year_treat>0 ~ year - first_year_treat,
      first_year_treat==0 ~ -1),
    treat_state=ifelse(first_year_treat>0, 1, 0),
    treat_time=ifelse(year>=first_year_treat & first_year_treat>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat3 <- est.dat %>% filter(BDTOT<75, distance>20) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_treat), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    first_year_treat=first(first_year_treat)) %>%
  mutate(
    event_time=case_when(
      first_year_treat>0 ~ year - first_year_treat,
      first_year_treat==0 ~ -1),
    treat_state=ifelse(first_year_treat>0, 1, 0),
    treat_time=ifelse(year>=first_year_treat & first_year_treat>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))


# Source analysis code files -----------------------------------------------

source('analysis/1-sum-stats.R')
source('analysis/2-margins.R')
source('analysis/3-closures-mergers.R')