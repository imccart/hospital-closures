# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   8/12/2025
## Description:   Run Analysis Files


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR, did, fixest, panelView, did2s, dotwhisker, mlogit,
               haven, sf, igraph, plotly, synthdid, BMisc, nnet, glmnet, zoo, purrr, grid, rlang, survival)

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


# Aggregate to state level --------------------------------------------------

## state.dat1 is the raw count of hospitals, closures, and mergers
state.dat1 <- est.dat %>% 
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    state_treat_year=first(state_treat_year)) %>%
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
    rate_changes=any_changes/hospitals_lag)

## state.dat2 is the same but with distance and size restrictions
state.dat2 <- est.dat %>%
  group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE), max_distance=max(distance, na.rm=TRUE)) %>% ungroup() %>%  
  filter(min_bedsize<=50) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    state_treat_year=first(state_treat_year)) %>%
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
    rate_changes=any_changes/hospitals_lag)


# Stacked data (treatment at hospital level) ----------------------------

## loop over all possible values of eff_year (year of CAH designation)
post.period <- 5
pre.period <- 5
stack.hosp <- tibble()
for (i in unique(est.dat$eff_year)) {

  ## define treated group relative to eff_year i within event window
  treat.dat <- est.dat %>%
    filter(eff_year==i,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="treated")

  ## define control group relative to eff_year i within event window
  ## never treated (+ never a CAH designation in the state)
  control.dat.never <- est.dat %>% 
    filter(is.na(eff_year), state_treat_year==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="never")

  ## not yet treated (+ no CAH designation *yet* in the state)
  control.dat.notyet <- est.dat %>% 
    filter( (eff_year>(i+post.period) | is.na(eff_year)),
            state_treat_year>(i+post.period),
            year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="notyet")
  
  ## inner join back to est.dat
  stack.dat.group <- est.dat %>% 
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    mutate(state=as.numeric(factor(MSTATE)),
      stacked_event_time=year-i,
      stack_group=i)

  stack.hosp <- bind_rows(stack.hosp, stack.dat.group)
 
}

## drop the eff_year==0 phantom treatment group
stack.hosp <- stack.hosp %>% ungroup() %>%
  filter(stack_group>0) %>%
  mutate(treated=ifelse(treat_type=="treated",1,0),
         control_any=ifelse(treat_type!="treated",1,0),
         control_notyet=ifelse(treat_type=="notyet",1,0),
         control_never=ifelse(treat_type=="never",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated,
         all_treated=sum(treated),
         all_control=sum(control_any),
         all_control_notyet=sum(control_notyet)) %>%
  group_by(stack_group) %>%
  mutate(group_treated=sum(treated),
         group_control_any=sum(control_any),
         group_control_notyet=sum(control_notyet)) %>%
  ungroup() %>%
  mutate(weight_any=case_when(
           treat_type=="treated" ~ 1,
           treat_type!="treated" ~ (group_treated/all_treated)/(group_control_any/all_control)),
         weight_notyet=case_when(
           treat_type=="treated" ~ 1,
           treat_type=="notyet" ~ (group_treated/all_treated)/(group_control_notyet/all_control_notyet))
  )


# Stacked data (treatment at state level) --------------------------------

## loop over all possible values of state_treat_year (year of CAH availability)
post.period <- 5
pre.period <- 5
stack.state1 <- tibble()
stack.state2 <- tibble()
state.count <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n()) %>%
  group_by(MSTATE) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(hospitals_lag = lag(hospitals)) %>%
  ungroup() %>%
  select(MSTATE, year, hospitals, hospitals_lag)

for (i in unique(est.dat$state_treat_year)) {
  ## define treated group relative to state_treat_year i within event window
  treat.dat <- est.dat %>%
    filter(state_treat_year==i,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="treated")

  ## define control group relative to state_treat_year i within event window
  ## never treated
  control.dat.never <- est.dat %>% 
    filter(state_treat_year==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="never")

  ## not yet treated
  control.dat.notyet <- est.dat %>% 
    filter( state_treat_year>(i+post.period),
            year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="notyet")

  ## inner join back to est.dat
  stack.dat.group <- est.dat %>% 
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    left_join(state.count, by=c("MSTATE","year")) %>%    
    mutate(state=as.numeric(factor(MSTATE)),
      stacked_event_time=year-i,
      stack_group=i)

  ## collapse to state-year level
  stack.dat1 <- stack.dat.group %>%
    group_by(MSTATE, year) %>%
    summarize(changes=sum(closed_merged), sum_cah=sum(cah, na.rm=TRUE),
              group_type=first(treat_type), closures=sum(closed), mergers=sum(merged),
              hospitals=first(hospitals), hospitals_lag=first(hospitals_lag)) %>%
    arrange(MSTATE, year) %>%
    group_by(MSTATE) %>%
    mutate(cumul_changes=cumsum(changes), rate_changes=changes/hospitals_lag, rate_closed=closures/hospitals_lag, rate_merged=mergers/hospitals_lag) %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i) %>%
    ungroup()

  ## collapse to state-year level with restrictions on hospital types
  stack.dat2 <- stack.dat.group %>% 
    group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE), max_distance=max(distance, na.rm=TRUE)) %>% ungroup() %>%  
    filter(min_bedsize<=50) %>%
    group_by(MSTATE, year) %>%
    summarize(changes=sum(closed_merged), sum_cah=sum(cah, na.rm=TRUE),
              group_type=first(treat_type), closures=sum(closed), mergers=sum(merged),
              hospitals=first(hospitals), hospitals_lag=first(hospitals_lag)) %>%
    arrange(MSTATE, year) %>%
    group_by(MSTATE) %>%
    mutate(cumul_changes=cumsum(changes), rate_changes=changes/hospitals_lag, rate_closed=closures/hospitals_lag, rate_merged=mergers/hospitals_lag) %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i) %>%
    ungroup()

  stack.state1 <- bind_rows(stack.state1, stack.dat1)
  stack.state2 <- bind_rows(stack.state2, stack.dat2)
 
}

## drop the first_year_treat==0 phantom treatment group
stack.state1 <- stack.state1 %>% ungroup() %>%
  filter(stack_group>0) %>%
  mutate(treated=ifelse(group_type=="treated",1,0),
         control_any=ifelse(group_type!="treated",1,0),
         control_notyet=ifelse(group_type=="notyet",1,0),
         control_never=ifelse(group_type=="never",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated,
         all_treated=sum(treated),
         all_control=sum(control_any),
         all_control_notyet=sum(control_notyet)) %>%
  group_by(stack_group) %>%
  mutate(group_treated=sum(treated),
         group_control_any=sum(control_any),
         group_control_notyet=sum(control_notyet)) %>%
  ungroup() %>%
  mutate(weight_any=case_when(
           group_type=="treated" ~ 1,
           group_type!="treated" ~ (group_treated/all_treated)/(group_control_any/all_control)),
         weight_notyet=case_when(
           group_type=="treated" ~ 1,
           group_type=="notyet" ~ (group_treated/all_treated)/(group_control_notyet/all_control_notyet))
  )

stack.state2 <- stack.state2 %>% filter(stack_group>0) %>%
  mutate(treated=ifelse(group_type=="treated",1,0),
         control_any=ifelse(group_type!="treated",1,0),
         control_notyet=ifelse(group_type=="notyet",1,0),
         control_never=ifelse(group_type=="never",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated,
         all_treated=sum(treated),
         all_control=sum(control_any),
         all_control_notyet=sum(control_notyet)) %>%
  group_by(stack_group) %>%
  mutate(group_treated=sum(treated),
         group_control_any=sum(control_any),
         group_control_notyet=sum(control_notyet)) %>%
  ungroup() %>%
  mutate(weight_any=case_when(
           group_type=="treated" ~ 1,
           group_type!="treated" ~ (group_treated/all_treated)/(group_control_any/all_control)),
         weight_notyet=case_when(
           group_type=="treated" ~ 1,
           group_type=="notyet" ~ (group_treated/all_treated)/(group_control_notyet/all_control_notyet))
  )

# Source analysis code files -----------------------------------------------

source('analysis/1-sum-stats.R')
source('analysis/2-margins.R')
source('analysis/3-closures-mergers.R')