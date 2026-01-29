# Build stacked dataset at hospital level -----------------------------------------

stack_hosp <- function(pre.period, post.period, state.period) {

  stack.hosp <- tibble()
  for (i in unique(est.dat$eff_year)) {

    ## define treated group relative to eff_year i within event window
    treat.dat <- est.dat %>%
      filter(eff_year == i,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "treated")

    ## define control group relative to eff_year i within event window
    ## never treated (+ never a CAH designation in the state)
    control.dat.never <- est.dat %>%
      filter(is.na(eff_year), state_treat_year == 0,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "never")

    ## not yet treated (+ no CAH designation *yet* in the state)
    control.dat.notyet <- est.dat %>%
      filter((eff_year > (i + post.period) | is.na(eff_year)),
             state_treat_year > (i + state.period),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "notyet")

    ## inner join back to est.dat
    stack.dat.group <- est.dat %>%
      inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet),
                 by = c("ID", "year")) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i)

    stack.hosp <- bind_rows(stack.hosp, stack.dat.group)
  }

  ## drop the eff_year==0 phantom treatment group
  stack.hosp <- stack.hosp %>% ungroup() %>%
    filter(stack_group > 0) %>%
    mutate(treated = ifelse(treat_type == "treated", 1, 0),
           control_any = ifelse(treat_type != "treated", 1, 0),
           control_notyet = ifelse(treat_type == "notyet", 1, 0),
           control_never = ifelse(treat_type == "never", 1, 0),
           post = ifelse(stacked_event_time >= 0, 1, 0),
           post_treat = post * treated,
           stacked_event_time_treat = stacked_event_time * treated,
           all_treated = sum(treated),
           all_control = sum(control_any),
           all_control_notyet = sum(control_notyet)) %>%
    group_by(stack_group) %>%
    mutate(group_treated = sum(treated),
           group_control_any = sum(control_any),
           group_control_notyet = sum(control_notyet)) %>%
    ungroup() %>%
    mutate(weight_any = case_when(
             treat_type == "treated" ~ 1,
             treat_type != "treated" ~ (group_treated / all_treated) / (group_control_any / all_control)),
           weight_notyet = case_when(
             treat_type == "treated" ~ 1,
             treat_type == "notyet" ~ (group_treated / all_treated) / (group_control_notyet / all_control_notyet))
    )

  stack.hosp
}


# Build stacked dataset at state level ---------------------------------------------

stack_state <- function(pre.period, post.period, state.period) {

  ## loop over all possible values of state_treat_year (year of CAH availability)
  stack.state <- tibble()
  state.count <- est.dat %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals = n(), .groups = "drop_last") %>%
    group_by(MSTATE) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(hospitals_lag = lag(hospitals)) %>%
    ungroup() %>%
    select(MSTATE, year, hospitals, hospitals_lag)

  for (i in unique(est.dat$state_treat_year)) {

    ## define treated group relative to state_treat_year i within event window
    treat.dat <- est.dat %>%
      filter(state_treat_year == i,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "treated")

    ## define control group relative to state_treat_year i within event window
    ## never treated
    control.dat.never <- est.dat %>%
      filter(state_treat_year == 0,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "never")

    ## not yet treated
    control.dat.notyet <- est.dat %>%
      filter(state_treat_year > (i + state.period),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "notyet")

    ## inner join back to est.dat
    stack.dat.group <- est.dat %>%
      inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by = c("ID", "year")) %>%
      left_join(state.count, by = c("MSTATE", "year")) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i)

    ## collapse to state-year level
    stack.dat <- stack.dat.group %>%
      group_by(MSTATE, year) %>%
      summarize(changes = sum(closed_merged),
                sum_cah = sum(cah, na.rm = TRUE),
                group_type = first(treat_type),
                closures = sum(closed),
                mergers = sum(merged),
                hospitals = first(hospitals),
                hospitals_lag = first(hospitals_lag),
                .groups = "drop") %>%
      arrange(MSTATE, year) %>%
      group_by(MSTATE) %>%
      mutate(cumul_changes = cumsum(changes),
             rate_changes = changes / hospitals_lag,
             rate_closed  = closures / hospitals_lag,
             rate_merged  = mergers / hospitals_lag) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i) %>%
      ungroup()

    stack.state <- bind_rows(stack.state, stack.dat)
  }

  ## drop the first_year_treat==0 phantom treatment group
  stack.state <- stack.state %>% ungroup() %>%
    filter(stack_group > 0) %>%
    mutate(treated = ifelse(group_type == "treated", 1, 0),
           control_any = ifelse(group_type != "treated", 1, 0),
           control_notyet = ifelse(group_type == "notyet", 1, 0),
           control_never = ifelse(group_type == "never", 1, 0),
           post = ifelse(stacked_event_time >= 0, 1, 0),
           post_treat = post * treated,
           stacked_event_time_treat = stacked_event_time * treated,
           all_treated = sum(treated),
           all_control = sum(control_any),
           all_control_notyet = sum(control_notyet)) %>%
    group_by(stack_group) %>%
    mutate(group_treated = sum(treated),
           group_control_any = sum(control_any),
           group_control_notyet = sum(control_notyet)) %>%
    ungroup() %>%
    mutate(weight_any = case_when(
             group_type == "treated" ~ 1,
             group_type != "treated" ~ (group_treated / all_treated) / (group_control_any / all_control)),
           weight_notyet = case_when(
             group_type == "treated" ~ 1,
             group_type == "notyet" ~ (group_treated / all_treated) / (group_control_notyet / all_control_notyet))
    )

  stack.state
}


