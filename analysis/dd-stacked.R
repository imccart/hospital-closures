## Stacked DID ----------------------------------------------------------------

## loop over all possible values of first_year_obs

time.window <- 3
stack.dat1 <- tibble()
stack.dat2 <- tibble()
for (i in unique(est.dat$first_year_obs)) {

  ## define treated group relative to first_year_obs i within event window
  treat.dat <- est.dat %>%
    filter(first_year_obs==i, 
           year>=(i-time.window), year<=(i+time.window)) %>%
    select(ID, year)

  ## define control group relative to first_year_obs i within event window  
  ## never treated
  control.dat.never <- est.dat %>%
    filter(first_year_obs==0,
           year>=(i-time.window), year<=(i+time.window)) %>%
    select(ID, year)

  ## not yet treated
  control.dat.notyet <- est.dat %>%
    filter( first_year_obs>(i+time.window), 
            year>=(i-time.window), year<=(i+time.window)) %>%
    select(ID, year)
  
  ## inner join back to est.dat
  stack.dat.group <- est.dat %>%
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    mutate(treat_type=case_when(
      first_year_obs==i ~ "treated",
      first_year_obs==0 ~ "never",
      first_year_obs>(i+time.window) ~ "notyet"))


  ## year of closure among hospitals that closed during event window
  year.closed <- stack.dat.group %>%
    filter(closed==1) %>%
    group_by(ID) %>%
    summarize(year_closed=min(year))

  ## collapse to state level
  stack.dat.group1 <- stack.dat.group %>%
    left_join(year.closed, by="ID") %>%
    mutate(closed_window=
            case_when(
              year_closed>= i ~ 1,
              TRUE ~ 0)) %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals=n(), closures=sum(closed_window), sum_cah=sum(cah, na.rm=TRUE),
              group_type=first(treat_type)) %>%
    ungroup() %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i)

  ## limit to similar hospital types and collapse to state level
  stack.dat.group2 <- stack.dat.group %>%
    left_join(year.closed, by="ID") %>%
    mutate(closed_window=
            case_when(
              year_closed>= i ~ 1,
              TRUE ~ 0)) %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals=sum(ipw), closures=sum(closed_window*ipw), sum_cah=sum(cah*ipw, na.rm=TRUE),
              group_type=first(treat_type)) %>%
    ungroup() %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i)          

  stack.dat1 <- bind_rows(stack.dat1, stack.dat.group1)
  stack.dat2 <- bind_rows(stack.dat2, stack.dat.group2)
 
}

## drop the first_year_obs==0 phantom treatment group
stack.dat1 <- stack.dat1 %>% filter(stack_group>0) %>%
  mutate(treated=ifelse(group_type=="treated",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated)
stack.dat2 <- stack.dat2 %>% filter(stack_group>0) %>%
  mutate(treated=ifelse(group_type=="treated",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated)


## run stacked DD
stack.mod1 <- feols(closures ~ post_treat | year^stack_group + state^stack_group,
                   data=stack.dat1,
                   cluster="state")
stack.mod2 <- feols(closures ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat1,
                   cluster="state")
stack.mod3 <- feols(closures ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat1 %>% filter(group_type!="never"),
                   cluster="state")
png("results/stacked-dd-all.png")                 
iplot(stack.mod2, i.select=2, 
      xlab = 'Time to treatment',
      main = 'Event study')
dev.off()

## re-run with some restrictions on hospital size and distance when aggregating to state level              
stack.mod1 <- feols(closures ~ post_treat | year^stack_group + state^stack_group,
                   data=stack.dat2,
                   cluster="state")
stack.mod2 <- feols(closures ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat2,
                   cluster="state")
stack.mod3 <- feols(closures ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat2 %>% filter(group_type!="never"),
                   cluster="state")                   
png("results/stacked-dd-weighted.png")                 
iplot(stack.mod2, i.select=2, 
      xlab = 'Time to treatment',
      main = 'Event study')
dev.off()
