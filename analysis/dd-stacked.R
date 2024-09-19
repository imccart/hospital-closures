## Stacked DD by State Treatment ----------------------------------------------

## loop over all possible values of first_year_treat
post.period <- 5
pre.period <- 5
stack.dat1 <- tibble()
stack.dat2 <- tibble()
for (i in unique(est.dat$first_year_treat)) {

  ## define treated group relative to first_year_treat i within event window
  treat.dat <- est.dat %>%
    filter(first_year_treat==i,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year)

  ## define control group relative to first_year_treat i within event window  
  ## never treated
  control.dat.never <- est.dat %>% 
    filter(first_year_treat==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year)

  ## not yet treated
  control.dat.notyet <- est.dat %>% 
    filter( first_year_treat>(i+post.period), 
            year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year)
  
  ## inner join back to est.dat
  stack.dat.group <- est.dat %>% 
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    mutate(treat_type=case_when(
      first_year_treat==i ~ "treated",
      first_year_treat==0 ~ "never",
      first_year_treat>(i+post.period) ~ "notyet"))

  ## collapse to state-year level
  stack.dat.group1 <- stack.dat.group %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals=n(), changes=sum(closed_merged), sum_cah=sum(cah, na.rm=TRUE),
              group_type=first(treat_type)) %>%
    arrange(MSTATE, year) %>%
    mutate(cumul_changes=cumsum(changes)) %>%
    ungroup() %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i)

  ## limit to similar hospital types and collapse to state level
  stack.dat.group2 <- stack.dat.group %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals=sum(ipw), changes=sum(closed_merged*ipw), sum_cah=sum(cah*ipw, na.rm=TRUE),
              group_type=first(treat_type)) %>%
    arrange(MSTATE, year) %>%
    mutate(cumul_changes=cumsum(changes)) %>%
    ungroup() %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i)          

  stack.dat1 <- bind_rows(stack.dat1, stack.dat.group1)
  stack.dat2 <- bind_rows(stack.dat2, stack.dat.group2)
 
}

## drop the first_year_treat==0 phantom treatment group
stack.dat1 <- stack.dat1 %>% ungroup() %>%
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
           group_type=="notyet" ~ (group_treated/all_treated)/(group_control_notyet/all_control_notyet)),
         share_changes=changes/hospitals
  )

stack.dat2 <- stack.dat2 %>% filter(stack_group>0) %>%
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


## run stacked DD
stack.mod1 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) + treated,
                   data=stack.dat1,
                   weights=stack.dat1$weight_any,
                   cluster="state")
stack.mod2 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) + treated,
                   data=stack.dat1 %>% filter(group_type!="never"),
                   weights=(stack.dat1 %>% filter(group_type!="never"))$weight_notyet,
                   cluster="state")
stack.mod3 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat1,
                   weights=stack.dat1$weight_any,
                   cluster="state")
stack.mod4 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat1 %>% filter(group_type!="never"),
                   weights=(stack.dat1 %>% filter(group_type!="never"))$weight_notyet,
                   cluster="state")
png("results/stacked-dd-all.png")             
iplot(stack.mod2, i.select=2, 
      xlab = 'Time to treatment',
      main = 'Event study')
dev.off()

## re-run with some restrictions on hospital size and distance when aggregating to state level              
stack.mod1 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | treated,
                   data=stack.dat2,
                   weights=stack.dat2$weight_any,
                   cluster="state")
stack.mod2 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | treated,
                   data=stack.dat2 %>% filter(group_type!="never"),
                   weights=(stack.dat2 %>% filter(group_type!="never"))$weight_notyet,
                   cluster="state")
stack.mod3 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat2,
                   cluster="state")
stack.mod4 <- feols(cumul_changes ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat2 %>% filter(group_type!="never"),
                   cluster="state")               
png("results/stacked-dd-weighted.png")                 
iplot(stack.mod2, i.select=2, 
      xlab = 'Time to treatment',
      main = 'Event study')
dev.off()




## Stacked DD by Hospital Treatment ----------------------------------------------

## loop over all possible values of eff_year (year of CAH designation)
post.period <- 5
pre.period <- 5
stack.hosp1 <- tibble()
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
    filter(BDTOT<75, is.na(eff_year), first_year_treat==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="never")

  ## not yet treated (+ no CAH designation *yet* in the state)
  control.dat.notyet <- est.dat %>% 
    filter( (eff_year>(i+post.period) | is.na(eff_year)), first_year_treat>(i+post.period),
            BDTOT<75,
            year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="notyet")
  
  ## inner join back to est.dat
  stack.dat.group <- est.dat %>% 
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    mutate(state=as.numeric(factor(MSTATE)),
      stacked_event_time=year-i,
      stack_group=i)

  stack.hosp1 <- bind_rows(stack.hosp1, stack.dat.group)
 
}

## drop the eff_year==0 phantom treatment group
stack.hosp1 <- stack.hosp1 %>% ungroup() %>%
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
