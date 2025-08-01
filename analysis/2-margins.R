
## Stacked Data by Hospital (hospital-level CAH designation) --------

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
    filter(is.na(eff_year), first_year_treat==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="never")

  ## not yet treated (+ no CAH designation *yet* in the state)
  control.dat.notyet <- est.dat %>% 
    filter( (eff_year>(i+post.period) | is.na(eff_year)),
            first_year_treat>(i+post.period),
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


## Stacked Data by Hospital (state-level CAH availability) --------

## loop over all possible values of first_year_treat (year of CAH availability)
post.period <- 5
pre.period <- 5
stack.hosp2 <- tibble()
for (i in unique(est.dat$first_year_treat)) {

  ## define treated group relative to first_year_treat i within event window
  treat.dat <- est.dat %>%
    filter(first_year_treat==i,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="treated")

  ## define control group relative to first_year_treat i within event window
  ## never treated
  control.dat.never <- est.dat %>% 
    filter(BDTOT<75, first_year_treat==0,
           year>=(i-pre.period), year<=(i+post.period)) %>%
    select(ID, year) %>%
    mutate(treat_type="never")

  ## not yet treated
  control.dat.notyet <- est.dat %>% 
    filter( first_year_treat>(i+post.period),
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

  stack.hosp2 <- bind_rows(stack.hosp2, stack.dat.group)
 
}

## drop the first_year_treat==0 phantom treatment group
stack.hosp2 <- stack.hosp2 %>% ungroup() %>%
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



## Synthetic DD ---------------------------------------------------------------

## Notes on Synth DD. 
##   - Must be data frame (tibble forces errors) and must have more than 1 pre-treatment period

# use stacked DD data at hospital level
synth.1998 <- stack.hosp1 %>%
            filter(stack_group==1998, !is.na(margin)) %>%
            select(ID, year, margin, post_treat)

balance.1998  <- as_tibble(makeBalancedPanel(synth.1998, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.1998))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- synthdid_se(synth.est, method="bootstrap")
sprintf('point estimate: %1.2f', synth.est)
sprintf('95%% CI (%1.2f, %1.2f)', synth.est - 1.96 * se, synth.est + 1.96 * se)
plot(synth.est, se.method='placebo')
synthdid_units_plot(synth.est, se.method='bootstrap')

plot(synth.est, overlay=1, se.method='bootstrap')

synthdid_plot(synth.est, facet.vertical=FALSE,
              control.name='control', treated.name='treated',
              lambda.comparable=TRUE, se.method = 'none',
              trajectory.linetype = 1, line.width=.75, effect.curvature=-.4,
              trajectory.alpha=.7, effect.alpha=.7,
              diagram.alpha=1, onset.alpha=.7) +
    theme(legend.position=c(.26,.07), legend.direction='horizontal',
          legend.key=element_blank(), legend.background=element_blank(),
          strip.background=element_blank(), strip.text.x = element_blank())


## Callaway and Sant'Anna -----------------------------------------------------

# treatment at hospital-level
cs.dat1 <- est.dat %>%
      mutate(ID=as.numeric(factor(ID)), treat_group=if_else(!is.na(eff_year),eff_year,0)) %>%
      filter(!is.na(margin), !is.na(year), !is.na(BDTOT), !is.na(distance), 
          BDTOT<30, distance>10) %>%
      select(ID, treat_group, year, margin, BDTOT, distance, own_type, teach_major) %>%
      mutate(own_type=as.factor(own_type))

csa.margin1 <- att_gt(yname="margin",
                   gname="treat_group",
                   idname="ID",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = cs.dat1,
                   xformla = ~BDTOT+distance,
                   base_period="universal",
                   est_method="ipw")
csa.att1 <- aggte(csa.margin1, type="simple", na.rm=TRUE)
summary(csa.att1)
csa.es1 <- aggte(csa.margin1, type="dynamic", na.rm=TRUE, min_e=-10, max_e=10)
summary(csa.es1)
ggdid(csa.es1)

raw.att <- tibble(
    group = csa.margin1$group,
    time = csa.margin1$t,
    att = csa.margin1$att,
    se = csa.margin1$se
  ) %>%
  mutate(event_time = time - group)

## Standard TWFE -----------------------------------------------------

twfe.dat <- est.dat %>%
      mutate(ID=as.numeric(factor(ID)),treat_group=if_else(!is.na(eff_year),eff_year,-1)) %>%
      filter(!is.na(margin), !is.na(year), BDTOT<30, distance>10) %>%
      mutate(own_type=as.factor(own_type),
             hosp_event_time=case_when(
                hosp_event_time> 10 ~ 10,
                hosp_event_time< -10 ~ -10,
                TRUE ~ hosp_event_time),
             ) %>%
      select(ID, year, margin, BDTOT, distance, eff_year,
            own_type, teach_major, hosp_event_time, treat_group) 

feols(margin~i(hosp_event_time, ref=-1) + BDTOT + distance | ID + year, 
      data=twfe.dat, cluster="ID") 
      
feols(margin~sunab(treat_group, hosp_event_time) + BDTOT + distance | ID + year, 
      data=twfe.dat, cluster="ID") %>%
  iplot(
    main     = "fixest::sunab",
    xlab     = "Time to treatment",
    ref.line = 0
    )
