## Prepare data for estimation ----------------------------------------------

est.dat <- final.dat %>%
  mutate(
    event_time=case_when(
      first_year_obs>0 ~ year - first_year_obs,
      first_year_obs==0 ~ -1),
    aha_id=as.numeric(ID),
    treat_state=ifelse(first_year_obs>0, 1, 0),
    treat_time=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0),
    compare_hosp=case_when(
      year<1999 & BDTOT<12 & distance>30 ~ 1,
      year>=1999 & BDTOT<25 & distance>30 ~ 1,
      TRUE ~ 0))

## plot closures over time by treat_state
est.dat %>%
  group_by(year, treat_state) %>%
  summarize(closures=sum(closed)) %>%
  ggplot(aes(x=year, y=closures, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of closures", color="Treatment group") +
  theme(legend.position="bottom")





## Inverse probability weighting ---------------------------------------------

## create propensity weights
logit.dat <- final.dat %>%
  group_by(ID) %>% 
  summarize(ever_cah=max(ever_cah),
            mean_bed=mean(BDTOT, na.rm=TRUE),
            mean_dist=mean(distance, na.rm=TRUE)) %>%
  filter(!is.na(mean_bed), !is.na(mean_dist))


logit.reg <- glm(ever_cah ~ mean_bed + mean_dist, data=logit.dat, 
          family=binomial(link='logit'))
id.weights <- logit.dat %>%
  mutate(ps = predict(logit.reg, type = 'response')) %>%
  filter(ps>0 & ps<1) %>%
  mutate(ipw = case_when(
    ever_cah == 1 ~ 1/ps,
    ever_cah == 0 ~ 1/(1-ps)
  ))

## add to estimation data
est.dat <- est.dat %>%
  left_join(id.weights %>% select(ID, ipw), by="ID")


## Initial TWFE ---------------------------------------------------------------

## feols using ipw as weights
mod.twfe1 <- feols(closed~i(event_time, treat_state, ref=-1) | year + ID,
                   cluster="ID",
                   weights=est.dat$ipw,
                 data=est.dat)
iplot(mod.twfe1, 
      xlab = 'Time to treatment',
      main = 'Event study')


mod.twfe2 <- feols(closed~i(event_time, treat_state, ref=-1) | ID + year,
                  cluster=~ID,
                  data=est.dat)
iplot(mod.twfe2, 
      xlab = 'Time to treatment',
      main = 'Event study')


## Sun and Abraham -----------------------------------------------------------

mod.sa1 <- feols(closed ~  sunab(first_year_obs, year)
              | year + ID,
              weights=est.dat$ipw,
              cluster="ID", 
              data=est.dat)
iplot(mod.sa1, 
      xlab = 'Time to treatment',
      main = 'Event study')


mod.sa2 <- feols(closed ~  sunab(first_year_obs, year)
              | year + ID,
              cluster="ID", 
              data=est.dat %>% filter(BDTOT<25, distance>30))
iplot(mod.sa2, 
      xlab = 'Time to treatment',
      main = 'Event study')
summary(mod.sa2, agg="ATT")
iplot(mod.sa2, 
      xlab = 'Time to treatment',
      main = 'Event study')



iplot(list(mod.sa1, mod.sa2, mod.sa3), ref.line=-1, xlab="Time", main="",
      ylab = "Estimate and 95% CI", pt.pch=c(20,17,15), pt.col="black", ci.col="black")
legend("bottomleft", pch=c(20,17,15), bty="n",
       legend=c("All Hospitals","CH Only","NCH Only"))
dev.copy(png,"results/f5-sa-ch.png")
dev.off()


## Callaway and Sant'Anna ----------------------------------------------------

state.dat1 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_obs), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat2 <- est.dat %>% filter(BDTOT<75, distance>20) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_obs), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

csa.state1 <- att_gt(yname="closures",
                   gname="cah_treat",
                   idname="state",
                   tname="year",
                   panel=TRUE,
                   control_group="notyettreated",
                   data = state.dat1,
                   est_method="reg")
csa.att1 <- aggte(csa.state1, type="simple")
summary(csa.att1)

csa.es1 <- aggte(csa.state1, type="dynamic")
summary(csa.es1)
ggdid(csa.es1)

csa.state2 <- att_gt(yname="closures",
                   gname="cah_treat",
                   idname="state",
                   tname="year",
                   panel=TRUE,
                   control_group="notyettreated",
                   data = state.dat2,
                   est_method="reg")
csa.att2 <- aggte(csa.state2, type="simple", na.rm=TRUE)
summary(csa.att2)

csa.es2 <- aggte(csa.state2, type="dynamic", na.rm=TRUE)
summary(csa.es2)
ggdid(csa.es2)


csa.mod1 <- att_gt(yname="closed",
                   gname="first_year_obs",
                   idname="aha_id",
                   tname="year",
                   panel=FALSE,
                   control_group="notyettreated",
                   data = est.dat,
                   est_method="reg")
csa.att <- aggte(csa.mod1, type="simple", na.rm=TRUE)
summary(csa.att)

csa.es <- aggte(csa.mod1, type="dynamic", na.rm=TRUE)
summary(csa.es)
ggdid(csa.es)


## Two-Stage Difference-in-Differences ----------------------------------------

mod.did2s <- did2s(est.dat,
  yname = "closed",
  first_stage = ~ 0 | year,
  second_stage = ~ i(event_time, ref=-1), 
  treatment="treat_time",
  weights="ipw",
  cluster="ID")
iplot(mod.did2s, 
      xlab = 'Time to treatment',
      main = 'Event study')

## Stacked DID ----------------------------------------------------------------

## loop over all possible values of first_year_obs

time.window <- 3
outcome.window <- 1
stack.dat <- tibble()
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
    filter( first_year_obs>(i+time.window+outcome.window), 
            year>=(i-time.window), year<=(i+time.window)) %>%
    select(ID, year)
  
  ## inner join back to est.dat
  stack.dat.group <- est.dat %>%
    inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by=c("ID","year")) %>%
    mutate(treat_type=case_when(
      first_year_obs==i ~ "treated",
      first_year_obs==0 ~ "never",
      first_year_obs>(i+time.window+outcome.window) ~ "notyet"))


  ## year of closure among hospitals that closed during event window
  year.closed <- stack.dat.group %>%
    filter(closed==1) %>%
    group_by(ID) %>%
    summarize(year_closed=min(year))

  ## expland closure definition to any closure within time window and collapse to state level
  stack.dat.group <- stack.dat.group %>%
    left_join(year.closed, by="ID") %>%
    mutate(closed_window=
            case_when(
              year_closed>=year & year_closed<=(year+outcome.window) ~ 1,
              TRUE ~ 0)) %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals=n(), closures=sum(closed_window), sum_cah=sum(cah, na.rm=TRUE),
              group_type=first(treat_type)) %>%
    ungroup() %>%
    mutate(state=as.numeric(factor(MSTATE)),
          stacked_event_time=year-i,
          stack_group=i)

  stack.dat <- bind_rows(stack.dat, stack.dat.group)
 
}

## drop the first_year_obs==0 phantom treatment group
stack.dat <- stack.dat %>% filter(stack_group>0) %>%
  mutate(treated=ifelse(group_type=="treated",1,0),
         post=ifelse(stacked_event_time>=0,1,0),
         post_treat=post*treated,
         stacked_event_time_treat=stacked_event_time*treated)
stack.dat %>% group_by(stacked_event_time, group_type) %>% summarize(closures=mean(closures))

## run stacked DD
stack.mod1 <- feols(closures ~ post_treat | year^stack_group + state^stack_group,
                   data=stack.dat,
                   cluster="state")
stack.mod2 <- feols(closures ~ i(stacked_event_time, ref=-1) + i(stacked_event_time_treat, ref=-1) | state^stack_group,
                   data=stack.dat,
                   cluster="state")
iplot(stack.mod2, i.select=2, 
      xlab = 'Time to treatment',
      main = 'Event study')     

## re-run with some restrictions on hospital size and distance when aggregating to state level              
