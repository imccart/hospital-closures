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

## plot closures, mergers, and both over time by treat_state
changes.plot <- est.dat %>%
  group_by(year, treat_state) %>%
  summarize(closures=sum(closed),
            mergers=sum(merged),
            all_changes=sum(closed_merged)) 

closure.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=closures, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of closures", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/closures.png", closure.plot, width = 6, height = 10, dpi = 300)


merger.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=mergers, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of mergers", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/merger.png", merger.plot, width = 6, height = 10, dpi = 300)


anychange.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=all_changes, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of closures or mergers", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/all-changes.png", anychange.plot, width = 6, height = 10, dpi = 300)


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

## Aggregate to state level --------------------------------------------------

state.dat1 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_obs), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    first_year_obs=first(first_year_obs)) %>%
  mutate(
    event_time=case_when(
      first_year_obs>0 ~ year - first_year_obs,
      first_year_obs==0 ~ -1),
    treat_state=ifelse(first_year_obs>0, 1, 0),
    treat_time=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat2 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=sum(ipw), cah_treat=min(first_year_obs), 
    closures=sum(closed*ipw), sum_cah=sum(cah*ipw, na.rm=TRUE),
    mergers=sum(merged*ipw), any_changes=sum(closed_merged*ipw),
    first_year_obs=first(first_year_obs)) %>%
  mutate(
    event_time=case_when(
      first_year_obs>0 ~ year - first_year_obs,
      first_year_obs==0 ~ -1),
    treat_state=ifelse(first_year_obs>0, 1, 0),
    treat_time=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat3 <- est.dat %>% filter(BDTOT<75, distance>20) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_obs), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), any_changes=sum(closed_merged),
    first_year_obs=first(first_year_obs)) %>%
  mutate(
    event_time=case_when(
      first_year_obs>0 ~ year - first_year_obs,
      first_year_obs==0 ~ -1),
    treat_state=ifelse(first_year_obs>0, 1, 0),
    treat_time=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))


## Initial TWFE ---------------------------------------------------------------


## feols without weights, hospital level
mod.twfe1 <- feols(closed~i(event_time, treat_state, ref=-1) | ID + year,
                  cluster="ID",
                  data=est.dat)
iplot(mod.twfe1, 
      xlab = 'Time to treatment',
      main = 'Event study')

# feols without weights, state level
mod.twfe2 <- feols(closures~i(event_time, treat_state, ref=-1) | state + year,
                  cluster="state",
                  data=state.dat1)
iplot(mod.twfe2, 
      xlab = 'Time to treatment',
      main = 'Event study')


## feols using ipw as weights, state level
mod.twfe3 <- feols(closures~i(event_time, treat_state, ref=-1) | year + state,
                   cluster="state",
                 data=state.dat2)
iplot(mod.twfe3, 
      xlab = 'Time to treatment',
      main = 'Event study')


## feols with distance and bed size restrictions, state level
mod.twfe4 <- feols(closures~i(event_time, treat_state, ref=-1) | year + state,
                   cluster="state",
                 data=state.dat3)
iplot(mod.twfe4, 
      xlab = 'Time to treatment',
      main = 'Event study')




## Sun and Abraham -----------------------------------------------------------

## hospital level, no weights
mod.sa1 <- feols(closed ~  sunab(first_year_obs, year)
              | year + ID,
              cluster="ID", 
              data=est.dat)
iplot(mod.sa1, 
      xlab = 'Time to treatment',
      main = 'Event study')

## state level, no weights
mod.sa2 <- feols(closures ~  sunab(first_year_obs, year)
              | year + state,
              cluster="state", 
              data=state.dat1)
iplot(mod.sa2, 
      xlab = 'Time to treatment',
      main = 'Event study')
summary(mod.sa2, agg="ATT")


## combine plots (outdated code for now)
# iplot(list(mod.sa1, mod.sa2, mod.sa3), ref.line=-1, xlab="Time", main="",
#     ylab = "Estimate and 95% CI", pt.pch=c(20,17,15), pt.col="black", ci.col="black")
#legend("bottomleft", pch=c(20,17,15), bty="n",
#       legend=c("All Hospitals","CH Only","NCH Only"))
#dev.copy(png,"results/f5-sa-ch.png")
#dev.off()


## Callaway and Sant'Anna ----------------------------------------------------

## hospital level
## won't run -- too many all-zero cells

## state level, no weights
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

## state level, with weights
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


## Two-Stage Difference-in-Differences ----------------------------------------

mod.did2s <- did2s(state.dat1,
  yname = "closures",
  first_stage = ~ 0 | year,
  second_stage = ~ i(event_time, ref=-1), 
  treatment="treat_time",
  cluster="state")
iplot(mod.did2s, 
      xlab = 'Time to treatment',
      main = 'Event study')

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
