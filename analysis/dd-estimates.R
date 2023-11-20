## Prepare data for estimation ----------------------------------------------

est.dat <- final.dat %>%
  mutate(
    event_time=case_when(
      ever_cah==1 ~ year - first_year_obs,
      ever_cah==0 ~ -1),
    sunab_cohort=ifelse(first_year_obs==0,-1,first_year_obs),
    aha_id=as.numeric(ID),
    treat_state=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0))

## create propensity weights
logit.dat <- final.dat %>%
  group_by(ID) %>% 
  summarize(ever_cah=max(ever_cah),
            mean_bed=mean(BDTOT, na.rm=TRUE),
            mean_dist=mean(min_dist, na.rm=TRUE)) %>%
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
mod.twfe1 <- feols(closed~i(event_time, ever_cah, ref=-1) | year + ID,
                   cluster="ID",
                   weights=est.dat$ipw,
                 data=est.dat)
iplot(mod.twfe1, 
      xlab = 'Time to treatment',
      main = 'Event study')


mod.twfe2 <- feols(closed~i(event_time, ever_cah, ref=-1) | ID + year,
                  cluster=~ID,
                  data=est.dat %>% filter(BDTOT<30, min_dist>30))
iplot(mod.twfe2, 
      xlab = 'Time to treatment',
      main = 'Event study')


## sun and abraham approach
mod.sa1 <- feols(closed ~  sunab(sunab_cohort, event_time)
              | year + ID,
              weights=est.dat$ipw,
              cluster="ID", 
              data=est.dat)
iplot(mod.sa1, 
      xlab = 'Time to treatment',
      main = 'Event study')


mod.sa2 <- feols(closed ~  sunab(sunab_cohort, event_time)
              | year + ID,
              cluster="ID", 
              data=est.dat %>% filter(BDTOT<25, min_dist>30))
iplot(mod.sa2, 
      xlab = 'Time to treatment',
      main = 'Event study')
summary(mod.sa2, agg="ATT")


iplot(list(mod.sa1, mod.sa2, mod.sa3), ref.line=-1, xlab="Time", main="",
      ylab = "Estimate and 95% CI", pt.pch=c(20,17,15), pt.col="black", ci.col="black")
legend("bottomleft", pch=c(20,17,15), bty="n",
       legend=c("All Hospitals","CH Only","NCH Only"))
dev.copy(png,"results/f5-sa-ch.png")
dev.off()


## callaway and sant'anna approach
state.dat1 <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(first_year_obs), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(state=as.numeric(factor(MSTATE)))

state.dat2 <- est.dat %>% filter(BDTOT<50, min_dist>30) %>%
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
                   panel=TRUE,
                   xformla= ~ BDTOT + min_dist,
                   control_group="notyettreated",
                   data = est.dat %>% filter(BDTOT<25, min_dist>30),
                   est_method="dr")
csa.att <- aggte(csa.mod1, type="simple", na.rm=TRUE)
summary(csa.att)

csa.es <- aggte(csa.mod1, type="dynamic", na.rm=TRUE)
summary(csa.es)
ggdid(csa.es)


## did2s approach
mod.did2s <- did2s(est.dat,
  yname = "closed",
  first_stage = ~ 0 | ID + year,
  second_stage = ~ i(event_time, ref=-1), treatment="treat_state",
  cluster="ID")
iplot(mod.did2s, 
      xlab = 'Time to treatment',
      main = 'Event study')
