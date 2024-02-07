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
