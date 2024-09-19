## Synthetic DD ---------------------------------------------------------------
## Notes on Synth DD. 
##   - Must be data frame (tibble forces errors) and must have more than 1 pre-treatment period

# use stacked DD data at hospital level
synth.2000 <- stack.hosp1 %>%
            mutate(margin=case_when(
                  !is.na(margin_1) ~ margin_1,
                  is.na(margin_1) & !is.na(margin_2) ~ margin_2,
                  is.na(margin_1) & is.na(margin_2) & !is.na(margin_3) ~ margin_3)) %>%
            filter(stack_group==2000, year>=1998, year<=2004, !is.na(margin)) %>%
            select(ID, year, margin, post_treat)

balance.2000  <- as_tibble(makeBalancedPanel(synth.2000, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.2000))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- sqrt(vcov(synth.est, method='placebo'))
sprintf('point estimate: %1.2f', synth.est)
sprintf('95%% CI (%1.2f, %1.2f)', synth.est - 1.96 * se, synth.est + 1.96 * se)
plot(synth.est, se.method='placebo')
synthdid_units_plot(synth.est, se.method='placebo')

plot(synth.est, overlay=1, se.method='placebo')

synthdid_plot(synth.est, facet.vertical=FALSE,
              control.name='control', treated.name='treated',
              lambda.comparable=TRUE, se.method = 'none',
              trajectory.linetype = 1, line.width=.75, effect.curvature=-.4,
              trajectory.alpha=.7, effect.alpha=.7,
              diagram.alpha=1, onset.alpha=.7) +
    theme(legend.position=c(.26,.07), legend.direction='horizontal',
          legend.key=element_blank(), legend.background=element_blank(),
          strip.background=element_blank(), strip.text.x = element_blank())


cs.dat <- est.dat %>%
      mutate(ID=as.numeric(factor(ID)), treat_group=if_else(!is.na(eff_year),eff_year,0)) %>%
      filter(!is.na(margin_1), !is.na(year), !is.na(BDTOT)) %>%
      filter(eff_year!=1997) %>%
      select(ID, treat_group, year, margin_1, BDTOT)

csa.margin <- att_gt(yname="margin_1",
                   gname="treat_group",
                   idname="ID",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = cs.dat,
                   xformla = ~BDTOT,
                   base_period="universal",
                   est_method="reg")
csa.att1 <- aggte(csa.margin, type="simple", na.rm=TRUE)
summary(csa.att1)
csa.es1 <- aggte(csa.margin, type="dynamic", na.rm=TRUE)
summary(csa.es1)
ggdid(csa.es1)

## clean up margin variables...