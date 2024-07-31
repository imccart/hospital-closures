## Synthetic DD ---------------------------------------------------------------
## Notes on Synth DD. 
##   - Must be data frame (tibble forces errors) and must have more than 1 pre-treatment period

# use stacked DD data at hospital level
synth.1999 <- stack.hosp1 %>% 
            filter(stack_group==2004, year>=1998) %>%
            filter(!is.na(op_margin_cut)) %>%
            select(ID, year, op_margin_cut, post_treat)

synth.1999  <- as_tibble(makeBalancedPanel(synth.1999, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(synth.1999))
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


# use stacked DD data at state level
synth.1994 <- as.data.frame(stack.dat1 %>% 
                            filter(stack_group==1995) %>%
  select(MSTATE, year, cumul_changes, post_treat))

setup <- panel.matrices(synth.1994)
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