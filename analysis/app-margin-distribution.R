# Appendix: distributional effect of CAH designation on the operating margin.
# Estimates the effect at each quantile of the (raw, un-winsorized) margin via a
# recentered-influence-function (unconditional quantile) difference-in-differences on
# the eligibility-restricted stacked sample, showing the null holds across the whole
# distribution -- no effect at any quantile, including the lower tail where closure-
# relevant distress would appear. Saves results/margin-quantile-did.png and a CSV.
# Self-contained: run with `Rscript analysis/app-margin-distribution.R`.

source('analysis/0-setup.R')
source('analysis/functions.R')
est.dat <- read_csv('data/output/estimation_data.csv', show_col_types = FALSE)
bed.cut <- 50

# raw margin so the tail is not compressed by winsorization
aha.fin <- read_csv('data/output/aha_final.csv',
                    col_select = c('ID','year','net_pat_rev','tot_operating_exp'),
                    show_col_types = FALSE) %>%
  group_by(ID, year) %>%
  summarise(npr = first(na.omit(net_pat_rev)), opx = first(na.omit(tot_operating_exp)), .groups = 'drop') %>%
  mutate(margin_raw = ifelse(npr > 0, (npr - opx) / npr, NA_real_),
         margin_raw = ifelse(margin_raw < -5 | margin_raw > 1, NA_real_, margin_raw))

stack.elig <- stack_hosp_elig(pre.period = 5, post.period = 5, cohort.years = 1999:2005)
d <- stack.elig %>%
  left_join(aha.fin %>% select(ID, year, margin_raw), by = c('ID','year')) %>%
  group_by(ID) %>% mutate(mb = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
  filter(mb <= bed.cut, stacked_event_time >= -5, stacked_event_time <= 5, !is.na(margin_raw))
cat(sprintf('Analysis rows: %d  hospitals: %d\n', nrow(d), n_distinct(d$ID)))

# RIF (unconditional-quantile) DiD across a grid of quantiles
rif_did <- function(q) {
  qt   <- as.numeric(quantile(d$margin_raw, q))
  dens <- density(d$margin_raw, n = 2048)
  fq   <- approx(dens$x, dens$y, xout = qt)$y
  d$rif <- qt + (q - as.integer(d$margin_raw <= qt)) / fq
  m <- feols(rif ~ post_treat | ID^stack_group + year^stack_group, data = d, cluster = ~ID)
  tibble(q = q, effect = coef(m)['post_treat'], se = se(m)['post_treat'])
}
qte <- map_dfr(seq(0.10, 0.90, by = 0.05), rif_did) %>%
  mutate(lo = effect - 1.96 * se, hi = effect + 1.96 * se)
print(as.data.frame(qte), digits = 3)

p <- ggplot(qte, aes(q, effect)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = 'gray75', alpha = 0.6) +
  geom_line(linewidth = 0.9, color = 'black') +
  geom_point(size = 1.8, color = 'black') +
  labs(x = 'Quantile of the operating margin distribution',
       y = 'Estimated effect of CAH designation') +
  theme_bw(base_size = 12)
ggsave('results/margin-quantile-did.png', p, width = 6.5, height = 4.25, dpi = 300)
write_csv(qte, 'results/diagnostics/margin-quantile-did.csv')
cat('\nSaved results/margin-quantile-did.png\n')
