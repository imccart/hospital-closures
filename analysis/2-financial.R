# Synthetic DD --------------------------------------------------------
## Notes on Synth DD. 
##   - Must be data frame (tibble forces errors) and must have more than 1 pre-treatment period

## first build stacked data at hospital level
stack.hosp <- stack_hosp(pre.period=5, post.period=5, state.period=0)

## select outcome variable ("margin" or "net_fixed")
outcome_var   <- "net_fixed"
outcome_sym   <- sym(outcome_var)

outcome_label <- case_when(
  outcome_var == "margin"    ~ "Operating margin",
  outcome_var == "net_fixed" ~ "Net fixed assets",
  TRUE                       ~ outcome_var
)

file_stub <- case_when(
  outcome_var == "margin"    ~ "margin",
  outcome_var == "net_fixed" ~ "netfixed",
  TRUE                       ~ outcome_var
)

# SYNTH DD for single cohort ------------------------------------------------

cohort.year <- 2000
synth.year <- stack.hosp %>% group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE)) %>% ungroup() %>%
            filter(stack_group==cohort.year, !is.na(!!outcome_sym), min_bedsize<=50) %>%
            select(ID, year, outcome=!!outcome_sym, post_treat)

balance.year  <- as_tibble(makeBalancedPanel(synth.year, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.year))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- synthdid_se(synth.est, method="jackknife")

w      <- attr(synth.est, "weights")$omega    # unit weights on controls
Y      <- setup$Y
N0     <- setup$N0
N      <- nrow(Y)
years  <- as.numeric(colnames(Y))

treated_path   <- colMeans(Y[(N0+1):N, , drop = FALSE])
synthetic_path <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))

plot_df.year <- bind_rows(
  tibble(year = years, group = "Synthetic control", value = synthetic_path),
  tibble(year = years, group = "Treated",           value = treated_path)
) %>%
  mutate(group = factor(group, levels = c("Synthetic control","Treated")))

# 95% CI label (top-right), no box/border
att      <- as.numeric(synth.est)
se_num   <- as.numeric(se)
ci_low   <- att - 1.96 * se_num
ci_high  <- att + 1.96 * se_num
att_lab <- sprintf(
  "ATT and 95%% CI for %d: %.2f, [%.2f, %.2f]",
  cohort.year, att, ci_low, ci_high
)

plot.year <- ggplot(plot_df.year, aes(x = year, y = value, linetype = group)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_vline(xintercept = cohort.year, linewidth = 1) +
  scale_linetype_manual(values = c("dashed","solid")) +
  scale_x_continuous(breaks = seq(min(plot_df.year$year), max(plot_df.year$year), by = 1)) +
  labs(x = "Year", y = outcome_label, linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.12),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm")
  ) +
  annotate("text", x = max(plot_df.year$year), y = Inf,
           label = att_lab, hjust = 1, vjust = 1.5, size = 3.5)

ggsave(
  sprintf("results/%s-sdid-yearly.png", file_stub),
  plot.year, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
)


# SYNTH DD across cohorts ------------------------------------------------
cohorts <- 1999:2002

run_sdid_cohort <- function(c){
  # Build cohort panel (re-using your style/objects)
  synth.c <- stack.hosp %>%
    group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
    filter(stack_group == c, !is.na(!!outcome_sym), min_bedsize <= 50) %>%
    select(ID, year, outcome=!!outcome_sym, post_treat)

  bal.c   <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  setup   <- panel.matrices(as.data.frame(bal.c))
  est     <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se_c    <- synthdid_se(est, method = "jackknife")

  # SDID unit weights (omega) on controls → synthetic-control path over time
  w       <- attr(est, "weights")$omega         # length = N0 (controls)
  Y       <- setup$Y
  N0      <- setup$N0
  N       <- nrow(Y)
  years   <- as.numeric(colnames(Y))            # calendar time used by panel.matrices
  tau     <- years - c                          # event time

  syn_path   <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))
  treat_path <- colMeans(Y[(N0 + 1):N, , drop = FALSE])
  Ntr        <- N - N0

  list(
    cohort = c,
    paths  = tibble(cohort = c, tau = tau,
                    treated = treat_path, synthetic = syn_path,
                    Ntr = Ntr),
    att_se = tibble(cohort = c,
                    att = as.numeric(est),
                    se  = as.numeric(se_c),
                    Ntr = Ntr)
  )
}

# Run all cohorts, collect paths and estimates
out        <- map(cohorts, run_sdid_cohort)
paths_all  <- bind_rows(map(out, "paths"))
atts_all   <- bind_rows(map(out, "att_se"))

# Weighted-average paths across cohorts at each event time tau (weights = # treated units)
agg_paths <- paths_all %>%
  group_by(tau) %>%
  summarise(
    treated   = weighted.mean(treated, w = Ntr, na.rm = TRUE),
    synthetic = weighted.mean(synthetic, w = Ntr, na.rm = TRUE),
    cohorts_n = n(),
    .groups = "drop"
  )

# Combined ATT and 95% CI (weighted by Ntr; SE uses independence approximation)
att_w   <- with(atts_all, sum(Ntr * att) / sum(Ntr))
se_w    <- with(atts_all, sqrt(sum( (Ntr / sum(Ntr))^2 * se^2 )))
ci_low  <- att_w - 1.96 * se_w
ci_high <- att_w + 1.96 * se_w

# Per-cohort CIs and formatted lines (ordered 1999→2003)
atts_table <- atts_all %>%
  mutate(ci_lo = att - 1.96 * se,
         ci_hi = att + 1.96 * se) %>%
  arrange(match(cohort, cohorts)) %>%
  transmute(line = sprintf("%d: %.2f [%.2f, %.2f]", cohort, att, ci_lo, ci_hi))

att_lab <- paste(
  "ATT and 95%CI",
  sprintf("Overall: %.2f [%.2f, %.2f]", att_w, ci_low, ci_high),
  paste(atts_table$line, collapse = "\n"),
  sep = "\n"
)

# Multi-line annotation text
att_header <- "ATT and 95%CI"
att_body   <- sub("^ATT and 95%CI\\n", "", att_lab)

# Plot: treated vs synthetic, event time, solid vline at 0, wider legend key
y_top  <- max(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE)
x_pos  <- max(agg_paths$tau) - 2 
yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))

p_evt <- ggplot(agg_paths, aes(x = tau)) +
  geom_line(aes(y = treated,  linetype = "Treated"),  linewidth = 0.8, color = "black") +
  geom_line(aes(y = synthetic,linetype = "Synthetic control"), linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1) +
  scale_linetype_manual(values = c("Treated" = "solid", "Synthetic control" = "dashed")) +
  scale_x_continuous(breaks = seq(min(agg_paths$tau), max(agg_paths$tau), by = 1)) +
  labs(x = "Event time", y = outcome_label, linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.12),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm")
  ) +
  annotate("text", x = x_pos, y = y_top,
           label = att_header, hjust = 0, vjust = 1,
           size = 3.6, fontface = "bold") +
  annotate("text", x = x_pos, y = y_top - 0.03*yrange,  # small line break below header
           label = att_body, hjust = 0, vjust = 1,
           size = 3.4, lineheight = 1.05)

ggsave(
  sprintf("results/%s-sdid.png", file_stub),
  p_evt, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
)


# Callaway and Sant'Anna -----------------------------------------------------

min.es <- -5
max.es <- 5

cs.dat <- est.dat %>%
      group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE), max_distance=max(distance, na.rm=TRUE)) %>% ungroup() %>%  
      filter(min_bedsize<=50) %>%
      mutate(y=!!outcome_sym,
            ID2=as.numeric(factor(ID)), 
            treat_group=case_when(
                !is.na(eff_year) ~ eff_year,
                is.na(eff_year) & state_treat_year > year + max.es ~ 0,
                is.na(eff_year) & state_treat_year==0 ~ 0,
                TRUE ~ NA )) %>%
      filter(!is.na(y), !is.na(year), !is.na(BDTOT), !is.na(distance),
         treat_group %in% c(0, 1999, 2000, 2001, 2002)) %>%
      select(ID, ID2, MSTATE, treat_group, year, y, BDTOT, distance, 
             own_type, teach_major, min_bedsize, max_distance) %>%
      mutate(own_type=as.factor(own_type))

csa.raw <- att_gt(yname="y",
                   gname="treat_group",
                   idname="ID2",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = cs.dat,
                   xformla = ~min_bedsize+max_distance,
                   base_period="universal",
                   est_method="ipw")
csa.att <- aggte(csa.raw, type="simple", na.rm=TRUE)
summary(csa.att)
csa.es <- aggte(csa.raw, type="dynamic", na.rm=TRUE, min_e=min.es, max_e=max.es)
summary(csa.es)


## Collect CS estimates
est.cs <- tibble(
  event_time = csa.es$egt,
  estimate   = csa.es$att,
  se         = csa.es$se
) %>%
  mutate(
    conf.low  = if_else(event_time!= -1, estimate - 1.96 * se, 0),
    conf.high = if_else(event_time!= -1, estimate + 1.96 * se, 0),
    estimator = "CS"
  ) %>%
  select(-se)

## Plot
plot.cs <- ggplot(est.cs, aes(x = event_time, y = estimate)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2.5, color = "black", stroke = 0.1, fill="white") +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + #               
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
  labs(x = "Event time", y = "Estimated Effects") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggsave(
  sprintf("results/%s-cs.png", file_stub),
  plot.cs, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
)