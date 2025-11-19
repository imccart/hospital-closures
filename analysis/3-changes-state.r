# =========================================================
# Select outcome (one of: "closures", "mergers", "changes")
# =========================================================
outcome <- "changes"

# Map outcome -> variable names, axis labels, filename slugs
vars <- list(
  closures = list(y_count="closures", y_rate="rate_closed",
                  axis_label="Closures", slug="closure-rate"),
  mergers  = list(y_count="mergers", y_rate="rate_merged",
                  axis_label="Mergers", slug="merger-rate"),
  changes  = list(y_count="changes", y_rate="rate_changes",
                  axis_label="Closures or mergers", slug="change-rate")
)
o <- vars[[outcome]]

# =====================================================================
# Synthetic DD ---------------------------------------------------------
# =====================================================================

## Single year
denom.2001 <- stack.state1 %>% 
  filter(stacked_event_time <= -1, treated==1, stack_group==2001) %>% 
  summarize(mean_hosp=mean(hospitals, na.rm=TRUE)) %>%
  pull(mean_hosp)

synth.2001 <- stack.state1 %>%     
  filter(stack_group == 2001) %>%
  transmute(ID = as.numeric(factor(MSTATE)),
            year,
            rate = .data[[o$y_count]],        # outcome = counts (closures/mergers/changes)
            post_treat)

balance.2001  <- as_tibble(makeBalancedPanel(synth.2001, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.2001))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- synthdid_se(synth.est, method="jackknife")

w      <- attr(synth.est, "weights")$omega
Y      <- setup$Y
N0     <- setup$N0
N      <- nrow(Y)
years  <- as.numeric(colnames(Y))

treated_path   <- colMeans(Y[(N0+1):N, , drop = FALSE])
synthetic_path <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))

plot_df.2001 <- bind_rows(
  tibble(year = years, group="Synthetic control",
         rate_closed = (synthetic_path/denom.2001)*100),
  tibble(year = years, group="Treated",
         rate_closed = (treated_path/denom.2001)*100)
)

# 95% CI label (top-right), no box/border
att      <- (as.numeric(synth.est)/denom.2001)*100
se_num   <- (as.numeric(se)/denom.2001)*100
ci_low   <- att - 1.96 * se_num
ci_high  <- att + 1.96 * se_num

y_top <- max(plot_df.2001$rate_closed, na.rm = TRUE)
x_pos <- max(plot_df.2001$year)

lab_txt <- sprintf("ATT and 95%% CI: %.2f [%.2f, %.2f]", att, ci_low, ci_high)

plot.sdid2001 <- ggplot(plot_df.2001, aes(x = year, y = rate_closed, linetype = group)) +
  geom_line(linewidth = 0.8, color = "black")   +
  geom_vline(xintercept = 2001, linewidth = 1)  +
  scale_linetype_manual(values = c("Synthetic control" = "dashed", "Treated" = "solid")) +
  labs(x = "Year", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))  +
  annotate("text", x = x_pos, y = y_top, label = lab_txt, hjust = 1, vjust = 1, size = 3.4)

ggsave(paste0("results/", o$slug, "-sdid-1999.png"), plot.sdid2001,
       width = 6.5, height = 4.25, dpi = 300, scale=1.5)

## All cohorts
denom.lag <- stack.state1 %>% 
  filter(stacked_event_time <= -1, treated==1) %>% 
  group_by(stack_group) %>% 
  summarize(mean_hosp=mean(hospitals, na.rm=TRUE))

cohorts <- c(1999, 2000, 2001, 2002)

run_sdid_state <- function(c){
  dat <- stack.state1 %>%
    filter(stack_group == c) %>%
    transmute(ID = as.numeric(factor(MSTATE)),
              year,
              rate = .data[[o$y_count]],     # outcome = counts (closures/mergers/changes)
              post_treat)

  bal <- as_tibble(makeBalancedPanel(dat, idname = "ID", tname = "year"))
  setup <- panel.matrices(as.data.frame(bal))

  est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se  <- synthdid_se(est, method = "jackknife")

  # paths
  w       <- attr(est, "weights")$omega
  Y       <- setup$Y
  N0      <- setup$N0
  N       <- nrow(Y)
  years   <- as.numeric(colnames(Y))
  tau     <- years - c

  syn_path   <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))
  treat_path <- colMeans(Y[(N0 + 1):N, , drop = FALSE])
  Ntr        <- N - N0

  denom_c <- denom.lag %>% filter(stack_group == c) %>% pull(mean_hosp)
  list(
    cohort = c,
    paths  = tibble(cohort = c, tau = tau,
                    treated = (treat_path/denom_c)*100,
                    synthetic = (syn_path/denom_c)*100,
                    Ntr = Ntr),
    att_se = tibble(cohort = c,
                    att = (as.numeric(est)/denom_c)*100,
                    se  = (as.numeric(se)/denom_c)*100,
                    Ntr = Ntr)
  )
}

# Run all cohorts, collect paths and estimates
out        <- map(cohorts, run_sdid_state)
paths_all  <- bind_rows(map(out, "paths"))
atts_all   <- bind_rows(map(out, "att_se"))

# Weighted-average paths across cohorts at each event time tau
agg_paths <- paths_all %>%
  group_by(tau) %>%
  summarise(
    treated   = weighted.mean(treated, w = Ntr, na.rm = TRUE),
    synthetic = weighted.mean(synthetic, w = Ntr, na.rm = TRUE),
    cohorts_n = n(),
    .groups = "drop"
  )

# Pooled ATT using inverse-variance weights (approximate independence across cohorts)
atts_all <- atts_all %>% mutate(var = se^2, w = 1/var)
att_w  <- with(atts_all, sum(w * att) / sum(w))
se_w   <- sqrt(1 / sum(atts_all$w))
ci_low <- att_w - 1.96 * se_w
ci_high<- att_w + 1.96 * se_w

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

att_header <- "ATT and 95%CI"
att_body   <- sub("^ATT and 95%CI\\n", "", att_lab)

# Plot SDID paths
y_top  <- max(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE)
x_pos  <- max(agg_paths$tau) - 2
yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))

plot.sdid <- ggplot(agg_paths, aes(x = tau)) +
  geom_line(aes(y = treated,   linetype = "Treated"),   linewidth = 0.8, color = "black") +
  geom_line(aes(y = synthetic, linetype = "Synthetic"), linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1) +
  scale_linetype_manual(values = c("Treated" = "solid", "Synthetic" = "dashed")) +
  scale_x_continuous(breaks = seq(min(agg_paths$tau), max(agg_paths$tau), by = 1)) +
  labs(x = "Event time", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
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
  annotate("text", x = x_pos, y = y_top - 0.03*yrange,
           label = att_body, hjust = 0, vjust = 1,
           size = 3.4, lineheight = 1.05)

ggsave(paste0("results/", o$slug, "-sdid.png"), plot.sdid,
       width = 6.5, height = 4.25, dpi = 300, scale = 1.5)

# =====================================================================
# Standard TWFE / SA / CS ---------------------------------------------
# =====================================================================

min.es <- -5
max.es <- 5

twfe.mod1 <- fepois(y ~ i(event_time, ref=-1) | MSTATE + year, 
  offset = ~ log(hospitals_lag),
  cluster=~MSTATE,
  data=state.dat1 %>% filter(state_treat_year == 0 | state_treat_year > 1995) %>%
    mutate(event_time=case_when(
      event_time > max.es ~ max.es,
      event_time < min.es ~ min.es,
      TRUE ~ event_time),
      y=.data[[o$y_count]])
) 

# Sun and Abraham ------------------------------------------------------

min.es <- -10
max.es <- 10

sa.mod1 <- fepois(
  y ~ sunab(state_treat_year, event_time, ref.p=c(-2,-1), bin.rel = "bin::2") | MSTATE + year,
  offset   = ~ log(hospitals_lag),
  cluster  = ~ MSTATE,
  data=state.dat1 %>% filter(state_treat_year == 0 | state_treat_year > 1995) %>%
    mutate(event_time=case_when(
      event_time > max.es ~ max.es,
      event_time < min.es ~ min.es,
      TRUE ~ event_time),
      y=.data[[o$y_count]])      
) 

# Callaway and Sant'Anna ------------------------------------------------

csa.mod1 <- att_gt(yname=o$y_rate,
                   gname="state_treat_year",
                   idname="state",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = state.dat1 %>% filter(state_treat_year == 0 | state_treat_year > 1995),
                   base_period="universal",
                   est_method="ipw",
                   clustervars="state",
                   bstrap=TRUE)
csa.att1 <- aggte(csa.mod1, type="simple", na.rm=TRUE)
summary(csa.att1)
csa.es1 <- aggte(csa.mod1, type="dynamic", na.rm=TRUE, min_e=-5, max_e=5)
summary(csa.es1)

# Graph of TWFE, CS, SA estimators -------------------------------------

base.rate <- state.dat1 %>%
  filter(event_time <= -1, treat==1,
         state_treat_year == 0 | state_treat_year > 1995) %>%
  summarise(rate = sum(.data[[o$y_count]], na.rm=TRUE) / sum(hospitals_lag, na.rm=TRUE)) %>%
  pull(rate)  

new.row <- tibble(
  estimate=0,
  conf.low=0,
  conf.high=0,
  event_time=-1
)

## TWFE
point.twfe1 <- as_tibble(twfe.mod1$coefficients, rownames="term") %>%
  filter(str_starts(term,"event_time::")) %>%
  rename(estimate=value) %>%
  mutate(event_time = as.integer(str_remove(term, "event_time::"))) %>%
  select(-term)

ci.twfe1 <- as_tibble(confint(twfe.mod1), rownames="term") %>%
  filter(str_starts(term,"event_time::")) %>%
  rename(conf.low = `2.5 %`, conf.high = `97.5 %`) %>%
  mutate(event_time = as.integer(str_remove(term, "event_time::"))) %>%
  select(-term)

est.twfe1 <- point.twfe1 %>%
  left_join(ci.twfe1, by="event_time") %>%
  bind_rows(new.row) %>% 
  mutate(estimator="TWFE") %>%
  arrange(event_time) %>%
  select(event_time, estimate, conf.low, conf.high, estimator) %>%
  mutate(estimate = (exp(estimate) - 1) * base.rate,
         conf.low = (exp(conf.low) - 1) * base.rate,
         conf.high= (exp(conf.high) - 1) * base.rate)

## SA
est.sa1 <- tidy(sa.mod1, conf.int = TRUE) %>%
  filter(str_detect(term, "^event_time::")) %>%
  mutate(event_time = as.integer(str_remove(term, "event_time::"))) %>%
  bind_rows(new.row) %>% 
  mutate(estimator = "SA",
         event_bin = floor(event_time / 2)) %>%
  select(event_time=event_bin, estimate, conf.low, conf.high, estimator) %>%
  mutate(estimate = (exp(estimate) - 1) * base.rate,
         conf.low = (exp(conf.low) - 1) * base.rate,
         conf.high= (exp(conf.high) - 1) * base.rate)

## CS
est.cs1 <- tibble(
  event_time = csa.es1$egt,
  estimate   = csa.es1$att,
  se         = csa.es1$se
) %>%
  mutate(
    conf.low  = if_else(event_time!= -1, estimate - 1.96 * se, 0),
    conf.high = if_else(event_time!= -1, estimate + 1.96 * se, 0),
    estimator = "CS"
  ) %>%
  select(-se)

## Combine and plot
est.all <- bind_rows(est.twfe1, est.cs1, est.sa1) %>%
  arrange(estimator, event_time) %>%
  mutate(
    Estimator = factor(estimator, levels = c("TWFE", "CS", "SA")),
    estimate = 100*estimate,
    conf.low = 100*conf.low,
    conf.high= 100*conf.high
  )  

plot.estimators <- ggplot(est.all, aes(x = event_time, y = estimate, shape = Estimator)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.2, alpha = 0.3, color = "black") +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2.5, color = "black", stroke = 0.1, fill="white") +
  scale_shape_manual(values = c(
    "TWFE" = 24,    # hollow triangle
    "CS" = 21,      # hollow circle
    "SA" = 22       # hollow square
  )) +
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
  labs(x = "Event time", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

plot.estimators2 <- ggplot(est.all %>% filter(estimator=="CS"), aes(x = event_time, y = estimate, shape = Estimator)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.2, alpha = 0.3, color = "black") +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2.5, color = "black", stroke = 0.1, fill="white") +
  scale_shape_manual(values = c(
    "CS" = 21      # hollow circle
  )) +
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
  labs(x = "Event time", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


ggsave(paste0("results/", o$slug, "-other.png"), plot.estimators2,
       width = 6.5, height = 4.25, dpi = 300, scale=1.5)
