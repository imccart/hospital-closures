# Select outcome (one of: "closures", "mergers", "changes")
outcome <- "closures"

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

# Build state-level dataset for analysis
stack.state <- stack_state(pre.period=5, post.period=5, state.period=0)


# Synthetic DD ---------------------------------------------------------

## Single year
cohort.year <- 2001
denom.year <- stack.state %>% 
  filter(stacked_event_time <= -1, treated==1, stack_group==cohort.year) %>% 
  summarize(mean_hosp=mean(hospitals, na.rm=TRUE)) %>%
  pull(mean_hosp)

synth.year <- stack.state %>%     
  filter(stack_group == cohort.year) %>%
  transmute(ID = as.numeric(factor(MSTATE)),
            year,
            rate = .data[[o$y_count]],        # outcome = counts (closures/mergers/changes)
            post_treat)

balance.year  <- as_tibble(makeBalancedPanel(synth.year, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.year))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- synthdid_se(synth.est, method="jackknife")

w      <- attr(synth.est, "weights")$omega
Y      <- setup$Y
N0     <- setup$N0
N      <- nrow(Y)
years  <- as.numeric(colnames(Y))

treated_path   <- colMeans(Y[(N0+1):N, , drop = FALSE])
synthetic_path <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))

plot_df.year <- bind_rows(
  tibble(year = years, group="Synthetic control",
         rate_out = (synthetic_path/denom.year)*100),
  tibble(year = years, group="Treated",
         rate_out = (treated_path/denom.year)*100)
)

# 95% CI label (top-right), no box/border
att      <- (as.numeric(synth.est)/denom.year)*100
se_num   <- (as.numeric(se)/denom.year)*100
ci_low   <- att - 1.96 * se_num
ci_high  <- att + 1.96 * se_num

y_top <- max(plot_df.year$rate_out, na.rm = TRUE)
x_pos <- max(plot_df.year$year)

lab_txt <- sprintf("ATT and 95%% CI: %.2f [%.2f, %.2f]", att, ci_low, ci_high)

plot.sdid.year <- ggplot(plot_df.year, aes(x = year, y = rate_out, linetype = group)) +
  geom_line(linewidth = 0.8, color = "black")   +
  geom_vline(xintercept = cohort.year, linewidth = 1)  +
  scale_linetype_manual(values = c("Synthetic control" = "dashed", "Treated" = "solid")) +
  labs(x = "Year", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))  +
  annotate("text", x = x_pos, y = y_top, label = lab_txt, hjust = 1, vjust = 1, size = 3.4)

ggsave(paste0("results/", o$slug, "-sdid-year.png"), plot.sdid.year,
       width = 6.5, height = 4.25, dpi = 300, scale=1.5)

## All cohorts
denom.lag <- stack.state %>% 
  filter(stacked_event_time <= -1, treated==1) %>% 
  group_by(stack_group) %>% 
  summarize(mean_hosp=mean(hospitals, na.rm=TRUE))

cohorts <- c(1999, 2000, 2001)

run_sdid_state <- function(c){
  dat <- stack.state %>%
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


# Callaway and Sant'Anna ------------------------------------------------

min.es <- -5
max.es <- 5

denom_cs <- state.dat %>% 
  filter(year < state_treat_year,
         state_treat_year %in% c(1999, 2000, 2001)) %>%
  summarise(mean_hosp = mean(hospitals, na.rm = TRUE)) %>%
  pull(mean_hosp)

scale_cs <- 100 / denom_cs

csa.mod1 <- att_gt(yname=o$y_count,
                   gname="state_treat_year",
                   idname="state",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = state.dat %>% filter(state_treat_year %in% c(0, 1999,2000,2001,2002)),
                   base_period="universal",
                   est_method="ipw",
                   clustervars="state",
                   bstrap=TRUE)
csa.att1 <- aggte(csa.mod1, type="simple", na.rm=TRUE)
summary(csa.att1)
csa.es1 <- aggte(csa.mod1, type="dynamic", na.rm=TRUE, min_e=-5, max_e=5)
summary(csa.es1)

est.cs1 <- tibble(
  event_time = csa.es1$egt,
  estimate   = scale_cs*csa.es1$att,
  se         = scale_cs*csa.es1$se
) %>%
  mutate(
    conf.low  = if_else(event_time!= -1, estimate - 1.96 * se, 0),
    conf.high = if_else(event_time!= -1, estimate + 1.96 * se, 0),
    estimator = "CS"
  ) %>%
  select(-se)


plot.cs <- ggplot(est.cs1, aes(x = event_time, y = estimate)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2.5, color = "black", stroke = 0.1, fill="white") +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + #             
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
  labs(x = "Event time", y = paste0(o$axis_label, " per 100 hospitals"), linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


ggsave(paste0("results/", o$slug, "-cs.png"), plot.cs,
       width = 6.5, height = 4.25, dpi = 300, scale=1.5)
