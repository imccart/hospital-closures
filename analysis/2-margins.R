
# Synthetic DD ---------------------------------------------------------------

## Notes on Synth DD. 
##   - Must be data frame (tibble forces errors) and must have more than 1 pre-treatment period

## use stacked DD data at hospital level for single cohort
synth.2001 <- stack.hosp1 %>% group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE)) %>% ungroup() %>%
            filter(stack_group==2001, !is.na(margin), min_bedsize<=75) %>%
            select(ID, year, margin, post_treat)

balance.2001  <- as_tibble(makeBalancedPanel(synth.2001, idname="ID", tname="year"))

setup <- panel.matrices(as.data.frame(balance.2001))
synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se  <- synthdid_se(synth.est, method="bootstrap")

w      <- attr(synth.est, "weights")$omega    # unit weights on controls
Y      <- setup$Y
N0     <- setup$N0
N      <- nrow(Y)
years  <- as.numeric(colnames(Y))

treated_path   <- colMeans(Y[(N0+1):N, , drop = FALSE])
synthetic_path <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))

plot_df.2001 <- dplyr::bind_rows(
  tibble::tibble(year = years, group = "Synthetic control", margin = synthetic_path),
  tibble::tibble(year = years, group = "Treated",            margin = treated_path)
) %>%
  dplyr::mutate(group = factor(group, levels = c("Synthetic control","Treated")))

# 95% CI label (top-right), no box/border
att      <- as.numeric(synth.est)
se_num   <- as.numeric(se)
ci_low   <- att - 1.96 * se_num
ci_high  <- att + 1.96 * se_num
att_lab  <- sprintf("ATT and 95%% CI: %.2f, [%.2f, %.2f]", att, ci_low, ci_high)

x_tr <- max(plot_df.2001$year)
y_tr <- max(plot_df.2001$margin, na.rm = TRUE)

plot.2001 <- ggplot(plot_df.2001, aes(x = year, y = margin, linetype = group)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 2001, linewidth = 1) +
  scale_linetype_manual(values = c("dashed","solid")) +  # synthetic dashed, treated solid
  scale_x_continuous(breaks = seq(min(plot_df.2001$year), max(plot_df.2001$year), by = 1)) +
  labs(x = "Year", y = "Margin", linetype = NULL) +
  theme_bw() +
  theme(legend.key.width = grid::unit(2, "cm")) +
  annotate("text", x = x_tr, y = y_tr, label = att_lab, hjust = 1, vjust = 1, size = 3.5)

ggsave("results/sdid-margin-2001.png", plot.2001, width = 6.5, height = 4.25, dpi = 300, scale=1.5)


## use stacked DD data at hospital level for all cohorts and combine
cohorts <- 1999:2003

run_sdid_cohort <- function(c){
  # Build cohort panel (re-using your style/objects)
  synth.c <- stack.hosp1 %>%
    group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
    filter(stack_group == c, !is.na(margin), min_bedsize <= 75) %>%
    select(ID, year, margin, post_treat)

  bal.c   <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  setup   <- panel.matrices(as.data.frame(bal.c))
  est     <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se_c    <- synthdid_se(est, method = "bootstrap")

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
x_pos  <- max(agg_paths$tau) - 2   # or your preferred x
yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))

p_evt <- ggplot(agg_paths, aes(x = tau)) +
  geom_line(aes(y = treated,  linetype = "Treated"),  linewidth = 0.8, color = "black") +
  geom_line(aes(y = synthetic,linetype = "Synthetic control"), linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1) +
  scale_linetype_manual(values = c("Treated" = "solid", "Synthetic control" = "dashed")) +
  scale_x_continuous(breaks = seq(min(agg_paths$tau), max(agg_paths$tau), by = 1)) +
  labs(x = "Event time (τ)", y = "Margin", linetype = NULL) +
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

ggsave("results/sdid-margin-eventtime.png", p_evt, width = 6.5, height = 4.25, dpi = 300, scale=1.5)


## Callaway and Sant'Anna -----------------------------------------------------

# treatment at hospital-level
cs.dat1 <- est.dat %>%
      mutate(ID2=as.numeric(factor(ID)), treat_group=if_else(!is.na(eff_year),eff_year,0)) %>%
      filter(!is.na(margin), !is.na(year), !is.na(BDTOT), !is.na(distance), 
          BDTOT<30, distance>10, ID!=6740050) %>%
      select(ID, ID2, treat_group, year, margin, BDTOT, distance, own_type, teach_major) %>%
      mutate(own_type=as.factor(own_type))

csa.margin1 <- att_gt(yname="margin",
                   gname="treat_group",
                   idname="ID2",
                   tname="year",
                   control_group="notyettreated",
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   data = cs.dat1,
                   xformla = ~BDTOT+distance,
                   base_period="universal",
                   est_method="ipw")
csa.att1 <- aggte(csa.margin1, type="simple", na.rm=TRUE)
summary(csa.att1)
csa.es1 <- aggte(csa.margin1, type="dynamic", na.rm=TRUE, min_e=-10, max_e=10)
summary(csa.es1)
png("results/cs-margin.png")
ggdid(csa.es1)
dev.off()            

## Identifying problematic IDs via spikes in ATT
raw.att <- tibble(
    group = csa.margin1$group,
    time = csa.margin1$t,
    att = csa.margin1$att,
    se = csa.margin1$se
  ) %>%
  mutate(event_time = time - group)

raw.att %>% filter(event_time %in% c(-10,-9, -8, -5,-2,-3)) %>% 
  group_by(event_time) %>% 
  summarize(matt=mean(att, na.rm=TRUE), maxatt=max(att, na.rm=TRUE), mina= min(att, na.rm=TRUE), n=n(),
            se=mean(se, na.rm=TRUE))

## Problematic IDs
#  6740050

## Standard TWFE -----------------------------------------------------

twfe.dat <- est.dat %>%
      mutate(ID=as.numeric(factor(ID)),treat_group=if_else(!is.na(eff_year),eff_year,-1)) %>%
      filter(!is.na(margin), !is.na(year), BDTOT<30, distance>10, ID!=6740050) %>%
      mutate(own_type=as.factor(own_type),
             hosp_event_time=case_when(
                hosp_event_time> 10 ~ 10,
                hosp_event_time< -10 ~ -10,
                TRUE ~ hosp_event_time),
             ) %>%
      select(ID, year, margin, BDTOT, distance, eff_year,
            own_type, teach_major, hosp_event_time, treat_group) 

feols(margin~i(hosp_event_time, ref=-1) + BDTOT + distance | ID + year, 
      data=twfe.dat, cluster="ID") 
      
feols(margin~sunab(treat_group, hosp_event_time) + BDTOT + distance | ID + year, 
      data=twfe.dat, cluster="ID") %>%
  iplot(
    main     = "fixest::sunab",
    xlab     = "Time to treatment",
    ref.line = 0
    )
