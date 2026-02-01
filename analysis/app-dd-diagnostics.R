# Diagnostics: SDID vs CS Estimate Differences
#
# Run after _run-analysis.r (requires est.dat, stack.hosp)
#
# Diagnostics:
#   1. Sample size comparison (balanced panel attrition)
#   2. Pre-trend analysis for divergent outcomes
#   3. Synthetic weight concentration

cohorts <- 1999:2001
outcomes_to_check <- c("BDTOT", "FTERN", "IPDTOT", "margin")

# =============================================================================
# DIAGNOSTIC 1: Sample Size Comparison
# =============================================================================

sample_comparison <- tibble(
  outcome = character(),
  cohort = integer(),
  sdid_n_treated = integer(),
  sdid_n_control = integer(),
  sdid_n_total = integer(),
  cs_n_treated = integer(),
  cs_n_control = integer(),
  cs_n_total = integer(),
  pct_lost_to_balance = numeric()
)

for (outcome_var in outcomes_to_check) {
  outcome_sym <- sym(outcome_var)

  for (c in cohorts) {
    # SDID sample (balanced panel)
    synth.c <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut) %>%
      select(ID, year, outcome = !!outcome_sym, treated = post_treat)

    # Before balancing
    pre_balance <- synth.c %>%
      group_by(treated) %>%
      summarize(n_hospitals = n_distinct(ID), .groups = "drop")

    # After balancing
    bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    post_balance <- bal.c %>%
      group_by(treated) %>%
      summarize(n_hospitals = n_distinct(ID), .groups = "drop")

    sdid_treated <- post_balance %>% filter(treated == 1) %>% pull(n_hospitals)
    sdid_control <- post_balance %>% filter(treated == 0) %>% pull(n_hospitals)
    if (length(sdid_treated) == 0) sdid_treated <- 0
    if (length(sdid_control) == 0) sdid_control <- 0

    # CS sample (unbalanced allowed)
    cs.c <- est.dat %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(min_bedsize <= bed.cut) %>%
      mutate(
        y = !!outcome_sym,
        treat_group = case_when(
          !is.na(eff_year) ~ eff_year,
          is.na(eff_year) & state_treat_year > year ~ 0,
          is.na(eff_year) & state_treat_year == 0 ~ 0,
          TRUE ~ NA
        )
      ) %>%
      filter(!is.na(y), !is.na(year), !is.na(BDTOT), !is.na(distance),
             treat_group %in% c(0, c))

    cs_treated <- cs.c %>% filter(treat_group == c) %>% summarize(n = n_distinct(ID)) %>% pull(n)
    cs_control <- cs.c %>% filter(treat_group == 0) %>% summarize(n = n_distinct(ID)) %>% pull(n)

    # Calculate attrition
    pre_bal_total <- sum(pre_balance$n_hospitals)
    post_bal_total <- sdid_treated + sdid_control
    pct_lost <- ifelse(pre_bal_total > 0, 100 * (1 - post_bal_total / pre_bal_total), NA)

    sample_comparison <- bind_rows(sample_comparison, tibble(
      outcome = outcome_var,
      cohort = c,
      sdid_n_treated = sdid_treated,
      sdid_n_control = sdid_control,
      sdid_n_total = sdid_treated + sdid_control,
      cs_n_treated = cs_treated,
      cs_n_control = cs_control,
      cs_n_total = cs_treated + cs_control,
      pct_lost_to_balance = round(pct_lost, 1)
    ))
  }
}

sample_summary <- sample_comparison %>%
  group_by(outcome) %>%
  summarize(
    sdid_treated_total = sum(sdid_n_treated),
    sdid_control_total = sum(sdid_n_control),
    cs_treated_total = sum(cs_n_treated),
    cs_control_total = sum(cs_n_control),
    avg_pct_lost = mean(pct_lost_to_balance, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# DIAGNOSTIC 2: Pre-Trend Analysis
# =============================================================================

pretrend_results <- tibble(
  outcome = character(),
  cohort = integer(),
  trend_coef = numeric(),
  trend_se = numeric(),
  trend_pval = numeric(),
  has_pretrend = logical()
)

for (outcome_var in outcomes_to_check) {
  outcome_sym <- sym(outcome_var)
  all_cohort_data <- tibble()

  for (c in cohorts) {
    cohort_dat <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut) %>%
      mutate(outcome = !!outcome_sym)

    pre_dat <- cohort_dat %>%
      filter(stacked_event_time < 0)

    if (nrow(pre_dat) > 100) {
      trend_reg <- tryCatch({
        feols(outcome ~ stacked_event_time * treated, data = pre_dat, cluster = ~ID)
      }, error = function(e) NULL)

      if (!is.null(trend_reg)) {
        coefs <- as.data.frame(coeftable(trend_reg))
        interaction_row <- grep("stacked_event_time:treated", rownames(coefs))

        if (length(interaction_row) > 0) {
          trend_coef <- coefs[interaction_row, "Estimate"]
          trend_se <- coefs[interaction_row, "Std. Error"]
          trend_pval <- coefs[interaction_row, "Pr(>|t|)"]

          pretrend_results <- bind_rows(pretrend_results, tibble(
            outcome = outcome_var,
            cohort = c,
            trend_coef = trend_coef,
            trend_se = trend_se,
            trend_pval = trend_pval,
            has_pretrend = trend_pval < 0.10
          ))
        }
      }
    }

    means_dat <- cohort_dat %>%
      group_by(stacked_event_time, treated) %>%
      summarize(mean_outcome = mean(outcome, na.rm = TRUE),
                se_outcome = sd(outcome, na.rm = TRUE) / sqrt(n()),
                n = n(),
                .groups = "drop") %>%
      mutate(cohort = c)

    all_cohort_data <- bind_rows(all_cohort_data, means_dat)
  }

  agg_means <- all_cohort_data %>%
    group_by(stacked_event_time, treated) %>%
    summarize(mean_outcome = mean(mean_outcome, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(group = ifelse(treated == 1, "Treated", "Control"))

  p <- ggplot(agg_means, aes(x = stacked_event_time, y = mean_outcome,
                              color = group, linetype = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Control" = "gray50", "Treated" = "black")) +
    scale_linetype_manual(values = c("Control" = "dashed", "Treated" = "solid")) +
    labs(
      title = paste0("Pre/Post Trends: ", outcome_var),
      x = "Event Time",
      y = "Mean Outcome",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    sprintf("results/diagnostics/diag-pretrend-%s.png", tolower(outcome_var)),
    p, width = 6, height = 4, dpi = 300
  )
}

# =============================================================================
# DIAGNOSTIC 3: Synthetic Weight Concentration
# =============================================================================

weight_results <- tibble(
  outcome = character(),
  cohort = integer(),
  n_controls = integer(),
  n_nonzero_weights = integer(),
  max_weight = numeric(),
  top5_weight_share = numeric(),
  herfindahl = numeric()
)

for (outcome_var in outcomes_to_check) {
  outcome_sym <- sym(outcome_var)

  for (c in cohorts) {
    synth.c <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut) %>%
      select(ID, year, outcome = !!outcome_sym, treated = post_treat)

    bal.c <- tryCatch({
      as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    }, error = function(e) NULL)

    if (!is.null(bal.c) && nrow(bal.c) > 0) {
      setup <- tryCatch({
        panel.matrices(as.data.frame(bal.c))
      }, error = function(e) NULL)

      if (!is.null(setup) && setup$N0 > 0 && setup$T0 > 0) {
        est <- tryCatch({
          synthdid_estimate(setup$Y, setup$N0, setup$T0)
        }, error = function(e) NULL)

        if (!is.null(est)) {
          w <- attr(est, "weights")$omega

          n_controls <- length(w)
          n_nonzero <- sum(w > 0.001)
          max_w <- max(w)
          top5_share <- sum(sort(w, decreasing = TRUE)[1:min(5, length(w))])
          hhi <- sum(w^2)

          weight_results <- bind_rows(weight_results, tibble(
            outcome = outcome_var,
            cohort = c,
            n_controls = n_controls,
            n_nonzero_weights = n_nonzero,
            max_weight = round(max_w, 3),
            top5_weight_share = round(top5_share, 3),
            herfindahl = round(hhi, 4)
          ))
        }
      }
    }
  }
}

# Save outputs ----------------------------------------------------------------

write_csv(sample_comparison, "results/diagnostics/diag-sample-comparison.csv")
write_csv(sample_summary, "results/diagnostics/diag-sample-summary.csv")
write_csv(pretrend_results, "results/diagnostics/diag-pretrend-results.csv")
write_csv(weight_results, "results/diagnostics/diag-weight-results.csv")
