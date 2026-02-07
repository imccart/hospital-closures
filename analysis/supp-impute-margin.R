## Impute margins within stack.hosp (post-stacking)
## Assumes stack.hosp already exists in the environment.
## Produces stack.hosp with a new column: margin_imp
stack.hosp_base <- stack.hosp %>%
  distinct(ID, year, .keep_all = TRUE) %>%
  group_by(ID) %>%
  mutate(
    mcare_share = mcare_discharges / tot_discharges,
    min_bedsize = min(BDTOT, na.rm = TRUE)
  ) %>%
  ungroup()

train_97 <- stack.hosp_base %>%
  filter(year == 1997,
         min_bedsize <= bed.cut,
         !is.na(net_pat_rev), !is.na(tot_operating_exp), !is.na(tot_charges),
         !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
         net_pat_rev > 0, tot_charges > 0, tot_operating_exp > 0,
         (is.na(cah) | cah == 0)) %>%
  mutate(
    margin_true = (net_pat_rev - tot_operating_exp) / net_pat_rev,
    margin_charges = (tot_charges - tot_operating_exp) / tot_charges,
    charge_to_rev = tot_charges / net_pat_rev
  ) %>%
  filter(is.finite(margin_true), is.finite(margin_charges), is.finite(charge_to_rev),
         margin_true > -1, margin_true < 1,
         charge_to_rev > 0.5, charge_to_rev < 10)

ctr_avg <- stack.hosp_base %>%
  filter(year %in% 1997:1998,
         min_bedsize <= bed.cut,
         !is.na(tot_charges), !is.na(net_pat_rev),
         tot_charges > 0, net_pat_rev > 0) %>%
  mutate(charge_to_rev = tot_charges / net_pat_rev) %>%
  filter(charge_to_rev > 0.5, charge_to_rev < 10) %>%
  group_by(ID) %>%
  summarize(charge_to_rev_avg = mean(charge_to_rev), .groups = "drop")

ols_fit <- lm(margin_true ~ margin_charges + mcare_share + charge_to_rev, data = train_97)

train_ranges <- train_97 %>%
  summarize(
    mc_min = min(margin_charges),
    mc_max = max(margin_charges),
    ms_min = min(mcare_share),
    ms_max = max(mcare_share),
    ctr_min = min(charge_to_rev),
    ctr_max = max(charge_to_rev)
  )

pre_range <- stack.hosp_base %>%
  filter(year %in% 1994:1996,
         min_bedsize <= bed.cut,
         !is.na(tot_charges), !is.na(tot_operating_exp),
         !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
         tot_charges > 0, tot_operating_exp > 0) %>%
  left_join(ctr_avg, by = "ID") %>%
  mutate(
    est_total_opex = tot_operating_exp / mcare_share,
    margin_charges = (tot_charges - est_total_opex) / tot_charges,
    charge_to_rev = charge_to_rev_avg
  ) %>%
  filter(!is.na(charge_to_rev),
         is.finite(margin_charges), is.finite(charge_to_rev),
         charge_to_rev > 0.5, charge_to_rev < 10) %>%
  filter(
    margin_charges >= train_ranges$mc_min, margin_charges <= train_ranges$mc_max,
    mcare_share >= train_ranges$ms_min, mcare_share <= train_ranges$ms_max,
    charge_to_rev >= train_ranges$ctr_min, charge_to_rev <= train_ranges$ctr_max
  )

pre_range <- pre_range %>%
  mutate(margin_hat = as.numeric(predict(ols_fit, newdata = pre_range)))

stack.hosp <- stack.hosp %>%
  left_join(pre_range %>% select(ID, year, margin_hat), by = c("ID", "year")) %>%
  mutate(margin_imp = ifelse(is.na(margin), margin_hat, margin))

  ## Diagnostics ------------------------------------------------------------
  stack.hosp_diag <- stack.hosp %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup()

  ## 1) 1998 out-of-sample accuracy on stack.hosp sample
  test_98 <- stack.hosp_base %>%
    filter(year == 1998,
           min_bedsize <= bed.cut,
           !is.na(net_pat_rev), !is.na(tot_operating_exp), !is.na(tot_charges),
           !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
           net_pat_rev > 0, tot_charges > 0, tot_operating_exp > 0,
           (is.na(cah) | cah == 0)) %>%
    mutate(
      margin_true = (net_pat_rev - tot_operating_exp) / net_pat_rev,
      margin_charges = (tot_charges - tot_operating_exp) / tot_charges,
      charge_to_rev = tot_charges / net_pat_rev
    ) %>%
    filter(is.finite(margin_true), is.finite(margin_charges), is.finite(charge_to_rev),
           margin_true > -1, margin_true < 1,
           charge_to_rev > 0.5, charge_to_rev < 10)

  if (nrow(test_98) > 0) {
    test_98 <- test_98 %>%
      mutate(margin_hat = as.numeric(predict(ols_fit, newdata = test_98)))

    diag_metrics <- tibble(
      sample = "1998_holdout",
      n = nrow(test_98),
      rmse = round(sqrt(mean((test_98$margin_hat - test_98$margin_true)^2)), 4),
      mae = round(mean(abs(test_98$margin_hat - test_98$margin_true)), 4),
      corr = round(cor(test_98$margin_hat, test_98$margin_true), 4),
      bias = round(mean(test_98$margin_hat - test_98$margin_true), 4)
    )
  } else {
    diag_metrics <- tibble(
      sample = "1998_holdout",
      n = 0, rmse = NA_real_, mae = NA_real_, corr = NA_real_, bias = NA_real_
    )
  }

  ## 2) Imputation coverage by year (stack.hosp sample)
  diag_coverage <- stack.hosp_diag %>%
    mutate(imputed = is.na(margin) & !is.na(margin_imp)) %>%
    group_by(year) %>%
    summarize(
      n = n(),
      missing_margin = sum(is.na(margin)),
      imputed_n = sum(imputed),
      imputed_share = round(imputed_n / n, 4),
      .groups = "drop"
    ) %>%
    arrange(year)

  ## 2b) Imputation coverage by year (bed-cut sample only)
  diag_coverage_bedcut <- stack.hosp_diag %>%
    filter(min_bedsize <= bed.cut) %>%
    mutate(imputed = is.na(margin) & !is.na(margin_imp)) %>%
    group_by(year) %>%
    summarize(
      n = n(),
      missing_margin = sum(is.na(margin)),
      imputed_n = sum(imputed),
      imputed_share = round(ifelse(missing_margin == 0, 0, imputed_n / missing_margin), 4),
      .groups = "drop"
    ) %>%
    arrange(year)

  ## 3) Distribution of imputed margins by year (where imputed)
  diag_dist_pre <- stack.hosp_diag %>%
    filter(is.na(margin) & !is.na(margin_imp)) %>%
    group_by(year) %>%
    summarize(
      n = n(),
      mean = round(mean(margin_imp), 4),
      median = round(median(margin_imp), 4),
      sd = round(sd(margin_imp), 4),
      p10 = round(quantile(margin_imp, 0.10), 4),
      p90 = round(quantile(margin_imp, 0.90), 4),
      .groups = "drop"
    ) %>%
    arrange(year)

  ## 4) Compare imputed vs observed means (no recentering applied)
  obs_97_98 <- stack.hosp_diag %>%
    filter(year %in% 1997:1998,
           min_bedsize <= bed.cut,
           !is.na(margin),
           (is.na(cah) | cah == 0)) %>%
    summarize(mean_obs = mean(margin), .groups = "drop")

  imp_pre <- stack.hosp_diag %>%
    filter(year %in% 1994:1996,
           is.na(margin) & !is.na(margin_imp)) %>%
    group_by(year) %>%
    summarize(mean_imp = mean(margin_imp), .groups = "drop")

  level_comparison <- imp_pre %>%
    mutate(mean_obs_97_98 = ifelse(nrow(obs_97_98) == 0, NA_real_, obs_97_98$mean_obs),
           gap = mean_imp - mean_obs_97_98)
cat("Imputation complete. New column: margin_imp\n")
cat("Missing margin before:", sum(is.na(stack.hosp$margin)), "\n")
cat("Missing margin after:", sum(is.na(stack.hosp$margin_imp)), "\n")
cat("\nDiagnostics: OLS coefficients\n")
print(summary(ols_fit)$coefficients)
cat("\nDiagnostics: 1998 holdout performance\n")
print(diag_metrics, n = Inf, width = Inf)
cat("\nDiagnostics: Imputation coverage by year\n")
print(diag_coverage, n = Inf, width = Inf)
cat("\nDiagnostics: Imputation coverage by year (bed-cut sample)\n")
print(diag_coverage_bedcut, n = Inf, width = Inf)
cat("\nDiagnostics: Imputed margin distribution by year\n")
print(diag_dist_pre, n = Inf, width = Inf)
cat("\nDiagnostics: Imputed vs observed mean comparison\n")
print(level_comparison, n = Inf, width = Inf)
