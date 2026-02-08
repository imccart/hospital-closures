## Impute margins within stack.hosp (post-stacking)
## Assumes stack.hosp already exists in the environment.
## Produces stack.hosp with a new column: margin_imp
##
## Two-stage approach:
## Stage 1 (shift model): In 1998 overlap data where both PPS-equivalent and HCRIS
##   variables are observed, model the translation from PPS-measured margin_charges
##   to HCRIS-measured margin_charges. This bridges the measurement systems.
## Stage 2 (margin model): Predict true margin from HCRIS-measured predictors
##   (margin_charges, mcare_share, charge_to_rev). Trained on 1997 HCRIS data.
## Application: For 1994-1996, construct PPS-measured predictors, apply Stage 1
##   to get HCRIS-equivalent margin_charges, then apply Stage 2 for true margin.

## De-duplicate for model fitting (one row per hospital-year)
stack.hosp_base <- stack.hosp %>%
  distinct(ID, year, .keep_all = TRUE) %>%
  group_by(ID) %>%
  mutate(
    mcare_share = mcare_discharges / tot_discharges,
    min_bedsize = min(BDTOT, na.rm = TRUE)
  ) %>%
  ungroup()

## ── Stage 1: Shift model ─────────────────────────────────────────────────────
## Train on 1998 v1996 rows where both measurement systems are observed
shift_train <- stack.hosp_base %>%
  filter(year == 1998,
         min_bedsize <= bed.cut,
         !is.na(tot_charges), !is.na(tot_operating_exp),
         tot_charges > 0, tot_operating_exp > 0,
         !is.na(pps_mcare_cost), !is.na(pps_pgm_cost),
         !is.na(pps_ip_charges), !is.na(pps_op_charges),
         !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
         (is.na(cah) | cah == 0)) %>%
  mutate(
    pps_total_charges = pps_ip_charges + pps_op_charges,
    margin_charges_hcris = (tot_charges - tot_operating_exp) / tot_charges,
    margin_charges_pps = (pps_total_charges - pps_mcare_cost) / pps_total_charges,
    pgm_cost_ratio = pps_pgm_cost / pps_total_charges
  ) %>%
  filter(pps_total_charges > 0,
         is.finite(margin_charges_hcris), is.finite(margin_charges_pps),
         is.finite(pgm_cost_ratio))

shift_fit <- lm(margin_charges_hcris ~ margin_charges_pps + pgm_cost_ratio + mcare_share,
                data = shift_train)

## ── Stage 2: Margin model ────────────────────────────────────────────────────
## Train on 1997 HCRIS data (HCRIS-measured predictors → true margin)
margin_train <- stack.hosp_base %>%
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

margin_fit <- lm(margin_true ~ margin_charges + mcare_share + charge_to_rev,
                 data = margin_train)

## Hospital-specific charge_to_rev from 1997-1998 (for prediction step)
ctr_avg <- stack.hosp_base %>%
  filter(year %in% 1997:1998,
         min_bedsize <= bed.cut,
         !is.na(tot_charges), !is.na(net_pat_rev),
         tot_charges > 0, net_pat_rev > 0) %>%
  mutate(charge_to_rev = tot_charges / net_pat_rev) %>%
  filter(charge_to_rev > 0.5, charge_to_rev < 10) %>%
  group_by(ID) %>%
  summarize(charge_to_rev_avg = mean(charge_to_rev), .groups = "drop")

## Training-support bounds (Stage 2 predictors)
margin_ranges <- margin_train %>%
  summarize(
    mc_min = min(margin_charges), mc_max = max(margin_charges),
    ms_min = min(mcare_share), ms_max = max(mcare_share),
    ctr_min = min(charge_to_rev), ctr_max = max(charge_to_rev)
  )

## ── Apply to 1994-1996 ──────────────────────────────────────────────────────
pre_range <- stack.hosp_base %>%
  filter(year %in% 1994:1996,
         min_bedsize <= bed.cut,
         !is.na(tot_charges), !is.na(tot_operating_exp),
         !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
         tot_charges > 0, tot_operating_exp > 0) %>%
  mutate(
    ## PPS-measured predictors (opertots/totalg definitions)
    margin_charges_pps = (tot_charges - tot_operating_exp) / tot_charges,
    pgm_cost_ratio = tot_operating_exp / tot_charges
  ) %>%
  filter(is.finite(margin_charges_pps), is.finite(pgm_cost_ratio))

## Stage 1: translate PPS margin_charges → HCRIS-equivalent
pre_range <- pre_range %>%
  mutate(margin_charges = as.numeric(predict(shift_fit, newdata = pre_range)))

## Add charge_to_rev from 1997-1998 hospital-specific average
pre_range <- pre_range %>%
  left_join(ctr_avg, by = "ID") %>%
  mutate(charge_to_rev = charge_to_rev_avg) %>%
  filter(!is.na(charge_to_rev))

## Enforce Stage 2 training-support bounds
pre_range <- pre_range %>%
  filter(margin_charges >= margin_ranges$mc_min, margin_charges <= margin_ranges$mc_max,
         mcare_share >= margin_ranges$ms_min, mcare_share <= margin_ranges$ms_max,
         charge_to_rev >= margin_ranges$ctr_min, charge_to_rev <= margin_ranges$ctr_max)

## Stage 2: predict true margin
pre_range <- pre_range %>%
  mutate(margin_hat = as.numeric(predict(margin_fit, newdata = pre_range)))

## Join back to stack.hosp
stack.hosp <- stack.hosp %>%
  left_join(pre_range %>% select(ID, year, margin_hat), by = c("ID", "year")) %>%
  mutate(margin_imp = ifelse(is.na(margin), margin_hat, margin))

## ── Diagnostics ──────────────────────────────────────────────────────────────
stack.hosp_diag <- stack.hosp %>%
  group_by(ID) %>%
  mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
  ungroup()

## 1) Full-pipeline validation on 1999 holdout
##    Apply Stage 1 (pps → hcris margin_charges) then Stage 2 (→ true margin)
test_99 <- stack.hosp_base %>%
  filter(year == 1999,
         min_bedsize <= bed.cut,
         !is.na(net_pat_rev), !is.na(tot_operating_exp), !is.na(tot_charges),
         net_pat_rev > 0, tot_charges > 0, tot_operating_exp > 0,
         !is.na(pps_mcare_cost), !is.na(pps_pgm_cost),
         !is.na(pps_ip_charges), !is.na(pps_op_charges),
         !is.na(mcare_share), mcare_share > 0, mcare_share <= 1,
         (is.na(cah) | cah == 0)) %>%
  mutate(
    margin_true = (net_pat_rev - tot_operating_exp) / net_pat_rev,
    pps_total_charges = pps_ip_charges + pps_op_charges,
    margin_charges_pps = (pps_total_charges - pps_mcare_cost) / pps_total_charges,
    pgm_cost_ratio = pps_pgm_cost / pps_total_charges,
    charge_to_rev = tot_charges / net_pat_rev
  ) %>%
  filter(pps_total_charges > 0,
         is.finite(margin_true), is.finite(margin_charges_pps), is.finite(pgm_cost_ratio),
         margin_true > -1, margin_true < 1,
         charge_to_rev > 0.5, charge_to_rev < 10)

if (nrow(test_99) > 0) {
  ## Apply both stages (sequential — Stage 2 needs Stage 1 output)
  test_99 <- test_99 %>%
    mutate(margin_charges = as.numeric(predict(shift_fit, newdata = test_99)))
  test_99 <- test_99 %>%
    mutate(margin_hat = as.numeric(predict(margin_fit, newdata = test_99)))

  diag_metrics <- tibble(
    sample = "1999_holdout_full_pipeline",
    n = nrow(test_99),
    rmse = round(sqrt(mean((test_99$margin_hat - test_99$margin_true)^2)), 4),
    mae = round(mean(abs(test_99$margin_hat - test_99$margin_true)), 4),
    corr = round(cor(test_99$margin_hat, test_99$margin_true), 4),
    bias = round(mean(test_99$margin_hat - test_99$margin_true), 4)
  )

  ## Also report Stage 1 accuracy separately
  test_99 <- test_99 %>%
    mutate(margin_charges_hcris = (tot_charges - tot_operating_exp) / tot_charges)

  diag_shift <- tibble(
    sample = "1999_shift_model_only",
    n = nrow(test_99),
    rmse = round(sqrt(mean((test_99$margin_charges - test_99$margin_charges_hcris)^2)), 4),
    mae = round(mean(abs(test_99$margin_charges - test_99$margin_charges_hcris)), 4),
    corr = round(cor(test_99$margin_charges, test_99$margin_charges_hcris), 4),
    bias = round(mean(test_99$margin_charges - test_99$margin_charges_hcris), 4)
  )
} else {
  diag_metrics <- tibble(
    sample = "1999_holdout_full_pipeline",
    n = 0, rmse = NA_real_, mae = NA_real_, corr = NA_real_, bias = NA_real_)
  diag_shift <- tibble(
    sample = "1999_shift_model_only",
    n = 0, rmse = NA_real_, mae = NA_real_, corr = NA_real_, bias = NA_real_)
}

## 2) Imputation coverage by year (bed-cut sample)
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

## 3) Distribution of imputed margins by year
diag_dist <- stack.hosp_diag %>%
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

## 4) Compare imputed vs observed means
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
cat("\nDiagnostics: Stage 1 (shift model) coefficients\n")
print(summary(shift_fit)$coefficients)
cat("\nDiagnostics: Stage 2 (margin model) coefficients\n")
print(summary(margin_fit)$coefficients)
cat("\nDiagnostics: Stage 1 accuracy (1999 holdout)\n")
print(diag_shift, n = Inf, width = Inf)
cat("\nDiagnostics: Full pipeline accuracy (1999 holdout)\n")
print(diag_metrics, n = Inf, width = Inf)
cat("\nDiagnostics: Imputation coverage by year (bed-cut sample)\n")
print(diag_coverage_bedcut, n = Inf, width = Inf)
cat("\nDiagnostics: Imputed margin distribution by year\n")
print(diag_dist, n = Inf, width = Inf)
cat("\nDiagnostics: Imputed vs observed mean comparison\n")
print(level_comparison, n = Inf, width = Inf)
