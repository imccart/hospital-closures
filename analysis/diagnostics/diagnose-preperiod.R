# Diagnose Pre-Period Data Availability
#
# Run after _run-analysis.r (needs est.dat and stack.hosp)
#
# Investigates why 5 pre-periods causes NA CIs for financial outcomes
# in the 1999 cohort

library(tidyverse)
library(BMisc)

# =============================================================================
# 1. Margin data availability by year for the 1999 cohort
# =============================================================================

cat("\n=== Margin Data Availability by Year (1999 Cohort) ===\n")
stack.hosp %>%
  filter(stack_group == 1999) %>%
  group_by(year) %>%
  summarize(
    n_total = n_distinct(ID),
    n_with_margin = sum(!is.na(margin)),
    pct_margin = round(mean(!is.na(margin)) * 100, 1)
  ) %>%
  print(n = 20)

# =============================================================================
# 2. Treated vs control with complete data for 5 pre-periods
# =============================================================================

cat("\n=== Complete Data by Treatment Status (5 Pre-Periods) ===\n")
completeness_5 <- stack.hosp %>%
  filter(stack_group == 1999, stacked_event_time >= -5) %>%
  group_by(ID) %>%
  summarize(
    treated = first(treated),
    years_with_margin = sum(!is.na(margin)),
    total_years = n()
  ) %>%
  mutate(complete = years_with_margin == total_years) %>%
  group_by(treated) %>%
  summarize(
    n_hospitals = n(),
    n_complete = sum(complete),
    pct_complete = round(mean(complete) * 100, 1)
  )
print(completeness_5)

# =============================================================================
# 3. Compare to 4 and 3 pre-periods
# =============================================================================

cat("\n=== Complete Data by Treatment Status (4 Pre-Periods) ===\n")
completeness_4 <- stack.hosp %>%
  filter(stack_group == 1999, stacked_event_time >= -4) %>%
  group_by(ID) %>%
  summarize(
    treated = first(treated),
    years_with_margin = sum(!is.na(margin)),
    total_years = n()
  ) %>%
  mutate(complete = years_with_margin == total_years) %>%
  group_by(treated) %>%
  summarize(
    n_hospitals = n(),
    n_complete = sum(complete),
    pct_complete = round(mean(complete) * 100, 1)
  )
print(completeness_4)

cat("\n=== Complete Data by Treatment Status (3 Pre-Periods) ===\n")
completeness_3 <- stack.hosp %>%
  filter(stack_group == 1999, stacked_event_time >= -3) %>%
  group_by(ID) %>%
  summarize(
    treated = first(treated),
    years_with_margin = sum(!is.na(margin)),
    total_years = n()
  ) %>%
  mutate(complete = years_with_margin == total_years) %>%
  group_by(treated) %>%
  summarize(
    n_hospitals = n(),
    n_complete = sum(complete),
    pct_complete = round(mean(complete) * 100, 1)
  )
print(completeness_3)

# =============================================================================
# 4. Balanced panel sizes
# =============================================================================

cat("\n=== Balanced Panel Sizes ===\n")

for (pre_pd in c(3, 4, 5)) {
  synth.test <- stack.hosp %>%
    filter(stack_group == 1999, !is.na(margin), stacked_event_time >= -pre_pd) %>%
    select(ID, year, margin, treated)

  balanced.test <- as_tibble(makeBalancedPanel(synth.test, idname = "ID", tname = "year"))

  bal_summary <- balanced.test %>%
    group_by(treated) %>%
    summarize(n = n_distinct(ID), .groups = "drop")

  cat(sprintf("\nPre-period = %d:\n", pre_pd))
  cat(sprintf("  Before balance: %d hospitals\n", n_distinct(synth.test$ID)))
  cat(sprintf("  After balance:  %d hospitals\n", n_distinct(balanced.test$ID)))
  cat(sprintf("  Treated: %d, Control: %d\n",
              bal_summary$n[bal_summary$treated == 1],
              bal_summary$n[bal_summary$treated == 0]))
}

# =============================================================================
# 5. Which years are causing the drop?
# =============================================================================

cat("\n=== Hospitals Lost at Each Pre-Period Year ===\n")

# Get IDs with complete data at each pre-period length
ids_by_preperiod <- map(3:5, function(pre_pd) {
  stack.hosp %>%
    filter(stack_group == 1999, !is.na(margin), stacked_event_time >= -pre_pd) %>%
    select(ID, year, margin, treated) %>%
    makeBalancedPanel(idname = "ID", tname = "year") %>%
    as_tibble() %>%
    distinct(ID) %>%
    pull(ID)
}) %>%
  set_names(c("pre3", "pre4", "pre5"))

cat(sprintf("IDs with 3 pre-periods: %d\n", length(ids_by_preperiod$pre3)))
cat(sprintf("IDs with 4 pre-periods: %d\n", length(ids_by_preperiod$pre4)))
cat(sprintf("IDs with 5 pre-periods: %d\n", length(ids_by_preperiod$pre5)))

# Which IDs are lost going from 4 to 5?
lost_4to5 <- setdiff(ids_by_preperiod$pre4, ids_by_preperiod$pre5)
cat(sprintf("\nIDs lost going from 4 to 5 pre-periods: %d\n", length(lost_4to5)))

# Why are they lost? Check which year is missing
if (length(lost_4to5) > 0) {
  cat("\nMissing years for hospitals lost at 5 pre-periods:\n")
  stack.hosp %>%
    filter(stack_group == 1999, ID %in% lost_4to5, stacked_event_time >= -5) %>%
    group_by(ID) %>%
    summarize(
      treated = first(treated),
      missing_years = paste(year[is.na(margin)], collapse = ", ")
    ) %>%
    count(treated, missing_years) %>%
    arrange(desc(n)) %>%
    print(n = 20)
}
