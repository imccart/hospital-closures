# extract-manual-990.R
# Generate list of nonprofit hospitals needing manual Form 990 data collection
# Requires: est.dat from analysis/_run-analysis.r (run through line 152+)
# Output: data/output/manual-990-ids.csv

library(tidyverse)

bed.cut <- 50  # match analysis parameter

## Get 1999 cohort treated hospitals after bed.cut filter
treated_1999 <- est.dat %>%
  filter(eff_year == 1999) %>%
  group_by(ID) %>%
  mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(min_bedsize <= bed.cut)

## Identify nonprofits (CNTRL 21 = church, 23 = other nonprofit)
nonprofits_1999 <- treated_1999 %>%
  group_by(ID) %>%
  summarize(
    MNAME = first(na.omit(MNAME)),
    MSTATE = first(na.omit(MSTATE)),
    MLOCZIP = first(na.omit(MLOCZIP)),
    MLOCCITY = first(na.omit(MLOCCITY)),
    CNTRL = first(na.omit(CNTRL)),
    eff_year = first(eff_year),
    min_bedsize = first(min_bedsize),
    # Form 990 coverage
    ever_margin_990 = any(!is.na(margin_990)),
    first_margin_990 = ifelse(any(!is.na(margin_990)), min(year[!is.na(margin_990)]), NA),
    last_margin_990 = ifelse(any(!is.na(margin_990)), max(year[!is.na(margin_990)]), NA),
    n_years_990 = sum(!is.na(margin_990)),
    # Years we need (1994-1998 for pre-period coverage with 1999 treatment)
    has_1994 = any(year == 1994 & !is.na(margin_990)),
    has_1995 = any(year == 1995 & !is.na(margin_990)),
    has_1996 = any(year == 1996 & !is.na(margin_990)),
    has_1997 = any(year == 1997 & !is.na(margin_990)),
    has_1998 = any(year == 1998 & !is.na(margin_990)),
    .groups = "drop"
  ) %>%
  filter(CNTRL %in% c(21, 23))  # nonprofits only

cat("=== Nonprofit hospitals in 1999 treated cohort ===\n")
cat("Total nonprofits:", nrow(nonprofits_1999), "\n")
cat("With any 990 data:", sum(nonprofits_1999$ever_margin_990), "\n")
cat("Without any 990 data:", sum(!nonprofits_1999$ever_margin_990), "\n")
cat("With pre-1999 990 data:", sum(!is.na(nonprofits_1999$first_margin_990) &
                                     nonprofits_1999$first_margin_990 < 1999), "\n")

## Identify which hospitals need manual lookup
## Priority: hospitals missing 1994-1998 data (the pre-period years we need)
manual_lookup <- nonprofits_1999 %>%
  mutate(
    needs_1994 = !has_1994,
    needs_1995 = !has_1995,
    needs_1996 = !has_1996,
    needs_1997 = !has_1997,
    needs_1998 = !has_1998,
    n_years_needed = needs_1994 + needs_1995 + needs_1996 + needs_1997 + needs_1998,
    ownership_type = case_when(
      CNTRL == 21 ~ "Nonprofit - Church",
      CNTRL == 23 ~ "Nonprofit - Other",
      TRUE ~ "Other"
    )
  ) %>%
  filter(n_years_needed > 0) %>%  # only hospitals missing at least one pre-period year
  select(
    ID, MNAME, MSTATE, MLOCCITY, MLOCZIP, ownership_type, min_bedsize, eff_year,
    ever_margin_990, first_margin_990, n_years_990,
    needs_1994, needs_1995, needs_1996, needs_1997, needs_1998, n_years_needed
  ) %>%
  arrange(MSTATE, MNAME)

cat("\n=== Hospitals needing manual 990 lookup ===\n")
cat("Total:", nrow(manual_lookup), "\n")
cat("Missing all 5 years (1994-1998):", sum(manual_lookup$n_years_needed == 5), "\n")
cat("Missing 1-4 years:", sum(manual_lookup$n_years_needed < 5), "\n")

cat("\n--- By state ---\n")
print(manual_lookup %>% count(MSTATE, name = "n_hospitals"), n = Inf)

## Write output
write_csv(manual_lookup, "data/output/manual-990-ids.csv")
cat("\nOutput written to data/output/manual-990-ids.csv\n")
cat("Total data points to collect:", sum(manual_lookup$n_years_needed),
    "(hospitals × missing years)\n")
