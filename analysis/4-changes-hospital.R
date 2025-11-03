# Discrete-time hazards for closure, merger, and either -------------------

survival.dat <- est.dat %>%
  group_by(ID) %>% 
  mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% 
  ungroup() %>%
  filter(min_bedsize <= 50, distance>10) %>%
  mutate(ID = as.numeric(factor(ID)),
         treat_group = case_when(
           !is.na(eff_year) ~ eff_year,
           is.na(eff_year) & state_treat_year > year + max.es ~ 0,
           is.na(eff_year) & state_treat_year == 0 ~ 0,
           TRUE ~ NA 
         )) %>%
  filter(!is.na(year), !is.na(closed), !is.na(margin), 
         treat_group %in% c(0, 1999, 2000, 2001, 2002)) %>%
  group_by(ID) %>%
  mutate(panel_start_year = min(year, na.rm = TRUE),
         duration_years = year - panel_start_year + 1) %>% # Year 1 = first year observed
  ungroup() %>%
  mutate(
    treated = ifelse(treat_group > 0, 1, 0), 
    post_treat = ifelse(year >= treat_group & treat_group > 0, 1, 0)
   ) %>%  
  select(ID, year, closed, duration_years, treated, post_treat,
         BDTOT, distance, own_type, MSTATE)


cox_closure_did_mod <- coxph(
  Surv(duration_years, closed) ~ post_treat + treated + factor(year) +
                                   factor(own_type) + factor(MSTATE),
  data = survival.dat,
  cluster = survival.dat$ID
)

summary(cox_closure_did_mod)         






# packages
library(dplyr); library(tidyr); library(stringr)
library(fixest); library(sandwich); library(lmtest)
library(ggplot2)

# 0) Data: one row per hospital-year, y = closure in year (0/1), at-risk rows only
risk.close <- risk.close %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  mutate(event = as.integer(y == 1),
         at_risk = as.integer(cumsum(event) == 0)) %>%
  ungroup() %>%
  filter(at_risk == 1)

# 1) Define first-treatment (cohort) year; NA for never-treated
risk.close <- risk.close %>%
  group_by(ID) %>%
  mutate(cohort = ifelse(any(treated == 1L), min(year[treated == 1L]), NA_integer_)) %>%
  ungroup()

# 2) Fit nonlinear ETWFE (cloglog): year FE (baseline hazard) + cohort FE + post-only G×T
#    i(cohort) = cohort FE; i(cohort, year, drop = year < cohort) builds 1{G=r}×1{t=τ} for τ≥r
m_etwfe <- feglm(
  y ~ i(cohort) + i(cohort, year, drop = year < cohort) + BDTOT + distance | year,
  data   = risk.close,
  family = binomial(link = "cloglog"),
  cluster = "ID"
)

# 3) Pull clustered VCOV and coefficient vector
V   <- vcov(m_etwfe, cluster = "ID")         # cluster-robust VCOV
bet <- coef(m_etwfe)

# 4) Identify the post-only cohort×time coefficients and parse (cohort, year)
#    fixest names them like "i(cohort,year,drop = ...)#cohort::2008#year::2011"
cn   <- names(bet)
is_gt <- grepl("^i\\(cohort,year", cn)

gt_tbl <- tibble(term = cn[is_gt], b = bet[is_gt]) %>%
  mutate(cohort = as.integer(str_match(term, "cohort::(\\d+)")[,2]),
         year   = as.integer(str_match(term,   "year::(\\d+)")[,2])) %>%
  mutate(event_time = year - cohort) %>%
  arrange(event_time, cohort, year)

# 5) Build weights for aggregating to event time k = year - cohort.
#    Here: weight by number at risk in that (cohort, year) cell (simple, transparent).
atrisk_counts <- risk.close %>%
  filter(!is.na(cohort), year >= cohort) %>%
  count(cohort, year, name = "n_atrisk")

gt_tbl <- gt_tbl %>namely%
  left_join(atrisk_counts, by = c("cohort","year")) %>%
  mutate(n_atrisk = ifelse(is.na(n_atrisk), 0L, n_atrisk))

# 6) For each event time k, form linear combination w'β and its clustered SE √(w' V w)
#    Build a weight vector aligned with coefficient order in V.
k_vals <- sort(unique(gt_tbl$event_time))          # typically k = 0,1,2,...
agg_list <- lapply(k_vals, function(k){
  rows_k <- gt_tbl$event_time == k & gt_tbl$n_atrisk > 0
  idx    <- match(gt_tbl$term[rows_k], colnames(V))  # positions in VCOV
  w_full <- rep(0, length(bet)); names(w_full) <- names(bet)
  # normalize weights to sum to 1 within k (can choose other schemes)
  w_k    <- gt_tbl$n_atrisk[rows_k]
  w_k    <- w_k / sum(w_k)
  w_full[idx] <- w_k
  b_k    <- sum(w_full * bet)
  se_k   <- sqrt( as.numeric(t(w_full) %*% V %*% w_full) )
  tibble(event_time = k, b = b_k, se = se_k)
})
es_sum <- bind_rows(agg_list) %>%
  mutate(ci_lo = b - 1.96*se,
         ci_hi = b + 1.96*se,
         HR    = exp(b),
         HR_lo = exp(ci_lo),
         HR_hi = exp(ci_hi))

# 7) Plot as a standard event-study (log-hazard ratio and/or hazard ratio)
p1 <- ggplot(es_sum, aes(x = event_time, y = b)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  labs(x = "Event time (years since first treatment, t=0 = first treated year)",
       y = "log hazard ratio (vs t = -1)",
       title = "Event-study: discrete-time hazard (cloglog ETWFE)") +
  theme_minimal()

p2 <- ggplot(es_sum, aes(x = event_time, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point() +
  geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.2) +
  labs(x = "Event time (years since first treatment, t=0 = first treated year)",
       y = "Hazard ratio (vs t = -1)") +
  theme_minimal()

print(p1); print(p2)
