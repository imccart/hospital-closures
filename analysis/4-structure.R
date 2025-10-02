# Structural Estimation: Hospital Exit and Merger Decisions -----------------------

# does CAH designation improve hospital survival through continuation, or does it block efficiency-enhancing mergers, ultimately increasing closures?

# Step 1: Prepare State Space -----------------------------------------------------
# d_it ∈ {0, 1, 2} for {continue, close, merge}; add year FE; construct COPA; keep required vars

vars_used <- c("cah", "distance", "margin_1", "BDTOT", "year", "year_f", "copa", "d_it")

struc.dat <- est.dat %>%
  mutate(
    d_it = case_when(
      closed == 1 ~ 1L,
      merged == 1 ~ 2L,
      TRUE ~ 0L
    ),
    year_f = factor(year)   # year FE (used throughout)
  ) %>%
  left_join(copa.dat, by = "MSTATE") %>%
  mutate(
    copa = as.integer(!is.na(copa_start) & year >= copa_start & year <= copa_end)
  ) %>%
  filter(complete.cases(across(all_of(vars_used))))  

# Step 2: Estimate Choice Probabilities (Nested Logit to relax IIA) ---------------

# Alternatives and long data with a choice-situation id (chid)
alts <- c("continue","close","merge")

long <- struc.dat %>%
  mutate(
    choice = factor(c("continue","close","merge")[d_it + 1L], levels = alts),
    chid   = paste0(ID, "_", year)   # one row per choice situation
  ) %>%
  select(all_of(vars_used), choice, chid, ID) %>%
  expand_grid(alt = alts) %>%
  mutate(chosen = as.integer(choice == alt))

mnl.data <- mlogit.data(
  long,
  choice  = "chosen",
  shape   = "long",
  alt.var = "alt",
  chid.var= "chid"
)

mnl.model <- mlogit(
  chosen ~ 1 | cah + BDTOT + distance + margin_1,
  data     = mnl.data,
  reflevel = "continue"
)

# Map fitted probabilities back to (ID, year), pivot wide, and trim
key <- long %>% distinct(chid, ID, year)

prob_long <- cbind(index(ml_data), p = fitted(mnl_model)) %>%
  as_tibble() %>%
  select(chid, alt, p) %>%
  left_join(key, by = "chid")

P_hat <- prob_long %>%
  pivot_wider(names_from = alt, values_from = p) %>%
  transmute(
    ID, year,
    p0 = pmax(pmin(continue, 1 - 1e-6), 1e-6),  # continue
    p1 = pmax(pmin(close,    1 - 1e-6), 1e-6),  # close
    p2 = pmax(pmin(merge,    1 - 1e-6), 1e-6)   # merge
  )

struct_data <- struct_data %>% left_join(P_hat, by = c("ID","year"))

# Step 3–5: Inversion (full sample), Emax, and Transitions ------------------------

# Hotz–Miller log-odds on FULL sample (no conditioning on d_it)
log_odds <- struct_data %>%
  mutate(
    v1_minus_v0 = log(p1) - log(p0),   # close vs continue
    v2_minus_v0 = log(p2) - log(p0)    # merge vs continue
  )

# Inversion regressions with year & market FE; cluster at market
close_model_base <- feols(v1_minus_v0 ~ distance + beds + copa | year_f + MSTATE,
                          data = log_odds, cluster = ~ MSTATE)

# Baseline merge model WITHOUT 'cah' (so CF2 can "turn it on")
merge_model_base <- feols(v2_minus_v0 ~ distance + beds + copa | year_f + MSTATE,
                          data = log_odds, cluster = ~ MSTATE)

# Extended merge model WITH 'cah' (for CF2)
merge_model_ext  <- feols(v2_minus_v0 ~ distance + beds + copa + cah | year_f + MSTATE,
                          data = log_odds, cluster = ~ MSTATE)

# Deterministic utilities up to v0=0
log_odds <- log_odds %>%
  mutate(
    v1_base = as.numeric(predict(close_model_base, newdata = log_odds)),
    v2_base = as.numeric(predict(merge_model_base, newdata = log_odds)),
    v2_ext  = as.numeric(predict(merge_model_ext,  newdata = log_odds)),
    v0 = 0
  )

# Flow payoff π(s): fit with FE, use fitted values (or set to 0 if you prefer)
pi_model <- feols(margin_1 ~ cah + beds + distance | year_f + MSTATE,
                  data = struct_data, cluster = ~ MSTATE)
struct_data <- struct_data %>% mutate(pi_hat = as.numeric(fitted(pi_model)))

# Helper: log-sum-exp
logsumexp3 <- function(a,b,c){ m <- pmax(a, pmax(b,c)); m + log(exp(a-m)+exp(b-m)+exp(c-m)) }

# Baseline one-shot Emax from base utilities (for initializing transitions)
init <- log_odds %>%
  transmute(ID, year, v0 = 0, v1 = v1_base, v2 = v2_base,
            Emax = logsumexp3(v0, v1_base, v2_base))

# Estimate transitions for Emax by action
# d=0 (continuation): E[Emax_{t+1} | s_t, d=0]
trans0 <- struct_data %>%
  select(ID, year, cah, distance, beds, year_f, d_it) %>%
  left_join(init, by = c("ID","year")) %>%
  group_by(ID) %>% arrange(year) %>%
  mutate(lead_Emax = lead(Emax)) %>%
  ungroup() %>%
  filter(d_it == 0, !is.na(lead_Emax))

T0 <- feols(lead_Emax ~ cah + distance + beds | year_f, data = trans0, cluster = ~ ID)

# d=2 (merger): try to estimate; if too few, calibrate relative to T0
trans2 <- struct_data %>%
  select(ID, year, cah, distance, beds, year_f, d_it) %>%
  left_join(init, by = c("ID","year")) %>%
  group_by(ID) %>% arrange(year) %>%
  mutate(lead_Emax = lead(Emax)) %>%
  ungroup() %>%
  filter(d_it == 2, !is.na(lead_Emax))

if(nrow(trans2) >= 50){
  T2 <- feols(lead_Emax ~ cah + distance + beds | year_f, data = trans2, cluster = ~ ID)
  E_Emax_next_d2_fun <- function(df) as.numeric(predict(T2, newdata = df))
} else {
  merger_delta <- -0.10  # calibration if few mergers observed
  E_Emax_next_d2_fun <- function(df) as.numeric(predict(T0, newdata = df)) + merger_delta
}

# Closure (d=1) absorbing
E_Emax_next_d1_const <- 0



# CCP fixed-point solver and counterfactuals --------------------------------------
beta <- 0.95

# Utility builders (baseline vs. extended merge)
build_v_base <- function(df){
  v1 <- as.numeric(predict(close_model_base, newdata = df))
  v2 <- as.numeric(predict(merge_model_base, newdata = df))
  list(v0 = rep(0, nrow(df)), v1 = v1, v2 = v2)
}
build_v_ext  <- function(df){
  v1 <- as.numeric(predict(close_model_base, newdata = df))
  v2 <- as.numeric(predict(merge_model_ext,  newdata = df))
  list(v0 = rep(0, nrow(df)), v1 = v1, v2 = v2)
}

# Expectation under d=0 and d=2
E_Emax_next_d0_fun <- function(df) as.numeric(predict(T0, newdata = df))

iterate_ccp <- function(df, build_v_fun, maxit = 100, tol = 1e-6){
  # df must contain: ID, year, year_f, cah, distance, beds, MSTATE, pi_hat
  v <- build_v_fun(df)
  Emax <- logsumexp3(v$v0, v$v1, v$v2)
  for(it in 1:maxit){
    # expected Emax next period under actions
    E_Emax_next_d0 <- E_Emax_next_d0_fun(df)
    E_Emax_next_d1 <- E_Emax_next_d1_const
    E_Emax_next_d2 <- E_Emax_next_d2_fun(df)

    v0_new <- df$pi_hat + beta * E_Emax_next_d0
    v1_new <- v$v1 + beta * E_Emax_next_d1         # closure: no future value
    v2_new <- v$v2 + beta * E_Emax_next_d2

    Emax_new <- logsumexp3(v0_new, v1_new, v2_new)
    if(max(abs(Emax_new - Emax)) < tol) break

    v$v0 <- v0_new; v$v1 <- v1_new; v$v2 <- v2_new; Emax <- Emax_new
  }
  denom <- exp(v$v0) + exp(v$v1) + exp(v$v2)
  tibble(ID = df$ID, year = df$year,
         p0 = as.numeric(exp(v$v0) / denom),
         p1 = as.numeric(exp(v$v1) / denom),
         p2 = as.numeric(exp(v$v2) / denom))
}

# BASELINE probabilities (observed regime: base merge model)
base_probs <- iterate_ccp(
  df = struct_data %>% select(ID, year, year_f, MSTATE, copa, cah, distance, beds, pi_hat),
  build_v_fun = build_v_base
) %>% rename(p_continue_base = p0, p_close_base = p1, p_merge_base = p2)

# CF1: No CAH designation (set cah = 0 throughout)
cf1_df <- struct_data %>% mutate(cah = 0) %>% select(ID, year, year_f, MSTATE, cah, copa, distance, beds, pi_hat)
cf1_probs <- iterate_ccp(
  df = cf1_df,
  build_v_fun = build_v_base
) %>% rename(p_continue_cf1 = p0, p_close_cf1 = p1, p_merge_cf1 = p2)

# CF2: CAH allowed to merge (use extended merge utility with 'cah')
cf2_df <- struct_data %>% select(ID, year, year_f, MSTATE, cah, copa, distance, beds, pi_hat)
cf2_probs <- iterate_ccp(
  df = cf2_df,
  build_v_fun = build_v_ext
) %>% rename(p_continue_cf2 = p0, p_close_cf2 = p1, p_merge_cf2 = p2)




# Join baseline and counterfactual probabilities ----------------------------------
compare_df <- base_probs %>%
  left_join(cf1_probs, by = c("ID","year")) %>%
  left_join(cf2_probs, by = c("ID","year"))

# Aggregate by year
compare_summary <- compare_df %>%
  group_by(year) %>%
  summarize(
    base_close = mean(p_close_base, na.rm = TRUE),
    cf1_close  = mean(p_close_cf1,  na.rm = TRUE),
    cf2_close  = mean(p_close_cf2,  na.rm = TRUE),
    base_merge = mean(p_merge_base, na.rm = TRUE),
    cf1_merge  = mean(p_merge_cf1,  na.rm = TRUE),
    cf2_merge  = mean(p_merge_cf2,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(
    cum_base_close = cumsum(base_close),
    cum_cf1_close  = cumsum(cf1_close),
    cum_cf2_close  = cumsum(cf2_close),
    cum_base_merge = cumsum(base_merge),
    cum_cf1_merge  = cumsum(cf1_merge),
    cum_cf2_merge  = cumsum(cf2_merge)
  )

# Plots (unchanged styling)
library(ggplot2)
library(tidyr)

# Annual Closure
compare_summary %>%
  pivot_longer(cols = c(base_close, cf1_close, cf2_close),
               names_to = "scenario", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line() + geom_point() +
  labs(title = "Predicted Closure Probabilities Over Time",
       y = "Probability of Closure", x = "Year", color = "Scenario") +
  theme_minimal()

# Annual Merger
compare_summary %>%
  pivot_longer(cols = c(base_merge, cf1_merge, cf2_merge),
               names_to = "scenario", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line() + geom_point() +
  labs(title = "Predicted Merger Probabilities Over Time",
       y = "Probability of Merger", x = "Year", color = "Scenario") +
  theme_minimal()

# Cumulative Closure
p3 <- compare_summary %>%
  pivot_longer(cols = c(cum_base_close, cum_cf1_close, cum_cf2_close),
               names_to = "scenario", values_to = "value") %>%
  mutate(scenario = recode(scenario,
                           cum_base_close = "Baseline",
                           cum_cf1_close  = "No CAH",
                           cum_cf2_close  = "CAH Allows Merger")) %>%
  ggplot(aes(x = year, y = value, linetype = scenario, group = scenario)) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(size = 1.5, shape = 16, color = "black") +
  scale_linetype_manual(values = c("Baseline" = "solid", "No CAH" = "dashed", "CAH Allows Merger" = "dotted")) +
  scale_linewidth_manual(values = c("Baseline" = 0.8, "No CAH" = 0.8, "CAH Allows Merger" = 0.8)) +
  labs(y = "Cumulative Probability of Closure", x = "Year", color = "Scenario") +
  theme_minimal()

ggsave("results/cuml_closures.png", plot = p3, width = 8, height = 5, dpi = 300)

# Cumulative Merger
p4 <- compare_summary %>%
  pivot_longer(cols = c(cum_base_merge, cum_cf1_merge, cum_cf2_merge),
               names_to = "scenario", values_to = "value") %>%
  mutate(scenario = recode(scenario,
                           cum_base_merge = "Baseline",
                           cum_cf1_merge  = "No CAH",
                           cum_cf2_merge  = "CAH Allows Merger")) %>%
  ggplot(aes(x = year, y = value, linetype = scenario, group = scenario)) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(size = 1.5, shape = 16, color = "black") +
  scale_linetype_manual(values = c("Baseline" = "solid", "No CAH" = "dashed", "CAH Allows Merger" = "dotted")) +
  scale_linewidth_manual(values = c("Baseline" = 0.8, "No CAH" = 0.8, "CAH Allows Merger" = 0.8)) +
  labs(y = "Cumulative Probability of Merger", x = "Year", color = "Scenario") +
  theme_minimal()

ggsave("results/cuml_mergers.png", plot = p4, width = 8, height = 5, dpi = 300)
