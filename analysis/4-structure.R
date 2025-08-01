# Structural Estimation: Hospital Exit and Merger Decisions -----------------------

# does CAH designation improve hospital survival through continuation, or does it block efficiency-enhancing mergers, ultimately increasing closures?

# Step 1: Prepare State Space -----------------------------------------------------
# We assume est.dat includes all state variables s_it and discrete choices
# d_it ∈ {0, 1, 2} for {continue, close, merge}

vars_used <- c("d_it", "cah", "distance", "margin_1", "beds", "year")

struct_data <- est.dat %>%
  mutate(
    d_it = case_when(
      closed == 1 ~ 1,
      merged == 1 ~ 2,
      TRUE ~ 0
    )
  ) %>%
  filter(complete.cases(across(all_of(vars_used)))) %>%
  left_join(copa.dat, by = c("MSTATE")) %>%
  mutate(copa = case_when(
    !is.na(copa_start) & year >= copa_start & year <= copa_end ~ 1,
    TRUE ~ 0
  ))

# Step 2: Estimate Choice Probabilities -------------------------------------------

# Fit multinomial logit model for P(d_it | s_it)
choice_model <- multinom(
  d_it ~ cah + BDTOT + distance + margin_1 + beds + year,
  data = struct_data
)

# Predicted choice probabilities
struct_data <- struct_data %>%
  mutate(
    p_continue = predict(choice_model, type = "probs")[,1],
    p_close    = predict(choice_model, type = "probs")[,2],
    p_merge    = predict(choice_model, type = "probs")[,3]
  )

# Step 3: Estimate Transition Probabilities ---------------------------------------

# Estimate transitions for next-period state variables (for d_it = 0 only)
trans_data <- struct_data %>%
  group_by(ID) %>%
  arrange(year) %>%
  mutate(lead_cah = lead(cah), lead_profit = lead(margin_1)) %>%
  filter(d_it == 0 & !is.na(lead_profit)) %>%
  ungroup()

# Use random forest or logistic/glmnet if high-dimensional
trans_model <- lm(lead_profit ~ cah + distance + beds + year, data = trans_data)

# Predict expected future values
trans_data <- trans_data %>%
  mutate(
    E_profit_next = predict(trans_model)
  )

# Step 4: Construct Value Function Components -------------------------------------

# Flow payoff from continuation
trans_data <- trans_data %>%
  mutate(
    pi = margin_1,  # or a linear index: pi = X * beta1
    V_continue = pi + 0.95 * E_profit_next  # β = 0.95
  )

# Step 5: Estimate Structural Parameters by Inversion ------------------------------

# Terminal values (closure and merger), function of observed covariates
# Use logit inversion:
# v_d = log(P_d) - log(P_0) = v_d - v_0

log_odds <- trans_data %>%
  mutate(
    v_close_diff = log(p_close / p_continue),
    v_merge_diff = log(p_merge / p_continue)
  )

# Estimate structural parameters by regressing inverted values
close_model <- lm(v_close_diff ~ distance + beds + year, data = log_odds)
merge_model <- lm(v_merge_diff ~ distance + beds + year + copa, data = log_odds)

# Output parameter estimates
summary(close_model)
summary(merge_model)


# Counterfactual 1: No CAH Designation ---------------------------------------------
# Create counterfactual dataset with CAH set to 0
cf1_data <- struct_data %>%
  mutate(cah = 0)

# Recompute expected future profit using transition model
cf1_data <- cf1_data %>%
  mutate(E_profit_next = predict(trans_model, newdata = cf1_data))

# Recompute continuation value
cf1_data <- cf1_data %>%
  mutate(
    pi = margin_1,
    V_continue = pi + 0.95 * E_profit_next
  )

# Predict choice probabilities under no CAH
cf1_data <- cf1_data %>%
  mutate(
    p_hat = predict(choice_model, newdata = cf1_data, type = "probs"),
    p_continue_cf1 = p_hat[,1],
    p_close_cf1    = p_hat[,2],
    p_merge_cf1    = p_hat[,3]
  ) %>%
  select(ID, year, p_continue_cf1, p_close_cf1, p_merge_cf1)


# Counterfactual 2: CAH allowed to merge --------------------------------------
# Step 1: Estimate extended merger model that includes CAH effect
cf2_data <- log_odds
merge_model_ext <- lm(v_merge_diff ~ distance + beds + year + copa + cah, data = cf2_data)

# Step 2: Extract gamma_cah coefficient
gamma_cah <- coef(merge_model_ext)["cah"]

# Step 3: Create counterfactual merger value differences with CAH included
cf2_data <- cf2_data %>%
  mutate(
    v_merge_diff_cf2 = predict(merge_model_ext, newdata = trans_data)
  )

# Step 4: Predict closure value difference as usual
cf2_data <- cf2_data %>%
  mutate(
    v_close_diff_pred = predict(close_model, newdata = cf2_data)
  )

# Step 5: Convert log-odds to choice probabilities
cf2_data <- cf2_data %>%
  mutate(
    exp_continue = 1,
    exp_close = exp(v_close_diff_pred),
    exp_merge = exp(v_merge_diff_cf2),
    denom = exp_continue + exp_close + exp_merge,
    p_continue_cf2 = exp_continue / denom,
    p_close_cf2    = exp_close / denom,
    p_merge_cf2    = exp_merge / denom
  ) %>%
  select(ID, year, p_continue_cf2, p_close_cf2, p_merge_cf2)


# Join back to baseline for comparison
compare_df <- struct_data %>%
  select(ID, year, p_close_base = p_close, p_merge_base = p_merge) %>%
  left_join(cf1_data, by = c("ID", "year")) %>%
  left_join(cf2_data, by = c("ID", "year"))

# Aggregate by year to visualize effects
compare_summary <- compare_df %>%
  group_by(year) %>%
  summarize(
    base_close = mean(p_close_base, na.rm=TRUE),
    cf1_close  = mean(p_close_cf1, na.rm=TRUE),
    cf2_close  = mean(p_close_cf2, na.rm=TRUE),
    base_merge = mean(p_merge_base, na.rm=TRUE),
    cf1_merge  = mean(p_merge_cf1, na.rm=TRUE),
    cf2_merge  = mean(p_merge_cf2, na.rm=TRUE)
  )

# Calculate cumulative probabilities
compare_summary <- compare_summary %>%
  mutate(
    cum_base_close = cumsum(base_close),
    cum_cf1_close  = cumsum(cf1_close),
    cum_cf2_close  = cumsum(cf2_close),
    cum_base_merge = cumsum(base_merge),
    cum_cf1_merge  = cumsum(cf1_merge),
    cum_cf2_merge  = cumsum(cf2_merge)
  )

# Annual Closure Plot
compare_summary %>%
  pivot_longer(cols = starts_with("base_close"):starts_with("cf2_close"),
               names_to = "scenario", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line() +
  geom_point() +
  labs(title = "Predicted Closure Probabilities Over Time",
       y = "Probability of Closure", x = "Year", color = "Scenario") +
  theme_minimal()

# Annual Merger Plot
compare_summary %>%
  pivot_longer(cols = starts_with("base_merge"):starts_with("cf2_merge"),
               names_to = "scenario", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line() +
  geom_point() +
  labs(title = "Predicted Merger Probabilities Over Time",
       y = "Probability of Merger", x = "Year", color = "Scenario") +
  theme_minimal()

# Cumulative Closure Plot
p3 <- compare_summary %>%
  pivot_longer(cols = starts_with("cum_base_close"):starts_with("cum_cf2_close"),
               names_to = "scenario", values_to = "value") %>%
  mutate(scenario = case_when(
    scenario == "cum_base_close" ~ "Baseline",
    scenario == "cum_cf1_close"  ~ "No CAH",
    scenario == "cum_cf2_close"  ~ "CAH Allows Merger"
  )) %>%               
  ggplot(aes(x = year, y = value, linetype = scenario, group = scenario)) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(size = 1.5, shape = 16, color = "black") +
  scale_linetype_manual(values = c(
    "Baseline" = "solid",
    "No CAH" = "dashed",
    "CAH Allows Merger" = "dotted"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline" = 0.8,
    "No CAH" = 0.8,
    "CAH Allows Merger" = 0.8
  )) +
  labs(y = "Cumulative Probability of Closure", x = "Year", color = "Scenario") +
  theme_minimal()

ggsave("results/cuml_closures.png", plot = p3, width = 8, height = 5, dpi = 300)

# Cumulative Merger Plot
p4 <- compare_summary %>%
  pivot_longer(cols = starts_with("cum_base_merge"):starts_with("cum_cf2_merge"),
               names_to = "scenario", values_to = "value") %>%
  mutate(scenario = case_when(
    scenario == "cum_base_merge" ~ "Baseline",
    scenario == "cum_cf1_merge"  ~ "No CAH",
    scenario == "cum_cf2_merge"  ~ "CAH Allows Merger"
  )) %>%
  ggplot(aes(x = year, y = value, linetype = scenario, group = scenario)) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(size = 1.5, shape = 16, color = "black") +
  scale_linetype_manual(values = c(
    "Baseline" = "solid",
    "No CAH" = "dashed",
    "CAH Allows Merger" = "dotted"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline" = 0.8,
    "No CAH" = 0.8,
    "CAH Allows Merger" = 0.8
  )) +
  labs(y = "Cumulative Probability of Merger", x = "Year", color = "Scenario") +
  theme_minimal()

ggsave("results/cuml_mergers.png", plot = p4, width = 8, height = 5, dpi = 300)
