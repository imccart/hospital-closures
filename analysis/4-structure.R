# Structural Estimation: Hospital Exit and Merger Decisions -----------------------

# Step 1: Prepare State Space -----------------------------------------------------
# We assume est.dat includes all state variables s_it and discrete choices
# d_it ∈ {0, 1, 2} for {continue, close, merge}

struct_data <- est.dat %>%
  mutate(
    d_it = case_when(
      closed == 1 ~ 1,
      merged == 1 ~ 2,
      TRUE ~ 0
    )
  ) 
  
  %>%
  left_join(copa, by = c("MSTATE", "year"))  # placeholder merge

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
  mutate(lead_cah = lead(cah), lead_profit = lead(profit_margin)) %>%
  filter(d_it == 0 & !is.na(lead_profit)) %>%
  ungroup()

# Use random forest or logistic/glmnet if high-dimensional
trans_model <- lm(lead_profit ~ cah + BDTOT + distance + beds + region + year, data = trans_data)

# Predict expected future values
trans_data <- trans_data %>%
  mutate(
    E_profit_next = predict(trans_model)
  )

# Step 4: Construct Value Function Components -------------------------------------

# Flow payoff from continuation
trans_data <- trans_data %>%
  mutate(
    pi = profit_margin,  # or a linear index: pi = X * beta1
    V_continue = pi + 0.95 * E_profit_next  # β = 0.95
  )

# Terminal values (closure and merger), function of observed covariates
trans_data <- trans_data %>%
  mutate(
    V_close = alpha0 + alpha1 * rural + alpha2 * year,   # placeholder
    V_merge = rho0 + rho1 * copa + rho2 * urban          # placeholder
  )

# Step 5: Estimate Structural Parameters by Inversion ------------------------------

# Use logit inversion:
# v_d = log(P_d) - log(P_0) = v_d - v_0

log_odds <- trans_data %>%
  mutate(
    v_close_diff = log(p_close / p_continue),
    v_merge_diff = log(p_merge / p_continue),
    v_close_rhs = V_close - V_continue,
    v_merge_rhs = V_merge - V_continue
  )

# Estimate structural parameters by regressing inverted values
close_model <- lm(v_close_diff ~ v_close_rhs, data = log_odds)
merge_model <- lm(v_merge_diff ~ v_merge_rhs, data = log_odds)

# Output parameter estimates
summary(close_model)
summary(merge_model)
