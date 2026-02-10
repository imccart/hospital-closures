# Discrete-time hazards for closure, merger, and either -------------------
survival.dat <- est.dat %>%
  group_by(ID) %>% 
  mutate(min_bedsize = min(BDTOT, na.rm = TRUE), max_distance=max(distance, na.rm=TRUE)) %>%
  ungroup() %>%
  filter(min_bedsize <= 50) %>%
  mutate(ID = as.numeric(factor(ID)),
         treat_group = case_when(
           !is.na(eff_year) ~ eff_year,
           is.na(eff_year) & state_treat_year > year ~ 0,
           is.na(eff_year) & state_treat_year == 0 ~ 0,
           TRUE ~ NA 
         )) %>%
  filter(!is.na(year), !is.na(closed),
         treat_group %in% c(0, 1999, 2000, 2001, 2002)) %>%
  group_by(ID) %>%
  mutate(panel_start_year = min(year, na.rm = TRUE),
         duration_years = year - panel_start_year + 1) %>% # Year 1 = first year observed
  ungroup() %>%
  mutate(
    treated = ifelse(treat_group > 0, 1, 0), 
    post_treat = ifelse(year >= treat_group & treat_group > 0, 1, 0),
    time_from_treat = ifelse(treat_group > 0, year - treat_group, -1),
    event_time = case_when(
      treated == 0 ~ -1,
      time_from_treat < - max.es ~ - max.es,
      time_from_treat > max.es ~ max.es,
      TRUE ~ time_from_treat
    ),
    event_m5 = ifelse(event_time == -5 & treated==1, 1, 0),
    event_m4 = ifelse(event_time == -4 & treated==1, 1, 0),
    event_m3 = ifelse(event_time == -3 & treated==1, 1, 0),
    event_m2 = ifelse(event_time == -2 & treated==1, 1, 0),
    event_m1 = ifelse(event_time == -1 & treated==1, 1, 0),
    event_0 = ifelse(event_time == 0 & treated==1, 1, 0),
    event_p1 = ifelse(event_time == 1 & treated==1, 1, 0),
    event_p2 = ifelse(event_time == 2 & treated==1, 1, 0),
    event_p3 = ifelse(event_time == 3 & treated==1, 1, 0),
    event_p4 = ifelse(event_time == 4 & treated==1, 1 , 0),
    event_p5 = ifelse(event_time == 5 & treated==1, 1 , 0 )) %>%  
  select(ID, year, closed, merged, closed_merged, duration_years, treated, post_treat, treat_group, state_treat_year,
         BDTOT, distance, own_type, MSTATE, min_bedsize, max_distance, time_from_treat, starts_with("event"))

clog.closed1 <- feglm( closed ~ post_treat + treated  | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5)
)

clog.closed2 <- feglm( closed ~ post_treat + treated + min_bedsize + max_distance | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5)
)


clog.merged1 <- feglm( merged ~ post_treat + treated  | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5)
)

clog.merged2 <- feglm( merged ~ post_treat + treated + min_bedsize + max_distance | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5),
)


clog.change1 <- feglm( closed_merged ~ post_treat + treated | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5)
)

clog.change2 <- feglm( closed_merged ~ post_treat + treated + min_bedsize + max_distance | MSTATE + year + duration_years,
  family = binomial(link = "cloglog"),
  data   = survival.dat %>% filter(time_from_treat>=-5, time_from_treat<=5)
)


# Export results table ------------------------------------------------------------

models <- list(
  clog.closed1,
  clog.closed2,
  clog.merged1,
  clog.merged2,
  clog.change1,
  clog.change2
)

## Helper: extract coef and clustered SE for a given coefficient name
get_coef_se <- function(m, coef_name) {
  V  <- vcov(m, cluster = ~ MSTATE)  # cluster-robust vcov at state level
  b  <- coef(m)[coef_name]
  se <- sqrt(diag(V))[coef_name]
  c(b = b, se = se)
}

## Extract for post_treat and treated across all 6 models
post_mat  <- t(sapply(models, get_coef_se, coef_name = "post_treat"))
treat_mat <- t(sapply(models, get_coef_se, coef_name = "treated"))

# Convert to hazard ratios and SEs on HR scale via delta method
post_hr  <- exp(post_mat[, "b.post_treat"])
post_se_hr <- post_hr * post_mat[, "se.post_treat"]

treat_hr  <- exp(treat_mat[, "b.treated"])
treat_se_hr <- treat_hr * treat_mat[, "se.treated"]

## Formatters
fmt_coef <- function(x) sprintf("%.3f", x)
fmt_se   <- function(x) sprintf("(%.3f)", x)

## Bottom panel entries
hospital_controls <- c("No", "Yes", "No", "Yes", "No", "Yes")
state_fe          <- rep("Yes", 6)
year_fe           <- rep("Yes", 6)
duration_fe       <- rep("Yes", 6)


## Build table data frame
table_df <- tibble(
  term = c(
    "Post-treatment",
    "",
    "Treated",
    "",
    "Hospital controls",
    "State FE",
    "Year FE",
    "Duration FE"
  )
)

for (j in 1:6) {
  col_vals <- c(
    fmt_coef(post_hr[j]),
    fmt_se(post_se_hr[j]),
    fmt_coef(treat_hr[j]),
    fmt_se(treat_se_hr[j]),
    hospital_controls[j],
    state_fe[j],
    year_fe[j],
    duration_fe[j]
  )
  table_df[[paste0("col", j)]] <- col_vals
}

## Column names and header groups
col_names <- c("", paste0("(", 1:6, ")"))
header_spanner <- c(
  " " = 1,
  "Closed" = 2,
  "Merged" = 2,
  "Closed or merged" = 2
)

## LaTeX table
latex_table_hazard <- table_df %>%
  kbl(
    format    = "latex",
    booktabs  = TRUE,
    caption   = "Discrete-time hazard models for closure, merger, and any change",
    col.names = col_names,
    align     = c("l", rep("c", 6)),
    escape    = TRUE
  ) %>%
  add_header_above(header_spanner) %>%
  pack_rows("Hazard ratios", 1, 4) %>%
  pack_rows("Specification", 5, 8) %>%
  kable_styling(
    font_size    = 9,
    latex_options = c("hold_position", "scale_down")
  )

save_kable(latex_table_hazard, file = "results/changes-hazard.tex")
