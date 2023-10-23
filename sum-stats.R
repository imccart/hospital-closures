# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   7/17/2023
## Description:   Initial Summary Statistics


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR)

# Read-in data ------------------------------------------------------------
aha.combine <- read_csv('data/input/aha_data')
cah.supplement <- read_csv('data/input/cah_data')



# Fuzzy match of CAH supplement -------------------------------------------

aha.small <- aha.combine %>%
  select(ID, SYSID, critical_access, name=MNAME, state=MSTATE, city=MLOCCITY, MLOCZIP, year) %>%
  mutate(zip=substr(MLOCZIP, 1, 5)) %>% 
  select(-MLOCZIP) %>%
  distinct(ID, name, state, city, zip) %>%
  group_by(ID, name, state, city, zip) %>%
  mutate(aha_id=cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(name, state, city), str_to_lower)
  

cah.small <- cah.supplement %>%
  select(name, city, state, zip, eff_date) %>%
  group_by(name, city, state, zip) %>%
  mutate(cah_id=cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(name, state, city), str_to_lower)

fuzzy.merge <- merge_plus(
  data1=aha.small,
  data2=cah.small,
  by=c("name","city", "state", "zip"),
  unique_key_1="aha_id",
  unique_key_2="cah_id",
  match_type="multivar",
  multivar_settings = build_multivar_settings(
    compare_type=c("stringdist","stringdist","indicator","indicator"),
    wgts=c(0.2, 0.2, 0.3, 0.3)
  )
)

fuzzy.match <- as_tibble(fuzzy.merge$matches) %>%
  select(ID, name_1, city_1, state_1, zip_1, name_2, city_2, state_2, zip_2, 
         name_compare, city_compare, state_compare, zip_compare, multivar_score, eff_date) %>%
  filter(!is.na(multivar_score)) %>%
  filter(city_compare>0.9, zip_compare==1, name_compare>0.69) %>%
  group_by(ID) %>%
  mutate(max_score=max(multivar_score),
         max_name_score=max(name_compare)) %>%
  filter(max_score==multivar_score) %>%
  filter(max_name_score==name_compare) %>%
  ungroup()
  



# Final data --------------------------------------------------------------

fuzzy.unique <- fuzzy.match %>% 
  select(ID, eff_date, aha_name=name_1, name_compare) %>% 
  mutate(cah_sup=1) %>%
  group_by(ID, aha_name) %>%
  mutate(first_date=min(eff_date)) %>%
  filter(n()==1)

aha.final <- aha.combine %>% 
  mutate(aha_name=str_to_lower(MNAME)) %>%
  left_join(fuzzy.unique,
            by=c('ID', 'aha_name')) %>%
  mutate(eff_year=year(first_date),
         critical_access = case_when(
           critical_access==0 & cah_sup==1 & year>=eff_year ~ 1,
           critical_access==1 ~ 1,
           TRUE ~ 0),
         hosp_type = case_when(
           SERV==10 ~ 'General',
           SERV==13 ~ 'Surgical',
           SERV==33 ~ 'Specialty - Respiratory',
           SERV==41 ~ 'Specialty - Cancer',
           SERV==42 ~ 'Specialty - Heart',
           SERV==44 ~ 'Specialty - OBGYN',
           SERV==45 ~ 'Specialty - ENT',
           SERV==47 ~ 'Specialty - Orthopedic',
           SERV==49 ~ 'Specialty - Other',
           SERV %in% c(50, 53, 55, 56, 57, 58, 59) ~ 'Specialty - Childrens',
           SERV==80 ~ 'Long-term Care',
           TRUE ~ 'Other'
         )) %>%
  filter(COMMTY=="Y",
         SERV==10)
  


# Quick summary -----------------------------------------------------------

## bed sizes of hospitals that closed
bed.closed <- aha.final %>% filter(change_type=="Closure") %>% 
  mutate(bed_bin=cut(BDTOT, breaks=c(0,5,10,15,20,25,30,40,50,75,100,150,200,250,300,400,500,999))) %>% 
  count(bed_bin) %>% arrange(bed_bin)

## CAH designation among hospitals that closed
cah.closed <- aha.final %>% mutate(all_cah=sum(critical_access)) %>% 
  group_by(change_type, critical_access) %>% summarize(n=n(), all_cah=first(all_cah), all_hosp=nrow(aha.combine))

## Count of CAH and Non-CAH over time
fig.hosp.type <- aha.final %>% group_by(year, critical_access) %>%
  summarize(hosp_count=n())  %>%
  ggplot(aes(x=year, y=hosp_count, group=critical_access)) + 
  geom_line() + geom_point() + theme_bw() +
  geom_text(data = aha.final %>% filter(year==2016) %>% group_by(critical_access, year) %>% summarize(hosp_count=n()), 
            aes(label = c("Non-CAH", "CAH"),
                x = year+1,
                y = hosp_count-100)) +
  scale_y_continuous(labels = comma,
                     breaks=seq(0, 8000, 500)) +
  scale_x_continuous(breaks=seq(1980, 2019, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(linetype="none") +
  labs(
    x="Year",
    y="Count of Hospital Types",
    caption="Limited to general, short-term, acute care community hospitals"
  )

## Share of CAH and Non-CAH closures over time
cah.change <- aha.final %>%
  filter(critical_access==1) %>%
  group_by(year, change_type) %>%
  summarize(cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_cah=sum(cah_change),
         prop_cah=cah_change/all_cah) %>%
  select(year, prop_cah, change_type) %>%
  ungroup()

non.cah.change <- aha.final %>%
  filter(critical_access==0) %>%
  group_by(year, change_type) %>%
  summarize(non_cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_non_cah=sum(non_cah_change),
         prop_non_cah=non_cah_change/all_non_cah) %>%
  select(year, prop_non_cah, change_type) %>%
  ungroup()

all.change <- aha.final %>%
  group_by(year, change_type) %>%
  summarize(any_change=n()) %>%
  group_by(year) %>%
  mutate(all_hosp=sum(any_change),
         prop_change=any_change/all_hosp) %>%
  select(year, change_type, prop_change) %>%
  ungroup()


fig.close <- all.change %>%
  left_join(non.cah.change,
            by=c("year","change_type")) %>%
  left_join(cah.change, 
            by=c("year", "change_type")) %>%
  filter(change_type %in% c("Closure", "Merger")) %>%
  select(year, change_type, prop_non_cah, prop_cah) %>%
  replace_na(list(prop_non_cah=0, prop_cah=0)) %>%
  pivot_longer(cols = starts_with("prop"), names_to = "type", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = type, linetype = change_type)) +
  geom_line() +
  ylim(0,0.05) + 
  scale_x_continuous(breaks=seq(1980, 2019, 1)) +
  labs(x = "Year", y = "Proportion", color = "Hospital Type", linetype = "Change Type") +
  theme_minimal() +
  scale_color_discrete(labels=c("CAH","Non-CAH")) + 
  guides(color = guide_legend(order = 2), linetype = guide_legend(order = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position=c(0.8,0.65))

merge.close <- aha.final %>% 
  group_by(year, change_type, critical_access) %>%
  summarize(hosp_count=n()) %>% ungroup() %>%
  group_by(year, critical_access) %>%
  mutate(hosp_type_count=sum(hosp_count)) %>% ungroup() %>%
  group_by(year, change_type) %>%
  mutate(change_type_count=sum(hosp_count)) %>% ungroup() %>%
  write_csv('data/output/merge_close.csv')



## identify each hospital's nearest neighbor and the distance to that neighbor using the haversine formula
aha.geo <- aha.final %>%
  mutate(zip=substr(MLOCZIP, 1, 5)) %>%
  select(ID, LAT, LONG, year, zip) %>%
  distinct(ID, LAT, LONG, year, zip) %>%
  mutate_at(vars(LAT, LONG), as.numeric)

list_results <- list()
for (yr in unique_years) {
  res <- aha.geo %>%
    filter(year == 1995) %>% # Filter data for the current year
    full_join(aha.geo %>% select(ID2=ID, LAT2=LAT, LONG2=LONG, zip2=zip, year), by = "year") %>%
    filter(ID != ID2) %>%
    mutate(dist = case_when(
      is.na(LAT) | is.na(LONG) | is.na(LAT2) | is.na(LONG2) ~ zip_distance(~zip, ~zip2),
      TRUE ~ haversine(cbind(LONG2, LAT2), cbind(LONG, LAT)))) %>%
    filter(!is.na(dist)) %>%
    group_by(ID, year) %>%
    mutate(min_dist=min(dist, na.rm=TRUE)) %>%
    filter(dist == min_dist) 
    
    %>%
    ungroup() %>%
    select(ID, ID2, dist, year)
  
  # Add results for the current year to the list
  list_results[[as.character(yr)]] <- res
}

# Combine results of all years into a single data frame
final_result <- bind_rows(list_results)

# Write the final result to CSV
write_csv(final_result, 'data/output/aha_geo_neighbors.csv')
