# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   2/08/2024
## Description:   Build Analytic Data


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR)

# Read-in data ------------------------------------------------------------
aha.combine <- read_csv('data/input/aha_data')
cah.supplement <- read_csv('data/input/cah_data')
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv')
col_names <- names(state.zip.xwalk)
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv', col_names=col_names, skip=2) %>%
  group_by(zcta5, stab) %>% mutate(state_zip=row_number()) %>% filter(state_zip==1) %>% ungroup() %>%
  group_by(zcta5) %>% mutate(zip_count=n()) %>% filter(zip_count==1) %>% ungroup() %>%
  select(zcta5, state_xwalk=stab)
source('data-code/functions.R')
source('data-code/api-keys.R')



# Fuzzy match of CAH supplement -------------------------------------------

aha.small <- aha.combine %>%
  select(ID, SYSID, critical_access, name=MNAME, state=MSTATE, city=MLOCCITY, MLOCZIP, year) %>%
  left_join(state.zip.xwalk, by=c("MLOCZIP"="zcta5")) %>%
  mutate(zip=substr(MLOCZIP, 1, 5)) %>% 
  select(-MLOCZIP) %>%
  mutate(state=case_when(
      !is.na(state) ~ state,
      is.na(state) & !is.na(state_xwalk) ~ state_xwalk
    )) %>%
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

fuzzy.match.full <- as_tibble(fuzzy.merge$matches) %>%
  select(ID, name_1, city_1, state_1, zip_1, name_2, city_2, state_2, zip_2, 
         name_compare, city_compare, state_compare, zip_compare, multivar_score, eff_date) 

fuzzy.match.full %>% distinct(name_2, eff_date) %>% nrow()

fuzzy.match.score <- fuzzy.match.full %>%
  filter(!is.na(multivar_score)) 

fuzzy.match.score %>% distinct(name_2, eff_date) %>% nrow()

fuzzy.match.highscore <- fuzzy.match.score %>%
  filter(city_compare>0.9, zip_compare==1, name_compare>0.69) 

fuzzy.match.highscore %>% distinct(name_2, eff_date) %>% nrow()

fuzzy.match <- fuzzy.match.highscore %>%
  group_by(ID) %>%
  mutate(max_score=max(multivar_score, na.rm=TRUE),
         max_name_score=max(name_compare, na.rm=TRUE)) %>%
  filter(max_score==multivar_score) %>%
  filter(max_name_score==name_compare) %>%
  ungroup()

fuzzy.match %>% distinct(name_2, eff_date) %>% nrow()


# Final data --------------------------------------------------------------

fuzzy.unique <- fuzzy.match %>% 
  group_by(ID) %>%
  summarize(first_date=min(eff_date, na.rm=TRUE)) %>%
  mutate(cah_sup=1) %>%
  ungroup()

aha.final <- aha.combine %>% 
  left_join(fuzzy.unique,
            by='ID') %>%
  left_join(state.zip.xwalk, by=c("MLOCZIP"="zcta5")) %>%
  mutate(MSTATE=case_when(
      !is.na(MSTATE) ~ MSTATE,
      is.na(MSTATE) & !is.na(state_xwalk) ~ state_xwalk
    )) %>%
  mutate(eff_year=year(first_date),
         cah = case_when(
           is.na(critical_access) & cah_sup==1 & year>=eff_year ~ 1,
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
    write_csv('data/output/aha_final.csv')         

#  filter(COMMTY=="Y",
#         SERV==10) %>%


## hospital closures and mergers
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
  mutate_at(vars(LAT, LONG), as.numeric) %>%
  mutate(LAT=ifelse(LAT==0, NA, LAT),
         LONG=ifelse(LONG==0, NA, LONG))


unique_years <- unique(aha.geo$year)
final.neighbors <- tibble()
for (yr in unique_years) {

  aha.pairs <- aha.geo %>%
    filter(year == yr) %>% # Filter data for the current year
    full_join(aha.geo %>% select(ID2=ID, LAT2=LAT, LONG2=LONG, zip2=zip, year), by = "year") %>%
    filter(ID != ID2)

  zip.pairs <- aha.pairs %>%
    filter(!is.na(zip), !is.na(zip2)) %>%
    distinct(zip, zip2)

  zip.dist <- zip_distance(zip.pairs$zip, zip.pairs$zip2) %>%
    rename(zip=zipcode_a, zip2=zipcode_b, dist_zip=distance) %>%
    filter(!is.na(dist_zip))

  aha.zipdist <- aha.pairs %>%
    left_join(zip.dist, by = c("zip", "zip2")) %>%
    select(ID, ID2, dist_zip, year) %>%
    filter(!is.na(dist_zip)) %>%
    group_by(ID, year) %>%
    summarize(min_dist_zip=min(dist_zip, na.rm=TRUE))

  aha.latlongdist <- aha.pairs %>%
    filter(!is.na(LAT), !is.na(LAT2)) %>%
    rowwise() %>%
    mutate(dist_latlong=haversine(c(LONG2, LAT2), c(LONG, LAT))) %>% ungroup() %>%
    group_by(ID, year) %>%
    summarize(min_dist_latlong=min(dist_latlong, na.rm=TRUE))

  aha.neighbors <- aha.geo %>% filter(year==yr) %>%
    left_join(aha.zipdist, by=c("ID", "year")) %>%
    left_join(aha.latlongdist, by=c("ID", "year")) %>%
    select(ID, min_dist_zip, min_dist_latlong, year) %>%
    mutate(distance=
      case_when(
        !is.na(min_dist_latlong) & !is.na(min_dist_zip) ~ pmin(min_dist_latlong, min_dist_zip),
        !is.na(min_dist_latlong) & is.na(min_dist_zip) ~ min_dist_latlong,
        is.na(min_dist_latlong) & !is.na(min_dist_zip) ~ min_dist_zip,
        TRUE ~ NA_real_))
  
  # Add results for the current year to the list
  final.neighbors <- bind_rows(final.neighbors, aha.neighbors)
}

nearest.neighbor <- final.neighbors %>%
    write_csv('data/output/aha_neighbors.csv')

## check means and count of hospitals with distances by year
check <- nearest.neighbor %>%
  group_by(year) %>%
  summarize(mean_dist=mean(distance, na.rm=TRUE),
            count_dist=sum(!is.na(distance)),
            count_all=n())

