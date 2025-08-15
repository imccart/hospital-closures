# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   1/29/2025
## Description:   Build Analytic Data


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               fedmatch, scales, zipcodeR)


# Read-in data ------------------------------------------------------------
aha.combine <- read_csv('data/input/aha_data.csv')
cah.supplement <- read_csv('data/input/cah_data.csv')
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv')
col_names <- names(state.zip.xwalk)
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv', col_names=col_names, skip=2) %>%
  group_by(zcta5, stab) %>% mutate(state_zip=row_number()) %>% filter(state_zip==1) %>% ungroup() %>%
  group_by(zcta5) %>% mutate(zip_count=n()) %>% filter(zip_count==1) %>% ungroup() %>%
  select(zcta5, state_xwalk=stab)

hcris.data <- read_tsv('data/input/hcris_data.txt') %>%
  rename(MCRNUM=provider_number)

form990.data <- read_tsv('data/input/form990_ahaid.txt') %>%
  mutate(total_revenue=abs(total_revenue),
         total_expenses=abs(total_expenses),
         total_assets=abs(total_assets),
         total_liabilities=abs(total_liabilities),
         margin=if_else(
              !is.na(total_revenue) & total_revenue>0 & !is.na(total_expenses) & total_expenses>0,
              (total_revenue-total_expenses)/total_revenue,
              NA),
         current_ratio=if_else(
              !is.na(total_assets) & total_assets>0 & !is.na(total_liabilities) & total_liabilities>0,
              total_assets/total_liabilities,
              NA)) %>%
  filter(margin>-1, margin<1)

source('data-code/functions.R')
source('data-code/api-keys.R')


# AHA data cleaning (manual) -----------------------------------------------
aha.combine <- aha.combine %>% 
  mutate(MSTATE=if_else(ID=="6540810", "MS", MSTATE))

# AHA ID to MCRNUM Crosswalk ----------------------------------------------

## Fuzzy match of names in AHA and HCRIS
aha.small <- aha.combine %>%
  select(ID, SYSID, critical_access, name=MNAME, state=MSTATE, city=MLOCCITY, MLOCZIP, year) %>%
  left_join(state.zip.xwalk, by=c("MLOCZIP"="zcta5")) %>%
  mutate(zip=substr(MLOCZIP, 1, 5)) %>% 
  select(-MLOCZIP) %>%
  mutate(state=case_when(
      !is.na(state) ~ state,
      is.na(state) & !is.na(state_xwalk) ~ state_xwalk
    )) %>%
  distinct(ID, name, state, city, zip, year) %>%
  group_by(ID, name, state, city, zip, year) %>%
  mutate(aha_id=cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(name, state, city), str_to_lower)
  

hcris.small <- hcris.data %>%
  mutate(zip=substr(zip, 1, 5)) %>%   
  select(name, city, state, zip, MCRNUM, year) %>%
  distinct() %>%
  group_by(name, city, state, zip, MCRNUM, year) %>%
  mutate(hcris_id=cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(name, state, city), str_to_lower)


fuzzy.merge.hcris <- merge_plus(
  data1=aha.small,
  data2=hcris.small,
  by=c("name","city", "state", "zip", "year"),
  unique_key_1="aha_id",
  unique_key_2="hcris_id",
  match_type="multivar",
  multivar_settings = build_multivar_settings(
    compare_type=c("stringdist","stringdist","indicator","indicator","indicator"),
    wgts=c(0.2, 0.2, 0.2, 0.2, 0.2)
  )
)


fuzzy.match.hcris <- as_tibble(fuzzy.merge.hcris$matches) %>%
  select(ID, name_1, city_1, state_1, zip_1, name_2, city_2, state_2, zip_2, MCRNUM, year_1,
         name_compare, city_compare, state_compare, zip_compare, year_compare, multivar_score) %>%
  filter(!is.na(multivar_score)) %>%
  filter(state_compare==1, year_compare==1) %>%
  group_by(MCRNUM, year_1) %>%
  mutate(max_score=max(multivar_score, na.rm=TRUE),
         max_name_score=max(name_compare, na.rm=TRUE)) %>%
  filter(max_score==multivar_score) %>%
  filter(max_name_score==name_compare) %>%
  ungroup() %>%
  group_by(ID, MCRNUM, year_1) %>%
  slice(1) %>%
  ungroup() %>%
  distinct(ID, MCRNUM, year=year_1) %>%
  arrange(MCRNUM, year, ID)

## Create crosswalk from AHA ID and MCRNUM
aha.crosswalk1 <- aha.combine %>%
  select(ID, MCRNUM, year) %>%
  filter(!is.na(MCRNUM)) %>%
  group_by(ID, MCRNUM, year) %>%
  mutate(count_pair=n()) %>%
  filter(count_pair==1) %>%
  select(ID, year, MCRNUM_xw1=MCRNUM) %>%
  ungroup()

aha.crosswalk2 <- aha.combine %>%
  distinct(ID, year) %>%
  left_join(fuzzy.match.hcris %>% select(ID, MCRNUM, year), by=c("ID","year")) %>%
  select(ID, year, MCRNUM_xw2=MCRNUM) %>%
  ungroup()

aha.crosswalk3 <- aha.combine %>%
  distinct(ID, year) %>%
  left_join(aha.crosswalk1, by=c("ID", "year")) %>%
  left_join(aha.crosswalk2, by=c("ID", "year")) %>%
  mutate(MCRNUM_xw3=if_else(!is.na(MCRNUM_xw1),MCRNUM_xw1, MCRNUM_xw2)) %>%
  filter(!is.na(MCRNUM_xw3)) %>%
  group_by(ID) %>%
  arrange(year) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, MCRNUM_xw3)


aha.crosswalk <- aha.combine %>%
  distinct(ID, year) %>%
  left_join(aha.crosswalk1, by=c("ID", "year")) %>%
  left_join(aha.crosswalk2, by=c("ID", "year")) %>%
  left_join(aha.crosswalk3, by="ID") %>%
  mutate(MCRNUM_xw=case_when(
    !is.na(MCRNUM_xw1) ~ MCRNUM_xw1,
    is.na(MCRNUM_xw1) & !is.na(MCRNUM_xw2) ~ MCRNUM_xw2,
    TRUE ~ MCRNUM_xw3
  )) %>%
  filter(!is.na(MCRNUM_xw)) %>%
  select(ID, MCRNUM_xw, year) %>%
  arrange(ID, MCRNUM_xw, year)

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

fuzzy.merge.cah <- merge_plus(
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

fuzzy.match.cah <- as_tibble(fuzzy.merge.cah$matches) %>%
  select(ID, name_1, city_1, state_1, zip_1, name_2, city_2, state_2, zip_2, 
         name_compare, city_compare, state_compare, zip_compare, multivar_score, eff_date) %>%
  filter(!is.na(multivar_score), state_compare==1, zip_compare==1, name_compare>0.8) %>%
  group_by(ID) %>%
  mutate(max_score=max(multivar_score, na.rm=TRUE),
         max_name_score=max(name_compare, na.rm=TRUE)) %>%
  filter(max_score==multivar_score) %>%
  filter(max_name_score==name_compare) %>%
  ungroup()

# Fuzzy match of EINs from Form 990 ---------------------------------------

form990.small <- form990.data %>%
  filter(is.na(ID_hospital_1)) %>%
  filter(str_detect(str_to_lower(name), "hospital|medical|health"))
  select(ein, name, state, zip, year) %>%
  distinct(ein, name, state, zip) %>%
  group_by(ein, name, state, zip) %>%
  mutate(irs_id=cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(name, state), str_to_lower)

fuzzy.merge.990 <- merge_plus(
  data1=aha.small,
  data2=form990.small,
  by=c("name", "state", "zip"),
  unique_key_1="aha_id",
  unique_key_2="irs_id",
  match_type="multivar",
  multivar_settings = build_multivar_settings(
    compare_type=c("stringdist","indicator","indicator"),
    wgts=c(0.2, 0.4, 0.4)
  )
)

fuzzy.match.990 <- as_tibble(fuzzy.merge.990$matches) %>%
  select(aha_id, ID, ein, irs_id, name_1, state_1, zip_1, name_2, state_2, zip_2, 
         name_compare, state_compare, zip_compare, multivar_score) %>%
  filter(!is.na(multivar_score), state_compare==1, zip_compare==1, name_compare>0.9) %>%
  group_by(ID) %>%
  mutate(max_score=max(multivar_score, na.rm=TRUE),
         max_name_score=max(name_compare, na.rm=TRUE)) %>%
  filter(max_score==multivar_score) %>%
  filter(max_name_score==name_compare) %>%
  ungroup()

# Final data --------------------------------------------------------------

fuzzy.unique.cah <- fuzzy.match.cah %>% 
  group_by(ID) %>%
  summarize(first_date=min(eff_date, na.rm=TRUE)) %>%
  mutate(cah_sup=1) %>%
  ungroup()

form990.id <- form990.data %>%
  filter(!is.na(ID_hospital_1)) %>%
  select(ID_hospital_1, ID_hospital_2, ID_hospital_3, ein, year) %>%
  pivot_longer(cols=c(ID_hospital_1, ID_hospital_2, ID_hospital_3), 
               names_to="rcount", values_to="ID") %>%
  filter(!is.na(ID), !is.na(ein)) %>%
  mutate(ID=as.character(ID)) %>%
  distinct(ID, ein)  # no duplicates within year

fuzzy.unique.990 <- fuzzy.match.990 %>% 
  distinct(ID, ein, multivar_score, name_compare, aha_name=name_1, ein_name=name_2) %>%
  group_by(ID) %>% mutate(rcount=row_number()) %>%
  pivot_wider(values_from=c("ein","multivar_score","name_compare","aha_name","ein_name"), names_from="rcount") %>%
  ungroup()

unique.990 <- bind_rows(fuzzy.unique.990, form990.id %>% select(ein_1=ein, ID))

aha.final <- aha.combine %>% 
  left_join(fuzzy.unique.cah, by='ID') %>%
  left_join(unique.990, by='ID') %>%
  left_join(form990.data %>% select(ein_1=ein, year, margin_1=margin, current_ratio_1=current_ratio), by=c('ein_1','year')) %>%
  left_join(form990.data %>% select(ein_2=ein, year, margin_2=margin, current_ratio_2=current_ratio), by=c('ein_2','year')) %>%
  left_join(form990.data %>% select(ein_3=ein, year, margin_3=margin, current_ratio_3=current_ratio), by=c('ein_3','year')) %>%
  left_join(state.zip.xwalk, by=c("MLOCZIP"="zcta5")) %>%
  mutate(MSTATE=case_when(
      !is.na(MSTATE) ~ MSTATE,
      is.na(MSTATE) & !is.na(state_xwalk) ~ state_xwalk
    )) %>%
  left_join(aha.crosswalk, by=c("ID","year")) %>%
  mutate(MCRNUM=case_when(
      !is.na(MCRNUM) ~ MCRNUM,
      is.na(MCRNUM) & !is.na(MCRNUM_xw) ~ MCRNUM_xw
    )) %>%
  left_join(hcris.data, by=c('MCRNUM','year')) %>%
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

