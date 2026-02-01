# Fuzzy Matching Wrapper
#
# Loads packages and data, then runs all fuzzy matching scripts:
#   - fuzzy990.R (Form 990 to AHA)
#   - fuzzyhcris.R (HCRIS to AHA)
#   - fuzzycah.R (CAH supplement to AHA)
#
# Run this script before build-data.R

# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, fedmatch)

source('data-code/functions.R')

# Read Source Data --------------------------------------------------------

aha.combine <- read_csv('data/input/aha_data.csv', show_col_types = FALSE)
hcris.data <- read_tsv('data/input/hcris_data.txt', show_col_types = FALSE) %>%
  rename(MCRNUM = provider_number)
cah.supplement <- read_csv('data/input/cah_data.csv', show_col_types = FALSE)
form990.data <- read_tsv('data/input/form990_ahaid.txt', show_col_types = FALSE)

col_names <- names(read_csv('data/input/zcta-to-county.csv', n_max = 0, show_col_types = FALSE))
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv', col_names = col_names, skip = 2, show_col_types = FALSE) %>%
  group_by(zcta5, stab) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(zcta5) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(zcta5, state_xwalk = stab)

# Prepare common AHA data -------------------------------------------------

# Valid US state/territory codes
valid_states <- c("al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl",
                  "ga", "hi", "id", "il", "in", "ia", "ks", "ky", "la", "me",
                  "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh",
                  "nj", "nm", "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri",
                  "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi",
                  "wy", "as", "gu", "mp", "pr", "vi")

aha.small <- aha.combine %>%
  select(ID, SYSID, critical_access, name = MNAME, state = MSTATE, city = MLOCCITY, MLOCZIP, year) %>%
  left_join(state.zip.xwalk, by = c("MLOCZIP" = "zcta5")) %>%
  mutate(
    zip = str_pad(substr(as.character(MLOCZIP), 1, 5), width = 5, side = "left", pad = "0"),
    # Clean invalid zips (should be 5 digits, not all zeros)
    zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_),
    # Clean invalid state codes before coalesce
    state = str_to_lower(state),
    state = if_else(state %in% valid_states, state, NA_character_),
    state = coalesce(state, str_to_lower(state_xwalk))
  ) %>%
  # Fill missing state/zip from other years of the same hospital
  group_by(ID) %>%
  mutate(
    state = coalesce(state, modal_value(state)),
    zip = coalesce(zip, modal_value(zip))
  ) %>%
  ungroup() %>%
  distinct(ID, name, state, city, zip, year) %>%
  group_by(ID, name, state, city, zip, year) %>%
  mutate(aha_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(
    name_raw = str_to_lower(name),
    name_clean = preprocess_hospital_name(name),
    city = str_to_lower(city)
  ) %>%
  filter(!is.na(name_clean), name_clean != "")

# Collapse AHA to unique hospitals (used by fuzzyhcris.R, fuzzycah.R)
# Convert ID to character for consistent joins across all matching scripts
# Uses modal_value() from functions.R to handle groups with all-NA values
aha.collapsed <- aha.small %>%
  mutate(ID = as.character(ID)) %>%
  group_by(ID) %>%
  summarize(
    name_clean = modal_value(name_clean),
    state = first(state),
    city = modal_value(city),
    zip = modal_value(zip),
    .groups = "drop"
  ) %>%
  mutate(aha_id = row_number())

# All historical name/state/zip variants per hospital (used by fuzzy990.R)
# Filters to complete records only - variants with NA state/zip can't be matched anyway
aha.variants <- aha.small %>%
  mutate(ID = as.character(ID)) %>%
  distinct(ID, name_clean, state, zip) %>%
  filter(!is.na(state), !is.na(zip)) %>%
  mutate(variant_id = row_number())

# Prepare Form 990 data (used by fuzzy990.R) ------------------------------

# Keyword filter for likely hospitals
hospital_keywords <- "hospital|medical|health"

# Exclude non-hospital entities that file Form 990s with similar names
exclude_patterns <- c(
  "foundation", "fdn", "fndtn", "fdtn",
  "trust", " tr$", " tr ",
  "auxiliary", "guild", "volunteer",
  "self.ins", "selfins",
  "employee benefit", "emp bene", "bene tr", "benefit plan",
  "insurance plan", "insurance guaranty", "insurance guarnaty",
  "development fund", "endowment fund", "dues fund", "health fund",
  "retirement"
)

form990.small <- form990.data %>%
  filter(is.na(ID_hospital_1)) %>%  # Only unmatched records
  filter(str_detect(str_to_lower(name), hospital_keywords)) %>%
  filter(!str_detect(str_to_lower(name), paste(exclude_patterns, collapse = "|"))) %>%
  mutate(
    zip = str_pad(substr(as.character(zip), 1, 5), width = 5, side = "left", pad = "0"),
    zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_),
    state = str_to_lower(state),
    state = if_else(state %in% valid_states, state, NA_character_)
  ) %>%
  filter(!is.na(state), !is.na(zip)) %>%  # Filter to complete records
  distinct(ein, name, state, zip) %>%
  mutate(
    form990_row_id = row_number(),
    name_raw = str_to_lower(name),
    name_clean = preprocess_hospital_name(name)
  ) %>%
  filter(!is.na(name_clean), name_clean != "")

# Run fuzzy matching scripts ----------------------------------------------

source('data-code/fuzzy990.R')
source('data-code/fuzzyhcris.R')
source('data-code/fuzzycah.R')
