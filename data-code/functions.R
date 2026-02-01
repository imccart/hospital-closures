# Modal value (most frequent) with NA handling --------------------------------
# Returns the most frequent non-NA value, or NA if all values are NA

modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}


haversine <- function(coord1, coord2) {
  R <- 6371.01  # Earth's radius in kilometers

  lat1 <- as.numeric(coord1[2]) * pi / 180
  long1 <- as.numeric(coord1[1]) * pi / 180
  lat2 <- as.numeric(coord2[2]) * pi / 180
  long2 <- as.numeric(coord2[1]) * pi / 180

  dlat <- lat2 - lat1
  dlong <- long2 - long1

  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlong/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))

  d <- R * c
  return(d)
}


# Vectorized haversine for efficient distance calculations ----------------------
# Computes great-circle distance between coordinate vectors
# Returns distance in miles (to match zipcodeR::zip_distance default)

haversine_vec <- function(lon1, lat1, lon2, lat2, units = "miles") {
  R <- if (units == "miles") 3958.8 else 6371.01

  lat1_rad <- lat1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180

  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  R * c
}


# Word-level Jaccard similarity -------------------------------------------------
# Compares two strings by tokenizing into words and computing Jaccard index

jaccard_words <- function(a, b) {
  # Handle NA or empty strings
  if (is.na(a) || is.na(b) || a == "" || b == "") return(0)

  # Tokenize into words
  words_a <- unlist(str_split(str_squish(a), "\\s+"))
  words_b <- unlist(str_split(str_squish(b), "\\s+"))

  # Remove empty tokens
  words_a <- words_a[words_a != ""]
  words_b <- words_b[words_b != ""]

  # Handle edge case of empty token sets
  if (length(words_a) == 0 || length(words_b) == 0) return(0)

  # Compute Jaccard: |intersection| / |union|
  intersection <- length(intersect(words_a, words_b))
  union <- length(union(words_a, words_b))

  return(intersection / union)
}

# Vectorized version for applying to data frame columns
jaccard_words_vec <- function(a_vec, b_vec) {
  mapply(jaccard_words, a_vec, b_vec, USE.NAMES = FALSE)
}


# Word-level Jaccard matching with blocking -------------------------------------
# Creates candidate pairs based on blocking variables, computes word Jaccard

match_jaccard_words <- function(data1, data2,
                                 name_col1, name_col2,
                                 id_col1, id_col2,
                                 block_vars,
                                 threshold = 0.75) {

  # Rename columns for consistency
  d1 <- data1 %>%
    select(id1 = all_of(id_col1),
           name1 = all_of(name_col1),
           all_of(block_vars))

  d2 <- data2 %>%
    select(id2 = all_of(id_col2),
           name2 = all_of(name_col2),
           all_of(block_vars))

  # Filter out NA block vars to prevent NA-to-NA matching
  d1 <- d1 %>% filter(if_all(all_of(block_vars), ~ !is.na(.)))
  d2 <- d2 %>% filter(if_all(all_of(block_vars), ~ !is.na(.)))

  # Join on blocking variables to create candidate pairs
  candidates <- inner_join(d1, d2, by = block_vars, relationship = "many-to-many")

  # Compute word-level Jaccard for each pair
  candidates <- candidates %>%
    mutate(score = jaccard_words_vec(name1, name2)) %>%
    filter(score >= threshold)

  return(candidates)
}


# Hospital name preprocessing ----------------------------------------------------
# Standardizes hospital names for better fuzzy matching

preprocess_hospital_name <- function(x) {
  x %>%
    str_to_lower() %>%

    # Remove health system prefixes
    str_remove("^hshs\\s+") %>%
    str_remove("^ssm\\s+") %>%
    str_remove("^hca\\s+") %>%
    str_remove("^banner\\s+") %>%
    str_remove("^ascension\\s+") %>%
    str_remove("^adventist\\s+") %>%
    str_remove("^trinity\\s+") %>%
    str_remove("^dignity\\s+") %>%
    str_remove("^tenet\\s+") %>%
    str_remove("^lifepoint\\s+") %>%
    str_remove("^community\\s+health\\s+systems\\s+") %>%
    str_remove("^universal\\s+health\\s+services\\s+") %>%

    # Standardize common abbreviations
    str_replace_all("\\bst\\.?\\s", "saint ") %>%
    str_replace_all("\\bmt\\.?\\s", "mount ") %>%
    str_replace_all("\\bmem\\.?\\b", "memorial") %>%
    str_replace_all("\\bgen\\.?\\b", "general") %>%
    str_replace_all("\\bcomm\\.?\\b", "community") %>%
    str_replace_all("\\breg\\.?\\b", "regional") %>%
    str_replace_all("\\bpresby\\.?\\b", "presbyterian") %>%
    str_replace_all("\\bmeth\\.?\\b", "methodist") %>%
    str_replace_all("\\bhosp\\.?\\b", "hospital") %>%
    str_replace_all("\\bmed\\.?\\b", "medical") %>%
    str_replace_all("\\bctr\\.?\\b", "center") %>%
    str_replace_all("\\bhlth\\.?\\b", "health") %>%
    str_replace_all("\\bsvc\\.?\\b", "services") %>%
    str_replace_all("\\bsvcs\\.?\\b", "services") %>%

    # Normalize facility type terms (treat as equivalent)
    str_replace_all("\\bmedical\\s+center\\b", "hospital") %>%
    str_replace_all("\\bhealth\\s+center\\b", "hospital") %>%
    str_replace_all("\\bhealthcare\\s+center\\b", "hospital") %>%
    str_replace_all("\\bhealth\\s+system\\b", "hospital") %>%

    # Remove corporate suffixes
    str_remove("\\s+inc\\.?$") %>%
    str_remove("\\s+corp\\.?$") %>%
    str_remove("\\s+llc\\.?$") %>%
    str_remove("\\s+lp\\.?$") %>%
    str_remove("\\s+association$") %>%
    str_remove("\\s+assn\\.?$") %>%
    str_remove("\\s+incorporated$") %>%
    str_remove("\\s+corporation$") %>%
    str_remove("\\s+company$") %>%
    str_remove("\\s+co\\.?$") %>%

    # Clean up whitespace
    str_trim() %>%
    str_squish()
}
