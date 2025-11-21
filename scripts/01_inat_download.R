# 01_inat_download.R 
#
# 1) Download iNat data
# 2) Download taxon information (e.g. family, class, genus...)
# 3) Clean observation fields
# 4) Calculate summaries
# 5) Write to file


suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(data.table)
})

# Parameters ------------------------------------
PROJECT_ID   <- 143591        # Accidental epiphytes project ID 

OUT_DIR      <- "data/processed"
CSV_PATH     <- file.path(OUT_DIR, "inat_observations.csv")
JSON_PATH    <- file.path(OUT_DIR, "inat_observations.json")
UA_STRING    <- "AccidentalEpiphyteProject" #For identification by the API

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PER_PAGE <- 200

# 1) Download iNat -------------------------------
fetch_page <- function(page) {
  
  url <- modify_url(
    "https://api.inaturalist.org/v1/observations",
    query = list(
      project_id = PROJECT_ID,
      per_page   = PER_PAGE,
      page       = page
    )
  )
  
  resp <- GET(url, user_agent(UA_STRING), timeout(60))
  stop_for_status(resp)
  txt <- content(resp, as = "text", encoding = "UTF-8")
  fromJSON(txt, flatten = TRUE)
}

## Begin download ------
# This includes observation fields but not taxonomy information (e.g. family)
all_results <- list()
page <- 1L

repeat {
  dat <- fetch_page(page)
  res <- dat$results
  
  all_results[[page]] <- res
  message(sprintf("Fetched page %d with %d records", page, nrow(res)))
  
  n   <- if (!is.null(res)) nrow(res) else 0L
  # stop after reaching the last page
  if (is.null(n) || n < 200L) {
    message("Reached the last page. Stopping.")
    break
  }
  
  
  page <- page + 1L
  Sys.sleep(0.25)   #set a time-out to not strain the API too much
}

raw_df <- bind_rows(all_results) %>% as_tibble()

write_json(
  raw_df,
  JSON_PATH,
  dataframe = "rows",   
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

# Read the JSON file into a list of rows
#raw_df <- fromJSON(JSON_PATH, flatten = TRUE)  %>%  as_tibble()

#glimpse(raw_df)

# . -----------------
# 2) Add taxonomy information ------------------------------
# Parse ancestry into ordered integer vectors and append the focal id
tax_ids <- raw_df %>%
  reframe(
    taxon_id       = as.integer(taxon.id),
    taxon_rank     = taxon.rank,
    taxon_name     = taxon.name,
    taxon_ancestry = taxon.ancestry
  ) %>%
  distinct(taxon_id, .keep_all = TRUE) %>%
  mutate(taxon_ancestry = replace_na(taxon_ancestry, "")) %>%
  mutate(anc_ids = str_split(taxon_ancestry, "/", simplify = FALSE)) %>%
  mutate(anc_ids = map(anc_ids, ~ .x[.x != ""])) %>%
  mutate(anc_ids = map(anc_ids, ~ as.integer(.))) %>%
  mutate(anc_ids = map2(anc_ids, taxon_id, c))   # append focal species' taxon id

# 2) Build a clean integer vector of all unique IDs and fetch ranks/names 
all_ids_vec <- tax_ids %>%
  pull(anc_ids) %>%
  unlist(use.names = FALSE) %>%
  tibble(anc_id = .) %>%
  drop_na() %>%               
  distinct() %>%
  pull(anc_id) %>%
  as.integer()

rank_lookup <- all_ids_vec %>%
  split(ceiling(seq_along(.) / 100)) %>%   # chunk into bits not to overload the API
  map_dfr(~ {
    paste0(
      "https://api.inaturalist.org/v1/taxa?",
      paste(paste0("id=", .x), collapse = "&"),   
      "&is_active=any"                             
    ) %>%
      GET(user_agent(UA_STRING), timeout(60)) %>%
      content(as = "text", encoding = "UTF-8") %>%
      fromJSON(flatten = TRUE) %>%
      pluck("results")
  }) %>%
  transmute(
    anc_id   = as.integer(id),
    anc_rank = rank,
    anc_name = name
  ) %>%
  distinct(anc_id, .keep_all = TRUE)

# # Refetch any missing IDs
missing_ids <- tibble(anc_id = all_ids_vec) %>%
  anti_join(rank_lookup)
# 
# rank_lookup <- bind_rows(
#   rank_lookup,
#   missing_ids$anc_id %>%
#     split(ceiling(seq_along(.) / 100)) %>%
#     map_dfr(~ {
#       paste0(
#         "https://api.inaturalist.org/v1/taxa?",
#         paste(paste0("id=", .x), collapse = "&"),
#         "&is_active=any"
#       ) %>%
#         GET(user_agent(UA_STRING), timeout(60)) %>%
#         content(as = "text", encoding = "UTF-8") %>%
#         fromJSON(flatten = TRUE) %>%
#         pluck("results")
#     }) %>%
#     transmute(
#       anc_id   = as.integer(id),
#       anc_rank = rank,
#       anc_name = name
#     )
# ) %>%
#   distinct(anc_id, .keep_all = TRUE)

# 3) Use ordered ancestry to pick the nearest family / class / genus 
tax_lineage <- tax_ids %>%
  unnest_longer(anc_ids, values_to = "anc_id") %>%
  group_by(taxon_id) %>%
  mutate(anc_pos = row_number()) %>%        
  ungroup() %>%
  left_join(rank_lookup) %>%
  filter(anc_rank %in% c("family", "class", "genus")) %>%
  group_by(taxon_id, anc_rank) %>%
  slice_max(anc_pos, n = 1, with_ties = FALSE) %>%  # nearest of each rank
  ungroup() %>%
  select(taxon_id, anc_rank, anc_id, anc_name) %>%
  pivot_wider(
    names_from  = anc_rank,
    values_from = c(anc_id, anc_name),
    names_sep   = "_"
  ) %>%
  transmute(
    taxon_id,
    family_id   = as.integer(anc_id_family),
    family_name = anc_name_family,
    class_id    = as.integer(anc_id_class),
    class_name  = anc_name_class,
    genus_id    = as.integer(anc_id_genus),
    genus_name  = anc_name_genus
  )



# . ---------------
# 3) Tidy up ---------------------------------

# Core, keep non-list columns
obs_core <- raw_df %>%
  mutate(id = as.character(id)) %>%
  select(
    id,
    quality_grade,
    observed_on = observed_on_string,
    description,
    captive,
    uri,
    location,
    taxon.common = taxon.preferred_common_name,
    user_login   = user.login,
    user_observations_count = user.observations_count,
    user_activity_count     = user.activity_count,
    user_species_count      = user.species_count,
    taxon.name  = taxon.name,
    taxon.endemic,
    taxon.rank,
    taxon.introduced,
    taxon.native,
    taxon_id = taxon.id,
    taxon.observations_count,
    taxon.ancestor_ids,
    comments,
    photos
  ) %>%
  # split "lat,lon"
  separate_wider_delim(location, delim = ",", names = c("latitude", "longitude"), cols_remove = FALSE) %>%
  mutate(
    latitude  = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  )

## Observation fields ------------------------
# observation fields are stored in list-column (OFV):  
# from list -> long -> wide 
obs_ofv <- raw_df %>%
  reframe(id = as.character(id), ofvs) %>%
  unnest(ofvs, keep_empty = TRUE, names_sep = "") %>%
  reframe(
    id,
    name  = paste0("field:", tolower(ofvsname)),
    value = ofvsvalue
  ) %>%
  distinct() %>%
  pivot_wider(names_from = name, values_from = value)

## Join & clean fields --------------------------------

obs_joined <- obs_core %>%
  left_join(obs_ofv) %>%
  # Standardize observation fields
  rename(
    DBH       = `field:dbh of host tree (m) (diameter at breast height)`,
    EpiHeight = `field:height on tree etc (m)`,
    EpiSize   = `field:size of accidental (m)`,
    EpiCount  = `field:count of individuals observed`,
    TreeHeight= `field:estimated tree height (m)`,
    Setting   = `field:setting of area`,
    Moss      = `field:moss?`,
    ObsGroup  = `field:observation group`,
    GrowingSite = `field:growing site on the tree`,
    TreeSp    = `field:host tree species`
  ) %>%
  # Turn text to numbers for
  mutate(
    across(c(EpiCount, EpiSize, EpiHeight),
           ~ .x %>%
             str_replace_all("many", "10") %>%
             str_replace_all("∞", "10") %>%
             str_replace_all("&", "-") %>% 
             str_replace_all(",", "."))
  )

## Numerical conversions & units ------------------------------
# Expand EpiHeight/EpiSize given EpiCount (range like "1-3" -> seq length EpiCount)

# Convert DBH/TreeHeight text to numerics (extract first number); 
obs_numeric <- obs_joined %>%
  mutate(
    DBH_num       = as.numeric(str_extract(DBH, "\\d*\\.?\\d+")),
    TreeHeight_num= as.numeric(str_extract(TreeHeight, "\\d*\\.?\\d+"))
  ) %>%
  # If values are in cm, convert to m 
  mutate(
    DBH_num = if_else(DBH_num == 999 | DBH_num == 9999, NA, DBH_num),
    DBH_num       = if_else(!is.na(DBH) & str_detect(DBH, "cm"), DBH_num/100, DBH_num),
    TreeHeight_num= if_else(!is.na(TreeHeight) & str_detect(TreeHeight, "cm"), TreeHeight_num/100, TreeHeight_num)
  )



## Photos -------------------------------
# Extract all photo URLs per observation
extract_photo_urls <- function(photos_df) {
  if (!is.data.frame(photos_df)) return(character(0))
  col <- if ("url" %in% names(photos_df)) {
    "url"
  } else if ("large_url" %in% names(photos_df)) {
    "large_url"
  } else {
    return(character(0))
  }
  urls <- photos_df[[col]]
  urls[!is.na(urls)]
}

obs_photos <- obs_numeric %>%
  mutate(
    photo_urls = map(photos, extract_photo_urls),
    # Render all photos in a simple HTML gallery
    photo_gallery = map_chr(
      photo_urls,
      ~ if (length(.x) == 0) "" else
        paste0(
          '<div>',
          paste0('<img src="', .x, '"style="margin:20px;">', collapse = ""),
          '</div>'
        )
    ),
    # Pipe-separated URLs so they survive in the CSV output
    photo_urls_concat = map_chr(photo_urls, ~ if (length(.x) == 0) "" else paste(.x, collapse = "|"))
  )

## Expand measurements by EpiCount ---------------------
expand_measurements <- function(df) {
  df %>%
    mutate(
      EpiCount_num  = as.integer(str_extract(EpiCount, "\\d+")),
      EpiCount_noNA = if_else(is.na(EpiCount_num), 1, as.integer(EpiCount_num)),
      EpiHeight_chr = as.character(EpiHeight),
      EpiSize_chr   = as.character(EpiSize)
    ) %>%
    rowwise() %>%
    mutate(
      EpiHeight_num = list({
        if (!is.na(EpiHeight_chr) && str_detect(EpiHeight_chr, "-")) {
          rng <- suppressWarnings(as.numeric(str_split(EpiHeight_chr, "-", simplify = TRUE)))
          if (length(rng) == 2 && !any(is.na(rng)) && !is.na(EpiCount_noNA)) {
            seq(from = rng[1], to = rng[2], length.out = EpiCount_noNA)
          } else rep(NA_real_, EpiCount_noNA)
        } else {
          val <- as.numeric(EpiHeight_chr)
          if (!is.na(val) && !is.na(EpiCount_noNA)) rep(val, EpiCount_noNA) else rep(NA_real_, EpiCount_noNA)
        }
      }),
      EpiSize_num   = list({
        if (!is.na(EpiSize_chr) && str_detect(EpiSize_chr, "-")) {
          rng <- as.numeric(str_split(EpiSize_chr, "-", simplify = TRUE))
          if (length(rng) == 2 && !any(is.na(rng)) && !is.na(EpiCount_noNA)) {
            seq(from = rng[1], to = rng[2], length.out = EpiCount_noNA)
          } else rep(NA_real_, EpiCount_noNA)
        } else {
          val <- as.numeric(EpiSize_chr)
          if (!is.na(val) && !is.na(EpiCount_noNA)) rep(val, EpiCount_noNA) else rep(NA_real_, EpiCount_noNA)
        }
      })
    ) %>%
    ungroup() %>%
    unnest(c(EpiHeight_num, EpiSize_num))
}

obs_expanded <- obs_photos %>%
  expand_measurements() %>%
  # Additional unit harmonization for expanded fields
  mutate(
    EpiSize_num   = if_else(!is.na(EpiSize)   & str_detect(EpiSize,   "cm"), EpiSize_num/100, EpiSize_num),
    EpiHeight_num = if_else(!is.na(EpiHeight) & str_detect(EpiHeight, "cm"), EpiHeight_num/100, EpiHeight_num)
  )


## Clean date column-----------
# Force English months/weekdays, regardless of OS locale
with_lctime_C <- function(expr) {
  old <- Sys.getlocale("LC_TIME")
  on.exit(try(Sys.setlocale("LC_TIME", old), silent = TRUE))
  try(Sys.setlocale("LC_TIME", "C"), silent = TRUE)
  force(expr)
}

normalize_observed_on <- function(x) {
  x %>%
    # Remove parentheses like "(PDT)" or "(GMT+13)" or "(GMT-5)"
    str_remove("\\([^)]*\\)") %>%
    # ISO 'T' to space
    str_replace_all("(?<=\\d)T(?=\\d)", " ") %>%
    # Drop trailing 'Z' and literal 'UTC'
    str_replace("\\s*Z$", "") %>%
    str_replace("\\bUTC\\b", "") %>%
    # --- Handle all GMT variants, allowing optional whitespace after GMT ---
    str_replace_all("GMT\\s*([+-])(\\d{2}):(\\d{2})\\b", "\\1\\2\\3") %>%  # GMT +HH:MM -> +HHMM
    str_replace_all("GMT\\s*([+-])(\\d{2})(\\d{2})\\b", "\\1\\2\\3")   %>%  # GMT +HHMM  -> +HHMM
    str_replace_all("GMT\\s*([+-])(\\d{2})\\b", "\\1\\200")            %>%  # GMT +HH    -> +HH00
    # Drop leftover 'GMT'
    str_replace_all("\\bGMT\\b", "") %>%
    # Ensure there is a space before any trailing offset (with or without colon)
    str_replace("(\\d)([+-]\\d{2}:\\d{2})$", "\\1 \\2") %>%
    str_replace("(\\d)([+-]\\d{4})$", "\\1 \\2") %>%
    # Convert colon style offsets at end -> +HHMM
    str_replace("([+-])(\\d{2}):(\\d{2})$", "\\1\\2\\3") %>%
    # Normalize separators and whitespace
    str_replace_all("\\.", "/") %>%
    str_squish()
}

parse_inat_dt <- function(x, tz = "UTC") {
  with_lctime_C({
    x1 <- normalize_observed_on(x)
    
    # Detect trailing numeric offset +HHMM/-HHMM
    has_offset <- str_detect(x1, "([+-]\\d{4})$")
    
    # Detect textual names
    has_wday  <- str_detect(x1, "^(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\\b")
    has_mname <- str_detect(x1, "\\b(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|January|February|March|April|May|June|July|August|September|October|November|December)\\b")
    has_names <- has_wday | has_mname
    
    n   <- length(x1)
    res <- as.POSIXct(rep(NA_real_, n), origin = "1970-01-01", tz = tz)
    
    # 1) Numeric-only + offset (fast C parser)
    idx <- which(has_offset & !has_names)
    if (length(idx)) {
      res[idx] <- suppressWarnings(
        parse_date_time2(
          x1[idx],
          orders = c(
            "Ymd HMS z", "Ymd HM z", "Ymd z",
            "mdY HMS z", "mdY HM z", "mdY z",
            "dmy HMS z", "dmy HM z", "dmy z"
          ),
          tz = tz
        )
      )
    }
    
    # 2) Name-based + offset (weekday/month names) — exact = TRUE
    idx <- which(has_offset & has_names & is.na(res))
    if (length(idx)) {
      res[idx] <- suppressWarnings(
        parse_date_time(
          x1[idx],
          orders = c(
            # e.g. "Wed Apr 29 2020 19:11:58 +0200"
            "a b d Y HMS z", "a b d Y HM z", "a b d Y z",
            # e.g. "October 05, 2002 11:13 +0200" (rare but possible)
            "B d, Y HMS z", "B d, Y HM z",
            "b d, Y HMS z", "b d, Y HM z"
          ),
          tz = tz,
          exact = TRUE
        )
      )
    }
    
    # 3) Numeric, no offset
    idx <- which(!has_offset & !has_names & is.na(res))
    if (length(idx)) {
      res[idx] <- suppressWarnings(
        parse_date_time2(
          x1[idx],
          orders = c(
            "Ymd HMS", "Ymd HM", "Ymd",
            "mdY HMS", "mdY HM", "mdY",
            "dmy HMS", "dmy HM", "dmy"
          ),
          tz = tz
        )
      )
    }
    
    # 4) Name-based, no offset — add Month-with-comma variants
    idx <- which(!has_offset & has_names & is.na(res))
    if (length(idx)) {
      res[idx] <- suppressWarnings(
        parse_date_time(
          x1[idx],
          orders = c(
            "a b d Y HMS", "a b d Y HM", "a b d Y",
            "B d, Y HMS", "B d, Y HM", "B d, Y",
            "b d, Y HMS", "b d, Y HM", "b d, Y"
          ),
          tz = tz,
          exact = TRUE
        )
      )
    }
    
    # 5) LAST-RESORT fallbacks with base strptime (targeted)
    need <- which(is.na(res))
    
    if (length(need)) {
      y <- x1[need]
      
      # 5a) Weekday abbr + month abbr + seconds + numeric offset
      # "%a %b %d %Y %H:%M:%S %z"
      m1 <- str_detect(y, "^(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\\s+[A-Za-z]{3}\\s+\\d{1,2}\\s+\\d{4}\\s+\\d{2}:\\d{2}:\\d{2}\\s+[+-]\\d{4}$")
      if (any(m1)) {
        res[need[m1]] <- as.POSIXct(strptime(y[m1], "%a %b %d %Y %H:%M:%S %z", tz = tz))
      }
      
      # 5b) Month FULL name with comma, no offset
      # "%B %d, %Y %H:%M:%S" or "%B %d, %Y %H:%M"
      m2 <- str_detect(y, "^[A-Za-z]+\\s+\\d{1,2},\\s+\\d{4}\\s+\\d{2}:\\d{2}(:\\d{2})?$")
      if (any(m2)) {
        y2 <- y[m2]
        has_sec <- str_detect(y2, ":\\d{2}:\\d{2}$")
        res[need[m2 & has_sec]] <- as.POSIXct(strptime(y2[has_sec], "%B %d, %Y %H:%M:%S", tz = tz))
        res[need[m2 & !has_sec]] <- as.POSIXct(strptime(y2[!has_sec], "%B %d, %Y %H:%M",    tz = tz))
      }
      
      # 5c) Pure dmy date like "30/09/2012"
      m3 <- str_detect(y, "^\\d{1,2}/\\d{1,2}/\\d{4}$")
      if (any(m3)) {
        res[need[m3]] <- as.POSIXct(strptime(y[m3], "%d/%m/%Y", tz = tz))
      }
    }
    
    res
  })
}


obs_date_clean <- obs_expanded %>%
  drop_na(observed_on) %>% 
  mutate(
    obs_dt   = parse_inat_dt(observed_on, tz = "UTC"),
    obs_date = as_date(obs_dt)
  ) %>%
  select(-obs_dt)

# 
# obs_date_clean %>% distinct(observed_on,obs_date) %>% view()

#glimpse(obs_final)


## Clean dbh column ---------------

# Replace a variety of separators with standard decimal dot
normalize_decimal_marks <- function(s) {
  s %>%
    # lowercase, trim
    tolower() %>%
    trimws() %>%
    # normalize unicode dashes to simple dash for ranges
    str_replace_all("[\u2012\u2013\u2014\u2212]", "-") %>%  # figure/en/en dash/minus
    # convert decimal comma to dot
    str_replace_all(",", ".") %>%
    # convert apostrophes/primes commonly used as decimal separators to dot
    str_replace_all("['’`]", ".") %>%
    # fix '0 4' style spacing as decimal
    str_replace_all("^([0-9])\\s+([0-9])$", "\\1.\\2") %>%
    # remove repeated dots like '0..5' -> '0.5'
    str_replace_all("\\.{2,}", ".") %>%
    # remove odd prefixes
    str_replace_all("^~\\s*", "") %>%
    str_replace_all("^c\\.?\\s*", "") %>%     # c. 15 cm
    str_replace_all("^ca\\.?\\s*", "") %>%
    # standardize ' m' spacing
    str_replace_all("\\s*m\\b", "m") %>%
    str_replace_all("\\s*cm\\b", "cm") %>%
    str_replace_all("\\s*mm\\b", "mm")
}

# Extract up to two numeric tokens (for ranges) and unit tag if present
# Returns a list with numeric vector and unit tag ("m","cm","mm", or NA)
extract_numbers_and_unit <- function(s) {
  unit <- if (str_detect(s, "cm\\b")) "cm"
  else if (str_detect(s, "mm\\b")) "mm"
  else if (str_detect(s, "m\\b"))  "m"
  else NA_character_
  
  # Allow numbers like .85, 0.85, 1., and scientific is not expected here
  nums <- str_extract_all(s, "(?<![A-Za-z])\\d*\\.?\\d+(?![A-Za-z])")[[1]]
  nums <- suppressWarnings(as.numeric(nums))
  nums <- nums[is.finite(nums)]
  list(nums = nums, unit = unit)
}

# Convert numeric(s) & unit into meters; handle ranges by mean
# Apply heuristic: if >5 and unit missing -> assume cm
convert_to_meters <- function(nums, unit) {
  if (length(nums) == 0) return(NA_real_)
  val <- if (length(nums) >= 2) mean(nums[1:2]) else nums[1]
  
  if (!is.na(unit)) {
    if (unit == "cm") return(val / 100)
    if (unit == "mm") return(val / 1000)
    if (unit == "m")  return(val)
  }
  # No unit tag: heuristic
  if (val > 5) return(val / 100)  # likely centimeters
  val
}

# Main parser: returns cleaned meters + flags/notes
parse_dbh_to_m <- function(dbh_chr) {
  s <- normalize_decimal_marks(dbh_chr)
  
  # special cases that should be NA
  if (is.na(s) || s == "" || s %in% c("na", "nan", "none")) {
    return(list(val = NA_real_, flags = NA_character_, notes = NA_character_))
  }
  if (str_detect(s, "[^0-9cm\\.\\-m ]") && !str_detect(s, "cm|mm|m\\b")) {
    # weird tokens like "0'rp" will likely be NA; we still try numbers below
    weird <- TRUE
  } else weird <- FALSE
  
  ex <- extract_numbers_and_unit(s)
  val_m <- convert_to_meters(ex$nums, ex$unit)
  
  flags <- character(0); notes <- character(0)
  if (!is.na(ex$unit)) {
    flags <- c(flags, paste0("unit_", ex$unit))
  }
  if (length(ex$nums) >= 2) {
    flags <- c(flags, "range_mean")
    notes <- c(notes, "used mean of first two numbers")
  }
  if (isTRUE(weird)) {
    flags <- c(flags, "nonstandard_chars")
    notes <- c(notes, "non-numeric tokens present")
  }
  
  # Heuristic conversion flag
  if (is.na(ex$unit) && !is.na(val_m)) {
    # recover original numeric before heuristic
    orig_num <- if (length(ex$nums)) ex$nums[1] else NA_real_
    if (!is.na(orig_num) && orig_num > 5) {
      flags <- c(flags, "converted_cm_to_m")
      notes <- c(notes, "value > 5 without unit; assumed cm and /100")
    }
  }
  
  # Post rules: non-positive, tiny, plausibility window
  if (!is.na(val_m) && val_m <= 0) {
    flags <- c(flags, "non_positive")
    notes <- c(notes, "set NA")
    val_m <- NA_real_
  }
  
  if (!is.na(val_m) && val_m < 0.01) {
    flags <- c(flags, "too_small")
    notes <- c(notes, "value < 0.01 m set NA")
    val_m <- NA_real_
  }
  
  if (!is.na(val_m) && (val_m < 0.03 || val_m > 3.0)) {
    flags <- c(flags, "outside_plausible_range")
    notes <- c(notes, "outside [0.03, 3.0] m set NA")
    val_m <- NA_real_
  }
  
  list(
    val   = val_m,
    flags = if (length(flags)) paste(flags, collapse = ";") else NA_character_,
    notes = if (length(notes)) paste(notes, collapse = ";") else NA_character_
  )
}

# Apply to inat data ----

# inat$dbh is the raw column (character/numeric mixed)
parsed <- purrr::pmap(
  list(obs_date_clean$DBH),
  ~ parse_dbh_to_m(..1)
)

inat <- obs_date_clean %>%
  mutate(
    dbh_txt   = tolower(trimws(as.character(DBH))),
    dbh_clean = purrr::map_dbl(parsed, "val"),
    dbh_flag  = purrr::map_chr(parsed, "flags"),
    dbh_note  = purrr::map_chr(parsed, "notes")
  )

#######
# How many had DBH, how many changed
inat %>%
  summarise(
    n_total = n(),
    n_with_dbh = sum(!is.na(dbh_txt)),
    n_changed = sum(!is.na(dbh_flag)),
    pct_changed = 100 * n_changed / n_with_dbh
  )

# What kinds of changes
inat %>%
  filter(!is.na(dbh_flag)) %>%
  separate_rows(dbh_flag, sep = ";") %>%
  count(dbh_flag, sort = TRUE)

# Distribution of cleaned values
summary(inat$dbh_clean)
quantile(inat$dbh_clean, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)

# Spot-check rows we converted from cm
obs_date_clean <- inat %>%
  reframe(id,user_login,DBH_num=dbh_clean) %>% distinct() %>% 
  right_join(obs_date_clean %>% select(-DBH_num) %>% distinct())




##############

# . ----------------
# 4) Calculate taxon summaries --------------------------------

taxon_summary <- obs_date_clean %>%
  distinct() %>%
  group_by(taxon_id, taxon.name) %>%
  summarise(
    n_obs  = n_distinct(id),
    sum_individuals = sum(as.numeric(replace_na(as.numeric(str_extract(EpiCount, "\\d+")), 1L)), na.rm = TRUE),
    .groups = "drop"
  )

# 5) Write outputs ---------------------------------
# Merge all datafranes
obs_final <- taxon_summary %>%
  left_join(obs_date_clean) %>% 
  left_join(tax_lineage) 

glimpse(obs_final)

# Non-list columns to CSV 
obs_csv <- obs_final %>% select(where(~ !is.list(.)))
glimpse(obs_csv)
fwrite(obs_csv, CSV_PATH)

# Full JSON (list-columns preserved)
write_json(obs_final, "inat_observations_cleaned.json", pretty = TRUE, auto_unbox = TRUE, na = "null")


