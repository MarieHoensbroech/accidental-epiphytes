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

# 1) Download iNat -------------------------------
fetch_page <- function(page) {
  
  url <- modify_url(
    "https://api.inaturalist.org/v1/observations",
    query = list(
      project_id = PROJECT_ID,
      per_page   = 200,
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
  
  
  # stop after reaching the last page
  if (is.null(n) || n < per_page) {
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
  anti_join(rank_lookup, by = "anc_id")
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
  left_join(rank_lookup, by = "anc_id") %>%
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
obs_core <- inat_df %>%
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
    taxon_id,
    taxon.observations_count,
    taxon.ancestor_ids,
    comments,
    photos,
    family_name, 
    genus_name
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
obs_ofv <- inat_df %>%
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
             str_replace_all("&", "-"))
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
      EpiCount_num  = as.numeric(str_extract(EpiCount, "\\d*\\.?\\d+")),
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


obs_final <- obs_expanded %>%
  left_join(tax_lineage) %>%
  identity()

#glimpse(obs_final)

# . ----------------
# 4) Calculate taxon summaries --------------------------------

taxon_summary <- obs_expanded %>%
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
  left_join(obs_expanded, by = c("taxon_id", "taxon.name")) %>% 
  left_join(tax_lineage, by = "taxon_id") 

glimpse(obs_final)

# Non-list columns to CSV 
obs_csv <- obs_final %>% select(where(~ !is.list(.)))
fwrite(obs_csv, CSV_PATH)

# Full JSON (list-columns preserved)
write_json(obs_final, JSON_PATH, pretty = TRUE, auto_unbox = TRUE, na = "null")


