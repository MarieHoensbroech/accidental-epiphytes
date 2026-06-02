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
  library(lubridate)
})

# Parameters ------------------------------------
OUT_DIR         <- "data/processed"
CSV_PATH        <- file.path(OUT_DIR, "inat_observations.csv")
JSON_PATH       <- file.path(OUT_DIR, "inat_observations_raw.json")
CLEAN_JSON_PATH <- file.path(OUT_DIR, "inat_observations_cleaned.json")
UA_STRING       <- "AccidentalEpiphyteProject" #For identification by the API

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PER_PAGE <- 200

# PROJECT_IDS <- c(
#   143591, # Accidental epiphytes project ID 
#   190374, # epiphytic-plants-on-plants-in-europe-s-l-and-mediterranean
#   160915, # Epiphytes on Phoenix palm (New Zealand)
#   160824, # Accidental epiphytes in New Zealand
#   279492,  # epiphyte-diversity-and-host-tree-associations-on-southern-vancouver-island
#   15728, #Moreton bay fig epiphytes (NZ)
# )

#UMBRELLA PROJECT _---NOT PROPERLY CURATED
PROJECT_IDS <- c("accidental-epiphytes-united")

# 1) Download iNat -------------------------------

fetch_page <- function(project_id, page) {
  
  url <- modify_url(
    "https://api.inaturalist.org/v1/observations",
    query = list(
      project_id  = project_id,
      iconic_taxa = "Plantae",
      per_page    = PER_PAGE,
      page        = page
    )
  )
  
  resp <- GET(url, user_agent(UA_STRING), timeout(60))
  stop_for_status(resp)
  
  txt <- content(resp, as = "text", encoding = "UTF-8")
  fromJSON(txt, flatten = TRUE)
}

all_results <- list()

for (project_id in PROJECT_IDS) {
  
  message(sprintf("Starting download for project %s", project_id))
  
  page <- 1L
  
  repeat {
    
    dat <- fetch_page(project_id, page)
    res <- dat$results
    
    n <- if (!is.null(res)) nrow(res) else 0L
    
    if (!is.null(res) && n > 0L) {
      res$source_project_id <- project_id
      all_results[[length(all_results) + 1L]] <- res
    }
    
    message(sprintf(
      "Fetched project %s, page %d with %d records",
      project_id, page, n
    ))
    
    if (is.null(n) || n < PER_PAGE) {
      message(sprintf("Reached last page for project %s. Stopping.", project_id))
      break
    }
    
    page <- page + 1L
    Sys.sleep(0.25)
  }
}

raw_df <- bind_rows(all_results) %>%
  as_tibble()

write_json(
  raw_df,
  JSON_PATH,
  dataframe = "rows",
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

# . -----------------
# 2) Add taxonomy information ------------------------------

fetch_taxa_chunk <- function(ids) {
  
  ids <- ids[!is.na(ids)]
  ids <- unique(as.integer(ids))
  
  if (length(ids) == 0) {
    return(tibble())
  }
  
  url <- paste0(
    "https://api.inaturalist.org/v1/taxa?",
    paste(paste0("id=", ids), collapse = "&"),
    "&is_active=any"
  )
  
  resp <- GET(url, user_agent(UA_STRING), timeout(60))
  stop_for_status(resp)
  
  txt <- content(resp, as = "text", encoding = "UTF-8")
  out <- fromJSON(txt, flatten = TRUE)
  
  if (is.null(out$results) || length(out$results) == 0) {
    return(tibble())
  }
  
  as_tibble(out$results)
}

fetch_taxa_by_ids <- function(ids) {
  
  ids <- ids[!is.na(ids)]
  ids <- unique(as.integer(ids))
  
  if (length(ids) == 0) {
    return(tibble())
  }
  
  ids %>%
    split(ceiling(seq_along(.) / 100)) %>%
    map_dfr(fetch_taxa_chunk) %>%
    distinct(id, .keep_all = TRUE)
}

fetch_observations_chunk <- function(ids) {
  
  ids <- ids[!is.na(ids)]
  ids <- unique(as.integer(ids))
  
  if (length(ids) == 0) {
    return(tibble())
  }
  
  url <- paste0(
    "https://api.inaturalist.org/v1/observations?",
    paste(paste0("id=", ids), collapse = "&")
  )
  
  resp <- GET(url, user_agent(UA_STRING), timeout(60))
  stop_for_status(resp)
  
  txt <- content(resp, as = "text", encoding = "UTF-8")
  out <- fromJSON(txt, flatten = TRUE)
  
  if (is.null(out$results) || length(out$results) == 0) {
    return(tibble())
  }
  
  as_tibble(out$results)
}

fetch_observations_by_ids <- function(ids) {
  
  ids <- ids[!is.na(ids)]
  ids <- unique(as.integer(ids))
  
  if (length(ids) == 0) {
    return(tibble())
  }
  
  ids %>%
    split(ceiling(seq_along(.) / 100)) %>%
    map_dfr(fetch_observations_chunk) %>%
    distinct(id, .keep_all = TRUE)
}

tax_ids <- raw_df %>%
  reframe(
    taxon_id       = as.integer(taxon.id),
    taxon_rank     = taxon.rank,
    taxon_name     = taxon.name,
    taxon_ancestry = taxon.ancestry
  ) %>%
  filter(!is.na(taxon_id)) %>%
  distinct(taxon_id, .keep_all = TRUE) %>%
  mutate(taxon_ancestry = replace_na(taxon_ancestry, "")) %>%
  mutate(anc_ids = str_split(taxon_ancestry, "/", simplify = FALSE)) %>%
  mutate(anc_ids = map(anc_ids, ~ .x[.x != ""])) %>%
  mutate(anc_ids = map(anc_ids, ~ as.integer(.x))) %>%
  mutate(anc_ids = map2(anc_ids, taxon_id, c))

all_ids_vec <- tax_ids %>%
  pull(anc_ids) %>%
  unlist(use.names = FALSE) %>%
  tibble(anc_id = .) %>%
  drop_na() %>%
  distinct() %>%
  pull(anc_id) %>%
  as.integer()

rank_lookup <- fetch_taxa_by_ids(all_ids_vec) %>%
  transmute(
    anc_id   = as.integer(id),
    anc_rank = rank,
    anc_name = name
  ) %>%
  distinct(anc_id, .keep_all = TRUE)

missing_ids <- tibble(anc_id = all_ids_vec) %>%
  anti_join(rank_lookup, by = "anc_id") %>%
  pull(anc_id)

if (length(missing_ids) > 0) {
  
  rank_lookup <- bind_rows(
    rank_lookup,
    fetch_taxa_by_ids(missing_ids) %>%
      transmute(
        anc_id   = as.integer(id),
        anc_rank = rank,
        anc_name = name
      )
  ) %>%
    distinct(anc_id, .keep_all = TRUE)
}

tax_lineage <- tax_ids %>%
  unnest_longer(anc_ids, values_to = "anc_id") %>%
  group_by(taxon_id) %>%
  mutate(anc_pos = row_number()) %>%
  ungroup() %>%
  left_join(rank_lookup, by = "anc_id") %>%
  filter(anc_rank %in% c("family", "class", "genus")) %>%
  group_by(taxon_id, anc_rank) %>%
  slice_max(anc_pos, n = 1, with_ties = FALSE) %>%
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

source_projects <- raw_df %>%
  mutate(id = as.character(id)) %>%
  group_by(id) %>%
  summarise(
    source_project_ids = paste(sort(unique(source_project_id)), collapse = ","),
    .groups = "drop"
  )

obs_core <- raw_df %>%
  mutate(id = as.character(id)) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    taxon.status = case_when(
      taxon.endemic    ~ "endemic",
      taxon.introduced ~ "introduced",
      taxon.native     ~ "native",
      TRUE             ~ NA_character_
    )
  ) %>%
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
    taxon.name   = taxon.name,
    community_taxon_id,
    taxon.rank,
    taxon.status,
    taxon_id = taxon.id,
    taxon.ancestor_ids,
    comments,
    photos
  ) %>%
  left_join(source_projects, by = "id") %>%
  separate_wider_delim(
    location,
    delim = ",",
    names = c("latitude", "longitude"),
    cols_remove = FALSE,
    too_few = "align_start"
  ) %>%
  mutate(
    latitude  = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  )

## Observation fields ------------------------
# observation fields are stored in list-column (OFV):  
# from list -> long -> wide 
#clean host species to genus



clean_to_genus <- function(x) {
  # normalize first
  y <- x %>%
    str_to_lower() %>%
    str_replace_all("\\bcf\\b", "")  
  
  # map common names/typos to Latin genus (English + German hints)
  y <- case_when(
    
    is.na(y) ~ "unknown",
    str_detect(y, "\\bdead\\b|\\btot\\b")                                ~ "robinia",
    str_detect(y, "\\bpine\\b")                                          ~ "pinus",
    str_detect(y, "\\bplum\\b|\\bcherry\\b")                             ~ "prunus",
    str_detect(y, "\\bahorn\\b")                                         ~ "acer",
    str_detect(y, "\\bsvensk\\b")                                   ~ "sorbus",
    str_detect(y, "\\brubinia(s)?\\b|\\brobinie(n|s)?\\b|\\brubinie(n|s)?\\b") ~ "robinia",
    str_detect(y, "\\bvogelbeere\\b|\\beuropean\\s+ash\\b")              ~ "sorbus",
    str_detect(y, "\\bphoenix\\b|\\bcanariensis\\b")                     ~ "phoenix",
    str_detect(y, "\\bpalm(s)?\\b")                                      ~ "arecaceae",
    str_detect(y, "\\bpoplar(s)?\\b|\\bpappel(n|s)?\\b")                 ~ "populus",
    str_detect(y, "\\bbeech(es)?\\b|\\bbuche(n|s)?\\b")                  ~ "fagus",
    str_detect(y, "\\boak(s)?\\b|\\beiche(n|s)?\\b")                     ~ "quercus",
    str_detect(y, "\\bwillow(s)?\\b|\\bweide(n|s)?\\b")                  ~ "salix",
    str_detect(y, "\\bbirch(es)?\\b|\\bbirke(n)?\\b")                    ~ "betula",
    str_detect(y, "\\btila\\b|\\btil(l)?ia\\b|\\blinde(n)?\\b")          ~ "tilia",
    str_detect(y, "\\bsambuc(a|us)(es)?\\b|\\bholunder\\b")              ~ "sambucus",
    
    str_detect(y, "\\btree\\s*fern\\b") ~ "treefern",
    str_detect(y, "\\bseq[uio]{3}a\\b") ~ "sequoia",
    str_detect(
      y,
      "(^\\s*$|^\\s*\\?\\s*$|\\b(unknown|un|unclear|unsure|not\\s+sure|see\\s+obs|see|photo|picture|n/?a|na|cf)\\b)"
    ) ~ "unknown",
    
    
    TRUE ~ y
  )
  
  y %>%
    str_replace_all("[^\\p{L} ]+", " ") %>%  # keep only letters & spaces
    str_squish() %>%
    str_to_title() %>%
    word(1)                                   # keep the Genus only
}


obs_ofv <- raw_df %>%
  reframe(id = as.character(id), ofvs) %>%
  unnest(ofvs, keep_empty = TRUE, names_sep = "") %>%
  reframe(
    id,
    name  = paste0("field:", tolower(ofvsname)),
    value = ofvsvalue
  ) %>%
  distinct() %>%
  group_by(id, name) %>%
  summarise(
    value = {
      v <- value[!is.na(value) & value != ""]
      if (length(v) == 0) NA_character_ else as.character(v[1])
    },
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = name, values_from = value)

## 4) Join & clean observation fields --------------------------------

existing_cols <- function(df, cols) {
  cols[cols %in% names(df)]
}

coalesce_fields <- function(df, cols) {
  
  cols_present <- existing_cols(df, cols)
  
  if (length(cols_present) == 0) {
    rep(NA_character_, nrow(df))
  } else {
    df %>%
      select(all_of(cols_present)) %>%
      mutate(across(everything(), as.character)) %>%
      coalesce(!!!.)
  }
}

to_na <- function(x) {
  
  x <- as.character(x)
  x <- str_squish(x)
  x_lower <- str_to_lower(x)
  
  x[x_lower %in% c("", "na", "n/a", "nan", "none", "unknown", "unclear", "?", "-", "999", "99")] <- NA_character_
  
  x
}

looks_like_inat_id <- function(x) {
  str_detect(as.character(x), "^\\d+$")
}

normalize_measurement_text <- function(x) {
  
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_squish() %>%
    str_replace_all("[\u2012\u2013\u2014\u2212]", "-") %>%
    str_replace_all(",", ".") %>%
    str_replace_all("many", "10") %>%
    str_replace_all("∞", "10") %>%
    str_replace_all("&", "-") %>%
    str_replace_all("^~\\s*", "") %>%
    str_replace_all("^c\\.?\\s*", "") %>%
    str_replace_all("^ca\\.?\\s*", "") %>%
    str_replace_all("\\s+cm\\b", "cm") %>%
    str_replace_all("\\s+mm\\b", "mm") %>%
    str_replace_all("\\s+m\\b", "m")
}

detect_unit_one <- function(x) {
  
  if (length(x) == 0 || is.na(x) || str_squish(as.character(x)) == "") {
    return(NA_character_)
  }
  
  x_low <- normalize_measurement_text(x)
  
  if (length(x_low) == 0 || is.na(x_low) || x_low == "") {
    return(NA_character_)
  }
  
  if (str_detect(x_low, "\\bcm\\b|centimet")) return("cm")
  if (str_detect(x_low, "\\bmm\\b|millimet")) return("mm")
  if (str_detect(x_low, "\\bin\\b|inch|inches|\"")) return("in")
  if (str_detect(x_low, "\\bft\\b|feet|foot|'")) return("ft")
  if (str_detect(x_low, "\\bm\\b|metre|meter")) return("m")
  
  NA_character_
}

extract_numbers_one <- function(x) {
  
  if (length(x) == 0 || is.na(x) || str_squish(as.character(x)) == "") {
    return(numeric(0))
  }
  
  x_norm <- normalize_measurement_text(x)
  
  if (length(x_norm) == 0 || is.na(x_norm) || x_norm == "") {
    return(numeric(0))
  }
  
  nums <- str_extract_all(x_norm, "\\d*\\.?\\d+")[[1]]
  nums <- suppressWarnings(as.numeric(nums))
  nums[is.finite(nums)]
}

convert_one_value_to_m <- function(val, unit = NA_character_, assume_cm_if_gt = NULL) {
  
  if (length(val) == 0 || is.na(val) || !is.finite(val)) {
    return(NA_real_)
  }
  
  if (length(unit) == 0 || is.na(unit)) {
    unit <- NA_character_
  }
  
  if (val %in% c(999, 9999)) {
    return(NA_real_)
  }
  
  if (identical(unit, "cm")) return(val / 100)
  if (identical(unit, "mm")) return(val / 1000)
  if (identical(unit, "in")) return(val * 0.0254)
  if (identical(unit, "ft")) return(val * 0.3048)
  if (identical(unit, "m"))  return(val)
  
  if (!is.null(assume_cm_if_gt) && val > assume_cm_if_gt) {
    return(val / 100)
  }
  
  val
}

convert_measurement_one_to_m <- function(x, assume_cm_if_gt = NULL, use_range_mean = TRUE) {
  
  nums <- extract_numbers_one(x)
  unit <- detect_unit_one(x)
  
  if (length(nums) == 0) {
    return(NA_real_)
  }
  
  val <- if (use_range_mean && length(nums) >= 2) {
    mean(nums[1:2], na.rm = TRUE)
  } else {
    nums[1]
  }
  
  out <- convert_one_value_to_m(
    val,
    unit = unit,
    assume_cm_if_gt = assume_cm_if_gt
  )
  
  if (length(out) == 0 || is.na(out) || !is.finite(out)) {
    return(NA_real_)
  }
  
  as.numeric(out)
}

convert_measurement_to_m <- function(x, assume_cm_if_gt = NULL, use_range_mean = TRUE) {
  
  vapply(
    x,
    FUN = function(xx) {
      convert_measurement_one_to_m(
        xx,
        assume_cm_if_gt = assume_cm_if_gt,
        use_range_mean = use_range_mean
      )
    },
    FUN.VALUE = numeric(1),
    USE.NAMES = FALSE
  )
}

measurement_values_to_m <- function(x, n_out = 1L, assume_cm_if_gt = NULL) {
  
  if (length(n_out) == 0 || is.na(n_out) || n_out < 1L) {
    n_out <- 1L
  }
  
  n_out <- as.integer(n_out)
  
  if (length(x) == 0 || is.na(x) || str_squish(as.character(x)) == "") {
    return(rep(NA_real_, n_out))
  }
  
  x_norm <- normalize_measurement_text(x)
  unit   <- detect_unit_one(x)
  nums   <- extract_numbers_one(x)
  
  vals <- if (length(nums) >= 2 && str_detect(x_norm, "-")) {
    seq(from = nums[1], to = nums[2], length.out = n_out)
  } else if (length(nums) >= 1) {
    rep(nums[1], n_out)
  } else {
    rep(NA_real_, n_out)
  }
  
  vapply(
    vals,
    FUN = function(v) {
      convert_one_value_to_m(
        v,
        unit = unit,
        assume_cm_if_gt = assume_cm_if_gt
      )
    },
    FUN.VALUE = numeric(1),
    USE.NAMES = FALSE
  )
}

clean_setting <- function(x) {
  x %>%
    to_na() %>%
    str_squish() %>%
    str_to_title()
}

clean_yes_noish <- function(x) {
  x %>%
    to_na() %>%
    str_squish() %>%
    str_to_lower()
}

##Synonyms-----------

field_map <- list(
  
  DBH = c(
    "field:dbh of host tree (m) (diameter at breast height)"
  ),
  
  EpiHeight = c(
    "field:height on tree etc (m)",
    "field:height off ground",
    "field:height (m)"
  ),
  
  EpiSize = c(
    "field:size of accidental (m)",
    "field:estimated height",
    "field:estimate the height of the plant (in m)"
  ),
  
  EpiCount = c(
    "field:count of individuals observed",
    "field:count",
    "field:abundance category"
  ),
  
  TreeHeight = c(
    "field:estimated tree height (m)"
  ),
  
  Setting = c(
    "field:setting of area",
    "field:built environment type",
    "field:habitat",
    "field:hábitat",
    "field:habitat (s afr)",
    "field:local habitat description",
    "field:specifics of the site",
    "field:about the site",
    "field:madagascar vegetation type"
  ),
  
  Moss = c(
    "field:moss?"
  ),
  
  ObsGroup = c(
    "field:observation group",
    "field:similar observation set",
    "field:related observations",
    "field:associated observation",
    "field:host observation",
    "field:linked observation"
  ),
  
  GrowingSite = c(
    "field:growing site on the tree",
    "field:location on tree",
    "field:physical habitat/substrate",
    "field:substrate/soil type",
    "field:substrate_detail2"
  ),
  
  TreeSp = c(
    "field:epiphytic on",
    "field:host tree species",
    "field:host plant",
    "field:host",
    "field:epiphyte host",
    "field:name of associated plant",
    "field:associated species",
    "field:associated species with names lookup",
    "field:host species with names lookup",
    "field:interaction->host id",
    "field:attached to: (interaction)",
    "field:tree",
    "field:other organism",
    "field:associated taxa",
    "field:interaction->epiphyte host",
    "field:epiphytic on",
    "field:host plant id",
    "field:interaction->host id",
    "field:attached to: (interaction)",
    "field:associated with: (interaction)",
    "field:feeding on"
  )
)

obs_joined_raw <- tibble(
  id              = obs_ofv$id,
  DBH_raw         = coalesce_fields(obs_ofv, field_map$DBH),
  EpiHeight_raw   = coalesce_fields(obs_ofv, field_map$EpiHeight),
  EpiSize_raw     = coalesce_fields(obs_ofv, field_map$EpiSize),
  EpiCount_raw    = coalesce_fields(obs_ofv, field_map$EpiCount),
  TreeHeight_raw  = coalesce_fields(obs_ofv, field_map$TreeHeight),
  Setting_raw     = coalesce_fields(obs_ofv, field_map$Setting),
  Moss_raw        = coalesce_fields(obs_ofv, field_map$Moss),
  ObsGroup        = coalesce_fields(obs_ofv, field_map$ObsGroup),
  GrowingSite_raw = coalesce_fields(obs_ofv, field_map$GrowingSite),
  TreeSp_raw      = coalesce_fields(obs_ofv, field_map$TreeSp)
)

host_numeric_ids <- obs_joined_raw %>%
  transmute(host_id = suppressWarnings(as.integer(TreeSp_raw))) %>%
  filter(!is.na(host_id), str_detect(as.character(host_id), "^\\d+$")) %>%
  distinct() %>%
  pull(host_id)

host_observation_lookup <- if (length(host_numeric_ids) > 0) {
  fetch_observations_by_ids(host_numeric_ids) %>%
    transmute(
      host_id = as.integer(id),
      host_obs_taxon_name = as.character(taxon.name)
    ) %>%
    distinct(host_id, .keep_all = TRUE)
} else {
  tibble(
    host_id = integer(),
    host_obs_taxon_name = character()
  )
}

host_taxon_lookup <- if (length(host_numeric_ids) > 0) {
  fetch_taxa_by_ids(host_numeric_ids) %>%
    transmute(
      host_id = as.integer(id),
      host_taxon_name = as.character(name)
    ) %>%
    distinct(host_id, .keep_all = TRUE)
} else {
  tibble(
    host_id = integer(),
    host_taxon_name = character()
  )
}

obs_joined <- obs_joined_raw %>%
  mutate(
    host_id_tmp = if_else(
      looks_like_inat_id(TreeSp_raw),
      suppressWarnings(as.integer(TreeSp_raw)),
      NA_integer_
    )
  ) %>%
  left_join(host_observation_lookup, by = c("host_id_tmp" = "host_id")) %>%
  left_join(host_taxon_lookup, by = c("host_id_tmp" = "host_id")) %>%
  mutate(
    DBH        = to_na(DBH_raw),
    EpiHeight  = to_na(EpiHeight_raw),
    EpiSize    = to_na(EpiSize_raw),
    EpiCount   = to_na(EpiCount_raw),
    TreeHeight = to_na(TreeHeight_raw),
    
    Setting     = clean_setting(Setting_raw),
    Moss        = clean_yes_noish(Moss_raw),
    GrowingSite = clean_yes_noish(GrowingSite_raw),
    ObsGroup = if_else(is.na(ObsGroup),id,ObsGroup),
    
    TreeSp = case_when(
      !is.na(host_obs_taxon_name) ~ as.character(host_obs_taxon_name),
      !is.na(host_taxon_name) ~ as.character(host_taxon_name),
      !looks_like_inat_id(TreeSp_raw) ~ as.character(str_squish(TreeSp_raw)),
      TRUE ~ as.character(TreeSp_raw)
    ),
    
    TreeSp = TreeSp %>%
      str_remove_all(fixed("*")) %>%
      str_squish(),
    
    TreeSp = to_na(TreeSp),
    TreeGenus = clean_to_genus(TreeSp)
  ) %>%
  select(
    id,
    DBH,
    EpiHeight,
    EpiSize,
    EpiCount,
    TreeHeight,
    Setting,
    Moss,
    ObsGroup,
    GrowingSite,
    TreeSp,
    TreeGenus,
    everything(),
    -host_id_tmp,
    -host_obs_taxon_name,
    -host_taxon_name
  ) %>%
  left_join(obs_core, by = "id")

## Numerical conversions & units ------------------------------

obs_numeric <- obs_joined %>%
  mutate(
    DBH_num        = convert_measurement_to_m(DBH, assume_cm_if_gt = 5),
    TreeHeight_num = convert_measurement_to_m(TreeHeight),
    
    EpiCount_num = EpiCount %>%
      normalize_measurement_text() %>%
      str_extract("\\d+") %>%
      as.integer(),
    
    EpiCount_num = if_else(
      is.na(EpiCount_num) | EpiCount_num < 1L,
      1L,
      EpiCount_num
    )
  )

## Photos -------------------------------

extract_photo_urls <- function(photos_df) {
  
  if (!is.data.frame(photos_df)) {
    return(character(0))
  }
  
  col <- case_when(
    "url" %in% names(photos_df) ~ "url",
    "large_url" %in% names(photos_df) ~ "large_url",
    TRUE ~ NA_character_
  )
  
  if (is.na(col)) {
    return(character(0))
  }
  
  urls <- photos_df[[col]]
  urls[!is.na(urls)]
}

obs_photos <- obs_numeric %>%
  mutate(
    photo_urls = map(photos, extract_photo_urls),
    
    photo_gallery = map_chr(
      photo_urls,
      ~ if (length(.x) == 0) "" else
        paste0(
          '<div>',
          paste0('<img src="', .x, '"style="margin:20px;">', collapse = ""),
          '</div>'
        )
    ),
    
    photo_urls_concat = map_chr(
      photo_urls,
      ~ if (length(.x) == 0) "" else paste(.x, collapse = "|")
    )
  )

## Expand measurements by EpiCount ---------------------

expand_measurements <- function(df) {
  
  df %>%
    rowwise() %>%
    mutate(
      EpiHeight_num = list(
        measurement_values_to_m(EpiHeight, EpiCount_num)
      ),
      EpiSize_num = list(
        measurement_values_to_m(EpiSize, EpiCount_num)
      )
    ) %>%
    ungroup() %>%
    unnest(c(EpiHeight_num, EpiSize_num))
}

obs_expanded <- obs_photos %>%
  expand_measurements()

## Clean date column-----------

with_lctime_C <- function(expr) {
  
  old <- Sys.getlocale("LC_TIME")
  on.exit(try(Sys.setlocale("LC_TIME", old), silent = TRUE))
  
  try(Sys.setlocale("LC_TIME", "C"), silent = TRUE)
  force(expr)
}

normalize_observed_on <- function(x) {
  
  x %>%
    str_remove("\\([^)]*\\)") %>%
    str_replace_all("(?<=\\d)T(?=\\d)", " ") %>%
    str_replace("\\s*Z$", "") %>%
    str_replace("\\bUTC\\b", "") %>%
    str_replace_all("GMT\\s*([+-])(\\d{2}):(\\d{2})\\b", "\\1\\2\\3") %>%
    str_replace_all("GMT\\s*([+-])(\\d{2})(\\d{2})\\b", "\\1\\2\\3") %>%
    str_replace_all("GMT\\s*([+-])(\\d{2})\\b", "\\1\\200") %>%
    str_replace_all("\\bGMT\\b", "") %>%
    str_replace("(\\d)([+-]\\d{2}:\\d{2})$", "\\1 \\2") %>%
    str_replace("(\\d)([+-]\\d{4})$", "\\1 \\2") %>%
    str_replace("([+-])(\\d{2}):(\\d{2})$", "\\1\\2\\3") %>%
    str_replace_all("\\.", "/") %>%
    str_squish()
}

parse_inat_dt <- function(x, tz = "UTC") {
  
  with_lctime_C({
    
    x1 <- normalize_observed_on(x)
    
    suppressWarnings(
      parse_date_time(
        x1,
        orders = c(
          "Ymd HMS z", "Ymd HM z", "Ymd z",
          "mdY HMS z", "mdY HM z", "mdY z",
          "dmy HMS z", "dmy HM z", "dmy z",
          "Ymd HMS", "Ymd HM", "Ymd",
          "mdY HMS", "mdY HM", "mdY",
          "dmy HMS", "dmy HM", "dmy",
          "a b d Y HMS z", "a b d Y HM z", "a b d Y z",
          "a b d Y HMS", "a b d Y HM", "a b d Y",
          "B d, Y HMS z", "B d, Y HM z", "B d, Y z",
          "B d, Y HMS", "B d, Y HM", "B d, Y",
          "b d, Y HMS z", "b d, Y HM z", "b d, Y z",
          "b d, Y HMS", "b d, Y HM", "b d, Y"
        ),
        tz = tz,
        exact = FALSE
      )
    )
  })
}

obs_date_clean <- obs_expanded %>%
  drop_na(observed_on) %>%
  mutate(
    obs_dt   = parse_inat_dt(observed_on, tz = "UTC"),
    obs_date = as_date(obs_dt)
  ) %>%
  select(-obs_dt)

## Clean dbh column ---------------

append_flag <- function(old, new) {
  case_when(
    is.na(old) ~ new,
    TRUE ~ paste(old, new, sep = ";")
  )
}

obs_date_clean <- obs_date_clean %>%
  mutate(
    dbh_txt = tolower(trimws(as.character(DBH))),
    
    dbh_flag = case_when(
      is.na(dbh_txt) ~ NA_character_,
      str_detect(dbh_txt, "\\bcm\\b|centimet") ~ "unit_cm",
      str_detect(dbh_txt, "\\bmm\\b|millimet") ~ "unit_mm",
      str_detect(dbh_txt, "\\bin\\b|inch|inches|\"") ~ "unit_inches",
      str_detect(dbh_txt, "\\bft\\b|feet|foot|'") ~ "unit_feet",
      !str_detect(dbh_txt, "\\bm\\b|cm\\b|mm\\b|in\\b|inch|feet|foot|ft|\"|'") &
        suppressWarnings(as.numeric(str_extract(dbh_txt, "\\d*\\.?\\d+"))) > 5 ~ "assumed_cm",
      TRUE ~ NA_character_
    ),
    
    dbh_flag = if_else(
      !is.na(DBH_num) & (DBH_num < 0.03 | DBH_num > 3.0),
      append_flag(dbh_flag, "outside_plausible_range"),
      dbh_flag
    ),
    
    DBH_num = if_else(
      !is.na(DBH_num) & (DBH_num < 0.03 | DBH_num > 3.0),
      NA_real_,
      DBH_num
    ),
    
    dbh_note = case_when(
      str_detect(replace_na(dbh_flag, ""), "outside_plausible_range") ~ "outside [0.03, 3.0] m; set NA",
      str_detect(replace_na(dbh_flag, ""), "assumed_cm") ~ "value > 5 without unit; assumed cm and divided by 100",
      TRUE ~ NA_character_
    )
  )

obs_date_clean %>%
  summarise(
    n_total = n(),
    n_with_dbh = sum(!is.na(dbh_txt)),
    n_flagged = sum(!is.na(dbh_flag)),
    pct_flagged = 100 * n_flagged / n_with_dbh
  )

obs_date_clean %>%
  filter(!is.na(dbh_flag)) %>%
  separate_rows(dbh_flag, sep = ";") %>%
  count(dbh_flag, sort = TRUE)

summary(obs_date_clean$DBH_num)

quantile(
  obs_date_clean$DBH_num,
  probs = c(0.01, 0.5, 0.99),
  na.rm = TRUE
)

# . ----------------
# 4) Calculate taxon summaries --------------------------------

taxon_summary <- obs_date_clean %>%
  distinct(id, taxon_id, taxon.name, EpiCount_num) %>%
  group_by(taxon_id, taxon.name) %>%
  summarise(
    n_obs = n_distinct(id),
    sum_individuals = sum(EpiCount_num, na.rm = TRUE),
    .groups = "drop"
  )

# 5) Write outputs ---------------------------------

obs_final <- obs_date_clean %>%
  left_join(taxon_summary, by = c("taxon_id", "taxon.name")) %>%
  left_join(tax_lineage, by = "taxon_id")

glimpse(obs_final)

obs_final %>% 
  distinct(id)

obs_csv <- obs_final %>%
  select(where(~ !is.list(.)))

glimpse(obs_csv)

fwrite(obs_csv, CSV_PATH)

write_json(
  obs_final,
  CLEAN_JSON_PATH,
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)



# 6) Contributors -------------
# Contributors ------------------------------------

# Helper: safely extract a column
get_col <- function(df, col, default = NA) {
  if (col %in% names(df)) {
    df[[col]]
  } else {
    rep(default, nrow(df))
  }
}

# Helper: check whether a list/string value is non-empty
has_value <- function(x) {
  if (is.null(x)) return(FALSE)
  if (length(x) == 0) return(FALSE)
  if (all(is.na(x))) return(FALSE)
  if (is.character(x) && all(str_squish(x) == "")) return(FALSE)
  TRUE
}

# Helper: convert possible list-column values to character
as_chr_safe <- function(x) {
  if (is.list(x)) {
    map_chr(x, ~ {
      if (!has_value(.x)) {
        NA_character_
      } else {
        paste(as.character(.x), collapse = ",")
      }
    })
  } else {
    as.character(x)
  }
}

# Helper: extract users from nested list-columns such as identifications / annotations
extract_nested_users <- function(df, list_col, role) {
  
  if (!list_col %in% names(df)) {
    return(tibble(
      user_id = integer(),
      user_login = character(),
      user_name = character(),
      roles = character()
    ))
  }
  
  map_dfr(df[[list_col]], function(x) {
    
    if (!is.data.frame(x) || nrow(x) == 0) {
      return(tibble())
    }
    
    x <- as_tibble(x)
    
    tibble(
      user_id = suppressWarnings(as.integer(coalesce(
        get_col(x, "user.id"),
        get_col(x, "user_id"),
        get_col(x, "user.id.y")
      ))),
      user_login = as.character(coalesce(
        get_col(x, "user.login"),
        get_col(x, "user_login"),
        get_col(x, "user.login_exact"),
        get_col(x, "user.login_autocomplete")
      )),
      user_name = as.character(coalesce(
        get_col(x, "user.name"),
        get_col(x, "user_name"),
        get_col(x, "user.name_autocomplete")
      )),
      roles = role
    )
  }) %>%
    filter(!is.na(user_id) | !is.na(user_login)) %>%
    distinct()
}

# Observers ------------------------------------

contributors_observers <- raw_df %>%
  transmute(
    user_id = suppressWarnings(as.integer(user.id)),
    user_login = as.character(user.login),
    user_name = as.character(user.name),
    roles = "observer"
  ) %>%
  filter(!is.na(user_id) | !is.na(user_login)) %>%
  distinct()


# Identifiers ------------------------------------

contributors_identifiers <- extract_nested_users(
  raw_df,
  list_col = "identifications",
  role = "identifier"
)


# Annotators------------------------------------

contributors_annotators <- extract_nested_users(
  raw_df,
  list_col = "annotations",
  role = "annotator"
)

# Curators ------------------------------------
# This tries to identify curators from nested user role fields where present.
# It does NOT make additional API calls.

extract_curator_users_from_nested <- function(df, list_col) {
  
  if (!list_col %in% names(df)) {
    return(tibble(
      user_id = integer(),
      user_login = character(),
      user_name = character(),
      roles = character()
    ))
  }
  
  first_existing_chr <- function(x, cols) {
    cols_present <- cols[cols %in% names(x)]
    
    if (length(cols_present) == 0) {
      return(rep(NA_character_, nrow(x)))
    }
    
    vals <- rep(NA_character_, nrow(x))
    
    for (cc in cols_present) {
      candidate <- as.character(x[[cc]])
      candidate <- na_if(str_squish(candidate), "")
      
      vals <- if_else(
        is.na(vals) & !is.na(candidate),
        candidate,
        vals
      )
    }
    
    vals
  }
  
  first_existing_int <- function(x, cols) {
    cols_present <- cols[cols %in% names(x)]
    
    if (length(cols_present) == 0) {
      return(rep(NA_integer_, nrow(x)))
    }
    
    vals <- rep(NA_integer_, nrow(x))
    
    for (cc in cols_present) {
      candidate <- suppressWarnings(as.integer(x[[cc]]))
      
      vals <- if_else(
        is.na(vals) & !is.na(candidate),
        candidate,
        vals
      )
    }
    
    vals
  }
  
  role_text_from_cols <- function(x, role_cols) {
    
    if (length(role_cols) == 0) {
      return(rep("", nrow(x)))
    }
    
    x %>%
      select(all_of(role_cols)) %>%
      mutate(across(everything(), ~ as.character(.x))) %>%
      mutate(across(everything(), ~ replace_na(.x, ""))) %>%
      unite("role_text", everything(), sep = ",", remove = TRUE, na.rm = TRUE) %>%
      pull(role_text) %>%
      str_to_lower()
  }
  
  map_dfr(df[[list_col]], function(x) {
    
    if (!is.data.frame(x) || nrow(x) == 0) {
      return(tibble())
    }
    
    x <- as_tibble(x)
    
    role_cols <- names(x)[str_detect(names(x), regex("role|curator", ignore_case = TRUE))]
    
    if (length(role_cols) == 0) {
      return(tibble())
    }
    
    role_text <- role_text_from_cols(x, role_cols)
    
    out <- tibble(
      user_id = first_existing_int(
        x,
        c(
          "user.id",
          "user_id",
          "user.id.y",
          "user_id.y"
        )
      ),
      user_login = first_existing_chr(
        x,
        c(
          "user.login",
          "user_login",
          "user.login_exact",
          "user.login_autocomplete",
          "user.login.y",
          "login"
        )
      ),
      user_name = first_existing_chr(
        x,
        c(
          "user.name",
          "user_name",
          "user.name_autocomplete",
          "user.name.y",
          "name"
        )
      ),
      role_text = role_text
    )
    
    out %>%
      filter(str_detect(role_text, "curator")) %>%
      transmute(
        user_id,
        user_login,
        user_name,
        roles = "curator"
      )
  }) %>%
    filter(!is.na(user_id) | !is.na(user_login)) %>%
    distinct()
}

contributors_curators <- bind_rows(
  extract_curator_users_from_nested(raw_df, "identifications"),
  extract_curator_users_from_nested(raw_df, "annotations")
) %>%
  distinct()

# Combine contributors ------------------------------------

contributors <- bind_rows(
  contributors_observers,
  contributors_identifiers,
  contributors_annotators,
  contributors_curators
) %>%
  mutate(
    user_login = na_if(str_squish(user_login), ""),
    user_name  = na_if(str_squish(user_name), "")
  ) %>%
  filter(!is.na(user_id) | !is.na(user_login)) %>%
  group_by(user_id, user_login) %>%
  summarise(
    user_name = {
      nm <- na.omit(user_name)
      if (length(nm) == 0) NA_character_ else nm[1]
    },
    roles = paste(sort(unique(roles)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(user_login)

glimpse(contributors)

CONTRIBUTORS_CSV_PATH <- file.path("outputs/Appendix3_inat_contributors.csv")

fwrite(contributors, CONTRIBUTORS_CSV_PATH)
