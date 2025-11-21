# 11a_dispersal_analyses.R
# - obtain information on dispersal syndrome from 
# 1) BIEN 
# 2) EUDIS 
# 3) TRY 
# clean and collate information 
# write output 
# for use in: 
# 11b_dispersal_analyses.R 
# Assess what dispersal syndromes accidental epiphytes use 
# 1) using traditional stats 
# 2) using Bayesian 
# 3) Plotting


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(viridis)
  library(forcats)
  library(scales)
  library(maps)
  library(patchwork)
  library(stringr)
  library(data.table)
  library(cowplot)
  library(grid)
  library(ggimage)
  library(rnaturalearth)
  library(stringr)
  library(brms)
  library(tidybayes)
  library(emmeans)
  library(multcompView)
  library(rstatix)
  library(janitor)
  library(BIEN)
  library(purrr)
})

# Parameters -------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"
IN_EUDIS      <- "data/raw/eudis_v6.csv"
IN_TRY        <- "data/raw/44974_01112025180040/44974.txt"

OUT_DIR_DISP  <- "outputs/11_dispersal/"
dir.create(OUT_DIR_DISP, showWarnings = FALSE, recursive = TRUE)
update_geom_defaults("point", list(alpha = 0.8))

my_sites   <- fread(IN_MY_SITES)   %>% as_tibble() %>% unique() %>% clean_names()
inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% unique() %>% clean_names()
inat_obs   <- fread(IN_INAT_OBS)   %>% as_tibble() %>% 
  mutate(lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>% 
  drop_na(latitude,longitude) %>% clean_names()

# Standard ggplot theme 
my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top"
  )

# _________ ----------
# Prep data -----------
df_host <- my_sites %>%
  mutate(
    lon = as.numeric(site_lon),
    lat = as.numeric(site_lat),
    total_trees = as.numeric(number_trees),
    obs_count_dbh = suppressWarnings(as.numeric(n_epis_per_site_cat)),
    spp_count_dbh = suppressWarnings(as.numeric(n_sp_per_site_cat)),
    obs_count = suppressWarnings(as.numeric(n_epis_per_site)),
    spp_count  = suppressWarnings(as.numeric(n_sp_per_site)),
    dbh_cat=tree_cat
  ) %>%
  group_by(site,dbh_cat) %>%
  reframe(
    total_trees = sum(number_trees, na.rm = TRUE),
    obs_count_dbh,spp_count_dbh,obs_count,spp_count,
    lon = first(lon),
    lat = first(lat),
    setting,
    obs_per10 = 10 * obs_count_dbh / total_trees,
    spp_per10 = 10 * spp_count_dbh   / total_trees,
    dominanttree1 = word(dominanttree1, 1),
    dominanttree2 = word(dominanttree2, 1),
    dominanttree3 = word(dominanttree3, 1)
  ) %>%
  filter(is.finite(lat), is.finite(lon), !(lat == 0 & lon == 0),
         dbh_cat!="0 - 10") %>%
  unique()

inat_field <-
  inat_field %>%
  drop_na(latitude,longitude) %>%
  reframe(site,
          taxon_name, genus_name, tree_sp, 
          epi_count=epi_count_no_na, epi_lat=latitude, epi_lon=longitude,
          dbh_cat=tree_cat, taxon_rank) %>%
  right_join(df_host) 

inat_obs <- inat_obs %>% #glimpse()
  reframe(taxon_name, genus_name, tree_sp, 
          epi_count=epi_count_no_na, 
          epi_lat=latitude, 
          epi_lon=longitude, 
          dbh=dbh_num, taxon_rank,
          family = family_name,
          setting) %>%
  drop_na(epi_lat,epi_lon)

Mode_ <- function(x) {
  x <- x[is.finite(match(x, x)) & !is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}
collapse_values <- function(x, k = 3L) {
  ux <- unique(na.omit(x))
  if (length(ux) <= k) paste(ux, collapse = " | ") else Mode_(x)
}

# _________________________ ----
# 1) Dispersal by species ---------------------------------------------------------------
# ## A) BIEN ---------------
# species_vec <- inat_obs %>%
#   distinct(taxon_name) %>% pull() %>% na.omit() %>% unique()
# 
# bien_ver <- tryCatch(BIEN_metadata_database_version(), error=function(e) NA_character_)
# message("Using BIEN DB version: ", bien_ver)
# 
# trait_catalog <- tryCatch(BIEN_trait_list() %>% clean_names(), error=function(e) tibble())
# target_traits <- c(
#   "fruit type","longest whole plant longevity","maximum whole plant longevity",
#   "whole plant dispersal syndrome"
# )
# 
# species_vec <- inat_obs %>% filter(taxon_rank%in%c("complex","species")) %>%
#   drop_na(taxon_name) %>% pull(taxon_name) %>% unique()
# 
# retry <- function(expr, times = 3, wait = 3) {
#   for (i in seq_len(times)) {
#     out <- try(eval.parent(substitute(expr)), silent = TRUE)
#     if (!inherits(out, "try-error")) return(out)
#     Sys.sleep(wait)
#   }
#   stop(out)
# }
# bien_raw <- tryCatch(retry(BIEN_trait_species(
#   species = species_vec,
#   all.taxonomy = TRUE,
#   source.citation = TRUE
# )), error=function(e) NULL)
# 
# bien_clean <- if (!is.null(bien_raw)) {
#   bien_raw %>% as_tibble() %>%
#     clean_names() %>%
#     reframe(
#       species_bien = scrubbed_species_binomial,
#       author = scrubbed_taxon_name_with_author,
#       trait = tolower(trait_name),
#       value = as.character(trait_value),
#       unit = unit,
#       source = source_citation
#     ) %>%
#     filter(!is.na(species_bien), !is.na(trait), !is.na(value)) %>%
#     filter(trait %in% tolower(target_traits)) %>%
#     drop_na(species_bien)
# } else {
#   tibble(species_bien=character(), trait=character(), value=character())
# }
# 
# bien_wide <- bien_clean %>%
#   group_by(species_bien, trait) %>%
#   summarise(value = collapse_values(value), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = trait, values_from = value)
# 
# dispersal_bien <- bien_wide %>%
#   reframe(species=species_bien, dispersal_mode=`whole plant dispersal syndrome`) %>%
#   mutate(source="BIEN") %>%
#   distinct()

## B) EUDIS ---------------
eudis <- fread(IN_EUDIS) %>% 
  as_tibble() %>% clean_names()

syndrome_cols <- eudis %>% select(endozoochorous:ballochorous) %>% names()

eudis_with_combo <- eudis %>%
  mutate(across(all_of(syndrome_cols), ~ as.integer(replace_na(., 0L)))) %>%
  rowwise() %>%
  mutate(
    n_hits = sum(c_across(all_of(syndrome_cols))),                                
    .hits  = list(syndrome_cols[as.logical(c_across(all_of(syndrome_cols)))]),    
    combo_label = case_when(
      n_hits > 1  ~ paste(unlist(.hits), collapse = "_"),                         
      n_hits == 1 ~ unlist(.hits)[1],
      TRUE        ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(across(all_of(syndrome_cols), ~ if_else(n_hits > 1, 0L, .x))) %>%        
  select(-.hits)


long_single <- eudis_with_combo %>%
  pivot_longer(
    cols = all_of(syndrome_cols),
    names_to = "syndrome",
    values_to = "present"
  ) %>%
  filter(present == 1L) %>%
  transmute(
    species = species_flora_europaea, # or species_powo
    dispersal_mode = syndrome,
    source = "EuDiS"
  )

long_combo <- eudis_with_combo %>%
  filter(n_hits > 1) %>%
  transmute(
    species = species_flora_europaea,
    dispersal_mode = combo_label,
    source = "EuDiS"
  )

eudis_long <- bind_rows(long_single, long_combo) %>% distinct()

## C) TRY ---------------
#Kattge, J, Boenisch, G, Diaz, S, et al. TRY plant trait database -
# enhanced coverage and open access. Glob Change Biol. 2020; 26: 119-188. 
# https://doi.org/10.1111/gcb.14904
try_raw <- fread(IN_TRY) %>%
  as_tibble() %>%
  clean_names()

map_try_value_to_key <- function(v) {
  v <- stringr::str_to_lower(v)
  v <- stringr::str_replace_all(v, "[^a-z]", " ")
  v <- stringr::str_squish(v)
  dplyr::case_when(
    str_detect(v, "caching|hoarding")                                  ~ "vertebrate-hoarding",
    str_detect(v, "anemo|wind|meteo")                                  ~ "anemochorous",
    str_detect(v, "human|hemero|commerce|vehicle|speiro|ethelo")       ~ "human-mediated",
    str_detect(v, "ant|myrme|herpo|elaiosom")                          ~ "myrmecochorous",
    str_detect(v, "agochor|autochor|unspecialised|unassisted|sema|baro|gravity") ~ "unassisted",
    str_detect(v, "endozoo|eaten|intern|endo")                         ~ "endozoochorous",
    str_detect(v, "animal|epizoo|extern|feather|fur|dyso|zoochorous")         ~ "epizoochorous",
    str_detect(v, "thalasso|sea|ocean|current")                        ~ "thalassochorous",
    str_detect(v, "nauto|water|rain|fresh|hydrochor")                  ~ "freshwater-hydrochorous",
    str_detect(v, "boleochor|ballochor|ballistic|explos")              ~ "ballochorous",
    TRUE ~ NA_character_
  )
}

try_species <- 
  try_raw %>%
  reframe(species=str_to_sentence(acc_species_name), 
          data_name, 
          value=orig_value_str) %>%
  filter(str_detect(data_name,"[Dd]ispersal")) %>%
  mutate(key=map_try_value_to_key(value)) %>%
  filter(!is.na(key)) %>%
  distinct(species, key) %>%
  group_by(species) %>%
  summarise(dispersal_mode = paste(unique(key), collapse = "_"),
            .groups="drop") %>%
  mutate(source="TRY")

## D) Combine & check-----------
## Which species are still missing a dispersal syndrome? ---------------
species_modes_all <- bind_rows(
  #dispersal_bien %>% mutate(source="BIEN"),
  eudis_long, 
  try_species
) %>%
  filter(!is.na(dispersal_mode), dispersal_mode!="") %>%
  mutate(source = factor(source, levels=c("BIEN","EuDiS","TRY"))) %>%
  arrange(species, source) %>%
  group_by(species) %>% slice(1) %>% ungroup()

inat_species <- inat_obs %>% distinct(taxon_name) %>% rename(species=taxon_name)

still_missing_species <- inat_species %>%
  anti_join(species_modes_all %>% select(species), by="species")

# ______________________________ ------------

# 2) Interpolate GENUS dispersal syndromes ---------------------------
# (TRY-based GENUS majority vote)


threshold <- 0.33  # proportion threshold

try_genus_majority <- try_species %>%
  # keep a clean key column name
  rename(key = dispersal_mode) %>%
  # extract genus
  mutate(genus = stringr::word(species, 1)) %>%
  # expand "a_b_c" -> rows "a","b","c"
  tidyr::separate_rows(key, sep = "_") %>%
  # clean and keep valid tokens
  filter(!is.na(key), key != "") %>%
  # count occurrences per genus × mode
  count(genus, key, name = "n") %>%
  group_by(genus) %>%
  # share of each mode within the genus
  mutate(prop = n / sum(n)) %>%
  # keep all modes above threshold
  filter(prop >= threshold) %>%
  # order 
  arrange(genus, desc(n), key) %>%
  # collapse per genus: "mode1_mode2"
  summarise(
    genus_fill = paste0(key, collapse = "_"),
    .groups = "drop"
  )

## Combine & check ----------
species_modes_genus_fill <- still_missing_species %>%
  mutate(genus = word(species,1)) %>%
  left_join(try_genus_majority, by="genus") %>%
  transmute(
    species, dispersal_mode = genus_fill,
    source = if_else(is.na(genus_fill), NA_character_, "TRY_GenusFill"),
    genus
  ) %>%
  filter(!is.na(dispersal_mode))

species_modes_all <- bind_rows(
  species_modes_all,
  species_modes_genus_fill
) %>%
  distinct(species, .keep_all = TRUE)

#What's missing now?
still_missing_species <- inat_species %>%
  anti_join(species_modes_all %>% select(species), by="species")

# _________________ ----------------

# 3) Interpolate FAMILY dispersal syndromes ---------
# (TRY-based FAMILY majority vote)

species_family_lk <- inat_obs %>%
  transmute(species = taxon_name, family) %>%
  distinct()

try_family_majority <- try_species %>%
  mutate(key=dispersal_mode) %>%
  filter(!is.na(key)) %>%
  distinct(species, key) %>%
  left_join(species_family_lk, by="species") %>%
  filter(!is.na(family)) %>%
  count(family, key, name = "n") %>%
  group_by(family) %>%
  # share of each mode within the genus
  mutate(prop = n / sum(n)) %>%
  # keep all modes above threshold
  filter(prop >= threshold) %>%
  # order for deterministic concatenation
  arrange(family, desc(n), key) %>%
  # collapse per family: "mode1_mode2"
  summarise(
    family_fill = paste0(key, collapse = "_"),
    .groups = "drop"
  )

species_modes_family_fill <- still_missing_species %>%
  left_join(species_family_lk, by="species") %>%
  left_join(try_family_majority, by="family") %>%
  transmute(
    species, dispersal_mode = family_fill,
    source = if_else(is.na(family_fill), NA_character_, "TRY_FamilyFill"),
    family
  ) %>%
  filter(!is.na(dispersal_mode))

## Combine & check ---------
species_modes_final <- bind_rows(
  species_modes_all,
  species_modes_family_fill
) %>%
  distinct(species, .keep_all = TRUE)

species_still_na <- inat_species %>%
  anti_join(species_modes_final %>% select(species), by="species")

message("Species still missing after TRY genus+family fills: ", nrow(species_still_na))

# ______________________________ -------

# 4) Normalise dispersal modes ---------------------------------------
#
# Order of tokens (within-combo sorting)
order_ref <- c(
  "endozoochorous","epizoochorous","thalassochorous","anemochorous",
  "myrmecochorous","vertebrate-hoarding","freshwater-hydrochorous","ballochorous",
  "unassisted","human-mediated"
)

# Detect tokens present in a raw string and return canonical hyphen-preserving keys
detect_mode_keys <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  s <- tolower(x)
  
  keys <- character(0)
  if (str_detect(s, "endo"))           keys <- c(keys, "endozoochorous")
  if (str_detect(s, "epizoo"))      keys <- c(keys, "epizoochorous")
  if (str_detect(s, "thalasso"))          keys <- c(keys, "thalassochorous")
  if (str_detect(s, "anemo"))                    keys <- c(keys, "anemochorous")
  if (str_detect(s, "myrme"))            keys <- c(keys, "myrmecochorous")
  if (str_detect(s, "hoard"))                          keys <- c(keys, "vertebrate-hoarding")
  if (str_detect(s, "fresh"))    keys <- c(keys, "freshwater-hydrochorous")
  if (str_detect(s, "ballochor"))keys <- c(keys, "ballochorous")
  if (str_detect(s, "unassi"))    keys <- c(keys, "unassisted")
  if (str_detect(s, "human"))
    keys <- c(keys, "human-mediated")
  
  unique(keys)
}

# Canonicalise a single 'dispersal_mode' into:
# - mode_key  : underscore-joined keys, with hyphens inside specific categories preserved
# - mode_label: Sentence case, space-separated (hyphens preserved inside categories)
# - n_modes   : number of distinct tokens
canonise_mode <- function(dispersal_mode) {
  keys <- detect_mode_keys(dispersal_mode)
  keys <- intersect(order_ref, keys)  # sort & keep known
  n <- length(keys)
  if (n == 0) {
    return(list(mode_key = NA_character_, mode_label = NA_character_, n_modes = 0L))
  }
  if (n >= 2) {
    return(list(mode_key = paste(keys, collapse = "_"),
                mode_label = "Mixed",
                n_modes = as.integer(n)))
  }
  
  mode_key   <- paste(keys, collapse = "_")
  mode_label <- mode_key %>%
    str_replace_all("_", " ") %>%
    str_to_sentence() %>%
    # keep hyphenated words as-is
    identity()
  
  list(mode_key = mode_key, mode_label = mode_label, n_modes = as.integer(n))
}

canonised <- species_modes_final$dispersal_mode %>%
  lapply(canonise_mode)

species_modes_final <- species_modes_final %>%
  mutate(
    mode_key   = vapply(canonised, `[[`, character(1), "mode_key"),
    mode_label = vapply(canonised, `[[`, character(1), "mode_label"),
    n_modes    = vapply(canonised, `[[`, integer(1),   "n_modes")
  )



# Named mapping from canonical mode tokens -> agent
mode_to_agent <- c(
  "endozoochorous"          = "animal",
  "epizoochorous"           = "animal",
  "myrmecochorous"          = "animal",
  "vertebrate-hoarding"     = "animal",
  "anemochorous"            = "wind",
  "freshwater-hydrochorous" = "water",
  "thalassochorous"         = "water",
  "ballochorous"            = "ballistic",  # or "unassisted" 
  "unassisted"              = "unassisted",
  "human-mediated"          = "animal"
)

# helper: handles only a single scalar mode_key
mode_key_to_agents <- function(mode_key) {
  # Guarantee scalar behavior
  if (length(mode_key) != 1L) {
    # Pick first non-NA element if a vector is accidentally passed
    mode_key <- mode_key[1]
  }
  
  if (isTRUE(is.na(mode_key)) || identical(mode_key, "")) return(character(0))
  
  keys <- unlist(strsplit(mode_key, "_", fixed = TRUE), use.names = FALSE)
  agents <- unname(mode_to_agent[keys])
  unique(agents[!is.na(agents)])
}
agent_order <- c("animal","wind","water","ballistic","unassisted") 


# Canonical order for stable concatenation
agent_order    <- c("human","animal","wind","water","ballistic","unassisted")

species_modes_final <- species_modes_final %>%
  # Make sure mode_key is character, not factor/list
  mutate(mode_key = as.character(mode_key)) %>%
  rowwise() %>%
  mutate(
    agents_vec = list(mode_key_to_agents(mode_key))
  ) %>%
  ungroup() %>%
  # Deterministic ordering of agents
  mutate(
    agents_vec = purrr::map(agents_vec, ~ .x[order(match(.x, agent_order), na.last = TRUE)]),
    dispersal_agent_multi = ifelse(
      lengths(agents_vec) == 0L,
      NA_character_,
      purrr::map_chr(agents_vec, ~ paste(.x, collapse = "_"))
    )
  )

species_modes_final <- species_modes_final %>%
  mutate(
    dispersal_agent_summary = case_when(
      is.na(dispersal_agent_multi) ~ NA_character_,
      str_detect(dispersal_agent_multi, "_") ~ "mixed",
      TRUE ~ dispersal_agent_multi # from 1c
    )
  )


# ______________-----------

# 5) merge with inat.obs and inat.field --------------------------------
inat_obs_with_disp <- inat_obs %>%
  mutate(species = taxon_name) %>%
  left_join(species_modes_final %>% select(species, mode_key, mode_label, n_modes, source,
                                           dispersal_agent_summary), by = "species")

inat_field_with_disp <- inat_field %>%
  mutate(species = taxon_name) %>%
  left_join(species_modes_final %>% select(species, mode_key, mode_label, n_modes, source,
                                           dispersal_agent_summary), by = "species")
# _______________ --------

# 6) Write outputs ----
write_csv(species_modes_final, file.path(OUT_DIR_DISP, "species_dispersal_modes.csv"))
write_csv(species_still_na,   file.path(OUT_DIR_DISP, "species_missing_after_fill.csv"))
write_csv(inat_obs_with_disp, file.path(OUT_DIR_DISP, "inat_obs_with_dispersal.csv"))
write_csv(inat_field_with_disp, file.path(OUT_DIR_DISP, "inat_field_with_dispersal.csv"))

# _______________---------

# 7) Quick Summaries and plots----
sum_by_mode <- inat_obs_with_disp %>%
  mutate(mode_disp = coalesce(mode_label, "Unknown"),
         dispersal_agent =  coalesce(dispersal_agent_summary, "unknown")) %>%
  group_by(dispersal_agent) %>%
  summarise(n_species = n_distinct(taxon_name),
            n_obs     = sum(epi_count, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_species), desc(n_obs))

print(sum_by_mode, n = 60)

# inat_obs_with_disp %>% 
#   filter(!taxon_name%in%c("Magnoliopsida","Plantae","Angiospermae"),
#          is.na(mode_key)) %>% view()

# Example plot
ggplot(sum_by_mode, aes(x = fct_reorder(dispersal_agent, n_species), y = n_species)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  labs(x = NULL, y = "Number of species") +
  my_theme14


# quick plot
p_disp_sp <-
  ggplot(sum_by_mode, aes(x = fct_reorder(dispersal_agent, n_species), y = n_species)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  labs(x = NULL, y = "Number of species") +
  my_theme14

p_disp_obs <-
ggplot(sum_by_mode, aes(x = fct_reorder(dispersal_agent, n_obs), y = n_obs)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  labs(x = NULL, y = "Number of observations") +
  my_theme14

p_disp_obs + p_disp_sp

