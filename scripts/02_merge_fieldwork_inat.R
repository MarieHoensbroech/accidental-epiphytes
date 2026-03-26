# 03_merge_fieldwork_inat.R
# 
# 1) Fix some site ids (observation groups) that I messed up
# 2) Standardize IDs (site, TreeID, EpiID), 
# 3) Merge iNaturalist observations (Step 1) with fieldwork site & tree data,
# 4) Summarize counts per site & DBH class
# 5) Add site area from shapefile for standardisation/effort-adjustment
# 6) Save it all

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tidyverse)
  library(sf)
})

# ------------------------------- Parameters -----------------------------------
INAT_CSV   <- "data/processed/inat_observations.csv"
FIELD_DIR       <- "data/raw/fieldwork"
SITE_INFO_CSV   <- file.path(FIELD_DIR, "site_info.csv")
TREE_COUNT_CSV  <- file.path(FIELD_DIR, "tree_count.csv")
WILLOW_INFO_CSV <- file.path(FIELD_DIR, "willow_info.csv")  
ID_MAP_CSV      <- file.path(FIELD_DIR, "inat_epiid_map.csv") 

USER_FILTER <- "marie-ho"  # or NULL

OUT_DIR         <- "data/processed"
OUT_SITES_CSV   <- file.path(OUT_DIR, "my.sites.csv")
OUT_MERGED_CSV  <- file.path(OUT_DIR, "inat.merged.csv")

# site polygons -> area_ha 
SITES_SHP <- "data/raw/fieldsitesGIS/polygons/fieldsites.shp"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read iNat ------------------------------------
inat_obs <- fread(INAT_CSV) %>% 
  as_tibble()  %>% 
  mutate(id = as.character(id))

# . ------------------
# 1) Manual ID mapping ------------------------
#(I uploaded the first few observations without an ObsGroup; I put them in the id_map CSV:
# corresponding iNat obs IDs & EpiIDs from field notes,
# needs to be read here to backfill missing/incorrect ObsGroup/EpiID.
id_map <- fread(ID_MAP_CSV) %>% 
  as_tibble() %>% 
  reframe(id = as.character(id), EpiID = as.character(EpiID))

# 2.1) Standardise iNat data--------------------------
## 2.1.1) site/Tree/EpiID ---------
# ObsGroup should look like S{01,02}T{01,02} (e.g., "S03T01", "S12T07").
# Then EpiID can be generated per observation

# helper 
normalize_site_tree_df <- function(x) {
  site_raw <- stringr::str_extract(x, "S\\d{1,2}")
  tree_raw <- stringr::str_extract(x, "T\\d{1,2}")
  site <- ifelse(!is.na(site_raw),
                 paste0("S", stringr::str_pad(stringr::str_extract(site_raw, "\\d+"), 2, pad = "0")),
                 NA)
  tree <- ifelse(!is.na(tree_raw),
                 paste0("T", stringr::str_pad(stringr::str_extract(tree_raw, "\\d+"), 2, pad = "0")),
                 NA)
  tibble::tibble(
    site   = site,
    Tree   = tree,
    TreeID = ifelse(is.na(site) | is.na(tree), NA, paste0(site, tree))
  )
}


# Make the mapping table 
id_map2 <- id_map %>%
  mutate(id = as.character(id)) %>%
  rename(EpiID_map = EpiID) %>%
  distinct(id, .keep_all = TRUE)

# Prepare inat_for_merge up to BaseGroup 
inat_for_merge <- inat_obs %>%
  select(id,ObsGroup,DBH_num,Setting,EpiCount_noNA,taxon_id) %>% 
  mutate(id = as.character(id)) %>%
  left_join(id_map2, by = "id") %>%
  { if ("EpiID" %in% names(.)) mutate(., EpiID = as.character(coalesce(EpiID, EpiID_map)))
    else rename(., EpiID = EpiID_map) } %>%
  mutate(EpiID = as.character(EpiID)) %>%
  select(-any_of("EpiID_map")) %>%
  # Corrections
  mutate(
    ObsGroup = if_else(
      id != 312915038 & !is.na(ObsGroup) & str_detect(ObsGroup, "S32"),
      str_replace(ObsGroup, "S32", "S33"),
      ObsGroup
    ),
    ObsGroup = if_else(
      str_detect(ObsGroup, "S5") & Setting == "Urban",
      str_replace(ObsGroup, "S5", "S6"),
      ObsGroup
    )
  ) %>%
  mutate(BaseGroup = coalesce(str_extract(EpiID, "S\\d{1,2}T\\d{1,2}"), ObsGroup)) %>%
  filter(!is.na(BaseGroup) & str_detect(BaseGroup, "S\\d{1,2}T\\d{1,2}"))

# Normalize 
norm <- normalize_site_tree_df(inat_for_merge$BaseGroup)

# Bind normalized columns
inat_for_merge <- bind_cols(inat_for_merge, norm)

# Create standardized numerics and build padded TreeID
inat_for_merge <- inat_for_merge %>%
  mutate(
    SiteNum = as.integer(str_extract(site, "\\d+")),
    TreeNum = as.integer(str_extract(Tree, "\\d+")),
    TreeID  = sprintf("S%02dT%02d", SiteNum, TreeNum)) %>% 
  group_by(TreeID) %>% 
  arrange(TreeID) %>% 
  mutate(EpiNum    = sprintf("%02d", row_number()),
         EpiID = paste0(TreeID, "_E", EpiNum)) %>% 
  reframe(id,TreeID,EpiID,site,DBH_num,EpiCount_noNA,taxon_id)


# . ----------------
## 2.1.2 DBH categories & per-site sums ----------------------
dbh_breaks  <- c(0, 0.10, 0.30, 0.60, 0.90, Inf)
dbh_labels  <- c("0 - 10", "10 - 30","30 - 60","60 - 90","90+")
inat_for_merge <- inat_for_merge  %>% 
  mutate(
    TreeCat = cut(DBH_num, breaks = dbh_breaks, labels = dbh_labels, right = FALSE)
  )

# Per-site & per-size summaries (epiphyte counts and species per site)
Epi_sum <- inat_for_merge %>% 
  select(id, TreeCat, EpiCount_noNA, taxon_id, site)  %>% 
  distinct()  %>% 
  group_by(site)  %>% 
  mutate(
    n_epis_per_site = sum(EpiCount_noNA, na.rm = TRUE),
    n_sp_per_site   = n_distinct(taxon_id)
  )  %>% 
  group_by(site, TreeCat)  %>% 
  reframe(
    n_epis_per_site = first(n_epis_per_site),
    n_sp_per_site   = first(n_sp_per_site),
    n_epis_per_site_cat = sum(EpiCount_noNA, na.rm = TRUE),
    n_sp_per_site_cat   = n_distinct(taxon_id)
  )


# . ------------
# 2.2) Standardise my data-----------
## 2.2.1) Load fieldwork tables -------------------------
site_info  <- fread(SITE_INFO_CSV) %>% as_tibble()
tree_count <- fread(TREE_COUNT_CSV) %>% as_tibble() %>% rename(site=Site)
willow_info<- fread(WILLOW_INFO_CSV) %>% as_tibble()

## 2.2.2) Normalize site code in all field tables to "Sxx"
normalize_site_code <- function(s) {
  s <- str_replace(s, "site", "S")
  s <- paste0("S", str_pad(str_extract(s, "\\d+"), 2, pad = "0"))
  s
}

site_info <- site_info  %>% 
  select(site,lat,lon,setting,dominanttree1:undergrowth) %>% 
  mutate(site = normalize_site_code(site)) 


willow_info <- willow_info  %>% 
  mutate(site = normalize_site_code(site),
         TreeCat = cut(dbh,
                       breaks = c(0, 0.10, 0.30, 0.60, 0.90, Inf),
                       labels = c("0 - 10","10 - 30","30 - 60","60 - 90","90+"),
                       right = FALSE))  %>% 
  group_by(site, TreeCat)  %>% 
  reframe(NumberTrees = n_distinct(treeid)) %>% unique()


# Tree count table (site-level inventory). 
tree_count <- tree_count  %>% 
  mutate(site = normalize_site_code(site))  %>% 
  select(-any_of(c("NumberEpis"))) %>% 
  drop_na(NumberTrees) %>% 
  mutate(TreeCat=str_replace(TreeCat,"--"," - ")) 

# Combine field tree counts with willow-specific counts if provided
tree_inventory <- tree_count %>%
  bind_rows(willow_info %>% filter(site!="S15")) %>% 
  arrange(site)

# 2.3) Build my.sites --------------------------------
my_sites <- tree_inventory  %>% 
  right_join(site_info)  %>% 
  left_join(Epi_sum %>% rename(site=site))  %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
  mutate(TreeCat=str_replace(TreeCat,"--"," - ")) %>% 
  rename(siteLat=lat,
         siteLon=lon) 


# . -----------

# 3) Compute area_ha from polygons -------------------------------------
polys <- st_read(SITES_SHP, quiet = TRUE)

if (is.na(st_crs(polys))) {
  warning("site polygons have no CRS; assuming EPSG:4326.")
  st_crs(polys) <- 4326
}

polys <- polys %>%
  mutate(
    area_ha = dplyr::case_when(
      "Area_km2" %in% names(polys) & is.finite(Area_km2) ~ as.numeric(Area_km2) * 100,
      TRUE ~ as.numeric(st_area(st_transform(polys, 6933))) / 1e4
    )
  ) %>%
  st_transform(4326) %>%
  st_drop_geometry() %>%
  reframe(Name, area_ha) %>%
  mutate(
    site = normalize_site_code(Name)
  )

# Merge area_ha into my_sites
my_sites <- my_sites %>% left_join(polys) %>% unique()

# Write my.sites
fwrite(my_sites, OUT_SITES_CSV)

# 4) Merge to iNat ---------------------------------
inat_merged <- inat_for_merge  %>% 
  left_join(my_sites) %>% 
  left_join(inat_obs %>% distinct(id, .keep_all = T))

inat_merged %>% 
  group_by(id) %>% 
  filter(n()>1) %>% view

fwrite(inat_merged, OUT_MERGED_CSV)


