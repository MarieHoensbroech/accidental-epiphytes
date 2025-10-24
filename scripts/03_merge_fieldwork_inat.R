# 03_merge_fieldwork_inat.R
# 
# 1) Fix some site ids (observation groups) that I messed up
# 2) Standardize IDs (Site, TreeID, EpiID), 
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

# Site polygons -> area_ha 
SITES_SHP <- "data/Raw/fieldsites/fieldsites.shp"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read iNat ------------------------------------
inat_tbl <- fread(INAT_CSV) %>% 
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
## 2.1.1) Site/Tree/EpiID ---------
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
    Site   = site,
    Tree   = tree,
    TreeID = ifelse(is.na(site) | is.na(tree), NA, paste0(site, tree))
  )
}


# Make the mapping column
id_map2 <- id_map %>%
  mutate(id = as.integer(id)) %>%          
  rename(EpiID_map = EpiID) %>%
  distinct(id, .keep_all = TRUE)


inat_for_merge <- inat_tbl %>%
  mutate(id = as.integer(id)) %>%                          # align to map's id type
  left_join(id_map2, by = "id") %>%                   
  # Create a clean EpiID column
  { if ("EpiID" %in% names(.)) mutate(., EpiID = coalesce(EpiID, EpiID_map))
    else rename(.,  EpiID = EpiID_map) } %>%
  select(-any_of("EpiID_map")) %>%
  # (Optional) filter by user
  #{ if (!is.null(USER_FILTER)) filter(., user_login %in% USER_FILTER) else . } %>%
  # Correction: all sites with ObsGroup S32 need to be S33, except for the one observation at S32
  mutate(
    ObsGroup = if_else(
      id != 312915038 & !is.na(ObsGroup) & str_detect(ObsGroup, "S32"),
      str_replace(ObsGroup, "S32", "S33"),
      ObsGroup
    )
  ) %>%
  #Correction: S5 -> S6
  mutate(
    ObsGroup = if_else(
      str_detect(ObsGroup,"S5") & Setting=="Urban",
      str_replace(ObsGroup, "S5", "S6"),
      ObsGroup
    )
  ) %>%
  # Build a single base group from EpiID (strip suffix) or ObsGroup
  mutate(BaseGroup = coalesce(str_extract(EpiID, "S\\d{1,2}T\\d{1,2}"), ObsGroup)) %>%
  filter(!is.na(BaseGroup) & str_detect(BaseGroup, "S\\d{1,2}T\\d{1,2}")) %>%
  # Normalize with vectorized helper 
  bind_cols(normalize_site_tree_df(.$BaseGroup)) %>%
  mutate(Epi = str_extract(EpiID, "E\\d{1,3}")) %>%
  group_by(TreeID)


# Generate EpiID sequence per TreeID (E01, E02, …) and rename coords
inat_for_merge <- inat_for_merge  %>% 
  group_by(TreeID)  %>% 
  arrange(TreeID, id, .by_group = TRUE)  %>% 
  mutate(EpiNum = str_pad(row_number(), 2, pad = "0"),
         EpiID  = paste0(TreeID, "_E", EpiNum))  %>% 
  ungroup()  %>% 
  rename(Epi.lat = latitude, Epi.long = longitude) 



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
  select(id, TreeCat, EpiCount_noNA, taxon_id, Site)  %>% 
  distinct()  %>% 
  group_by(Site)  %>% 
  mutate(
    n_epis_per_site = sum(EpiCount_noNA, na.rm = TRUE),
    n_sp_per_site   = n_distinct(taxon_id)
  )  %>% 
  group_by(Site, TreeCat)  %>% 
  reframe(
    n_epis_per_site = first(n_epis_per_site),
    n_sp_per_site   = first(n_sp_per_site),
    n_epis_per_site_cat = sum(EpiCount_noNA, na.rm = TRUE),
    n_sp_per_site_cat   = n_distinct(taxon_id)
  )


# . ------------
# 2.2) Standardise my data-----------
## 2.2.1) Load fieldwork tables -------------------------
site_info  <- fread(SITE_INFO_CSV) 
tree_count <- fread(TREE_COUNT_CSV) %>% as_tibble()
willow_info<- fread(WILLOW_INFO_CSV) %>% as_tibble()

## 2.2.2) Normalize Site code in all field tables to "Sxx"
normalize_site_code <- function(s) {
  s <- str_replace(s, "Site", "S")
  s <- paste0("S", str_pad(str_extract(s, "\\d+"), 2, pad = "0"))
  s
}

site_info <- site_info  %>% 
  select(Site,lat,lon,dominanttrees) %>% 
  mutate(Site = normalize_site_code(Site)) 


willow_info <- willow_info  %>% 
  mutate(Site = normalize_site_code(site),
         TreeCat = cut(dbh,
                       breaks = c(0, 0.10, 0.30, 0.60, 0.90, Inf),
                       labels = c("0 - 10","10 - 30","30 - 60","60 - 90","90+"),
                       right = FALSE))  %>% 
  group_by(Site, TreeCat)  %>% 
  reframe(NumberTrees = n_distinct(treeid), Type = "Willow") %>% unique()


# Tree count table (site-level inventory). 
tree_count <- tree_count  %>% 
  mutate(Site = normalize_site_code(Site))  %>% 
  select(-any_of(c("NumberEpis"))) %>% 
  drop_na(NumberTrees) %>% 
  mutate(TreeCat=str_replace(TreeCat,"--"," - "))

# Combine field tree counts with willow-specific counts if provided
tree_inventory <- bind_rows(
  tree_count %||% tibble(),
  willow_info %||% tibble()
) %>% 
  arrange(Site) %>% 
  mutate(Type=if_else(Type=="Urban","Forest",Type))

# 2.3) Build my.sites --------------------------------
my_sites <- tree_inventory  %>% 
  full_join(site_info)  %>% 
  full_join(Epi_sum)  %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
  mutate(TreeCat=str_replace(TreeCat,"--"," - ")) %>% 
  rename(SiteLat=lat,
         SiteLon=lon) 


# . -----------

# 3) Compute area_ha from polygons -------------------------------------
polys <- st_read(SITES_SHP, quiet = TRUE)

if (is.na(st_crs(polys))) {
  warning("Site polygons have no CRS; assuming EPSG:4326.")
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
  select(Site, area_ha) %>%
  mutate(
    Site = paste0("S", str_pad(str_extract(Site, "(?<=^Site)\\d+"), 2, pad = "0", side = "left"))
  )

# Merge area_ha into my_sites
my_sites <- my_sites %>% left_join(polys, by = c("site" = "Site"))

# Write my.sites
fwrite(my_sites, OUT_SITES_CSV)

# 4) Merge to iNat ---------------------------------
inat_merged <- inat_for_merge  %>% 
  left_join(my_sites)

fwrite(inat_merged, OUT_MERGED_CSV)
