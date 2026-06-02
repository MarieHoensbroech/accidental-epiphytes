# 06a_environmental_variables_inat_presences.R

rm(list = ls())

# ----------------------------- 0) Packages ------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidyr)
  library(readr)
  library(purrr)
  library(forcats)
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(brms)
  library(posterior)
  library(data.table)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(lwgeom)
  library(archive)
  library(stars)
  library(ncdf4)
  library(tmap)
  library(tmaptools)
})

set.seed(42)  # Reproducibility


# -------------------------------- 0) Parameters --------------------------------
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"

DIR_ENV    <- "data/environmental"
DIR_CHELSA <- file.path(DIR_ENV, "CHELSA")

OUT_DIR <- "outputs/06_bayesian_model_relative_niche/environmental"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_FILE <- file.path(OUT_DIR, "inat_presences_with_4chelsa.csv")


# -------------------------------- 0) Load data --------------------------------
taxa<-fread(IN_INAT_FIELD) %>% 
  distinct(taxon_id) %>% deframe()

#inat_raw %>% distinct(taxon.rank)
inat_raw <- fread(IN_INAT_OBS) %>% 
  filter(taxon_id%in%taxa) %>% 
  distinct(id, .keep_all = T)

inat_raw <- inat_raw %>% 
  as_tibble() %>%
  filter(taxon_id%in%taxa,
         taxon.rank%in%c("species","complex","subspecies","hybrid"),
         quality_grade=="research"#| !is.na(community_taxon_id)
         ) %>% 
  mutate(
    taxon.id=taxon_id,
    lat = latitude, lon = longitude,
    lat_mid5 = floor(lat / 5) * 5 + 0.5
  ) %>% 
  drop_na(lat, lon) %>% 
  select(lat, lon, taxon.name, id, quality_grade,taxon.id) %>% 
  group_by(taxon.name) %>% 
  mutate(n=n()) %>% 
  filter(n>=20)


# --------------------------- 1) Clean & prepare data ---------------------------
inat_sf_global <- st_as_sf(inat_raw, coords = c("lon","lat"), crs = 4326, remove = FALSE) 
inat_coords_global <- inat_sf_global %>% 
  reframe(
    taxon.id,
    lat, lon,
    id = as.character(id),
    site = id,
    geometry,
    taxon.name,
    presence = 1
  ) %>% 
  distinct(id, .keep_all = TRUE)%>% 
  mutate(row_id__tmp = row_number())

pts_v <- terra::vect(inat_coords_global)

# 2) Find environmental rasters--------------
f_ai <- list.files(
  DIR_CHELSA,
  pattern = "chelsa.*ai.*1981.*2010.*v.*2.*1.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
) %>%
  .[1]

f_bio10 <- list.files(
  DIR_CHELSA,
  pattern = "chelsa.*bio10.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
) %>%
  .[1]

f_bio6 <- list.files(
  DIR_CHELSA,
  pattern = "chelsa.*bio6.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
) %>%
  .[1]

f_wind <- list.files(
  DIR_CHELSA,
  pattern = "chelsa.*sfc.*wind.*mean.*1981.*2010.*v.*2.*1.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
) %>%
  .[1]

message("Using CHELSA files:")
message("  AI   : ", f_ai)
message("  BIO10: ", f_bio10)
message("  BIO6 : ", f_bio6)
message("  WIND : ", f_wind)

# ----------------------------- 3) Load rasters --------------------------------

r_ai    <- terra::rast(f_ai)
r_bio10 <- terra::rast(f_bio10)
r_bio6  <- terra::rast(f_bio6)
r_wind  <- terra::rast(f_wind)

# ----------------------------- 4) Extract values ------------------------------

out_extract <- tibble(
  row_id__tmp = inat_coords_global$row_id__tmp,
  chelsa_ai_1981_2010_v_2_1 = terra::extract(r_ai, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_bio10 = terra::extract(r_bio10, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_bio6 = terra::extract(r_bio6, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_sfc_wind_mean_1981_2010_v_2_1 = terra::extract(r_wind, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric()
)

# ----------------------------- 5) Join back -----------------------------------
out <- inat_coords_global %>%
  left_join(out_extract, by = "row_id__tmp") %>%
  mutate(
    moisture = as.numeric(chelsa_ai_1981_2010_v_2_1),
    max_temp = as.numeric(chelsa_bio10),
    min_temp = as.numeric(chelsa_bio6),
    wind     = as.numeric(chelsa_sfc_wind_mean_1981_2010_v_2_1)
  ) %>%
  select(-row_id__tmp)

# ----------------------------- 6) Save ----------------------------------------

write_csv(out, OUT_FILE)


