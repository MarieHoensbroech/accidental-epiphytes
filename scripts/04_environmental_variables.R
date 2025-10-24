# 04_environmental_variables.R
#
#  
# 1) generate background points 
# 2) Combine with presence points
# 3) Extract environmental variables
# 4) Save to file

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(data.table)
  library(sf)
  library(terra)
  library(janitor)
  library(lubridate)
  library(stringr)
  library(glue)
})


# 0) Load presence points ------
my_sites_sf <- fread("data/processed/my.sites.csv") %>% 
  mutate(lon = as.numeric(SiteLon), lat = as.numeric(SiteLat)) %>%
  filter(is.finite(lon), is.finite(lat), !(lat == 0 & lon == 0)) %>% 
  select(site=Site, lat, lon) %>% 
  st_as_sf(coords=c("lon","lat"), crs=4326, remove = F)


inat_obs_sf <- fread("data/processed/inat_observations.csv") %>% as_tibble() %>% unique() %>% 
  drop_na(latitude, longitude) %>%
  mutate(site = sprintf("N%05d", dplyr::row_number())) %>%
  reframe(site, lat = latitude, lon = longitude) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)

pts_all <- bind_rows(my_sites_sf, inat_obs_sf)






# 1) Generate global background points ------------
set.seed(42)
n_points <- 10000

# Define global bounding box (adjust if needed)
global_bbox <- list(xmin = -180, xmax = 180, ymin = -90, ymax = 90)

generate_global_points <- function(n, bbox, presence_sf) {
  # Create random coordinates
  coords <- data.frame(
    lon = runif(n, bbox$xmin, bbox$xmax),
    lat = runif(n, bbox$ymin, bbox$ymax)
  )
  
  # Remove duplicates with presence points
  presence_coords <- st_coordinates(presence_sf)
  coords <- coords[!paste(coords$lon, coords$lat) %in% paste(presence_coords[,1], presence_coords[,2]), ]
  
  # Convert to sf object
  st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
}

bg_sf <- generate_global_points(n_points, global_bbox, pts_all)

# Create 10 km buffers around presence and background points ---
buffer_presence <- pts_all %>% st_transform(3857) %>% st_buffer(10000) %>% st_transform(4326)
buffer_background <- bg_sf %>% st_transform(3857) %>% st_buffer(10000) %>% st_transform(4326)

# Save background points ---
bg_dt <- data.table(
  site = paste0("bg", seq_len(nrow(bg_sf))),
  lon = st_coordinates(bg_sf)[,1],
  lat = st_coordinates(bg_sf)[,2]
)

fwrite(bg_dt, "data/processed/background_points.csv")

# 2) Combine with presence points ----------

bg_sf <- bg_dt %>% as_tibble() %>% st_as_sf(coords=c("lon","lat"), crs = 4326, remove=F)

rbind(bg_sf, pts_all) %>% unique() -> pts


# 3) Extract environmental variables ------------
## Prep -------------
# Environmental data root and files
DIR_ENV    <- "data/environmental"
PATH_ACCESS<- file.path(DIR_ENV, "traveltime/accessibility_to_cities_2015_v1.0.tif")
PATH_TD    <- file.path(DIR_ENV, "tree.density/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif")
PATH_FROST <- file.path(DIR_ENV, "frost.days/cru_ts4.08.1901.2023.frs.dat.nc")
PATH_AGE   <- file.path(DIR_ENV, "forest.age/2025417153414222_BGIForestAgeMPIBGC1.0.0.nc")
DIR_CANOPY <- file.path(DIR_ENV, "canopy.height/chm")
DIR_CHELSA <- file.path(DIR_ENV, "CHELSA/bio")
DIR_LC     <- file.path(DIR_ENV, "landcover")
PATH_WB    <- file.path(DIR_ENV, "woody.biomass.density/Aboveground_Live_Woody_Biomass_Density.geojson")

PATH_CHELSA <- list.files(
  DIR_CHELSA,
  pattern = glob2rx("CHELSA_*.tif"),
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# EarthEnv consensus classes
CONS_CLASSES <- c(1,2,3,4,5,6,7,8,9)


## Extract -------
occ_covs_tbl <- pts %>% st_drop_geometry() %>% select(site, lon, lat) %>% 
  unique()

### 3a) CHELSA ---------------
for (f in PATH_CHELSA) {
  r  <- terra::rast(f)
  nm <- stringr::str_extract(basename(f), "^CHELSA_bio\\d+")
  if (is.na(nm)) nm <- tools::file_path_sans_ext(basename(f))
  vals <- terra::extract(r, pts)
  occ_covs_tbl[[nm]] <- as.numeric(vals[,2])
}

### 3b) Accessibility ----------
r <- terra::rast(PATH_ACCESS)
occ_covs_tbl$access_min <- as.numeric(terra::extract(r, pts)[,2])

### 3c) EarthEnv land cover -----------
lc_files <- file.path(DIR_LC, sprintf("consensus_full_class_%d.tif", CONS_CLASSES))

for (k in seq_along(lc_files)) {
  r  <- terra::rast(lc_files[k])
  nm <- sprintf("clc_%d", CONS_CLASSES[k])
  occ_covs_tbl[[nm]] <- as.numeric(terra::extract(r, pts)[,2])
}

### 3d) Tree density -----------
r <- terra::rast(PATH_TD)
occ_covs_tbl$tree_density <- as.numeric(terra::extract(r, pts)[,2])

### 3e) Frost days (average across bands) --------
r <- terra::rast(PATH_FROST)
occ_covs_tbl$frost_days <- as.numeric(terra::extract(r, pts, fun = mean, na.rm = TRUE)[,2])

### 3f) Canopy height (VRT mosaic) ----------
chm_files <- list.files(DIR_CANOPY, pattern = "\\.tif$", full.names = TRUE)
chm_vrt <- tryCatch(terra::vrt(chm_files), error = function(e) NULL)
occ_covs_tbl$canopy_h_mean <- as.numeric(terra::extract(chm_vrt, pts)[,2])

### 3g) Forest age (BGI) ------------
r <- tryCatch(
  terra::rast(PATH_AGE, subds = "ForestAge_TC000"),
  error = function(e) terra::rast(PATH_AGE)
)
occ_covs_tbl$forest_age <- as.numeric(terra::extract(r, pts)[,2])





#4) Write to csv -----------

write_csv(occ_covs_tbl, "data/processed/10_occ_covs.csv")
