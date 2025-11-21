# 04_environmental_variables_prep.R 
#
# 1) Prep data (filter to Europe, create background points)
# 2) Extract
# 3) find out which environmental variables are most descriptive of accidental
#    epiphyte presence and are not collinear

rm(list=ls())

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
  library(terra)
  library(stars)
  library(ncdf4)
  library(tmap)
  library(tmaptools)
})

set.seed(42)  # Reproducibility for sampling background points, etc.


# 0) Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"

OUT_DIR_ENV  <- "outputs/04_environmental"
dir.create(OUT_DIR_ENV, showWarnings = FALSE, recursive = TRUE)

# restrict iNat points to within 100 km of the European coast?
use_coastal <- FALSE

# How many background points to sample for presence-only integration?
n_bg <- 10000

# Natural Earth land polygons (10m scale ~ detailed)
land <- rnaturalearth::ne_download(
  scale = 10,
  type  = "land",
  category = "physical",
  returnclass = "sf"
)

# Work in an equal-area projection for area calculation
# EPSG:6933 = NSIDC EASE-Grid 2.0 Global equal-area
land_eq <- st_transform(land, 6933) %>%
  st_cast("POLYGON") %>%
  mutate(poly_id = dplyr::row_number(),
         area_m2 = as.numeric(st_area(geometry)))

# 0) Load data --------
my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% 
  reframe(site, dbh_cat=TreeCat,
          setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon) %>% 
  drop_na(lat,lon)


inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% #glimpse()
  reframe(site, dbh_cat=TreeCat,
          setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          count=EpiCount_noNA,
          tree_id = TreeID,
          quality_grade) %>% 
  drop_na(lat,lon) %>% 
  unique()

inat_raw   <- fread(IN_INAT_OBS) %>% as_tibble() %>% #glimpse()
  filter(!id%in%c(inat_field$id)) %>% 
  mutate(lat=latitude, lon=longitude,
         lat_mid5 = floor(lat / 5) * 5 + 0.5) %>% 
  drop_na(lat,lon) %>% 
  select(lat,lon,taxon.name,id,quality_grade)





# 1) Clean & prepare data---------------------------
## 1A) Field data ------------
my_obs_coords <- 
  my_sites %>% 
  unique() %>% 
  full_join(inat_field) %>% 
  reframe(site,lat,lon,
          id=if_else(is.na(id),"no_id",as.character(id)),
          taxon.name=if_else(is.na(taxon.name),"no_taxon",taxon.name),
          quality_grade=if_else(is.na(quality_grade),"no_taxon",quality_grade),
          presence=if_else(taxon.name=="no_taxon",0,1)) %>% 
  distinct(id,site, .keep_all = T) %>% 
  st_as_sf(coords=c("lon","lat"),crs=4326, remove = F) 

## 1B) iNaturalist---------------
# Filter to European Atlantic coast 
# Choose a projected CRS for Europe (meters) for buffering
target_crs <- 3035  # ETRS89 / LAEA Europe

# Define a WGS84 bbox that covers the Atlantic-facing W Europe
bbox_ll <- st_as_sfc(
  st_bbox(c(xmin = -10, ymin = 35, xmax = 15, ymax = 72), crs = st_crs(4326)),
  crs = st_crs(4326)
)
# Transform bbox to target CRS
bbox_eu <- st_transform(bbox_ll, target_crs)

# load & prepare layers
land <- rnaturalearth::ne_download(scale = 10, type = "land", 
                                   category = "physical", returnclass = "sf") %>%
  st_make_valid() 

land_eu<- land %>% 
  st_transform(target_crs) %>%
  st_intersection(bbox_eu)


# points -> sf and project
inat_sf <- st_as_sf(inat_raw, coords = c("lon","lat"), crs = 4326, remove = FALSE) %>%
  st_transform(target_crs)

# Keep only the points inside bbox
inat_sf_atl <- st_crop(inat_sf, st_bbox(bbox_eu))
inat_coords <- inat_sf_atl %>% 
  st_transform(crs=4326) %>% 
  reframe(lat,lon,
          id=as.character(id),
          site=id,
          geometry,
          quality_grade,taxon.name,
          presence=1) %>% 
  distinct(id, .keep_all = T)

## 1C) Create EUROPE BACKGROUND points --------------------------
# Draw random points from the same region as filtered iNat,

# Load land polygons and bring them to 3035
land_10m <- rnaturalearth::ne_download(
  scale = 10, type = "land", category = "physical", returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_transform(3035)

# Restrict to bbox (keeps only land inside bbox)
land_clip <- st_intersection(land_10m, bbox_eu) %>%
  st_collection_extract("POLYGON") %>%
  filter(st_is_valid(geometry)) %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  filter(area_m2 > 1e5) %>%  # keep polygons > 0.1 km²; adjust as needed
  select(-area_m2)

# Sample points uniformly over land polygons
# exact = TRUE ensures exactly n_bg points overall (sf >= 1.0)
set.seed(42)
bg_pts_sfc <- st_sample(land_clip, size = n_bg, type = "random", exact = TRUE)

# Convert to sf and attach CRS 
bg_pts_sf <- st_as_sf(bg_pts_sfc) %>%
  st_set_crs(st_crs(land_clip)) %>%
  st_transform(4326) 

bg_pts_sf<-
  bg_pts_sf %>% 
  mutate(site = sprintf("B%05d", dplyr::row_number()),
         id = "no_id",
         taxon.name="no_taxon",
         quality_grade="no_grade",
         presence=0) %>% 
  rename(geometry = x) %>% 
  mutate(
    lon = st_coordinates(st_geometry(.))[, 1],
    lat = st_coordinates(st_geometry(.))[, 2]
  )

## 1D) JOIN All points for env extraction-------------

pts <- my_obs_coords %>% 
  bind_rows(inat_coords) %>% 
  bind_rows(bg_pts_sf) %>% 
  st_crop(bbox_ll) 

#save
pts %>% 
  st_drop_geometry() %>% 
  fwrite(file.path(OUT_DIR_ENV,"EUROPE_field+inat+bg_points.csv"))


## 1D) Quick plot to check extent------------------
countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
ggplot() +
  geom_sf(data = countries, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = pts, 
          size = 1.6, alpha = 0.8, show.legend = TRUE) +
  #scale_color_manual(values = c(island = "lightblue", mainland = "#d95f02"), na.translate = FALSE) +
  coord_sf(expand = FALSE)  +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    legend.position = "right",
    axis.title = element_blank()
  )

## 1E) Global LAI Merge (Run once)--------------------------
#DATA SOURCE: https://zenodo.org/records/14709655
nc_path <- "data/environmental/c_gls_LAI300-RT0_202507310000_GLOBE_OLCI_V1.1.1_nc/c_gls_LAI300-RT0_202507310000_GLOBE_OLCI_V1.1.1_nc/c_gls_LAI300-RT0_202507310000_GLOBE_OLCI_V1.1.1.nc"

# Quick open with terra 
r <- rast(nc_path)
r
nlyr(r)               # number of layers
names(r)              # layer names

# Read main variable as SpatRaster (first layer if single-var file)
lai_global <- rast(nc_path, subds = "LAI")  # if 'varname' is known; otherwise 'r' is already the layer

setMinMax(lai_global,force=TRUE)
# library(leaflet)
# library(terra)
# library(viridisLite)
# 
# # Example extent (Hondarribia), crop first
# hond <- ext(-1.9, -1.7, 43.34, 43.38)
# lai_hond <- crop(lai, hond)
# 
# pal <- colorNumeric(viridis(256, option = "G"), domain = c(0, 8), na.color = "transparent")
# 
# leaflet()  %>% 
#   addProviderTiles(providers$CartoDB.Positron) %>%
#   addRasterImage(lai, colors = pal, opacity = 0.6, project = TRUE) %>%
#   addLegend(pal = pal, values = c(0,8), title = "LAI (m²/m²)") %>%
#   addScaleBar()






### 1E.1) URBAN -------------------
#DATA SOURCE: https://zenodo.org/records/14709655

# zip_path <- "data/environmental/urban.tree.lai/global_15s.tar.gz"
# 
# archive::archive(zip_path)
# 
# 
# # Extract into a temp or project folder
# out_dir <- "data/environmental/urban.tree.lai"
# dir.create(out_dir, showWarnings = FALSE)
# archive::archive_extract(zip_path, dir = out_dir)
# 
# nc_files <- list.files(out_dir, pattern = "\\.nc$", full.names = TRUE)
# nc_files
# 
# # via terra (terra reads NetCDF too)
# r <- terra::rast(nc_files[23])
# r
# 
# 
# # Rename bands 
# names(r)[1:12]  <- paste0("LAI", sprintf("%02d", 1:12))  # Jan..Dec
# names(r)[13:24] <- paste0("SAI"  , sprintf("%02d", 1:12))  # Jan..Dec
# 
# # Create two stacks
# lai_urban <- r[[1:12]]
# sai   <- r[[13:24]]
# 
# 
# # Raster has 12 bands (Jan–Dec), I pick one month to plot: July
# lai_jul <- r[[7]]                # 7 = July
# names(lai_jul) <- "LAI (July)"
# setMinMax(lai_jul)
# 
# # Simple plot with a perceptually uniform palette
# pal <- hcl.colors(20, "YlGn")    # good for vegetation
# plot(lai_jul, col = pal, main = "Urban Tree LAI (July)")
# 
#  vals <- terra::spatSample(lai_jul, size = 2000, method = "random", na.rm = FALSE, values = TRUE)
# summary(vals)
# 
# pal <- hcl.colors(32, "YlGn")
# plot(lai_jul, col = pal, main = "Urban Tree LAI (July)",
#      range = c(0, 8),   # adjust if min/max suggest a different span
#      mar = c(3,3,2,6))  # wider right margin for legend
# 
# 
# paris <- ext(2.0, 2.6, 48.7, 49.1)
# plot(crop(lai_jul, paris), col = pal, range = c(0,8),
#      main = "Urban Tree LAI (July) — Paris")
# 
# 
# # Hendaye
# hendaye_ext <- ext(-1.82, -1.75, 43.35, 43.38)
# plot(crop(lai_jul, hendaye_ext), col = hcl.colors(32, "YlGn"),
#      range = c(0, 8), main = "Urban Tree LAI — Hendaye")
# 
# # Hondarribia
# hondarribia_ext <- ext(-1.82, -1.78, 43.34, 43.38)
# plot(crop(lai_jul, hondarribia_ext), col = hcl.colors(32, "YlGn"),
#      range = c(0, 8), main = "Urban Tree LAI — Hondarribia")
# #
## 
# 
# lai_hond <- crop(lai_jul, hendaye_ext)
# 
# # tmap in "plot" mode for static figures
# tmap_mode("plot")
# 
# # Interactive mode uses web tiles such as OpenStreetMap
# tmap_mode("view")
# 
# tm_basemap("OpenStreetMap") +      # basemap underlay (works only in 'view' mode)
#   tm_shape(lai_jul) +
#   tm_raster(
#     palette = "YlGn",
#     alpha   = 0.6,                 
#     title   = "Urban Tree LAI (m²/m²)",
#     style   = "cont",
#     breaks  = seq(0, 8, by = 1)
#   ) +
#   tm_scale_bar(position = c("left","bottom")) +
#   tm_compass(position = c("right","top"))
# 


#### Combine global and urban ---------------
# Step 1: Align lai_global to lai_urban's resolution and extent
# x<-lai_global
# y <- lai_jul
# 
# # Convert both rasters to float (drop scale/offset)
# x_float <- writeRaster(x, "x_float.tif", overwrite=TRUE, datatype="FLT4S")
# y_float <- writeRaster(y, "y_float.tif", overwrite=TRUE, datatype="FLT4S")
# 
# # Resample using float versions
# r <- resample(x_float, y_float, method="bilinear",threads=T)
# summary(r)  # should stay ~0..7
# 
# # Step 4: Fill NA values in lai_jul
# r_filled <- cover(r, y)
# 
# # Step 34 Save the result
# writeRaster(r_filled, "data/environmental/lai/lai_july_filled.tif", overwrite = TRUE)
# 
# 
# # Example extent (Hondarribia), crop first
# hond <- ext(-3, -1.7, 43.34, 44)
# lai_hond <- crop(r_filled, hond)
# 
# pal <- colorNumeric(viridis(256, option = "G"), domain = c(0, 7), na.color = "transparent")
# 
# leaflet() %>%
#   addProviderTiles(providers$CartoDB.Positron) %>%
#   addRasterImage(r_filled, colors = pal, opacity = 0.6, project = TRUE) %>%
#   addLegend(pal = pal, values = c(0,7), title = "LAI (m²/m²)") %>%
#   addScaleBar()
# 



# ## 1F) Landcover winner (run once)----------------------------------------------------
# DIR_ENV <- "data/environmental"
# DIR_LC  <- file.path(DIR_ENV, "landcover")
# 
# # 12 classes:
# CONS_CLASSES <- 1:12
# lc_files <- file.path(DIR_LC, sprintf("consensus_full_class_%d.tif", CONS_CLASSES))
# 
# # Example AOI extent (Hondarribia)
# hond <- ext(-3, -1.7, 43.34, 44)
# 
# # Stack
# lc_stack <- rast(lc_files)            # one band per class
# 
# # Crop to AOI (faster plotting/testing)
# # lcc <- crop(lc_stack, hond)
# lcc<-lc_stack
# 
# # Winner-take-all collapse 
# # Pick the class (band index) with the maximum value per pixel
# lc_idx <- which.max(lcc)                      # integer 1..N
# # Mask out pixels where the max score is 0 or NA
# lc_max <- app(lcc, max, na.rm = TRUE)
# lc_idx[lc_max <= 0 | is.na(lc_max)] <- NA
# 
# # Factor labels (sentence case) 
# # Match the number of layers in  (should be 12)
# class_labels <- c(
#   "Evergreen/deciduous needleleaf trees",
#   "Evergreen broadleaf trees",
#   "Deciduous broadleaf trees",
#   "Mixed/other trees",
#   "Shrubs",
#   "Herbaceous vegetation",
#   "Cultivated and managed vegetation",
#   "Regularly flooded vegetation",
#   "Urban/built-up",
#   "Snow/ice",
#   "Barren",
#   "Open water"
# )
# 
# # Trim labels if stack has fewer than 12 layers
# n_classes <- nlyr(lcc)
# class_labels <- class_labels[1:n_classes]
# 
# # Attach factor levels
# lc_idx <- as.factor(lc_idx)
# levels(lc_idx) <- data.frame(ID = 1:n_classes, class = class_labels)
# 
# # Define output path in the same folder as input rasters
# out_file <- file.path(DIR_LC, "consensus_landcover_winner.tif")
# 
# # Write the raster
# writeRaster(lc_idx,
#             filename   = out_file,
#             overwrite  = TRUE,
#             datatype   = "INT1U")  # integer, unsigned (good for class IDs)

# ---- Leaflet map 
# Convert to 'raster' for addRasterImage()
#lc_idx_r <- raster::raster(lc_idx)

# # Build a stable color palette (length == n_classes)
# class_cols <- c(
#   "#1b9e77", # Evergreen/deciduous needleleaf trees
#   "#66a61e", # Evergreen broadleaf trees
#   "#7570b3", # Deciduous broadleaf trees
#   "#d95f02", # Mixed/other trees
#   "#e6ab02", # Shrubs
#   "#a6761d", # Herbaceous vegetation
#   "#e7298a", # Cultivated & managed
#   "#1f78b4", # Regularly flooded vegetation
#   "#b2df8a", # Urban/built-up
#   "#a6cee3", # Snow/ice
#   "#666666", # Barren
#   "#3892E0"  # Open water
# )[1:n_classes]
# 
# # Use numeric class IDs 1..n_classes as the domain
# pal <- colorFactor(palette = class_cols,
#                    domain  = 1:n_classes,
#                    na.color = "transparent")
# 
# # Pull labels from the factor levels to keep legend and raster in sync
# levs <- levels(lc_idx)[[1]]              # data.frame with columns ID, class
# legend_ids    <- levs$ID
# legend_labels <- levs$class
# 
# 
# leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
#   addProviderTiles(providers$CartoDB.Positron) %>%
#   addRasterImage(lc_idx_r, colors = pal, opacity = 0.85, project = TRUE) %>%
#   # Manual legend: provide colors + labels explicitly
#   addLegend(position = "bottomleft",
#             colors  = class_cols,
#             labels  = class_labels,
#             title   = "Land cover",
#             opacity = 1) %>%
#   addScaleBar(position = "bottomright")



# ____________________________ -----------------------
# 2) Extract environmental -------------------
## 2.1) Prep -------------
# Environmental data root and files
DIR_ENV    <- "data/environmental"
PATH_ACCESS<- file.path(DIR_ENV, "traveltime/accessibility_to_cities_2015_v1.0.tif")
PATH_TD    <- file.path(DIR_ENV, "tree.density/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif")
PATH_FROST <- file.path(DIR_ENV, "frost.days/cru_ts4.08.1901.2023.frs.dat.nc")
PATH_AGE   <- file.path(DIR_ENV, "forest.age/2025417153414222_BGIForestAgeMPIBGC1.0.0.nc")
DIR_CANOPY <- file.path(DIR_ENV, "canopy.height/chm")
DIR_CHELSA <- file.path(DIR_ENV, "CHELSA/bio")
DIR_LC     <- file.path(DIR_ENV, "landcover/")
PATH_WB    <- file.path(DIR_ENV, "woody.biomass.density/Aboveground_Live_Woody_Biomass_Density.geojson")
PATH_LAI<- file.path(DIR_ENV, "lai/lai_july_filled.tif")

PATH_CHELSA <- list.files(
  DIR_CHELSA,
  pattern = glob2rx("CHELSA_*.tif"),
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# Filter only .tif files for classes 1–9
PATH_LC <- list.files(DIR_LC,
                      pattern = "^(consensus_full_class_[1-9]\\.tif|consensus_landcover_winner\\.tif)$",
                      full.names = TRUE)


## 2.2) Extract------------
occ_covs_tbl <- pts %>% st_drop_geometry %>% 
  distinct(id, site, .keep_all = T)



{
  ### a) CHELSA ---------------
  for (f in PATH_CHELSA) {
    r  <- terra::rast(f)
    nm <- stringr::str_extract(basename(f), "^CHELSA_bio\\d+")
    if (is.na(nm)) nm <- tools::file_path_sans_ext(basename(f))
    vals <- terra::extract(r, pts)
    occ_covs_tbl[[nm]] <- as.numeric(vals[,2])
  }
  
  ### b) Accessibility ----------
  #This layer has a lot of NAs (-9999) on small islands (on which I have field sites). 
  #I'm filling those NAs with the nearest neighbour.
  # Read raster & make NAs
  r <- rast(PATH_ACCESS)
  r <- ifel(r < 0, NA, r)
  #writeRaster(r, file.path(DIR_ENV,"traveltime/access.min.NA.tif"))  #ONCE ONLY
  #r<-rast(file.path(DIR_ENV,"traveltime/access.min.NA.tif"))
  # Helper: great-circle (haversine) distance in meters for lon/lat
  .hav <- function(lon1, lat1, lon2, lat2) {
    R <- 6371000
    toRad <- pi/180
    dlon <- (lon2 - lon1) * toRad
    dlat <- (lat2 - lat1) * toRad
    a <- sin(dlat/2)^2 + cos(lat1*toRad)*cos(lat2*toRad)*sin(dlon/2)^2
    2*R*asin(pmin(1, sqrt(a)))
  }
  
  # Helper: nearest NON-NA raster value per point 
  # r      : SpatRaster 
  # pts    : sf or SpatVector points
  # win    : starting half-window in CELLS (e.g., 5 => 11x11); expands if needed
  # winmax : maximum half-window in CELLS (failsafe)
  extract_nearest_nonNA <- function(r, pts, win = 5, winmax = 60) {
    if (!inherits(pts, "SpatVector")) pts <- terra::vect(pts)
    pts <- terra::project(pts, terra::crs(r))
    xy  <- terra::crds(pts, df = TRUE)
    
    # direct extract first
    v0  <- terra::extract(r, pts, method = "simple", ID = FALSE)[[1]]
    out <- v0
    todo <- which(is.na(out))
    if (!length(todo)) return(out)
    
    rx <- terra::res(r)[1]; ry <- terra::res(r)[2]
    
    .hav <- function(lon1, lat1, lon2, lat2) {
      R <- 6371000
      toRad <- pi/180
      dlon <- (lon2 - lon1) * toRad
      dlat <- (lat2 - lat1) * toRad
      a <- sin(dlat/2)^2 + cos(lat1*toRad)*cos(lat2*toRad)*sin(dlon/2)^2
      2 * R * asin(pmin(1, sqrt(a)))
    }
    
    for (i in todo) {
      w <- win; found <- FALSE
      while (!found && w <= winmax) {
        e <- terra::ext(xy$x[i] - w*rx, xy$x[i] + w*rx,
                        xy$y[i] - w*ry, xy$y[i] + w*ry)
        rc <- terra::crop(r, e, snap = "out")
        vals <- terra::values(rc, mat = FALSE)
        if (all(is.na(vals))) { w <- w + 5; next }
        
        ok <- which(!is.na(vals))
        if (length(ok)) {
          cc <- terra::xyFromCell(rc, ok)
          d  <- .hav(xy$x[i], xy$y[i], cc[,1], cc[,2])
          out[i] <- vals[ ok[ which.min(d) ] ]
          found <- TRUE
        } else {
          w <- w + 5
        }
      }
    }
    out
  }
  
  
  # Extract with nearest NON-NA fallback
  vals_nn <- extract_nearest_nonNA(r, pts, win = 5, winmax = 80)
  
  # Assign
  occ_covs_tbl$access_min <- as.numeric(vals_nn)
  
  # env_vars_sel$access_min <- as.numeric(vals_nn)
  # 
  # (env_vars_sel) %>% filter(is.na(access_min)) %>%select(site,lat,lon,access_min) %>%  view()
  # (Optional) quick QA
  # sum(is.na(values(r_filled_land)))     # should be << masked version
  # plot(r_filled_land); points(pts_v)
  ### c) Land cover -----------
  for (f in PATH_LC) {
    r <- terra::rast(f)
    
    # Extract numeric class ID from filename
    class_id <- str_extract(basename(f), "\\d+")
    
    # Create name LC1, LC2, ..., LC9
    nm <- paste0("LC", class_id)
    
    nm <- if (!is.na(class_id)) {
      paste0("LC", class_id)  # LC1, LC2, ..., LC9
    } else {
      "LC_winner"              # For consensus_landcover_winner.tif
    }
    
    
    # Extract values at points 
    vals <- terra::extract(r, pts)
    
    # Store extracted values as numeric
    occ_covs_tbl[[nm]] <- as.numeric(vals[, 2])
  }
  
  
  ### d) Tree density -----------
  r <- terra::rast(PATH_TD)
  occ_covs_tbl$tree_density <- as.numeric(terra::extract(r, pts)[,2])
  
  ### e) Frost days (average across bands) --------
  r <- terra::rast(PATH_FROST)
  occ_covs_tbl$frost_days <- as.numeric(terra::extract(r, pts, fun = mean, na.rm = TRUE)[,2])
  
  ### f) Canopy height (VRT mosaic) ----------
  # chm_files <- list.files(DIR_CANOPY, pattern = "\\.tif$", full.names = TRUE)
  # chm_vrt <- tryCatch(terra::vrt(chm_files), error = function(e) NULL)
  # occ_covs_tbl$canopy_h_mean <- as.numeric(terra::extract(chm_vrt, pts)[,2])
  # 
  ### g) Forest age (BGI) ------------
  r <- tryCatch(
    terra::rast(PATH_AGE, subds = "ForestAge_TC000"),
    error = function(e) terra::rast(PATH_AGE)
  )
  occ_covs_tbl$forest_age <- as.numeric(terra::extract(r, pts)[,2])
  
  
  ### h) LAI----------
  
  r <- terra::rast(PATH_LAI)
  occ_covs_tbl$LAI <- as.numeric(terra::extract(r, pts)[,2])
  }




## 2.3) CENTRE ALL --------------
field_lat_mean <- mean(unique(my_sites$lat), na.rm = TRUE)

occ_covs_tbl_z <- occ_covs_tbl %>%
  mutate(
    # center latitude to the field mean first
    lat_cdeg = lat - field_lat_mean
  ) %>%
  mutate(
    # z-score all numeric columns except lon, lat, and lat_cdeg
    across(
      where(is.numeric) & !any_of(c("lon", "lat", "lat_cdeg")),
      ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>% 
  select(!matches("CHELSA_kg"))



## 2.4) Write to csv----
write_csv(occ_covs_tbl, file.path(OUT_DIR_ENV,"EUROPE_occ_covs_tbl.csv"))
write_csv(occ_covs_tbl_z, file.path(OUT_DIR_ENV,"EUROPE_occ_covs_tbl_z.csv"))
#occ_covs_tbl_z<-fread("data/08_pca/occ_covs_tbl_z") %>% as_tibble()



#________________________________---------------------
# 3) Classify points as island / mainland ---------------------------

## 3a) Prep------------
# Natural Earth paths 
NE_BASE <- "C:/Users/Marie/Documents/LatGradientsAccEpis/github/data/environmental/island"
NE_LAND_DIR <- file.path(NE_BASE, "ne_10m_land")                 # folder containing the land shapefile/GeoPackage
NE_LAND_FILE <- list.files(NE_LAND_DIR, pattern = "\\.(gpkg|shp)$", full.names = TRUE)[1]
NE_MINOR_FILE <- file.path(NE_BASE, "ne_10m_minor_islands", "ne_10m_minor_islands.shp")


# Load NE layers
stopifnot(file.exists(NE_LAND_FILE), file.exists(NE_MINOR_FILE))
ne_land  <- st_read(NE_LAND_FILE, quiet = TRUE)
ne_minor <- st_read(NE_MINOR_FILE, quiet = TRUE)

# Split land by scalerank
land_major   <- ne_land  %>% dplyr::filter(!is.na(scalerank) & scalerank < 3)   # mainland / major landmasses
land_islands <- ne_land  %>% dplyr::filter(!is.na(scalerank) & scalerank >= 3)  # islands 


# Use an equal-area or leave in geographic; for intersects, CRS must match
common_crs <- st_crs(pts)
if (is.na(common_crs)) common_crs <- st_crs(4326)

pts   <- st_transform(pts, common_crs)
major  <- st_transform(land_major, common_crs)
islands  <- st_transform(land_islands, common_crs)
minor <- st_transform(ne_minor, common_crs)

# Fix invalid geometries if needed
suppressWarnings({
  if (any(!st_is_valid(land)))  major  <- st_make_valid(major)
  if (any(!st_is_valid(land)))  islands  <- st_make_valid(islands)
  if (any(!st_is_valid(minor))) minor <- st_make_valid(minor)
})


## 3b) CLASSIFY ------------
# Add buffer to land and minor island polygons (e.g., 500 meters)
buffer_dist <- units::set_units(5000, "m")  # adjust as needed

# Ensure CRS is projected (not geographic) for buffering
pts_crs <- st_crs(pts)
if (pts_crs$units_gdal == "degree") {
  pts <- st_transform(pts, 3857)  # Web Mercator for buffering
  major <- st_transform(major, 3857)
  islands <- st_transform(islands, 3857)
  minor <- st_transform(minor, 3857)
}

# Apply buffer
land_buf  <- st_buffer(major, dist = buffer_dist)
island_buf <- st_buffer(islands, dist = buffer_dist)
minor_buf <- st_buffer(minor, dist = buffer_dist)

# Intersect with buffered layers
minor_hit <- lengths(st_intersects(pts, minor_buf)) > 0
island_hit <- lengths(st_intersects(pts, island_buf)) > 0
land_hit <- lengths(st_intersects(pts, land_buf)) > 0

# Compose label
island_class <- ifelse(minor_hit|island_hit, "island",
                       ifelse(land_hit, "land", "island"))

## 3d) JOIN-----------
env_island <- st_drop_geometry(pts) %>%
  left_join(out) 
fwrite(env_island, "points_with_islands.csv")


occ_covs_tbl %>% 
  left_join(out)->occ_covs_tbl
occ_covs_tbl_z %>% 
  left_join(out)->occ_covs_tbl_z

## 2.4) Write to csv----
write_csv(occ_covs_tbl, file.path(OUT_DIR_ENV,"EUROPE_occ_covs_tbl.csv"))
write_csv(occ_covs_tbl_z, file.path(OUT_DIR_ENV,"EUROPE_occ_covs_tbl_z.csv"))
#occ_covs_tbl_z<-fread("data/08_pca/occ_covs_tbl_z") %>% as_tibble()

# Visualise results
library(leaflet)

leaflet(data = out) %>%    addTiles() %>%
  addCircleMarkers(
    lng = ~lon, lat = ~lat,
    label = ~paste0("ID: ", id, "<br>Site: ", site, "<br>Class: ", island_class),
    color = ~ifelse(island_class == "island", "blue",
                    ifelse(island_class == "land", "green", "gray")),
    radius = 5, fillOpacity = 0.7
  )







#____________________________-----------------
# 4) Shortlist environmental variables--------------------------
## 4A) Correlation tests-----------
# Which ones correlate with latitude?
#occ_covs_tbl<-fread("data/08_pca/EUROPE_occ_covs_tbl_z.csv")
dat_env <-  occ_covs_tbl %>% 
  select(!matches(c("CHELSA_kg","LC"))) %>% 
  select(!matches("_z")) %>% 
  left_join(env_island) %>% 
  mutate(is_island=if_else(island_class=="island",1,0)) %>% 
  select(-island_class)


dat_env %>% 
  pivot_longer(-c(1:7), names_to = "var",values_to = "value") %>%   
  ggplot(aes(x = lat, y = value)) +
  facet_wrap(~var,scales="free_y")+
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm",
              se = FALSE) +
  labs(x = "Latitude", y = "Variable") ->p_lat_cor

ggsave(file.path(OUT_DIR_ENV,"env~lat.correlation.jpeg"),
       p_lat_cor,
       width=20,height=20, dpi=800)

#Which ones correlate with presence?
dat_env %>% 
  mutate(taxon.name=if_else(taxon.name%in%c("",NA), "no_taxon",taxon.name),
         presence=if_else(taxon.name=="no_taxon",0,1)) %>% 
  pivot_longer(-c(1:11,presence), names_to = "var",values_to = "value") %>%   
  ggplot(aes(y = presence,x = value)) +
  facet_wrap(~var,scales="free_x")+
  #geom_point(alpha = 0.4) +
  geom_smooth(method = "glm", se = FALSE, 
              method.args = list(family="binomial")) +
  labs(x = "Latitude", y = "Variable") ->p_pres_cor

ggsave(file.path(OUT_DIR_ENV,"pres.correlation.jpeg"),
       p_pres_cor,
       width=20,height=20, dpi=800)


cor_rank <- dat_env %>%
  select(-matches("_z")) %>%
  mutate(
    taxon.name = if_else(is.na(taxon.name) | taxon.name == "", "no_taxon", taxon.name),
    presence   = as.integer(taxon.name != "no_taxon")
  ) %>%
  # keep env variables: drop first 11 meta columns + presence
  pivot_longer(cols = -c(1:11, presence),
               names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(
    r = cor(value, presence, use = "pairwise.complete.obs"),
    n = sum(is.finite(value) & is.finite(presence)),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(r)))



#Which ones correlate with each other?
# Presence 
dat_env2 <- dat_env %>%
  select(-matches("_z")) %>%
  mutate(
    taxon.name = if_else(is.na(taxon.name) | taxon.name == "", "no_taxon", taxon.name),
    presence   = as.integer(taxon.name != "no_taxon")
  ) %>% as_tibble() %>% 
  select(1:6,presence,everything())

meta_cols <- names(dat_env2)[1:7]  
# Numeric env matrix (no meta, no presence)
env_num <- dat_env2 %>%
  select(-all_of(meta_cols), -presence) %>%
  select(where(is.numeric)) %>%
  as.data.frame()


# Choose how many top-ranked candidates to consider (e.g., top 40 by |r|)
TOP_M  <- 40
thresh <- 0.60   # max allowed |pairwise correlation|

cands <- cor_rank %>%
  arrange(desc(abs(r))) %>%
  slice_head(n = TOP_M) %>%
  pull(var)

# Keep only candidates that are in env_num
cands <- intersect(cands, colnames(env_num))

# Correlation matrix among candidates
cm <- cor(env_num[, cands, drop = FALSE], use = "pairwise.complete.obs")

# Walk down the ranking, enforce max |r| <= thresh
keep <- character(0)
for (v in cands) {
  if (length(keep) == 0) {
    keep <- c(keep, v)
  } else {
    if (all(abs(cm[v, keep]) <= thresh, na.rm = TRUE)) {
      keep <- c(keep, v)
    }
  }
  if (length(keep) == 12) break
}

vars_10_select <- keep
print(vars_10_select)

## 4B) Create lookup table------
var_lookup <- tribble(
  ~var,                                   ~label,
  "CHELSA_bio7",                          "Temperature annual range (°C)",
  "CHELSA_npp_1981-2010_V.2.1",           "Net primary productivity (gC m⁻² yr⁻¹)",
  "CHELSA_bio3",                          "Isothermality (%)",
  "access_min",                           "Travel time to nearest city (mins)",
  "CHELSA_bio6",                          "Minimum temperature of coldest month (°C)",
  "CHELSA_vpd_max_1981-2010_V.2.1",       "Max vapour pressure deficit (kPa)",
  "CHELSA_fgd_1981-2010_V.2.1",           "First day of growing season",
  "CHELSA_sfcWind_max_1981-2010_V.2.1",   "Max near-surface wind speed (ms⁻¹)",
  "CHELSA_bio15",                         "Precipitation seasonality (CV, %)",
  "forest_age",                           "Forest age (years)",
  "LAI", "LAI (m² m-²)",
  "landcover", "Land cover class",
  "lat_cdeg", "Centred latitude (°)",
  "is_island","Island"
  
) 

##4C) filter chosen variables -----------
lat_mean_field <- mean(unique(my_sites$lat), na.rm = TRUE)

env_vars_sel <- dat_env %>%
  mutate(lat_cdeg = lat - lat_mean_field) %>%
  left_join(env_island) %>% 
  reframe(site,id,lon,lat,
          lat_cdeg,taxon.name,quality_grade,
          presence,across(all_of(vars_10_select)))

##4D) Write output----------
env_vars_sel %>% fwrite(file.path(OUT_DIR_ENV,"env_vars_sel.csv"))
var_lookup %>% fwrite(file.path(OUT_DIR_ENV,"var_lookup.csv"))

