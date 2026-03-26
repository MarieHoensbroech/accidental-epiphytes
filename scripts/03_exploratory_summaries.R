# 03b_exploratory+summaries.R
# 1) Calculate track length + non-overlap checks for MULTILINESTRING (WGS84) ----
# 1b) Plot our track + field sites
# 2) Summary stats
# 

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(sf)
  library(units)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
  library(stringr)
  library(lubridate)
  
})

IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"
OUT_DIR_SUM <- "outputs/03_exploratory"
dir.create(OUT_DIR_SUM)

# Theme 
my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top"
  )

my_theme14_inside <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top",
    axis.ticks.length = unit(-0.2, "cm"),
    axis.text.y = element_text(margin = margin(l = 20, r = -32)),
    axis.text.x = element_text(margin = margin(t = -18, b = 30))
  )

#Load data ---------
my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% #glimpse()
  reframe(site, dbh_cat=TreeCat, setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon,area_ha)

inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% 
  reframe(site, dbh_cat=TreeCat, setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          epi_height=EpiHeight_num,
          genus_name, obs_date, obs_group=ObsGroup) %>% 
  unique() 

all_inat   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  reframe(lat=latitude,
          lon=longitude,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          taxon.rank,
          epi_height = EpiHeight_num,
          genus_name,
          user_login,
          obs_date) %>% 
  unique() %>% 
  mutate(epi_height=if_else(epi_height>=15,epi_height/100,epi_height))

my_inat <- all_inat %>% filter(user_login=="marie-ho")

# 1) TRACK ---------

## 1a) Calculate non-overlapping sailing distance ---------
# Read data 
tracks <- st_read("data/raw/fieldsitesGIS/tracks/tracks3110hondarribia.shp", quiet = TRUE)


sf::sf_use_s2(TRUE)  # geodesic on lon/lat
tracks_ll <- st_transform(tracks, 4326)

# Geodesic (no projection) — count overlaps once
g_union_ll <- st_union(tracks_ll)           # dissolve to unique linework (GEOMETRYCOLLECTION)
g_merge_ll <- st_line_merge(g_union_ll)     # merge contiguous segments
len_union_geo_km <- set_units(st_length(g_merge_ll), "km")

cat(sprintf("Non-overlapping total length (geodesic, WGS84): %0.3f km\n",
            as.numeric(len_union_geo_km)))

# Planar version (project to UTM first), sometimes preferable for local datasets
centroid_ll <- st_transform(st_centroid(st_union(st_geometry(tracks_ll))), 4326)
lon <- st_coordinates(centroid_ll)[, "X"]; lat <- st_coordinates(centroid_ll)[, "Y"]
utm_zone <- floor((lon + 180)/6) + 1
epsg_utm <- (if (lat >= 0) 32600 else 32700) + utm_zone  

tracks_utm <- st_transform(tracks_ll, epsg_utm)

g_union_utm <- st_union(tracks_utm)
g_merge_utm <- st_line_merge(g_union_utm)
len_union_planar_km <- set_units(st_length(g_merge_utm), "km")

cat(sprintf("Non-overlapping total length (planar, UTM EPSG:%s): %0.3f km\n",
            epsg_utm, as.numeric(len_union_planar_km)))


#nautical miles
len_union_planar_km / 1.805


## 1b) PLOT track--------


# Compute UTM EPSG dynamically
centroid_ll <- st_transform(st_centroid(st_union(st_geometry(g_merge_utm))), 4326)
lon <- st_coordinates(centroid_ll)[, "X"]
lat <- st_coordinates(centroid_ll)[, "Y"]
utm_zone <- floor((lon + 180) / 6) + 1
epsg_utm <- (if (lat >= 0) 32600 else 32700) + utm_zone

# Transform track to UTM
tracks_utm <- st_transform(tracks_ll, epsg_utm)
g_union_utm <- st_union(tracks_utm)
g_merge_utm <- st_line_merge(g_union_utm)

# Compute bounding box and expand by ±100 km
bbox <- st_bbox(g_merge_utm)
buffer_m <- 100000  # 100 km
bbox_expanded <- bbox
bbox_expanded[c("xmin","ymin")] <- bbox_expanded[c("xmin","ymin")] - buffer_m
bbox_expanded[c("xmax","ymax")] <- bbox_expanded[c("xmax","ymax")] + buffer_m

# Get country polygons and transform to UTM
countries <- ne_countries(scale = "medium", returnclass = "sf")
countries_utm <- st_transform(countries, epsg_utm)

# Crop countries to expanded bbox
bbox_poly <- st_as_sfc(bbox_expanded, crs = epsg_utm)
countries_crop <- st_crop(countries_utm, bbox_poly)


# Convert my_sites to sf 
my_sites_sf <- st_as_sf(my_sites %>% distinct(site,lat,lon), coords = c("lon", "lat"), crs = 4326)  
my_sites_sf <- st_transform(my_sites_sf, crs = st_crs(countries_crop))  # Match map CRS

# Plot
p_track <-
  ggplot() +
  geom_sf(data = countries_crop, fill = "grey90", color = "grey70") +
  geom_sf(data = st_as_sf(g_merge_utm), aes(color = "Sailed track"), size = 1) +
  geom_sf(data = my_sites_sf, aes(fill = "Field sites"), size = 5, alpha = 0.6, 
          shape = 21, stroke = NA) +
  
  coord_sf(xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
           ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
           expand = FALSE) +
  scale_color_manual(name = NULL, values = c("Sailed track" = "#F4C542")) +  ##8C7A64, #7Ad151FF
  scale_fill_manual(name = NULL, values = c("Field sites" = "#8C7A64")) + ##E8E1D4, #2A7882AF
  
  annotation_scale(location = "tr", width_hint = 0.3) +
  scale_y_continuous(expand = TRUE) +
  my_theme14_inside +
  theme(legend.position = c(0.05, 0.95),  
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white",0.3), color = NA),
        legend.direction = "horizontal")

p_track

ggsave(file.path(OUT_DIR_SUM, "tracks+sites_map.png"), 
       p_track, width = 5, height = 7, dpi = 300)

ggsave(file.path(OUT_DIR_SUM, "tracks+sites_map.svg"), 
       p_track, width = 5, height = 7, dpi = 300)


# _______________________ ------------
# _______________________ ------------


# 2) Summaries-----------

SINCE_DATE   <- as.Date("2025-04-01")

# Normalize ranks considered as "species-level" (counts unique taxon_id)
species_rank <- c("species", "subspecies", "variety", "form", "infrahybrid", "hybrid", "complex")

## 2A) TREES & AREA ----
n_trees_sampled <- sum(my_sites$total_trees, na.rm = TRUE)
n_trees_with_epis <- inat_field %>% n_distinct(inat_field$obs_group)

# Area: prefer a per-site area column if present; else assume constant plot size
n_sites <- dplyr::n_distinct(my_sites$site)
area_ha <- sum(my_sites %>% distinct(site,area_ha) %>% pull())

## 2B) FIELD WORK: individuals & species (inat_field) -----
field_observations <- inat_field %>% 
  distinct(id) %>%
  summarise(n_distinct(id, na.rm=T)) %>% deframe()

field_individuals <- inat_field %>% 
  distinct(site,obs_count) %>%
  summarise(sum(obs_count, na.rm=T)) %>% deframe()


field_species <- inat_field %>%
  filter(taxon.rank%in%c(species_rank)) %>%
  summarise(n_species = n_distinct(taxon_id, na.rm=T)) %>%
  deframe()


## 2C) TOTALS SINCE APR 2025 ----
# My iNat
my_inat2 <- my_inat %>%
  anti_join(inat_field) %>%
  filter(is.na(obs_date) | obs_date >= SINCE_DATE) 

my_obs_since <- n_distinct(my_inat2$id)
my_individuals_since <- my_inat2 %>% distinct(id,count) %>% summarise(sum(count, na.rm = TRUE)) %>% pull()
my_species_since <- my_inat2 %>%
  filter(!is.na(taxon_id), taxon.rank%in%c(species_rank)) %>%
  summarise(n_species = n_distinct(taxon_id)) %>%
  pull(n_species) %||% 0L

# All iNat 
all_inat2 <- all_inat %>%
  filter(user_login!="marie-ho") %>% 
  distinct(id,obs_date,count,taxon_id,taxon.rank,user_login) %>% 
  filter(is.na(obs_date) | obs_date >= SINCE_DATE)

all_obs_since <- n_distinct(all_inat2$id)
all_individuals_since <- sum(all_inat2$count, na.rm = TRUE)

all_species_since <- all_inat2 %>%
  filter(!is.na(taxon_id), taxon.rank%in%c(species_rank)) %>%
  summarise(n_species = n_distinct(taxon_id)) %>%
  pull(n_species) %||% 0L


all_users<-  all_inat2 %>% summarise(n_distinct(user_login)) %>% deframe()

# All observations


## 2D) Formatted summary  ----
poster_stats <- tibble::tibble(
  metric = c(
    "Area sampled (ha)",
    "Trees sampled",
    "Trees with accidental epiphytes",
    "Fieldwork: observations",
    "Fieldwork: individuals",
    "Fieldwork: species",
    "My iNat since 2025-04-01: observations",
    "My iNat since 2025-04-01: individuals",
    "My iNat since 2025-04-01: species",
    "Other iNat since 2025-04-01: observations",
    "Other iNat since 2025-04-01: individuals",
    "Other iNat since 2025-04-01: species",
    "Number of citizen scientists that contributed"
  ),
  value = c(
    round(area_ha, 2),
    n_trees_sampled,
    n_trees_with_epis,
    field_observations,
    field_individuals,
    field_species,
    my_obs_since,
    my_individuals_since,
    my_species_since,
    all_obs_since,
    all_individuals_since,
    all_species_since,
    all_users
  )
)

print(poster_stats)
fwrite(poster_stats, file.path(OUT_DIR_SUM,"poster_summary.csv"))

## 2E) Poster one-liners-------------
{
cat(sprintf("Sampled %s trees &", format(n_trees_sampled, big.mark = "")))
cat(sprintf(" on %.2f ha\n", area_ha))
cat(sprintf("Observed %s accidental epiphyte indiciduals & %s species\n",
            format(field_individuals, big.mark = ","),
            format(field_species, big.mark = ",")))
cat(sprintf("My iNat since %s – individuals: %s | species: %s\n",
            format(SINCE_DATE), format(my_individuals_since, big.mark = ","),
            format(my_species_since, big.mark = ",")))
cat(sprintf("All iNat since %s – individuals: %s | species: %s\n",
            format(SINCE_DATE), format(all_individuals_since, big.mark = ","),
            format(all_species_since, big.mark = ",")))
}

