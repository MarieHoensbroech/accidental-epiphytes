# 11_koeppen_zones.R
#
# 1) Create per-taxon overlap maps of Observation-derived ranges vs iNat ranges
# 2) Create a combined overlap map for all taxa
# 3) Compute area-weighted composition of Köppen–Geiger climate classes
#      inside the union of Observation-derived ranges, and visualize as
#      a bar chart + map overlay
# 4) Plot
#

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(exactextractr)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(tidyterra)  
  library(ggnewscale)  
})

# Parameters ------------------------------------
GPKG_OBS   <- "data/spatial/obs_ranges.gpkg"
LYR_OBS    <- "obs_range"

GPKG_INAT  <- "data/spatial/inat_ranges.gpkg"
LYR_INAT   <- "inat_range"

IN_INAT_OBS <- "data/processed/inat_observations.csv"  

# NatCap COG (Köppen–Geiger, 1-km, 1991–2020)
KG_URL <- paste0(
  "https://storage.googleapis.com/natcap-data-cache/global/koppen_geiger_climatezones/",
  "koppen_geiger_climatezones_1991_2020_1km.tif"
)

OUT_DIR_RANGES  <- "outputs/figures/11_koeppen/ranges"
OUT_DIR_KG      <- "outputs/figures/11_koeppen"
OUT_DIR_TABLES  <- "outputs/tables/11_koeppen"
dir.create(OUT_DIR_RANGES, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_KG,     showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_TABLES, showWarnings = FALSE, recursive = TRUE)


# 0) Load core data --------------------------------
obs_ranges  <- st_read(GPKG_OBS,  layer = LYR_OBS,  quiet = TRUE) %>% st_make_valid()
inat_ranges <- st_read(GPKG_INAT, layer = LYR_INAT, quiet = TRUE) %>% st_make_valid()

# Normalize CRS (WGS84)
if (is.na(st_crs(obs_ranges)))  st_crs(obs_ranges)  <- 4326
if (is.na(st_crs(inat_ranges))) st_crs(inat_ranges) <- 4326
obs_ranges  <- st_transform(obs_ranges, 4326)
inat_ranges <- st_transform(inat_ranges, 4326)

# taxon_id -> taxon.name
taxon_names <- fread(IN_INAT_OBS) %>% 
  reframe(taxon_id = as.integer(taxon_id),
          taxon_name = as.character(`taxon.name`))  %>% 
  distinct()

# World basemap (single union geometry helps when plotting bbox)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Helpers (geometry & plotting) 
buffered_extent <- function(geom, buffer_km = 300) {
  # returns bbox (in lon/lat) around geom buffered by buffer_km
  g <- st_union(geom)
  g_m <- st_transform(g, 3857)
  g_b <- st_buffer(g_m, dist = buffer_km * 1000)
  st_bbox(st_transform(g_b, 4326))
}

as_multipolygon <- function(sfc) {
  sfc <- st_make_valid(sfc)
  suppressWarnings(st_collection_extract(sfc, "POLYGON")) %>%
    suppressWarnings(st_cast("MULTIPOLYGON", warn = FALSE))
}

# 1) Per-taxon overlap maps -----------------------------
# Build a taxon table present in either dataset
taxa <- sort(unique(c(obs_ranges$taxon_id, inat_ranges$taxon_id)))
taxa_tbl <- tibble(taxon_id = taxa) %>%
  left_join(taxon_names, by = "taxon_id")

for (i in seq_along(taxa)) {
  tid <- taxa[i]
  msg("Per-taxon plot [%d/%d] taxon_id=%s", i, length(taxa), tid)
  
  obs_g  <- obs_ranges  %>%
    filter(taxon_id == tid) %>%
    st_geometry() %>%
    as_multipolygon()
  inat_g <- inat_ranges %>%
    filter(taxon_id == tid) %>%
    st_geometry() %>%
    as_multipolygon()
  
  if ((length(obs_g) == 0 || any(st_is_empty(obs_g))) &&
      (length(inat_g) == 0 || any(st_is_empty(inat_g)))) {
    next
  }
  
  # Build combined sf for plotting
  layers <- list()
  if (length(inat_g) > 0 && !any(st_is_empty(inat_g))) {
    layers[[length(layers) + 1]] <- st_sf(type = "iNaturalist range", geometry = inat_g)
  }
  if (length(obs_g) > 0 && !any(st_is_empty(obs_g))) {
    layers[[length(layers) + 1]] <- st_sf(type = "Observed range",   geometry = obs_g)
  }
  layers_sf <- do.call(rbind, layers)
  
  # Figure bbox
  bb <- buffered_extent(layers_sf, buffer_km = 300)
  
  # Title
  taxon_label <- taxa_tbl$taxon_name[taxa_tbl$taxon_id == tid]
  taxon_label <- if (length(taxon_label) && !is.na(taxon_label)) taxon_label[1] else paste0("taxon_id ", tid)
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80") +
    geom_sf(data = layers_sf, aes(fill = type), alpha = 0.5, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Observed range" = "#F7F056", "iNaturalist range" = "#90C987"),
                      name = "Range type") +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"])) +
    ggtitle(glue::glue("Range overlap for: *{taxon_label}*")) +
    theme_minimal() +
    theme(plot.title = ggtext::element_markdown(size = 13),
          legend.position = "bottom")
  
  outf <- file.path(OUT_DIR_RANGES,
                    sprintf("taxon_range_%s_%s.jpg",
                            tid,
                            str_replace_all(taxon_label, "[^A-Za-z0-9]+", "_")))
  ggsave(outf, p, width = 9, height = 5.5, dpi = 300, bg = "white")
}

# 2) Combined overlap map -------------------------------
obs_union  <- obs_ranges  %>%
  st_geometry() %>%
  as_multipolygon() %>%
  st_union()

inat_union <- inat_ranges %>%
  st_geometry() %>%
  as_multipolygon() %>%
  st_union()

combined_layers <- dplyr::bind_rows(
  st_sf(type = "iNaturalist range", geometry = inat_union),
  st_sf(type = "Observed range",   geometry = obs_union)
)

p_combined <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80") +
  geom_sf(data = combined_layers, aes(fill = type),
          alpha = 0.5, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Observed range" = "#F7F056", "iNaturalist range" = "#90C987"),
                    name = "Range type") +
  coord_sf() +
  labs(title = "Combined range overlap across all taxa") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR_RANGES, "combined_range_map.jpg"),
       p_combined, width = 12, height = 8, dpi = 300, bg = "white")

# 3) Köppen–Geiger composition --------------------------
kg <- terra::rast(KG_URL)  # 1..30 class ids, NAs over ocean
names(kg) <- "kg"

# Crop to speed up (to observed union bbox)
kg_crop <- tryCatch(terra::crop(kg, terra::vect(st_as_sfc(st_bbox(obs_union)))), error = function(e) kg)

# Area per cell (km^2), on the kg grid
area_km2 <- terra::cellSize(kg_crop, unit = "km")
names(area_km2) <- "cell_area_km2"

# Exact area-weighted extraction by class inside Observed union
stacked <- c(kg_crop, area_km2)

# Keep only polygons and make them valid
obs_union <- obs_union %>%
  sf::st_collection_extract("POLYGON") %>%
  sf::st_make_valid()

# Run exact_extract
res_tbl <- exactextractr::exact_extract(
  stacked,
  obs_union, 
  fun = function(df, ...) {
    df %>%
      dplyr::filter(!is.na(kg)) %>%
      dplyr::mutate(partial_area = cell_area_km2 * coverage_fraction) %>%
      dplyr::group_by(kg) %>%
      dplyr::summarise(area_km2 = sum(partial_area, na.rm = TRUE), .groups = "drop")
  },
  summarize_df = TRUE,
  progress = TRUE
)

total_area <- sum(res_tbl$area_km2, na.rm = TRUE)
kg_summary <- res_tbl %>%
  mutate(percent = 100 * area_km2 / total_area) %>%
  arrange(desc(percent))

# Legend/table for KG IDs -> symbols & colors
koppen_map <- tibble::tribble(
  ~kg, ~koppen, ~description, ~r, ~g, ~b,
  1, "Af","Tropical, rainforest", 0, 0, 255,
  2, "Am","Tropical, monsoon", 0,120,255,
  3, "Aw","Tropical, savannah", 70,170,250,
  4, "BWh","Arid, desert, hot", 255, 0, 0,
  5, "BWk","Arid, desert, cold", 255,150,150,
  6, "BSh","Arid, steppe, hot", 245,165,  0,
  7, "BSk","Arid, steppe, cold",255,220,100,
  8, "Csa","Temperate, dry summer, hot",255,255, 0,
  9, "Csb","Temperate, dry summer, warm",200,200, 0,
  10, "Csc","Temperate, dry summer, cold",150,150, 0,
  11, "Cwa","Temperate, dry winter, hot",150,255,150,
  12, "Cwb","Temperate, dry winter, warm",100,200,100,
  13, "Cwc","Temperate, dry winter, cold", 50,150, 50,
  14, "Cfa","Temperate, no dry, hot",200,255, 80,
  15, "Cfb","Temperate, no dry, warm",100,255, 80,
  16, "Cfc","Temperate, no dry, cold", 50,200,  0,
  17, "Dsa","Cold, dry summer, hot",255,  0,255,
  18, "Dsb","Cold, dry summer, warm",200,  0,200,
  19, "Dsc","Cold, dry summer, cold",150, 50,150,
  20, "Dsd","Cold, dry summer, very cold",150,100,150,
  21, "Dwa","Cold, dry winter, hot",170,175,255,
  22, "Dwb","Cold, dry winter, warm", 90,120,220,
  23, "Dwc","Cold, dry winter, cold", 75, 80,180,
  24, "Dwd","Cold, dry winter, very cold", 50,  0,135,
  25, "Dfa","Cold, no dry, hot",  0,255,255,
  26, "Dfb","Cold, no dry, warm", 55,200,255,
  27, "Dfc","Cold, no dry, cold",  0,125,125,
  28, "Dfd","Cold, no dry, very cold",  0, 70, 95,
  29, "ET","Polar, tundra", 178,178,178,
  30, "EF","Polar, frost", 102,102,102
) %>%
  mutate(hex = sprintf("#%02X%02X%02X", r, g, b))

kg_summary <- kg_summary %>%
  left_join(koppen_map, by = "kg") %>%
  arrange(desc(percent))

# Persist table
readr::write_csv(kg_summary, file.path(OUT_DIR_TABLES, "koeppen_composition.csv"))

# 4) Plot ----------
## 4.1) Bar plot (percent by class) ----
kg_cols <- setNames(kg_summary$hex, kg_summary$koppen)

p_bar <- ggplot(kg_summary,
                aes(x = reorder(koppen, percent), y = percent, fill = koppen)) +
  geom_col(width = 0.9) +
  coord_flip() +
  geom_text(aes(label = description),
            hjust = 1.02, size = 3.1, color = "black") +
  scale_fill_manual(values = kg_cols, guide = "none") +
  labs(x = "Köppen class", y = "% of observed range union",
       title = "Area-weighted Köppen–Geiger composition") +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR_KG, "koeppen_bar_percent.jpg"),
       p_bar, width = 7.5, height = 6.5, dpi = 300, bg = "white")

## 4.2) Map overlay (KG raster + observed union) ----
kg_fac <- terra::as.factor(kg_crop)
kg_hex <- setNames(koppen_map$hex, as.character(koppen_map$kg))
kg_lab <- setNames(koppen_map$koppen, as.character(koppen_map$kg))

p_kg_map <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
  tidyterra::geom_spatraster(data = kg_fac, alpha = 0.6) +
  scale_fill_manual(
    guide = "none",
    values = kg_hex,
    breaks = names(kg_lab),
    labels = kg_lab,
    name = "Köppen–Geiger",
    na.value = "transparent"
  ) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = st_sf(range_type = "Observed range", geometry = obs_union),
          fill = NA, color = "black", linewidth = 0.4) +
  coord_sf(ylim = c(-60, 80), expand = FALSE) +
  labs(title = "Observed ranges over Köppen–Geiger climate zones") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR_KG, "koeppen_map_overlay.jpg"),
       p_kg_map, width = 10, height = 7.2, dpi = 300, bg = "white")

# . --------------
# Check: Find accidental epiphytes in ET (tundra) occurrences check ---------------------
pts <- st_read(GPKG_PTS, layer = LYR_PTS, quiet = TRUE)
if (is.na(st_crs(pts))) st_crs(pts) <- 4326
pts <- st_transform(pts, 4326)

msg("Extracting KG ids for points…")
pts$kg_id <- terra::extract(kg_crop, terra::vect(pts))[,1]
occ_et <- pts %>% filter(kg_id == 29)  # 29 == ET

n_total <- nrow(pts); n_et <- nrow(occ_et)
msg("Occurrences in ET: %s / %s (%.1f%%)", n_et, n_total, 100 * n_et / max(n_total, 1))

# Create masked ET raster for display within observed range bbox
bb  <- st_bbox(obs_union)
ext <- terra::ext(bb$xmin, bb$xmax, bb$ymin, bb$ymax)
et_only <- terra::ifel(kg_crop == 29, 29, NA)
et_c <- terra::crop(et_only, ext)
et_fac <- terra::as.factor(et_c)

et_hex <- c(`29` = kg_hex[["29"]])

p_et <- ggplot() +
  tidyterra::geom_spatraster(data = et_fac, alpha = 0.4) +
  scale_fill_manual(values = et_hex, breaks = "29",
                    labels = "ET — Polar tundra", guide = "none",
                    na.value = "transparent") +
  ggnewscale::new_scale_color() +
  geom_sf(data = st_sf(geometry = obs_union), fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = occ_et, aes(color = "ET occurrences"), size = 1.1, alpha = 0.85, show.legend = TRUE) +
  scale_color_manual(values = c("ET occurrences" = "red"), name = NULL) +
  coord_sf(xlim = c(bb$xmin, bb$xmax), ylim = c(bb$ymin, bb$ymax), expand = FALSE) +
  labs(title = "Occurrences falling in Köppen ET (Polar tundra)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR_KG, "et_occurrences_map.jpg"),
       p_et, width = 8, height = 6, dpi = 300, bg = "white")
