# 08_epiphyte_quotients.R
# 
#   1) Compute true-epiphyte (TEpis) and accidental-epiphyte (AEpis) quotients per
#   2) TDWG Level-3 region using GIFT checklists, restricting AEpis to regions
#   3) intersecting observed ranges from Step 2, and export quotients + maps.
#

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(sf)
  library(GIFT)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ----------------------------- Parameters -------------------------------------
IN_INAT_OBS   <- "data/processed/inat_observations.csv"      # Step 1 output
IN_TEPIS_CSV  <- "data/raw/ZotzEpiList.csv"                  # 
GPKG_OBS      <- "data/spatial/obs_ranges.gpkg"              # Step 2 output
LYR_OBS       <- "obs_range"

OUT_DIR_TBL   <- "outputs/tables/08_epiphyte_quotient"
OUT_DIR_FIG   <- "outputs/figures/08_epiphyte_quotient"
dir.create(OUT_DIR_TBL, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# Helpers ----------------------------------------
# Chunked wrapper around GIFT_checklists_raw()
get_raw_checklists_by_list <- function(list_ids, taxon_name = "Tracheophyta",
                                       floristic_group = "all", chunk = 50) {
  chunks <- split(list_ids, ceiling(seq_along(list_ids) / chunk))
  map_dfr(chunks, ~ GIFT_checklists_raw(
    list_ID = .x,
    taxon_name = taxon_name,
    floristic_group = floristic_group,
    namesmatched = FALSE
  ))
}

# Safe multipolygon coercion
as_multipolygon <- function(sfc) {
  sfc %>%
    st_make_valid() %>%
    suppressWarnings(st_collection_extract("POLYGON")) %>%
    suppressWarnings(st_cast("MULTIPOLYGON", warn = FALSE))
}

# . -----------
# 0) Read inputs ---------------------------------
inat_obs <- fread(IN_INAT_OBS, show_col_types = FALSE)
obs_ranges <- st_read(GPKG_OBS, layer = LYR_OBS, quiet = TRUE)

# TEpis list 
tepis_csv <- fread(IN_TEPIS_CSV, header = T) 


# 1) Build TEpis / AEpis name sets ---------------------
## 1.1) TEpis from Zotz list -------
TEpis_names <- if (nrow(tepis_csv)) {
  tepis_csv %>% mutate(work_species = str_squish(str_to_lower(Species))) %>% pull(work_species)
} else character(0)

## 2.1) AEpis: species observed in iNat project -------
AEpis_names <- 
  inat_obs %>% 
  filter(taxon.rank%in%c("species", "variety", "subspecies", "hybrid")) %>% 
  filter(!is.na(`taxon.name`)) %>%
  mutate(work_species = str_squish(str_to_lower(taxon.name))) %>% 
  filter(!work_species%in%TEpis_names) %>% ###FILTERED OUT 
  pull(work_species) %>% 
  unique() 

# Map taxon_id -> work_species to align obs_ranges with AEpis_names
taxon_lookup <- inat_obs %>%
  filter(taxon.rank%in%c("species", "variety", "subspecies", "hybrid")) %>% 
  reframe(taxon_id = as.integer(taxon_id),
            work_species = str_squish(str_to_lower(`taxon.name`))) %>%
  distinct() %>%
  filter(!is.na(taxon_id), !is.na(work_species)) %>% 
  filter(!work_species%in%TEpis_names) ###FILTERED OUT 


# 2) GIFT: TDWG Level-3 country codes--------------------
rich_all <- GIFT_richness(taxon_name = "Tracheophyta")           # entity_ID, totals
regions  <- GIFT_regions()                                       # TDWG metadata
tdwg3    <- regions %>%
  filter(!is.na(TDWG_lvl3_ID), TDWG_lvl3_ID != 0) %>%
  distinct(entity_ID, geo_entity, TDWG_lvl3_ID)

# Join richness to TDWG L3, compute denominators
rich_tdwg3 <- rich_all %>%
  inner_join(tdwg3, by = "entity_ID") %>%
  transmute(
    entity_ID, geo_entity,
    total_all  = suppressWarnings(as.numeric(total)),
    native_all = suppressWarnings(as.numeric(native)),
    naturalized = suppressWarnings(as.numeric(naturalized))
  ) %>%
  mutate(
    total_all = dplyr::coalesce(
      total_all,
      { s <- rowSums(cbind(native_all, naturalized), na.rm = TRUE)
      s[is.na(native_all) & is.na(naturalized)] <- NA_real_
      s }
    )
  )

# TDWG shapes for intersections
shp_tdwg3 <- GIFT_shapes(entity_ID = tdwg3$entity_ID)

# 3) GIFT: species checklists -----------------------
lists_meta <- GIFT_lists() %>%
  filter(entity_ID %in% tdwg3$entity_ID, suit_geo == 1) %>%
  select(ref_ID, list_ID, entity_ID, geo_entity)

raw_chk <- get_raw_checklists_by_list(
  lists_meta$list_ID,
  taxon_name = "Tracheophyta",
  floristic_group = "all",
  chunk = 75
)

checklists_tdwg3 <- raw_chk %>%
  select(list_ID, species = work_species) %>%
  mutate(work_species = str_squish(str_to_lower(species))) %>%
  left_join(lists_meta, by = "list_ID") %>%
  distinct(entity_ID, work_species, .keep_all = FALSE) %>%
  select(entity_ID, work_species)

# 4) Range overlap -----------------
# compute where observed ranges of acc epis intersect botanical countries

# Prepare obs_ranges with species names and valid geometries
if (is.na(st_crs(obs_ranges))) st_crs(obs_ranges) <- 4326
obs_ranges <- st_transform(obs_ranges, 4326)

obs_ranges_named <- obs_ranges %>%
  right_join(taxon_lookup, by = "taxon_id") %>% ###CHANGED TO RIGHT JOIN TO EXCLUDE NON WORK SPECIES
  dplyr::filter(!is.na(work_species)) %>%
  dplyr::mutate(geometry = as_multipolygon(geom)) %>%
  dplyr::filter(!st_is_empty(geometry))

# Keep only AE species
obs_ae <- obs_ranges_named %>%
  filter(work_species %in% AEpis_names) %>%
  select(work_species, geometry)

# Intersect with TDWG polygons 
# old_s2 <- sf::sf_use_s2()
# on.exit(sf::sf_use_s2(old_s2), add = TRUE)
# sf::sf_use_s2(FALSE)

allowed_pairs <- st_join(
  obs_ae,
  shp_tdwg3 %>% select(entity_ID),
  join = st_intersects, left = FALSE
) %>%
  st_drop_geometry() %>%
  distinct(work_species, entity_ID)

# 5) Calculate counts & quotients per botanical country ---------------------
TEpis_per_region <- checklists_tdwg3 %>%
  filter(work_species %in% TEpis_names) %>%
  count(entity_ID, name = "n_TEpis")

AEpis_per_region <- checklists_tdwg3 %>%
  filter(work_species %in% AEpis_names) %>%
  inner_join(allowed_pairs, by = c("work_species", "entity_ID")) %>%
  count(entity_ID, name = "n_AEpis")

quotients <- rich_tdwg3 %>%
  left_join(TEpis_per_region, by = "entity_ID") %>%
  left_join(AEpis_per_region, by = "entity_ID") %>%
  mutate(
    n_TEpis = coalesce(n_TEpis, 0L),
    n_AEpis = coalesce(n_AEpis, 0L),
    epiphyte_quotient_all   = if_else(total_all  > 0, n_TEpis / total_all,  NA_real_),
    accidental_quotient_all = if_else(total_all  > 0, n_AEpis / total_all,  NA_real_),
    epiphyte_quotient_native   = if_else(native_all > 0, n_TEpis / native_all, NA_real_),
    accidental_quotient_native = if_else(native_all > 0, n_AEpis / native_all, NA_real_),
    divergence_all   = epiphyte_quotient_all - accidental_quotient_all,
    divergence_native= epiphyte_quotient_native - accidental_quotient_native
  )

# Save table
out_tbl <- quotients %>%
  arrange(desc(epiphyte_quotient_all)) %>%
  select(entity_ID, geo_entity, n_TEpis, n_AEpis,
         total_all, native_all,
         epiphyte_quotient_all, accidental_quotient_all,
         divergence_all)

readr::write_csv(out_tbl, file.path(OUT_DIR_TBL, "tdwg_quotients.csv"))



# . -----------
# 6) Maps (facets + divergence) ------------------------
map_df <- shp_tdwg3 %>%
  left_join(quotients) %>%
  select(geo_entity,
         epiphyte_quotient_all, accidental_quotient_all,
         divergence_all,
         geometry) %>%
  st_as_sf()

## 6.1) Faceted A+E quotient maps -----------
map_long <- map_df %>%
  tidyr::pivot_longer(
    cols = c(epiphyte_quotient_all, accidental_quotient_all),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = recode(
    metric,
    epiphyte_quotient_all   = "Epiphyte quotient (TEpis / total flora)",
    accidental_quotient_all = "Accidental epiphyte quotient (AEpis / total flora)"
  ))

p_quot <- ggplot(map_long) +
  geom_sf(aes(fill = value), color = NA) +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "grey90",
    labels = label_percent(accuracy = 0.1),
    name = "Quotient"
  ) +
  facet_wrap(~ metric, ncol = 2) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 9, face = "bold")
  )

ggsave(file.path(OUT_DIR_FIG, "FILTERED_tdwg_quotients_facets.jpg"),
       p_quot, width = 12, height = 6.5, dpi = 600, bg = "white")

## 6.2) Divergence map (TEpis - AEpis) ----------
max_abs <- max(abs(map_df$divergence_all), na.rm = TRUE)
p_div <- ggplot(map_df) +
  geom_sf(aes(fill = divergence_all), color = NA) +
  scale_fill_distiller(
    palette = "BrBG", direction = -1,
    limits = c(-max_abs, max_abs),
    labels = label_percent(accuracy = 0.1),
    na.value = "grey90",
    name = "\u0394 quotient (TEpis − AEpis)" #DELTA
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank())

ggsave(file.path(OUT_DIR_FIG, "FILTERED_tdwg_divergence.jpg"),
       p_div, width = 10, height = 6, dpi = 600, bg = "white")


## 6.3) Combine and save maps ----------
p_eq_combined <- (p_quot / p_div) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1),
              guides = "collect") & theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2))
  

ggsave(file.path(OUT_DIR_FIG, "FILTERED_tdwg_eq_combined.jpg"),
       p_eq_combined, width = 14, height = 8, dpi = 600, bg = "white")

# Save workspace to a file
save.image(file = paste0(".Rdata/08_epiphyte_quotient", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))



