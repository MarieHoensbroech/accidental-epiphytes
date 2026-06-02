# 06d_inat_ranges_background.R
# Generate 10,000 random background points within each iNaturalist open range polygon

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(purrr)
  library(stringr)
  library(tibble)
  library(magrittr)
})

set.seed(42)

IN_DIR  <- "outputs/06_bayesian_model_relative_niche/inat_ranges"
OUT_DIR <- "outputs/06_bayesian_model_relative_niche/inat_ranges"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_CSV  <- file.path(OUT_DIR, "background_points.csv")

N_POINTS <- 10000

range_files <- list.files(
  IN_DIR,
  pattern = "\\.gpkg$",
  full.names = TRUE
)


bg_sf <- range_files %>%
  map_dfr(
    ~ .x %>%
      st_read(quiet = TRUE) %>%
      st_make_valid() %>%
      summarise(geometry = st_union(geom)) %>%
      st_sample(size = N_POINTS, type = "random", exact = TRUE) %>%
      st_as_sf() %>%
      mutate(
        file_name = basename(.x),
        taxon_id = str_extract(file_name, "^[^_]+"),
        species_label = file_name %>%
          str_remove("^[^_]+_") %>%
          str_remove("\\.gpkg$"),
        point_id = row_number()
      )
  ) %>%
  select(taxon_id, species_label, point_id, geometry=x)


bg_sf <- bg_sf %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  relocate(taxon_id, species_label, point_id, lon, lat, geometry)

bg_csv <- bg_sf %>%
  st_drop_geometry() %>%
  mutate(
    bg_id = row_number()
  ) %>%
  select(bg_id, taxon_id, species_label, point_id, lon, lat)

bg_sf <- bg_sf %>%
  left_join(
    bg_csv %>% select(bg_id, taxon_id, species_label, point_id),
    by = c("taxon_id", "species_label", "point_id")
  ) %>%
  relocate(bg_id, taxon_id, species_label, point_id, lon, lat, geometry)

write_csv(bg_csv, OUT_CSV)

bg_csv %>%
  count(taxon_id, species_label, name = "n_background_points") %>%
  write_csv(file.path(OUT_DIR, "background_points_counts.csv"))

message("DONE. Saved:")
message("  ", OUT_CSV)
message("  ", file.path(OUT_DIR, "background_points_counts.csv"))

