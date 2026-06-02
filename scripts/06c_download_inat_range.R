# 06c_download_inat_ranges.R
# Download iNaturalist Open Range Map (Geomodel Expected Nearby Map) polygons
# for species selected in 06b (species_counts file)

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(purrr)
  library(stringr)
})

IN_SPECIES <- "outputs/06_bayesian_model_relative_niche/environmental/species_selected_counts_min20.csv"
OUT_DIR    <- "outputs/06_bayesian_model_relative_niche/inat_ranges"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Base URL described by iNat Open Range Map Dataset page:
# https://inaturalist-open-data.s3.us-east-1.amazonaws.com/geomodel/geojsons/latest/<taxon_id>.geojson
BASE_URL <- "https://inaturalist-open-data.s3.us-east-1.amazonaws.com/geomodel/geojsons/latest"

species <- read_csv(IN_SPECIES, show_col_types = FALSE) %>%
  mutate(taxon_id = as.character(taxon_id))

download_one <- function(taxon_id, species_label) {
  url <- sprintf("%s/%s.geojson", BASE_URL, taxon_id)
  out_file <- file.path(OUT_DIR, sprintf("%s_%s.gpkg",
                                         taxon_id,
                                         str_replace_all(species_label, "[^A-Za-z0-9]+", "_")))
  # # If already downloaded, skip
  # if (file.exists(out_file)) 
  #   {
  #   return(tibble(taxon_id, species_label, status="exists", file=out_file, url=url))
  # }
  
  # Try to read from URL
  res <- tryCatch({
    g <- suppressWarnings(st_read(url, quiet = TRUE))
    # Add identifiers if missing
    if (!("taxon_id" %in% names(g))) g$taxon_id <- taxon_id
    if (!("species_label" %in% names(g))) g$species_label <- species_label
    st_write(g, out_file, delete_dsn = TRUE, quiet = TRUE)
    tibble(taxon_id, species_label, status="downloaded", file=out_file, url=url)
  }, error = function(e) {
    tibble(taxon_id, species_label, status="missing_or_error",
           file=NA_character_, url=url, error=as.character(e$message))
  })
  
  res
}

log_tbl <- pmap_dfr(list(species$taxon_id, species$species_label), download_one)
write_csv(log_tbl, file.path(OUT_DIR, "download_log.csv"))

message("Done. See: ", file.path(OUT_DIR, "download_log.csv"))




#plot
# Plot downloaded iNaturalist open range polygons one by one
# with country background clipped to each range extent

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

IN_DIR  <- "outputs/06_bayesian_model_relative_niche/inat_ranges"
OUT_DIR <- file.path(IN_DIR,"range_plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# set TRUE to write PNGs to disk
SAVE_PNG <- TRUE

# small buffer around extent (in degrees; data are lon/lat)
PAD_X <- 1
PAD_Y <- 1

# read country background
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# list all downloaded range files
gpkg_files <- list.files(IN_DIR, pattern = "\\.gpkg$", full.names = TRUE)

for (f in gpkg_files) {
  
  message("Plotting: ", basename(f))
  
  # 1. read range file
  g <- st_read(f, quiet = TRUE) %>%
    st_make_valid()
  
  # 2. drop empty geometries
  g <- g %>%
    filter(!st_is_empty(geom))
  
  # skip if empty after filtering
  if (nrow(g) == 0) next
  
  # 3. make sure CRS exists
  if (is.na(st_crs(g))) {
    st_crs(g) <- 4326
  }
  
  # 4. match CRS to world layer
  g <- g %>%
    st_transform(st_crs(world))
  
  # 5. extract labels
  species_name <- if ("species_label" %in% names(g)) unique(g$species_label)[1] else tools::file_path_sans_ext(basename(f))
  taxon_id     <- if ("taxon_id" %in% names(g)) unique(g$taxon_id)[1] else NA_character_
  
  # 6. get bounding box and padded extent
  bb <- st_bbox(g)
  
  xlim <- c(bb["xmin"] - PAD_X, bb["xmax"] + PAD_X)
  ylim <- c(bb["ymin"] - PAD_Y, bb["ymax"] + PAD_Y)
  
  # 7. build bounding box polygon
  bbox_poly <- st_as_sfc(
    
    st_bbox(c(
      xmin = as.numeric(xlim[1]),
      xmax = as.numeric(xlim[2]),
      ymin = as.numeric(ylim[1]),
      ymax = as.numeric(ylim[2])
    ), crs = st_crs(g))
    
  )
  
  # 8. crop country background to that extent
  world_clip <- world %>%
    st_crop(st_bbox(bbox_poly))
  
  # 9. build plot
  p <- ggplot() +
    geom_sf(
      data = world_clip,
      fill = "grey92",
      colour = "grey55",
      linewidth = 0.2
    ) +
    geom_sf(
      data = g,
      fill = "#2C7FB8",
      colour = "#08306B",
      alpha = 0.45,
      linewidth = 0.4
    ) +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    labs(
      title = species_name,
      subtitle = if (!is.na(taxon_id)) paste("taxon_id:", taxon_id) else NULL,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold")
    )
  
  # 10. print to viewer
  print(p)
  
  # 11. optionally save
  if (SAVE_PNG) {
    out_png <- file.path(
      OUT_DIR,
      paste0(str_replace_all(species_name, "[^A-Za-z0-9]+", "_"), ".png")
    )
    
    ggsave(
      filename = out_png,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}
