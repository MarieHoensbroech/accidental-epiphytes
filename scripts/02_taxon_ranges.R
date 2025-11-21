# 02_ranges_from_inat_and_observations.R
# 
# 1) Download taxon ranges from iNat
# 2) Generate accidental epiphyte ranges from project observations
# 3) Calculate %overlap
# 4) Save to file

#rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(dplyr)
  library(purrr)
  library(httr)
  library(jsonlite)
  library(stringr)
  library(tidyr)
  library(units)
  library(igraph)
  library(data.table)
  library(future.apply)
  library(future)
  library(dbscan)
  library(progressr)
  library(lwgeom)

})

#Set up parallel plan
plan(multisession)  
options(future.globals.maxSize = 2e9) 
options(future.maxWorkers = parallel::detectCores() - 1)

# ------------------------------- Parameters -----------------------------------
INPUT_CSV        <- "data/processed/inat_observations.csv"
OUT_DIR_SPATIAL  <- "data/spatial"
GPKG_INAT        <- file.path(OUT_DIR_SPATIAL, "inat_ranges.gpkg")
GPKG_OBS         <- file.path(OUT_DIR_SPATIAL, "obs_ranges.gpkg")
GPKG_POINTS      <- file.path(OUT_DIR_SPATIAL, "obs_points.gpkg")
LOOKUP_CSV       <- file.path(OUT_DIR_SPATIAL, "range_lookup.csv")

# iNat geomodel export (public S3)
INAT_GEOMODEL_BASE <- "https://inaturalist-open-data.s3.us-east-1.amazonaws.com/geomodel/geojsons/latest"

# Observation-derived range builder parameters
CLUSTER_KM     <- getOption("ODR.CLUSTER_KM", 100)     # threshold for occurrence cluster distances
BUFFER_KM      <- getOption("ODR.BUFFER_KM", 100)      # radius around occurrences
LINE_WIDTH_KM  <- getOption("ODR.LINE_WIDTH_KM", 200)  # used for "linebuffer" strategy -- 100 km each side
HULL_METHOD    <- getOption("ODR.HULL_METHOD", "convex")   # "convex" or "concave" (concave is very computationally intensive)
CONCAVITY      <- getOption("ODR.CONCAVITY", 2)            # only if concave
GEOM_STRATEGY  <- getOption("ODR.GEOM_STRATEGY", "hull")   # "hull", "circles_union", or "linebuffer2"
#  - "hull":     for clusters with n>=3 observations -> convex/concave hull + buffer; 
#  - "circles_union": for n = 1 or close pairs -> circle 
#  - "linebuffer2": for n==2 within threshold -> line to connect points
CLIP_TO_LAND   <- getOption("ODR.CLIP_TO_LAND", TRUE)
OUT_CRS        <- getOption("ODR.OUT_CRS", 4326)
GPKG_OBS       <- getOption("ODR.GPKG_OBS", "data/spatial/obs_ranges.gpkg")
ODR_LOG        <- getOption("ODR.LOG", "data/spatial/odr_log.txt")

UA_STRING        <- "AccidentalEpiphyteProject"  #used by the iNat API for identifying devices 
RATE_LIMIT_S     <- 0.15   

dir.create(OUT_DIR_SPATIAL, showWarnings = FALSE, recursive = TRUE)

# Load occurrence data -----------------------------------

obs_tbl <- fread(INPUT_CSV)

# Minimal columns needed
required_cols <- c("id", "taxon_id", "taxon.name", "latitude", "longitude")

# Clean coordinates and build sf points (keep lon/lat columns)
obs_sf <- obs_tbl %>%
  filter(!is.na(longitude), !is.na(latitude),
         between(longitude, -180, 180), between(latitude, -90, 90)) %>%
  mutate(across(c(longitude, latitude), as.numeric)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

taxon_ids <- unique(obs_sf$taxon_id)

# 1) iNaturalist ranges -------------------------------
inat_geom_url <- function(taxon_id) {
  paste0(INAT_GEOMODEL_BASE, "/", taxon_id, ".geojson")
}

fetch_inat_range <- function(taxon_id) {
  url <- inat_geom_url(taxon_id)
  res <- tryCatch(
    GET(url, user_agent(UA_STRING), timeout(30)),
    error = function(e) NULL
  )
  if (!is.null(res) && identical(res$status_code, 200L)) {
    txt <- content(res, as = "text", encoding = "UTF-8")
    tf  <- tempfile(fileext = ".geojson")
    writeLines(txt, tf)
    sfobj <- tryCatch(st_read(tf, quiet = TRUE), error = function(e) NULL)
    if (!is.null(sfobj) && nrow(sfobj) > 0) {
      return(sf::st_geometry(sfobj)[[1]])
    }
  }
  st_geometrycollection()
}

taxon_ids <- sort(unique(na.omit(obs_sf$taxon_id)))

inat_sfg_list <- vector("list", length(taxon_ids))
names(inat_sfg_list) <- taxon_ids

for (i in seq_along(taxon_ids)) {
  id <- taxon_ids[i]
  message(sprintf("  [%d/%d] taxon_id=%s", i, length(taxon_ids), id))
  inat_sfg_list[[i]] <- fetch_inat_range(id)
  Sys.sleep(RATE_LIMIT_S)
}

# Convert list-of-sfg -> sfc
inat_sfc <- st_sfc(inat_sfg_list, crs = 4326)
inat_ranges_sf <- st_sf(taxon_id = taxon_ids, geometry = inat_sfc) %>%
  mutate(inat_range_wkt = st_as_text(geometry))

# 1b) Write iNat ranges ---------
if (file.exists(GPKG_INAT)) unlink(GPKG_INAT)
st_write(inat_ranges_sf, GPKG_INAT, layer = "inat_range", quiet = TRUE)

#inat_ranges_sf <- st_read(GPKG_INAT, layer = "inat_range", quiet = TRUE)

# . -------------------
# 2) Observation-derived ranges (ODR) ----------------------------------
sf::sf_use_s2(TRUE)

# Helpers

aeqd_from_centroid <- function(geom_ll) {
  cxy <- sf::st_coordinates(sf::st_centroid(sf::st_union(geom_ll)))
  lon0 <- cxy[1]; lat0 <- cxy[2]
  paste0("+proj=aeqd +lat_0=", lat0, " +lon_0=", lon0,
         " +datum=WGS84 +units=m +no_defs")
}

sanitize_polygons <- function(sfc) {
  if (length(sfc) == 0) return(sf::st_sfc())
  sfc <- suppressWarnings(sf::st_make_valid(sfc))
  sfc <- suppressWarnings(sf::st_collection_extract(sfc, "POLYGON"))
  sfc <- sfc[!sf::st_is_empty(sfc)]
  sfc
}

# Density-based clustering without dense distance matrix
cluster_by_distance_km <- function(pts_ll, cluster_km) {
  n <- nrow(pts_ll)
  if (n == 0) return(integer(0))
  
  aeqd  <- aeqd_from_centroid(pts_ll)
  xy    <- sf::st_coordinates(sf::st_transform(pts_ll, aeqd))
  cl    <- dbscan::dbscan(xy, eps = cluster_km * 1000, minPts = 1)
  labs  <- cl$cluster
  # dbscan uses 0 for noise; treat each as its own cluster id after max label
  if (any(labs == 0)) {
    labs0 <- which(labs == 0)
    labs[labs0] <- seq_len(length(labs0)) + max(labs)
  }
  as.integer(labs)
}

# Build a circle around a point 
# pick local equal-area for precise radius 
# return in out_crs
circle_around_point <- function(pt_ll, radius_km = BUFFER_KM, out_crs = OUT_CRS) {
  aeqd <- aeqd_from_centroid(pt_ll)
  pt_local <- sf::st_transform(pt_ll, aeqd)
  circ <- sf::st_buffer(sf::st_geometry(pt_local), dist = radius_km * 1000)
  circ <- sanitize_polygons(circ)
  sf::st_transform(circ, out_crs)
}

# Build geometry for a single cluster according to GEOM_STRATEGY.
# - If cluster size == 1: returns a circle around the point (always).
# - If any observation is "far" (> CLUSTER_KM) from all others, it becomes its own cluster -> circle.
build_cluster_geom <- function(cluster_pts_ll,
                               buffer_km      = BUFFER_KM,
                               line_width_km  = LINE_WIDTH_KM,
                               hull_method    = HULL_METHOD,
                               concavity      = CONCAVITY,
                               geom_strategy  = GEOM_STRATEGY,
                               out_crs        = OUT_CRS) {
  
  # Deduplicate coordinates to avoid zero-length artifacts
  xy   <- sf::st_coordinates(cluster_pts_ll)[, 1:2, drop = FALSE]
  keep <- !duplicated.data.frame(as.data.frame(round(xy, 7)))
  cluster_pts_ll <- cluster_pts_ll[keep, , drop = FALSE]
  
  n <- nrow(cluster_pts_ll)
  if (n == 0) return(sf::st_sfc())
  
  # Always draw circles for single observations
  if (n == 1) {
    return(circle_around_point(cluster_pts_ll, radius_km = buffer_km, out_crs = out_crs))
  }
  
  # For multi-point clusters, choose geometry strategy
  if (identical(geom_strategy, "circles_union")) {
    # Union of per-point circles (works for any n>=2)
    geoms <- lapply(seq_len(n), function(i) {
      circle_around_point(cluster_pts_ll[i, , drop = FALSE], radius_km = buffer_km, out_crs = out_crs)
    })
    geom <- do.call(c, geoms)
    return(sanitize_polygons(sf::st_make_valid(sf::st_union(geom))))
  }
  
  if (identical(geom_strategy, "linebuffer2") && n == 2) {
    # Two points within threshold -> line buffer between them
    aeqd      <- aeqd_from_centroid(cluster_pts_ll)
    pts_local <- sf::st_transform(cluster_pts_ll, aeqd)
    m    <- sf::st_coordinates(pts_local)[, 1:2, drop = FALSE]
    line <- sf::st_sfc(sf::st_linestring(m), crs = sf::st_crs(pts_local))
    geom <- sf::st_buffer(line, dist = (line_width_km * 1000) / 2,
                          endCapStyle = "ROUND", joinStyle = "ROUND")
    geom <- sanitize_polygons(geom)
    return(sf::st_transform(geom, out_crs))
  }
  
  # Default "hull" strategy:
  #   n == 2 (within threshold): line buffer
  #   n >= 3: convex/concave hull then buffer
  aeqd      <- aeqd_from_centroid(cluster_pts_ll)
  pts_local <- sf::st_transform(cluster_pts_ll, aeqd)
  
  if (n == 2) {
    m    <- sf::st_coordinates(pts_local)[, 1:2, drop = FALSE]
    line <- sf::st_sfc(sf::st_linestring(m), crs = sf::st_crs(pts_local))
    geom <- sf::st_buffer(line, dist = (line_width_km * 1000) / 2,
                          endCapStyle = "ROUND", joinStyle = "ROUND")
    geom <- sanitize_polygons(geom)
    return(sf::st_transform(geom, out_crs))
  }
  
  # n >= 3
  geom <- tryCatch({
    union_pts <- sf::st_union(sf::st_geometry(pts_local))
    hull <- if (identical(hull_method, "concave")) {
      if (requireNamespace("concaveman", quietly = TRUE)) {
        concaveman::concaveman(pts_local, concavity = concavity)
      } else if (requireNamespace("lwgeom", quietly = TRUE)) {
        lwgeom::st_concave_hull(union_pts, ratio = 1 / concavity)
      } else {
        sf::st_convex_hull(union_pts)
      }
    } else {
      sf::st_convex_hull(union_pts)  # default (fast)
    }
    sf::st_buffer(hull, dist = buffer_km * 1000)
  }, error = function(e) {
    # Fallback: union of per-point circles
    geoms <- lapply(seq_len(n), function(i) {
      circle_around_point(cluster_pts_ll[i, , drop = FALSE], radius_km = buffer_km, out_crs = out_crs)
    })
    do.call(c, geoms) |> sf::st_union()
  })
  
  geom <- sanitize_polygons(geom)
  sf::st_transform(geom, out_crs)
}

# Per‑taxon builder that expects a subset of points for that taxon only

make_obs_range_for_taxon_subset <- function(pts_taxon_ll,
                                            cluster_km     = CLUSTER_KM,
                                            buffer_km      = BUFFER_KM,
                                            line_width_km  = LINE_WIDTH_KM,
                                            hull_method    = HULL_METHOD,
                                            concavity      = CONCAVITY,
                                            geom_strategy  = GEOM_STRATEGY,
                                            clip_to_land   = CLIP_TO_LAND,
                                            land_union     = NULL,
                                            out_crs        = OUT_CRS) {
  if (nrow(pts_taxon_ll) == 0) return(sf::st_sfc())
  
  # Cluster by CLUSTER_KM
  comps <- cluster_by_distance_km(pts_taxon_ll, cluster_km = cluster_km)
  pts_taxon_ll$cluster_id <- as.integer(comps)
  
  geoms <- lapply(split(pts_taxon_ll, pts_taxon_ll$cluster_id), function(g) {
    build_cluster_geom(g,
                       buffer_km     = buffer_km,
                       line_width_km = line_width_km,
                       hull_method   = hull_method,
                       concavity     = concavity,
                       geom_strategy = geom_strategy,
                       out_crs       = out_crs)
  })
  rng <- do.call(c, geoms)
  if (length(rng) == 0) return(sf::st_sfc())
  
  # Clip to land 
  if (isTRUE(clip_to_land)) {
    if (is.null(land_union)) stop("clip_to_land=TRUE but land_union is NULL")
    land_u <- sf::st_transform(land_union, sf::st_crs(rng))
    rng    <- suppressWarnings(sf::st_intersection(rng, land_u))
    rng    <- sanitize_polygons(rng)
    if (length(rng) == 0) return(sf::st_sfc())
  }
  
  # One union + validity per taxon
  sf::st_make_valid(sf::st_union(rng))
}

# Main driver

build_odr <- function(obs_sf,
                      taxon_ids,
                      workers      = 6,
                      write_gpkg   = TRUE,
                      gpkg_path    = GPKG_OBS,
                      layer_name   = "obs_range",
                      clip_to_land = CLIP_TO_LAND,
                      land_union   = NULL,
                      log_path     = ODR_LOG,
                      geom_strategy= GEOM_STRATEGY) {
  
  stopifnot(inherits(obs_sf, "sf"))
  if (is.na(sf::st_crs(obs_sf))) stop("obs_sf must have a valid CRS.")
  if (!"taxon_id" %in% names(obs_sf)) stop("obs_sf must contain column 'taxon_id'.")
  
  # Keep only requested taxa
  obs_sf <- obs_sf[obs_sf$taxon_id %in% taxon_ids, , drop = FALSE]
  
  # Work in lon/lat inside the pipeline
  pts_ll <- if (sf::st_is_longlat(obs_sf)) obs_sf else sf::st_transform(obs_sf, 4326)
  
  # Preserve order of the provided taxon_ids
  fac <- factor(pts_ll$taxon_id, levels = taxon_ids)
  pts_by_taxon <- split(pts_ll, fac, drop = TRUE)
  
  # Prepare land union once if clipping is requested and none supplied
  if (isTRUE(clip_to_land) && is.null(land_union)) {
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
      stop("clip_to_land=TRUE requires 'rnaturalearth' or supply 'land_union'.")
    }
    land_union <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
      sf::st_make_valid() |>
      sf::st_transform(4326)
  }
  
  # Parallel plan
  future::plan(future::multisession, workers = workers)
  
  progressr::handlers(global = TRUE)
  p <- progressr::progressor(steps = length(taxon_ids))
  
  # Parallel-safe logging (collect in memory then write once)
  results <- future.apply::future_lapply(
    seq_along(taxon_ids),
    function(i) {
      id <- taxon_ids[i]
      pts_taxon <- pts_by_taxon[[as.character(id)]]
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      log_line <- sprintf("[%s] Processing taxon_id %s (%d of %d)", ts, id, i, length(taxon_ids))
      p(sprintf("taxon %s (%d/%d)", id, i, length(taxon_ids)))
      
      if (is.null(pts_taxon) || nrow(pts_taxon) == 0) {
        return(list(id = id, geom = sf::st_sfc(), log = log_line))
      }
      
      geom <- make_obs_range_for_taxon_subset(
        pts_taxon_ll  = pts_taxon,
        clip_to_land  = clip_to_land,
        land_union    = land_union,
        geom_strategy = geom_strategy
      )
      list(id = id, geom = geom, log = log_line)
    },
    future.seed     = TRUE,
    future.packages = c("sf","dplyr","units","dbscan","progressr","rnaturalearth","lwgeom","concaveman")
  )
  
  # Write a single log file once (no contention)
  log_lines <- vapply(results, function(x) x$log, character(1))
  if (!is.null(log_path)) writeLines(log_lines, con = log_path, sep = "\n")
  
  # Combine geometries (keep only non-empty; maintain taxon_ids order)
  geoms <- lapply(results, `[[`, "geom")
  keep  <- vapply(geoms, function(g) length(g) > 0, logical(1))
  if (!any(keep)) {
    return(sf::st_sf(taxon_id = taxon_ids[0], geometry = sf::st_sfc(crs = OUT_CRS)))
  }
  
  obs_sfc <- do.call(c, geoms[keep])
  sf::st_crs(obs_sfc) <- OUT_CRS
  taxon_kept <- taxon_ids[keep]
  
  obs_ranges_sf <- sf::st_sf(
    taxon_id = taxon_kept,
    geometry = obs_sfc
  ) |>
    dplyr::mutate(obs_range_wkt = sf::st_as_text(geometry))
  
  if (isTRUE(write_gpkg)) {
    if (file.exists(gpkg_path)) unlink(gpkg_path)
    sf::st_write(obs_ranges_sf, gpkg_path, layer = layer_name, quiet = TRUE)
  }
  
  message(sprintf("Wrote observation-derived ranges: %s (layer: %s)", gpkg_path, layer_name))
  invisible(obs_ranges_sf)
}

# Run
future::plan(future::multisession, workers = 6)
progressr::handlers("txtprogressbar")
obs_ranges_sf <- build_odr(
  obs_sf       = obs_sf,                  
  taxon_ids    = taxon_ids,             
  workers      = 6,
  write_gpkg   = TRUE,
  gpkg_path    = "data/spatial/obs_ranges.gpkg",
  clip_to_land = TRUE,
  geom_strategy= "hull"                   # or "circles_union"
)




# . ----------------------
# 3) Calculate % overlap ----------------------------------
WORK_CRS  <- 54009     # Mollweide (meters) 
GRID_SIZE <- 1         # desired snap grid in meters
HEAL_BUF  <- 0.25      # tiny healing buffer in meters for intersection (0 disables)

## Clean, transform, (snap), validate, union-by-taxon, area -----------
make_clean <- function(x, work_crs = WORK_CRS, grid_size = GRID_SIZE) {
  x %>%
    dplyr::filter(!sf::st_is_empty(geometry)) %>%
    sf::st_zm(drop = TRUE, what = "ZM") %>%
    sf::st_transform(work_crs) %>%
    st_snap_to_grid(size = grid_size) %>%     
    sf::st_make_valid() %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) %>%
    dplyr::select(taxon_id, geometry, area_m2)
}

## Intersection area: polygon-only  -----
robust_intersection_area <- function(a, b, tol = HEAL_BUF) {
  a <- a |> sf::st_zm(drop = TRUE, what = "ZM") |> sf::st_make_valid()
  b <- b |> sf::st_zm(drop = TRUE, what = "ZM") |> sf::st_make_valid()
  if (!is.null(tol) && tol > 0) {
    a <- sf::st_buffer(a, tol)
    b <- sf::st_buffer(b, tol)
  }
  inter <- tryCatch(sf::st_intersection(a, b), error = function(e) sf::st_sfc(crs = sf::st_crs(a)))
  if (length(inter) == 0 || sf::st_is_empty(inter) |> all(TRUE)) return(0)
  inter_poly <- tryCatch(sf::st_collection_extract(inter, "POLYGON"),
                         error = function(e) sf::st_sfc(crs = sf::st_crs(a)))
  if (length(inter_poly) == 0 || sf::st_is_empty(inter_poly) |> all(TRUE)) return(0)
  val <- suppressWarnings(sf::st_area(inter_poly) |> sum(na.rm = TRUE))
  as.numeric(val %||% 0)
}
`%||%` <- function(x, y) if (is.null(x) || anyNA(x)) y else x

## Quick ggplot for one taxon (iNat fill, Obs outline, Intersection) ---
plot_taxon_overlap <- function(taxon_id, inat_clean, obs_clean, title_prefix = "Taxon") {
  gi <- inat_clean %>% dplyr::filter(taxon_id == !!taxon_id)
  go <- obs_clean  %>% dplyr::filter(taxon_id == !!taxon_id)
  
  # Try intersection; may be empty or non-polygonal
  inter <- tryCatch(
    suppressWarnings(sf::st_intersection(gi, go)),
    error = function(e) sf::st_sfc(crs = sf::st_crs(gi))
  )
  
  # Extract polygonal parts if any; if not, keep empty sfc
  inter_poly <- tryCatch(
    suppressWarnings(sf::st_collection_extract(inter, "POLYGON")),
    error = function(e) sf::st_sfc(crs = sf::st_crs(gi))
  )
  
  # Robust emptiness checker that never returns NA
  has_polys <- FALSE
  if (length(inter_poly) > 0) {
    emp <- tryCatch(sf::st_is_empty(inter_poly), error = function(e) TRUE)
    # st_is_empty returns a logical vector; treat NA as TRUE (i.e. empty)
    has_polys <- any(isFALSE(emp))
  }
  
  ggplot() +
    geom_sf(data = gi, fill = "#3A86FF3A", color = "#3A86FF", linewidth = 0.4) +
    geom_sf(data = go, fill = NA,         color = "#FF006E", linewidth = 0.6) +
    { if (has_polys) geom_sf(data = inter_poly, fill = "#8338EC66", color = "#8338EC", linewidth = 0.4) else NULL } +
    theme_minimal() +
    labs(title = paste(title_prefix, taxon_id, "- iNat (blue) vs Observed (red), intersection (purple)"))
}

## Build cleaned layers ----------------------------------------------
inat_clean <- make_clean(inat_ranges_sf, work_crs = WORK_CRS, grid_size = GRID_SIZE)

# Keep only taxa that exist in iNat after cleaning
obs_clean  <- obs_ranges_sf %>%
  dplyr::filter(taxon_id %in% inat_clean$taxon_id) %>%
  make_clean(work_crs = WORK_CRS, grid_size = GRID_SIZE)

## Compute overlap -------------------------------------
range_overlap <- tibble::tibble(
  taxon_id         = inat_clean$taxon_id,
  area_inat_m2     = NA_real_,
  area_obs_m2      = NA_real_,
  area_overlap_m2  = NA_real_,
  pct_inat         = NA_real_,
  pct_obs          = NA_real_,
  status           = NA_character_
)

for (k in seq_len(nrow(range_overlap))) {
  tx <- range_overlap$taxon_id[k]
  gi_row <- inat_clean %>% dplyr::filter(taxon_id == tx)
  go_row <- obs_clean  %>% dplyr::filter(taxon_id == tx)
  
  have_inat <- nrow(gi_row) == 1
  have_obs  <- nrow(go_row) == 1
  
  ai <- if (have_inat) gi_row$area_m2[1] else 0
  ao <- if (have_obs)  go_row$area_m2[1] else 0
  
  range_overlap$area_inat_m2[k] <- ai
  range_overlap$area_obs_m2[k]  <- ao
  
  if (!have_inat && !have_obs) {
    range_overlap$status[k] <- "missing_both"
    range_overlap$area_overlap_m2[k] <- 0
    next
  }
  if (!have_inat) {
    range_overlap$status[k] <- "no_inat"
    range_overlap$area_overlap_m2[k] <- 0
    range_overlap$pct_inat[k] <- NA_real_
    range_overlap$pct_obs[k]  <- if (ao > 0) 0 else NA_real_
    next
  }
  if (!have_obs) {
    range_overlap$status[k] <- "no_obs"
    range_overlap$area_overlap_m2[k] <- 0
    range_overlap$pct_obs[k]  <- NA_real_
    range_overlap$pct_inat[k] <- if (ai > 0) 0 else NA_real_
    next
  }
  
  aol <- robust_intersection_area(gi_row$geometry[[1]], go_row$geometry[[1]], tol = HEAL_BUF)
  range_overlap$area_overlap_m2[k] <- aol
  range_overlap$pct_inat[k] <- if (ai > 0) 100 * aol / ai else NA_real_
  range_overlap$pct_obs[k]  <- if (ao > 0) 100 * aol / ao else NA_real_
  range_overlap$status[k]   <- if (aol > 0) "overlap" else "no_overlap"
}

## Inspect & Plot a random 'no_overlap' ------------------------------
missing_overlap <- range_overlap %>%
  dplyr::filter(status == "no_overlap") %>%
  dplyr::arrange(dplyr::desc(pct_inat))

if (nrow(missing_overlap) > 0) {
  taxon_to_plot <- missing_overlap %>% dplyr::slice_sample(n = 1) %>% dplyr::pull(taxon_id)
  message("Random 'no_overlap' taxon to inspect: ", taxon_to_plot)
  print(plot_taxon_overlap(taxon_to_plot, inat_clean, obs_clean))
} else {
  message("No 'no_overlap' taxa found; nothing to plot.")
}


# ______________------------
# 4) Join & save ----------------------------------
# Save the observation points
if (file.exists(GPKG_POINTS)) unlink(GPKG_POINTS)
st_write(obs_sf, GPKG_POINTS, layer = "observations", quiet = TRUE)

# Join WKT lookup 
lookup <- obs_ranges_sf %>%
  st_drop_geometry() %>%
  select(taxon_id, obs_range_wkt) %>%
  left_join(
    inat_ranges_sf %>%
      st_drop_geometry() %>%
      select(taxon_id, inat_range_wkt),
    by = "taxon_id"
  ) %>% as_tibble() %>% 
  left_join(range_overlap) %>% 
  mutate(status = if_else(is.na(status), "range_missing", status))

fwrite(lookup, LOOKUP_CSV)

# Attach WKT lookup to point table
obs_range_sf <- obs_sf %>%
  left_join(lookup)

st_write(obs_range_sf, "data/spatial/obs_range_sf.json", driver = "GeoJSON", delete_dsn = T)
