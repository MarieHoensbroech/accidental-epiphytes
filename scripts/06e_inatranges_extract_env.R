#06e_inatranges_extract_env.R
#extract environment from background points -> species' available/realised niche


rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(sf)
  library(terra)
  library(readr)
})

set.seed(42)

# ----------------------------- 1) Paths ---------------------------------------

IN_POINTS <-  "outputs/06_bayesian_model_relative_niche/inat_ranges/background_points.csv"

OUT_DIR  <- "outputs/06_bayesian_model_relative_niche/environmental"
OUT_FILE <- file.path(OUT_DIR, "background_points_with_4chelsa.csv")

DIR_ENV    <- "data/environmental"
DIR_CHELSA <- file.path(DIR_ENV, "CHELSA")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- 2) Load points ---------------------------------

pts_tbl <- read_csv(IN_POINTS, show_col_types = FALSE) %>%
  mutate(row_id__tmp = row_number())

pts_valid <- pts_tbl %>%
  drop_na(lon, lat)

pts_sf <- pts_valid %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

pts_v <- terra::vect(pts_sf)

# ----------------------------- 3) Find CHELSA files ---------------------------

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

# ----------------------------- 4) Load rasters --------------------------------

r_ai    <- terra::rast(f_ai)
r_bio10 <- terra::rast(f_bio10)
r_bio6  <- terra::rast(f_bio6)
r_wind  <- terra::rast(f_wind)

# ----------------------------- 5) Extract values ------------------------------

out_extract <- tibble(
  row_id__tmp = pts_valid$row_id__tmp,
  chelsa_ai_1981_2010_v_2_1 = terra::extract(r_ai, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_bio10 = terra::extract(r_bio10, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_bio6 = terra::extract(r_bio6, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric(),
  chelsa_sfc_wind_mean_1981_2010_v_2_1 = terra::extract(r_wind, pts_v, ID = FALSE) %>% .[[1]] %>% as.numeric()
)

# ----------------------------- 6) Join back -----------------------------------

out <- pts_tbl %>%
  left_join(out_extract, by = "row_id__tmp") %>%
  mutate(
    moisture = as.numeric(chelsa_ai_1981_2010_v_2_1),
    max_temp = as.numeric(chelsa_bio10),
    min_temp = as.numeric(chelsa_bio6),
    wind     = as.numeric(chelsa_sfc_wind_mean_1981_2010_v_2_1)
  ) %>%
  select(-row_id__tmp)

# ----------------------------- 7) Save ----------------------------------------

write_csv(out, OUT_FILE)
# 
# sanity <- tibble(
#   in_file    = IN_POINTS,
#   out_file   = OUT_FILE,
#   n_in       = nrow(pts_tbl),
#   n_out      = nrow(out),
#   n_na_lonlat = sum(is.na(pts_tbl$lon) | is.na(pts_tbl$lat)),
#   na_moist   = sum(!is.finite(out$moisture)),
#   na_maxT    = sum(!is.finite(out$max_temp)),
#   na_minT    = sum(!is.finite(out$min_temp)),
#   na_wind    = sum(!is.finite(out$wind))
# )
# 
# write_csv(sanity, file.path(OUT_DIR, "sanity.csv"))

message("DONE. Saved:")
message("  ", OUT_FILE)
#message("  ", file.path(OUT_DIR, "sanity.csv"))
