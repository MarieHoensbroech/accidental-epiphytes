#06f_build_use_availability_table.R
#merge presence and background points
#use: epiphyte presence=1
#availability: epiphyte presence=0 (pseudo-absence)


rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

IN_PRES <- "outputs/06_bayesian_model_relative_niche/environmental/inat_presences_with_4chelsa.csv"
IN_BG   <- "outputs/06_bayesian_model_relative_niche/environmental/background_points_with_4chelsa.csv"

OUT_DIR <- "outputs/06_bayesian_model_relative_niche/model_data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


bg <- read_csv(IN_BG, show_col_types = FALSE) %>%
  transmute(
    y = 0L,
    taxon_id = as.character(taxon_id),
    species_label,
    lon, lat,
    chelsa_ai_1981_2010_v_2_1,
    chelsa_bio10,
    chelsa_bio6,
    chelsa_sfc_wind_mean_1981_2010_v_2_1
  )


pres <- read_csv(IN_PRES, show_col_types = FALSE) %>%
  transmute(
    y = 1L,
    taxon_id = as.character(taxon.id),
    species_label=paste(word(taxon.name,1),word(taxon.name,2),sep="_"),
    lon, lat,
    chelsa_ai_1981_2010_v_2_1,
    chelsa_bio10,
    chelsa_bio6,
    chelsa_sfc_wind_mean_1981_2010_v_2_1
  ) %>% 
  filter(species_label%in%unique(bg$species_label))

dat <- bind_rows(pres, bg) %>%
  filter(
    is.finite(chelsa_ai_1981_2010_v_2_1),
    is.finite(chelsa_bio10),
    is.finite(chelsa_bio6),
    is.finite(chelsa_sfc_wind_mean_1981_2010_v_2_1)
  )

write_csv(dat, file.path(OUT_DIR, "use_availability_4chelsa.csv"))
message("DONE: ", file.path(OUT_DIR, "use_availability_4chelsa.csv"))
