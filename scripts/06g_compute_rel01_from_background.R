# 06g_compute_rel01_from_background.R
#calculate max and min from each environmental variable per species
#turn to relative values

# 0) Load---------
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

IN_DAT <- "outputs/06_bayesian_model_relative_niche/model_data/use_availability_4chelsa.csv"
OUT_DIR <- "outputs/06_bayesian_model_relative_niche/model_data/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- read_csv(IN_DAT, show_col_types = FALSE) %>%
  mutate(
    taxon_id = as.character(taxon_id),
    species_label = as.factor(species_label)
  )


#convert to relative range-----------
dat2 <- dat %>%
  group_by(species_label) %>%
  summarise(
    bg_n_moist = sum(y == 0 & is.finite(chelsa_ai_1981_2010_v_2_1)),
    bg_min_moist = if (bg_n_moist > 0) min(chelsa_ai_1981_2010_v_2_1[y == 0], na.rm = TRUE) else NA_real_,
    bg_max_moist = if (bg_n_moist > 0) max(chelsa_ai_1981_2010_v_2_1[y == 0], na.rm = TRUE) else NA_real_,
    
    bg_n_maxT = sum(y == 0 & is.finite(chelsa_bio10)),
    bg_min_maxT = if (bg_n_maxT > 0) min(chelsa_bio10[y == 0], na.rm = TRUE) else NA_real_,
    bg_max_maxT = if (bg_n_maxT > 0) max(chelsa_bio10[y == 0], na.rm = TRUE) else NA_real_,
    
    bg_n_minT = sum(y == 0 & is.finite(chelsa_bio6)),
    bg_min_minT = if (bg_n_minT > 0) min(chelsa_bio6[y == 0], na.rm = TRUE) else NA_real_,
    bg_max_minT = if (bg_n_minT > 0) max(chelsa_bio6[y == 0], na.rm = TRUE) else NA_real_,
    
    bg_n_wind = sum(y == 0 & is.finite(chelsa_sfc_wind_mean_1981_2010_v_2_1)),
    bg_min_wind = if (bg_n_wind > 0) min(chelsa_sfc_wind_mean_1981_2010_v_2_1[y == 0], na.rm = TRUE) else NA_real_,
    bg_max_wind = if (bg_n_wind > 0) max(chelsa_sfc_wind_mean_1981_2010_v_2_1[y == 0], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) 

dat3 <-
  dat2 %>% 
  right_join(dat, by = "species_label") %>%
  group_by(species_label,taxon_id) %>% 
  reframe(y,lon,lat,
    moist_rel01 = case_when(
      bg_n_moist == 0 ~ NA,
      !is.finite(bg_max_moist - bg_min_moist) ~ NA,
      (bg_max_moist - bg_min_moist) <= 0 ~ NA,
      TRUE ~ (chelsa_ai_1981_2010_v_2_1 - bg_min_moist) / (bg_max_moist - bg_min_moist)
    ),
    maxT_rel01 = case_when(
      bg_n_maxT == 0 ~ NA,
      !is.finite(bg_max_maxT - bg_min_maxT) ~ NA,
      (bg_max_maxT - bg_min_maxT) <= 0 ~ NA,
      TRUE ~ (chelsa_bio10 - bg_min_maxT) / (bg_max_maxT - bg_min_maxT)
    ),
    minT_rel01 = case_when(
      bg_n_minT == 0 ~ NA,
      !is.finite(bg_max_minT - bg_min_minT) ~ NA,
      (bg_max_minT - bg_min_minT) <= 0 ~ NA,
      TRUE ~ (chelsa_bio6 - bg_min_minT) / (bg_max_minT - bg_min_minT)
    ),
    wind_rel01 = case_when(
      bg_n_wind == 0 ~ NA,
      !is.finite(bg_max_wind - bg_min_wind) ~ NA,
      (bg_max_wind - bg_min_wind) <= 0 ~ NA,
      TRUE ~ (chelsa_sfc_wind_mean_1981_2010_v_2_1 - bg_min_wind) / (bg_max_wind - bg_min_wind)
    )
  )


any(is.na(dat3))

write_csv(dat3, file.path(OUT_DIR, "use_availability_rel01.csv"))
message("DONE: ", file.path(OUT_DIR, "use_availability_rel01.csv"))
