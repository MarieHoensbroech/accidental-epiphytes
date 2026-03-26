# 06b_species_env_relative_range.R
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table); library(janitor)
  library(stringr); library(readr)
})

# 0) Setup ---------------------------------------------------------------------
IN_ENV         <- "outputs/06_bayesian_model_relative_niche/environmental/inat_presences_with_4chelsa.csv"

OUT_DIR <- "outputs/06_bayesian_model_relative_niche/environmental/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

min_occ <- 20L

# Predictors
# Moisture-related (AI or CMI)
moisture_var <- "chelsa_ai_1981_2010_v_2_1"

# Temperature-related (two variables)
max_temp_var <- "chelsa_bio10"   # warmest quarter mean temp
min_temp_var <- "chelsa_bio6"    # coldest month min temp

# Wind-related
wind_var <- "chelsa_sfc_wind_mean_1981_2010_v_2_1"

# 1) Load ----------------------------------------------------------------------
env_tbl <- fread(IN_ENV) %>%
  as_tibble() %>%
  clean_names()


# 2) Species selection (min_occ, present in fieldwork) --------------------------

inat_species_pool <- env_tbl %>%
  #filter(quality_grade=="research") %>% 
  mutate(species_label=paste(word(taxon_name,1),word(taxon_name,2),sep="_")) %>% 
  filter(is.finite(moisture), is.finite(max_temp), is.finite(min_temp), is.finite(wind))

species_counts <- inat_species_pool %>%
  group_by(taxon_id, species_label) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  filter(n_obs >= min_occ) %>%
  arrange(desc(n_obs))

inat_species_sel <- inat_species_pool %>%
  semi_join(species_counts, by = c("taxon_id", "species_label"))

# 3) Compute within-species env scaling -----------------------------------------

inat_species_scaled <- inat_species_sel %>%
  group_by(taxon_id, species_label) %>%
  mutate(
    moist_min = min(moisture, na.rm = TRUE),
    moist_max = max(moisture, na.rm = TRUE),
    moist_rng = moist_max - moist_min,
    moist_rel01 = if_else(is.finite(moist_rng) & moist_rng > 0,
                          (moisture - moist_min) / moist_rng,
                          0.5),
    
    max_temp_min = min(max_temp, na.rm = TRUE),
    max_temp_max = max(max_temp, na.rm = TRUE),
    max_temp_rng = max_temp_max - max_temp_min,
    max_temp_rel01 = if_else(is.finite(max_temp_rng) & max_temp_rng > 0,
                             (max_temp - max_temp_min) / max_temp_rng,
                             0.5),
    
    min_temp_min = min(min_temp, na.rm = TRUE),
    min_temp_max = max(min_temp, na.rm = TRUE),
    min_temp_rng = min_temp_max - min_temp_min,
    min_temp_rel01 = if_else(is.finite(min_temp_rng) & min_temp_rng > 0,
                             (min_temp - min_temp_min) / min_temp_rng,
                             0.5),
    
    wind_min = min(wind, na.rm = TRUE),
    wind_max = max(wind, na.rm = TRUE),
    wind_rng = wind_max - wind_min,
    wind_rel01 = if_else(is.finite(wind_rng) & wind_rng > 0,
                         (wind - wind_min) / wind_rng,
                         0.5),
    
    moist_mean = mean(moisture, na.rm = TRUE),
    moist_sd   = sd(moisture, na.rm = TRUE),
    moist_z_within = if_else(is.finite(moist_sd) & moist_sd > 0,
                             (moisture - moist_mean) / moist_sd,
                             0),
    
    max_temp_mean = mean(max_temp, na.rm = TRUE),
    max_temp_sd   = sd(max_temp, na.rm = TRUE),
    max_temp_z_within = if_else(is.finite(max_temp_sd) & max_temp_sd > 0,
                                (max_temp - max_temp_mean) / max_temp_sd,
                                0),
    
    min_temp_mean = mean(min_temp, na.rm = TRUE),
    min_temp_sd   = sd(min_temp, na.rm = TRUE),
    min_temp_z_within = if_else(is.finite(min_temp_sd) & min_temp_sd > 0,
                                (min_temp - min_temp_mean) / min_temp_sd,
                                0),
    
    wind_mean = mean(wind, na.rm = TRUE),
    wind_sd   = sd(wind, na.rm = TRUE),
    wind_z_within = if_else(is.finite(wind_sd) & wind_sd > 0,
                            (wind - wind_mean) / wind_sd,
                            0)
  ) %>%
  ungroup()

species_summary <- inat_species_scaled %>%
  group_by(taxon_id, species_label) %>%
  summarise(
    n_obs = n(),
    moist_min = min(moisture, na.rm = TRUE),
    moist_max = max(moisture, na.rm = TRUE),
    max_temp_min = min(max_temp, na.rm = TRUE),
    max_temp_max = max(max_temp, na.rm = TRUE),
    min_temp_min = min(min_temp, na.rm = TRUE),
    min_temp_max = max(min_temp, na.rm = TRUE),
    wind_min  = min(wind, na.rm = TRUE),
    wind_max  = max(wind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_obs))


# 4) Export -------------------------------------------------------------------

write_csv(species_counts,  file.path(OUT_DIR, "species_selected_counts_min20.csv"))
write_csv(species_summary, file.path(OUT_DIR, "species_selected_summary_env_ranges.csv"))
write_csv(inat_species_scaled, file.path(OUT_DIR, "inat_species_selected_scaled_rows.csv"))


message("DONE. Outputs in: ", OUT_DIR)
