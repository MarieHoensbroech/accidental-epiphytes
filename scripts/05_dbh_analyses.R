# 05_dbh_analyses.R
#
# Assess DBH effects on accidental epiphyte observation count and richness
# 1) using traditional stats
# 2) using Bayesian
# 3) Plotting
# 4) DBH distribution of all iNat accidental epiphyte observations
rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)     
  library(viridis)
  library(forcats)
  library(scales)
  library(maps)
  library(patchwork)
  library(stringr)
  library(data.table)
  library(cowplot)
  library(grid)
  library(ggimage)
  library(rnaturalearth)
  library(stringr)
  library(brms)
  library(tidybayes)
  library(emmeans)
  library(multcompView)
  library(rstatix)
  library(janitor)
})

# ------------------------------- Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"

OUT_DIR_DBH <- "outputs/05_dbh_analyses"
dir.create(OUT_DIR_DBH, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 12))
update_geom_defaults("point", list(alpha = 0.8))

my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% unique() %>% clean_names()
inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% unique() %>%  clean_names()
inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>% 
  drop_na(latitude,longitude) %>%  clean_names()


# Standard ggplot theme 

my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top"
  )

# . ----------
# Prep data -----------
df_dbh <- my_sites %>%
  mutate(
    lon = as.numeric(site_lon),
    lat = as.numeric(site_lat),
    total_trees = as.numeric(number_trees),
    obs_count_dbh = suppressWarnings(as.numeric(n_epis_per_site_cat)),
    spp_count_dbh   = suppressWarnings(as.numeric(n_sp_per_site_cat)),
    obs_count = suppressWarnings(as.numeric(n_epis_per_site)),
    spp_count  = suppressWarnings(as.numeric(n_sp_per_site)),
    dbh_cat=tree_cat
  ) %>%
  group_by(site,dbh_cat) %>%
  reframe(
    total_trees = sum(number_trees, na.rm = TRUE),
    obs_count_dbh,spp_count_dbh,obs_count,spp_count,
    lon = first(lon),
    lat = first(lat),
    setting,
    obs_per10 = 10 * obs_count_dbh / total_trees,
    spp_per10 = 10 * spp_count_dbh   / total_trees
  ) %>%
  filter(is.finite(lat), is.finite(lon), !(lat == 0 & lon == 0),
         dbh_cat!="0 - 10") %>%
  unique()


inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>% 
  drop_na(latitude,longitude)


# 1) KW-test and Dunn's post-hoc -------
## 1a) DBH~Observation count -----------
df_dbh %>% 
  group_by(setting) %>% 
  kruskal_test(obs_per10~dbh_cat)
# dbh has significant influence on observation counts in forests, not in willows

df_dbh %>% 
  group_by(setting) %>% 
  dunn_test(obs_per10~dbh_cat)
# only between dbhs 10-30 and 90+


## 1b) DBH~species count
df_dbh %>% 
  group_by(setting) %>% 
  kruskal_test(spp_per10~dbh_cat)
# dbh has significant influence on species counts in forests, not in willows

df_dbh %>% 
  group_by(setting) %>% 
  dunn_test(spp_per10~dbh_cat)
# only between dbhs 10-30 and 90+


# . ------------
# 2) Bayesian ------------
df_bayes <- df_dbh %>%
  filter(dbh_cat!="0 - 10",
         total_trees > 0) %>% 
  reframe(site,
          obs_count,
          spp_count,
          
          dbh_cat,
          setting,
          
          # Exposure offset on the log scale (adjusts for how many trees I sampled at a site)
          log_exposure = log(total_trees)
  )

priors_nb <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  set_prior("exponential(1)", class = "shape")  # negbin overdispersion
)


## 2.1) Observation counts -------------
fit_obs <- brm(
  formula = obs_count ~ setting * dbh_cat + (1 | site) + offset(log_exposure),
  data    = df_bayes,
  family  = negbinomial(),
  prior   = priors_nb,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  seed    = 42
)

# Posterior predictive checks
ppc_obs <- pp_check(fit_obs)     # visualize 
ppc_obs

summ_obs <- summary(fit_obs)     # summary
summ_obs


as_draws_df(fit_obs) %>%
  dplyr::select(dplyr::starts_with("b_")) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "bterm", values_to = "draw") %>%
  dplyr::mutate(term = stringr::str_remove(bterm, "^b_")) %>%
  dplyr::group_by(term) %>%
  dplyr::summarise(
    mean    = mean(draw),
    ci_low  = quantile(draw, 0.025),
    ci_high = quantile(draw, 0.975),
    pd      = as.numeric(bayestestR::p_direction(draw)),
    rope_perc = bayestestR::rope(draw, range = c(-log(1.1), log(1.1)))$ROPE_Percentage[1],
    rate_ratio_mean = exp(mean),
    rate_ratio_ci_l = exp(ci_low),
    rate_ratio_ci_h = exp(ci_high),
    .groups = "drop"
  ) %>%
  dplyr::arrange(term) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3))) %>%
  fwrite(file.path(OUT_DIR_DBH,"obs.count~dbh+setting_fixed_effects_summary.csv"))

# Observed zero rate
obs_z <- df_bayes %>%
  summarise(observed_zero_rate = mean(obs_count == 0, na.rm = TRUE))

# Predicted zero rate (overall)
pred_draws <- posterior_predict(fit_obs, ndraws = 1000)
pred_z <- tibble(predicted_zero_rate = mean(rowMeans(pred_draws == 0)))

# Bind and write
bind_cols(obs_z, pred_z) %>%
  data.table::fwrite(file.path(OUT_DIR_DBH,"obs.count~dbh+setting_zero_rates.csv"))

## 2.2) Species counts -----------
fit_spp <- brm(
  formula = spp_count ~ setting * dbh_cat + (1 | site) + offset(log_exposure),
  data    = df_bayes,
  family  = zero_inflated_negbinomial(),
  prior   = priors_nb,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  seed    = 43
)

as_draws_df(fit_spp) %>%
  dplyr::select(dplyr::starts_with("b_")) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "bterm", values_to = "draw") %>%
  dplyr::mutate(term = stringr::str_remove(bterm, "^b_")) %>%
  dplyr::group_by(term) %>%
  dplyr::summarise(
    mean    = mean(draw),
    ci_low  = quantile(draw, 0.025),
    ci_high = quantile(draw, 0.975),
    pd      = as.numeric(bayestestR::p_direction(draw)),
    rope_perc = bayestestR::rope(draw, range = c(-log(1.1), log(1.1)))$ROPE_Percentage[1],
    rate_ratio_mean = exp(mean),
    rate_ratio_ci_l = exp(ci_low),
    rate_ratio_ci_h = exp(ci_high),
    .groups = "drop"
  ) %>%
  dplyr::arrange(term) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3))) %>%
  fwrite(file.path(OUT_DIR_DBH,"spp.count~dbh+setting_fixed_effects_summary.csv"))

# Observed zero rate
obs_z <- df_bayes %>%
  summarise(observed_zero_rate = mean(spp_count == 0, na.rm = TRUE))

# Predicted zero rate (overall)
pred_draws <- posterior_predict(fit_spp, ndraws = 1000)
pred_z <- tibble(predicted_zero_rate = mean(rowMeans(pred_draws == 0)))

# Bind and write
bind_cols(obs_z, pred_z) %>%
  data.table::fwrite(file.path(OUT_DIR_DBH,"spp.count~dbh+setting_zero_rates.csv"))



# . ----------
# 3) Plots ------
p_dbh_obs <- ggplot(df_dbh, aes(setting, obs_count_dbh/total_trees, fill = dbh_cat)) +
  geom_boxplot(
    width = 0.7, outlier.shape = NA, alpha = 0.9, color = "grey20",
    linewidth = 1
  ) +
  scale_fill_viridis_d("DBH category (cm)", option = "D", end = 0.9) +
  scale_y_continuous(
    name = "Observation count",
    expand = expansion(mult = c(0.02, 0.12))  # extra top space so letters don't clip
  ) +
  labs(x = "Setting") +
  my_theme14 +
  theme(legend.position = "top") 

p_dbh_sp <- ggplot(df_dbh, aes(setting, spp_count_dbh/total_trees, fill = dbh_cat)) +
  geom_boxplot(
    width = 0.7, outlier.shape = NA, alpha = 0.9, color = "grey20",
    linewidth = 1
  ) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  scale_y_continuous(
    name = "Species count",
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(x = "Setting") +
  my_theme14 

# Combined plot
p_dbh_combo <- p_dbh_obs / p_dbh_sp + patchwork::plot_layout(axes = "collect")
p_dbh_combo

ggsave(file.path(OUT_DIR_DBH, "obs_and_species_per_DBH.jpg"),
       p_dbh_combo, dpi = 500, width = 6, height = 6)
ggsave(file.path(OUT_DIR_DBH, "obs_and_species_per_DBH.svg"),
       p_dbh_combo, dpi = 500, width = 6, height = 6)

## DBH structure heatmap by site -------------------
site_levels <- sort(unique(df_dbh$site))
midpoints <- c("10 - 30" = 20, "30 - 60" = 45, "60 - 90" = 75, "90+" = 100)

size_struct <- df_dbh %>%
  group_by(site, dbh_cat) %>%
  summarise(total_trees) %>%
  group_by(site) %>%
  mutate(prop = total_trees / sum(total_trees, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mid = midpoints[as.character(dbh_cat)]) %>%
  group_by(site) %>%
  mutate(weighted_mid = sum(prop * mid, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Site = fct_reorder(site, weighted_mid)) %>%
  filter(!is.na(mid))

p_heat <- ggplot(size_struct, aes(x = site, y = dbh_cat, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(labels = percent_format(accuracy = 1),
                       name = "Proportion of trees") +
  labs(x = "Site (ordered by mean DBH)", y = "DBH class (cm)",
       title = "Size-class structure by site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(OUT_DIR_DBH, "dbh_size_classes_by_site.jpeg"),
       p_heat, dpi = 500, width = 10, height = 8)







# . -------------
# 4) All iNat ------------
## 4.1) Observation count ~ DBH  ----------------------

inat_obs_dbh <- inat_obs %>%
  select(id, DBH_num, EpiCount_num, taxon.name) %>%
  filter(is.finite(DBH_num), is.finite(EpiCount_num))

# Sum individuals per DBH bin (0..8 m by 0.3 m)
bin_edges <- seq(0, 10, by = 0.3)
bin_labs  <- paste0(bin_edges[-length(bin_edges)], "–", bin_edges[-1], " m")

p_indiv_dbh <- inat_obs_dbh %>%
  mutate(DBH_bin = cut(DBH_num, breaks = bin_edges, include.lowest = TRUE, labels = bin_labs)) %>%
  group_by(DBH_bin) %>%
  summarise(
    sum_indiv = sum(EpiCount_num, na.rm = TRUE),
    n_tax     = n_distinct(taxon.name),
    n_obs     = n_distinct(id),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = DBH_bin, y = sum_indiv)) +
  geom_point(aes(size = n_tax), color = "#264653") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              color = "black", se = TRUE) +
  labs(
    x = "DBH bin (m)",
    y = "Number of observed individuals",
    title = "Abundance of individuals across DBH"
  ) +
  scale_size("No. of taxa")+
  my_theme14+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file.path(OUT_DIR_DBH, "ALL_INAT_individuals_vs_dbh.jpg"),
       p_indiv_dbh, width = 7, height = 5, dpi = 300)

## 4.2) Mean individuals per observation ~ DBH ---------
p_mean_dbh <- inat_obs_dbh %>%
  mutate(DBH_bin = cut(DBH_num, breaks = bin_edges, include.lowest = TRUE, labels = bin_labs)) %>%
  group_by(DBH_bin) %>%
  summarise(
    mean_indiv = mean(EpiCount_num, na.rm = TRUE),
    n_tax      = n_distinct(taxon.name),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = DBH_bin, y = mean_indiv)) +
  geom_point(aes(size = n_tax), color = "#1D3557") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              color = "black", se = TRUE) +
  labs(
    x = "DBH bin (m)",
    y = "Mean individuals per observation",
    title = "Mean observed individuals across DBH"
  ) +
  my_theme14+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR_DBH, "ALL_INAT_mean_individuals_vs_dbh.jpg"),
       p_mean_dbh, width = 7, height = 5, dpi = 300)




