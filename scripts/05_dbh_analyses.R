# 05_dbh_analyses.R
#
# Assess DBH effects on accidental epiphyte observation count and richness
# 1) using traditional stats
# 2) using Bayesian
# 3) Plotting
# 4) DBH distribution of all iNat accidental epiphyte observations

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
})

# ------------------------------- Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"

OUT_DIR_DBH <- "outputs/figures/05_dbh_analyses"
dir.create(OUT_DIR_INAT,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_SITES, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 12))
update_geom_defaults("point", list(alpha = 0.8))

my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% unique()
inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% unique()
inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>% 
  drop_na(latitude,longitude)


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
  mutate(site=Site,
         type=Type,
    lon = as.numeric(SiteLon),
    lat = as.numeric(SiteLat),
    total_trees = as.numeric(NumberTrees),
    obs_count_dbh = suppressWarnings(as.numeric(n_epis_per_site_cat)),
    spp_count_dbh   = suppressWarnings(as.numeric(n_sp_per_site_cat)),
    obs_count = suppressWarnings(as.numeric(n_epis_per_site)),
    spp_count  = suppressWarnings(as.numeric(n_sp_per_site)),
    dbh_cat=TreeCat
  ) %>%
  group_by(site,dbh_cat) %>%
  reframe(
    total_trees = sum(NumberTrees, na.rm = TRUE),
    obs_count_dbh,spp_count_dbh,obs_count,spp_count,
    lon = first(lon),
    lat = first(lat),
    type,
    obs_per10 = 10 * obs_count_dbh / total_trees,
    spp_per10 = 10 * spp_count_dbh   / total_trees
  ) %>%
  filter(is.finite(lat), is.finite(lon), !(lat == 0 & lon == 0),
         dbh_cat!="0 - 10",
         !type==""|is.na(type)) %>%
  unique()


inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>% 
  drop_na(latitude,longitude)


# 1) KW-test and Dunn's post-hoc -------
## 1a) DBH~Observation count -----------
df_dbh %>% 
  group_by(type) %>% 
  kruskal_test(obs_per10~dbh_cat)
# dbh has significant influence on observation counts in forests, not in willows

df_dbh %>% 
  group_by(type) %>% 
  dunn_test(obs_per10~dbh_cat)
# only between dbhs 10-30 and 90+


## 1b) DBH~species count
df_dbh %>% 
  group_by(type) %>% 
  kruskal_test(spp_per10~dbh_cat)
# dbh has significant influence on species counts in forests, not in willows

df_dbh %>% 
  group_by(type) %>% 
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
          type,
          
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
  formula = obs_count ~ type * dbh_cat + (1 | site) + offset(log_exposure),
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

# Compare observed vs predicted zero inflation
obs_zero_rate <- df_bayes %>% summarise(zero_rate = mean(obs_count == 0))
pred_draws <- posterior_predict(fit_obs, ndraws = 1000)
pred_zero_rate <- mean(rowMeans(pred_draws == 0))  # overall predicted zero rate

obs_zero_rate
pred_zero_rate


### Summary -------
# Model: Negative Binomial GLMM with log link and offset(log_exposure)
# Formula: obs_count ~ type * dbh_cat + (1 | site) + offset(log_exposure)

# Convergence:
# - All Rhat = 1.00; high ESS -> model converged well
# - shape ~ 0.57 -> strong overdispersion -> negative binomial appropriate
# - sd(Intercept) ~ 2.33 -> large site-to-site variation -> site as random factor good

# Reference levels:
# - Baseline: Forest type, DBH category "10 - 30"

# Fixed effects (log-rate scale):
# - typeWillow = 3.21 -> Willow has ~25× higher observation count than Forest at baseline DBH
# - dbh_cat effects (Forest only):
#   - 30 - 60: 1.05 -> ~3× increase
#   - 60 - 90: 2.16 -> ~9× increase
#   - 90+:     3.34 -> ~28× increase
# -> monotonic increase in accidental epiphyte count with DBH

# Interactions (type × dbh_cat):
# - All CIs include 0 -> no strong evidence of differing DBH effect by type

# Zero rate:
# - Observed: 48.4%; Predicted: 42.4%
# -> Model slightly underpredicts zeros 

# Overall:
# - there is strong evidence that Willow sites and larger trees 
#   host more  accidental epiphytes
#
# - DBH effect is strong and consistent; no clear interaction with type
# - Site-level random effects are substantial and necessary



## 2.2) Species counts -----------
fit_spp <- brm(
  formula = spp_count ~ type * dbh_cat + (1 | site) + offset(log_exposure),
  data    = df_bayes,
  family  = zero_inflated_negbinomial(),
  prior   = priors_nb,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  seed    = 43
)

ppc_spp <- pp_check(fit_spp)   
ppc_spp
spp_summary <- summary(fit_spp)
spp_summary

# Compare observed vs predicted zero inflation
pred_draws_spp <- posterior_predict(fit_spp, ndraws = 1000)
pred_zero_rate_spp <- mean(rowMeans(pred_draws_spp == 0))  # overall predicted zero rate

obs_zero_rate
pred_zero_rate_spp


### Summary -------
# Model: Zero-Inflated Negative Binomial GLMM
# Formula: spp_count ~ type * dbh_cat + (1 | site) + offset(log_exposure)

# Convergence:
# - All Rhat = 1.00; high ESS → model converged well
# - shape ~ 2.29 -> moderate overdispersion; NB appropriate
# - sd(Intercept) ~ 1.19 -> substantial site-level variation -> site as random effect good

# Zero inflation:
# - zi ~ 0.13 (95% CI: 0.01–0.30) -> excess zeros?
# - Observed zero rate: 48.4%; Predicted: 48.8% -> model fits zeros well

# Reference levels:
# - Baseline: Forest type, DBH category "10 - 30"

# Fixed effects (log-rate scale):
# - typeWillow = 3.75 -> Willow has ~43× higher species count than Forest at baseline DBH
# - dbh_cat effects (Forest only):
#   - 30 - 60: 0.85 -> ~2× increase
#   - 60 - 90: 2.06 -> ~8× increase
#   - 90+:     3.29 -> ~27 increase
# → Strong positive effect of DBH on species count

# Interactions (type × dbh_cat):
# - All CIs include 0 -> no strong evidence of differing DBH effect by type

# Overall:
# - Willow sites and larger trees host significantly more accidental epiphyte species
# - DBH effect is strong and consistent; no clear interaction with type
# - Zero-inflation component improves fit slightly
# - Site-level random effects are important and well-estimated








# . ----------
# 3) Plots ------
p_dbh_obs <- ggplot(df_dbh, aes(type, obs_count_dbh/total_trees, fill = dbh_cat)) +
  geom_boxplot(
    width = 0.7, outlier.shape = NA, alpha = 0.9, color = "grey20",
    linewidth = 1
  ) +
  scale_fill_viridis_d("DBH category (cm)", option = "D", end = 0.9) +
  scale_y_continuous(
    name = "Observation count",
    expand = expansion(mult = c(0.02, 0.12))  # extra top space so letters don't clip
  ) +
  labs(x = "Sampling site type") +
  my_theme14 +
  theme(legend.position = "top") 

p_dbh_sp <- ggplot(df_dbh, aes(type, spp_count_dbh/total_trees, fill = dbh_cat)) +
  geom_boxplot(
    width = 0.7, outlier.shape = NA, alpha = 0.9, color = "grey20",
    linewidth = 1
  ) +
  scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
  scale_y_continuous(
    name = "Species count",
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(x = "Sampling site type") +
  my_theme14 

# Combined plot
p_dbh_combo <- p_dbh_obs / p_dbh_sp + patchwork::plot_layout(axes = "collect")
p_dbh_combo

ggsave(file.path(OUT_DIR_DBH, "obs_and_species_per_DBH.jpg"),
       p_dbh_combo, dpi = 500, width = 9, height = 7)

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

ggsave(file.path(OUT_DIR_DBH, "individuals_vs_dbh.jpg"),
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

ggsave(file.path(OUT_DIR_DBH, "mean_individuals_vs_dbh.jpg"),
       p_mean_dbh, width = 7, height = 5, dpi = 300)




