#11_niche_model_simple.R
#Find out where accidental epiphytes are likely to occur 
#based on environmental predictors identified in step 4

# Impute sfcWind_max by nearest (lat, lon) neighbor (great-circle distance)


#Poster Caption:
  #Posterior estimates (mean ± 95% credible intervals) of environmental 
#predictors and their influence on accidental epiphyte presence. 
#Posterior estimates summarize the most likely effect sizes given 
#the data and model assumptions. The x-axis shows effect size on 
#the log-odds scale: values > 0 increase the probability of presence, 
#values < 0 decrease it, and intervals overlapping zero indicate weak or 
#uncertain effects. Effect sizes can be converted to odds ratios by 
#exponentiation (e.g., +2 ≈ 7× higher odds of presence; −2 ≈ 7× lower odds).
#Ecological significance:
  #Accidental epiphytes are favored by mild winter temperatures and 
#dense vegetation, while dry air, strong winds, and island isolation 
#reduce suitability—indicating reliance on sheltered, humid microhabitats. 
#In contrast, true epiphytes occupy a more specialized niche, adapted to 
#canopy conditions with traits for water and nutrient capture. These 
#patterns suggest that accidental epiphytes may represent an evolutionary 
#stepping stone toward obligate epiphytism under persistent selection 
#for arboreal life.

#Species occurrence is strongly favored by mild winter temperatures and dense 
#vegetation (high LAI), while dry air (high vapour pressure deficit), 
#strong winds, and island isolation reduce suitability. 
#This suggests the species thrives in sheltered, humid, and vegetated habitats.
rm(list=ls())

# 0) Load packages ---------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(brms)           # Bayesian GLMM
  library(bayesplot)      # Posterior visualization
  library(loo)            # Model comparison and LOO
  library(terra)          # Raster handling
  library(sf)             # Spatial operations
  library(future.apply)   # Parallel predictions
  library(viridisLite)
  library(patchwork)
  library(ggrepel)
  library(loo)
})

# Define theme
my_theme14 <- theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top"
  )

lookup_table <-fread("outputs/04_environmental/var_lookup.csv")

set.seed(2025)

# Load and prepare data ---------------------
## Load site-level effort data
my_sites <- fread("data/processed/my.sites.csv") %>% as_tibble() %>% 
  mutate(area_ha=if_else(area_ha==0,mean(area_ha,na.rm=T),area_ha))


## Load environmental variables 
env_vars_sel <- fread("outputs/04_environmental/env_vars_sel.csv") %>%
  as_tibble() %>%
  mutate(taxon.name = if_else(taxon.name %in% c("", NA), "no_taxon", taxon.name),
         LAI = if_else(is.na(LAI), 0, LAI),
         is_forest=if_else(is.na(forest_age),0,1)
  )



inat <- fread("data/processed/inat_observations.csv") %>% as_tibble() %>% 
  reframe(id=as.character(id),user_login, dbh=DBH) %>% 
  distinct()

host<-fread("data/processed/tree_genus_clean.csv") %>% 
  reframe(id=as.character(id),
          host_tree=tree_genus) 

#########___________###########


## Join with site-level effort data and inat user_login, dbh, and host genus---------
env_joined <- env_vars_sel %>%
  left_join(my_sites %>% group_by(site) %>% 
              reframe(NumberTrees=sum(NumberTrees),
                      area_ha) %>% 
              unique) %>%
  mutate(source=case_when(str_detect(site,"S")~"field",
                          str_detect(site,"B")~"background",
                          .default="iNat"),
         id=as.character(id)) %>%
  mutate(effort = if_else(source == "field", NumberTrees / area_ha, 0),
         log_effort = if_else(source == "field", log(effort + 1), 0)) %>% 
  left_join(inat) %>% 
  left_join(host) %>% 
  mutate(site=if_else(!str_detect(site,"S"),NA,site))

make_chelsa_short <- function(nms) {
  # 1) Remove leading "CHELSA_"
  out <- sub("^CHELSA_", "", nms)
  
  # 2) Remove suffixes like "_1981-2010_V.2.1" (including the preceding underscore)
  #    This pattern removes any trailing underscore + digits-digits + anything
  #    e.g., "_1981-2010_V.2.1", "_1981-2010", "_1979-2018_V3", etc.
  out <- sub("_[0-9]{4}-[0-9]{4}.*$", "", out)
  
  # Return cleaned names
  out
}

# Apply to a data.frame in place (safe: only alters names)
rename_chelsa_cols <- function(df) {
  old <- names(df)
  new <- make_chelsa_short(old)
  names(df) <- new
  df
}

lookup_table<- lookup_table %>% 
  pivot_wider(names_from = var,values_from = label) %>% 
  rename_chelsa_cols() %>% 
  pivot_longer(everything(),names_to = "var",values_to="label")


env_joined <- env_joined %>% rename_chelsa_cols() 

# Explanation:
# - effort = NumberTrees / area_ha as proxy for sampling intensity
# - log_effort ensures offset is on log scale (required for GLMM)
# - For iNat and background points, offset = 0


# Bayesian 1.0--------------


# Given the data (iNat + some field sites), these additions could be included:

# 1) Keep access_min as a fixed effect (proxy for iNat sampling bias). 
# Do not put it into an offset together with log_effort. Offsets are for known 
# multiplicative exposure. access_min should be estimated as a coefficient.

# 2) Random intercept for site 
# To avoid the giant "no_site" level overwhelming others, the cleanest tactic 
# is to exclude it from the RE via a by-smooth random effect (s(..., bs="re",
# by=...)) so that "no_site" gets 0 contribution, while field sites get partial 
# pooling.

# 3) Observer random intercept (1 | user_login) (optional).
# Useful if many observers have >1 contribution; if thousands of singletons, 
# it increases cost with limited gain. Lump the tail to "other_user".

# 4) Spatial residual term: start with thin‑plate smooth s(lon, lat, k = 80–150) 
# for speed and stability.
# 
# Consider gp(lon, lat) later for a GP kernel; it’s heavier but 
# principled. Thin‑plate splines are usually the best first step in brms 
# for large binary datasets.
# 
# Combine (2) + (4): s(lon,lat) captures broad spatial structure; 
# (re site) captures site‑specific deviations beyond that.


## 1.1) Data prep --------
dat <- env_joined %>%
  # --- Select relevant columns ---
  select(presence, site, lat, lon, log_effort, source, user_login,
         is_island, is_forest, host_tree, dbh, lat_cdeg,
         bio7:bio15, LAI, npp, vpd_max, sfcWind_max, access_min, forest_age) %>%
  
  # --- Factor cleanup and lumping ---
  mutate(
    site       = fct_na_value_to_level(as.factor(site), level = "no_site") |> fct_drop(),
    user_login = fct_na_value_to_level(as.factor(user_login), level = "no_user") |> fct_drop(),
    host_tree  = fct_na_value_to_level(as.factor(host_tree), level = "Unknown"),
    host_tree_lumped = fct_lump_min(host_tree, min = 5, other_level = "other_host") |> fct_drop(),
    user_login_lumped = fct_lump_min(user_login, min = 3, other_level = "other_user") |> fct_drop(),
    is_island  = fct_relevel(as.factor(is_island), "0"),
    
    # --- Indicators ---
    is_field   = as.numeric(site != "no_site"),
    is_forest  = as.numeric(is.finite(forest_age)),
    
    # --- Source type ---
    source_type = case_when(
      site != "no_site"                         ~ "field",
      as.character(user_login) == "no_user"     ~ "background",
      TRUE                                      ~ "inat"
    )  %>%  factor(levels = c("field","inat","background")),
    
    is_inat        = as.numeric(source_type == "inat"),
    is_background  = as.numeric(source_type == "background"),
    
    # --- Top users ---
    user_top = case_when(
      user_login_lumped %in% c("gz_uol","marie-ho","nadjakrebs","spencer_") ~ as.character(user_login_lumped),
      TRUE ~ "other"
    )  %>%  factor(levels = c("other","gz_uol","marie-ho","nadjakrebs","spencer_")),
    
    is_inat_nontop = as.numeric(source_type == "inat" & user_top == "other"),
    
    # --- Continuous covariates ---
    across(c(bio7, npp, bio3, bio6, vpd_max, sfcWind_max, bio15, LAI, lat_cdeg,
             access_min, forest_age, lon, lat), as.numeric),
    
    # --- Effort offset ---
    log_effort = if_else(is.finite(log_effort), log_effort, 0),
    
    # --- Forest age fill for non-forest rows ---
    forest_age_fill = if_else(is_forest == 1, forest_age, 0)
  ) %>%
  
  # --- Scale numeric covariates (except lon/lat if used for spatial smooth) ---
  mutate(
    across(c(dbh, forest_age, bio7, npp, bio3, bio6, vpd_max, sfcWind_max, bio15, LAI, lat_cdeg),
           ~ ifelse(is.finite(.x), as.numeric(scale(.x)), .x))
  )


# Sanity checks
# Check the model frame actually used (to verify row count)
nrow(model.frame(fit_s_re))

# Check for remaining NAs in crucial terms
anyNA(dat$site); anyNA(dat$user_login_lumped); anyNA(dat$forest_age_fill)
anyNA(dat$lon); anyNA(dat$lat); anyNA(dat$log_effort)

#Write
fwrite(dat,"outputs/11_niche_model/model_data.csv")

# Priors: remove 'sd' (no grouped REs left), keep 'sds' for s(..., bs="re") and other splines
pri <- c(
  prior(normal(0, 1), class = "b"),          # fixed effects
  prior(normal(0, 5), class = "Intercept"),  # intercept
  prior(exponential(1.5), class = "sds")     # s(site,...), s(lon,lat), s(user...), s(forest_age,...)
  #prior(exponential(1), class = "sds", coef = "saccess_min*")  # tighten smooths
  )

### Formula---------
#use forest_age_fill, add source_type, top-user fixed effects, and spline RE for non-top iNat users
form <- bf(
  presence ~
    # Fixed effects (standardized)
    bio7 + npp + bio3  + bio6 + vpd_max + sfcWind_max +
    bio15 + LAI + lat_cdeg + is_island + access_min +
    
    # Source type (baseline = "field")
    #source_type +
    
    # Forest-only smooth on forest_age + level shift
    s(forest_age_fill, by = is_forest, bs = "tp") +
    is_forest +
    
    
    # # Access effect: baseline smooth + iNat-specific deviation
    # s(access_min, bs = "tp", k = 10) +                          # baseline smooth (all sources)
    # s(access_min, by = I(source == "inat"), bs = "tp", k = 10) +# extra wiggle only for iNat
    
   # source +
    
    
    # Site RE (field rows only)
    s(site, bs = "re", by = is_field) +
    
    # Observer bias:
    # - Fixed effects for top users (gz_uol, marie-ho, nadjakrebs, spencer_), baseline = "other"
    #user_top +
    # - Random-effect spline for long tail of iNat observers only
    #s(user_login_lumped, bs = "re", by = is_inat_nontop) +
    
    # Spatial residual (thin-plate regression spline)
    #s(lon, lat, bs = "tp", k = 60) +
    
    # Known effort
    offset(log_effort)
)

# 
# library(dplyr)
# library(geosphere)   # distHaversine
# 
# lat_col <- "lat"; lon_col <- "lon"; val_col <- "sfcWind_max"
# df<-dat
# 
# # Split known vs missing (base R to avoid tidy-eval hiccups)
# known_idx <- which(!is.na(df[[val_col]]))
# miss_idx  <- which(is.na(df[[val_col]]))
# 
# if (length(miss_idx) > 0 && length(known_idx) > 0) {
#   # Build matrices in the REQUIRED order: (lon, lat)
#   known_xy <- cbind(lon = df[[lon_col]][known_idx],
#                     lat = df[[lat_col]][known_idx])
#   miss_xy  <- cbind(lon = df[[lon_col]][miss_idx],
#                     lat = df[[lat_col]][miss_idx])
#   
#   # Great-circle distances (meters)
#   D <- distm(miss_xy, known_xy, fun = distHaversine)
#   
#   # Nearest neighbor index per missing row
#   nn_idx_in_known <- max.col(-D)  # argmin by flipping sign
#   
#   # Impute
#   df[[val_col]][miss_idx] <- df[[val_col]][known_idx][nn_idx_in_known]
# }
# 
# # df$sfcWind_max now filled by nearest lat–lon neighbor
# 
# df->dat

## 1.2) Run ----------------
options(mc.cores = parallel::detectCores())
fit_s_re <- brm(
  formula = form,
  data    = dat,
  family  = bernoulli(link = "logit"),
  prior   = pri,
  chains  = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 13),
  seed    = 1
)

saveRDS(fit_s_re, "outputs/11_niche_model/fit_s_re.rds")

summary(fit_s_re)




# Spatial blocking involves dividing the study area into spatial units (blocks) 
# and using them to structure model fitting or evaluation. It helps:
#   
# Prevent spatial autocorrelation bias: Nearby observations often share similar 
# environments, which can inflate model performance if not accounted for.
#
# Improve generalization: By evaluating model performance across spatial blocks, 
# you test how well the model predicts in new areas, not just where it was trained.
#
# Support spatial random effects: Blocking allows to include terms like 
# s(lon, lat) or (1 | site) to capture unmeasured spatial variation.
# 
# Effect on model & interpretation:
#   
# The model becomes more robust to spatial clustering.
# Coefficients reflect broader ecological patterns, not just local sampling density.
# Random effects (e.g., site-level) absorb local deviations, improving fixed-effect estimates.

# 1. Spatial smooth (s(lon, lat, bs = "tp", k = 100))
# 
# This term captures fine-scale spatial autocorrelation and local patterns that 
# arise from geography—such as microclimate, habitat fragmentation, or dispersal limitations.
# It accounts for spatial structure that isn’t explained by the environmental 
# predictors, reducing residual spatial dependence and improving inference validity.
# 
# 2. Centered latitude (lat_cdeg)
# 
# Latitude is a broad-scale ecological gradient strongly correlated with climate, 
# seasonality, and species pools.
# Including it explicitly allows to model large-scale latitudinal trends 
# in accidental epiphyte presence, which might not be fully captured by the 
# smooth (especially if the smooth is penalized or for more interpretable 
# coefficients for latitude).
# 
# Why both?
#   
# Only use the spatial smooth -> captures patterns but loses 
# interpretability for latitude effects.
# Only use latitude -> ignores local spatial structure and risks 
# inflated Type I errors due to spatial autocorrelation.
# Together, they let you:
#   
# Quantify the effect of latitude (important for ecological hypotheses).
# Control for spatial autocorrelation and local variation.

## 2) Posterior summaries ----------

### 2.1) LOO ----------

# LOO stands for Leave-One-Out Cross-Validation, a method to estimate a model’s 
# predictive accuracy by systematically leaving out each observation and 
# evaluating how well the model predicts it.
#
# In Bayesian models, loo() uses Pareto-smoothed importance sampling (PSIS) 
# to approximate LOO efficiently.
#
# What LOO gives you:
#   
# elpd_loo: Expected log predictive density (higher = better).
# looic: LOO Information Criterion (lower = better).
# p_loo: Effective number of parameters (model complexity).
# Pareto k diagnostics: Flags influential observations.
# 
# Comparing Models with LOO using loo_compare() to see which has 
#better predictive performance.

# Interpretation:
#   
# If fit_s_re has a higher elpd_loo and lower looic, it’s a better model.
#
# If the difference in ELPD is large and the standard error is small, 
# the improvement is statistically meaningful.
# 
# If Pareto k > 0.7, consider using reloo = TRUE to refit problematic points.
# 

# PSIS-LOO with moment matching improves stability when Pareto-k is high

# Check which priors apply
get_prior(formula = form, data = dat, family = bernoulli())

# LOO with moment matching (since your previous fit had some k>0.7)
loo_s_re <- loo(fit_s_re)
plot(loo_s_re)
write_csv(as_tibble(loo_s_re$estimates), "outputs/11_niche_model/loo_estimates.csv")

# Inspect site effects
#ranef(fit_s_re, robust = TRUE)$ssite %>% head()

# Check the forest-age smooth
plot(conditional_smooths(fit_s_re), ask = FALSE)


# If any Pareto k > 0.7, run loo(fit, reloo = TRUE) to robustify the estimate. 

### 2.2) Posterior summaries -------
# Partial dependence plots (marginal effects)
ce <- conditional_effects(fit_s_re, "bio7:npp")
ce_plot<-
  plot(ce, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey()
ggsave("outputs/11_niche_model/conditional_effects_bio7-npp.png", ce_plot, 
       width = 10, height = 8)


pp_check(fit_s_re)
### 2.3) Check collinearity---------

posterior_samples <- as_draws_df(fit_s_re)
# Pairwise posterior correlations among b_ terms
# Or compute prior VIFs on standardized X before fitting:
library(performance)
check_collinearity(lm(scale(presence) ~ bio15+bio6 + bio7 + bio3 , data = dat))
check_collinearity(lm(scale(presence) ~ forest_age + LAI , data = dat))
check_collinearity(lm(scale(presence) ~ is_island + sfcWind_max , data = dat))
check_collinearity(lm(scale(presence) ~ bio3 + npp, data = dat))

### 2.4) variable importance plot ---------------------
# Explanation:
# - Rank predictors by absolute posterior mean effect size
# - Include 95% credible intervals
# - This helps visualize which predictors have the strongest influence on occurrence


### 2.5) Fixed-effect summary -------------------------------------
posterior_summary <- as_tibble(summary(fit_s_re)$fixed, rownames = "variable") %>%
  filter(!str_detect(variable, "Intercept")) %>%
  rename(
    Q2.5  = `l-95% CI`,
    Q97.5 = `u-95% CI`
  )

###2.6) Build a single spatial "s(lon, lat)" magnitude from posterior draws ----
# brms names population-level coefs with "b_"
draws <- as_draws_df(fit_s_re)

# # Find the exact names brms used for the two linear spatial components
# sl_terms <- names(draws)[grepl("^b_slonlat_\\d+$", names(draws))]
# # Safety check: proceed only if both present
# if (length(sl_terms) >= 2) {
#   b1 <- draws[[sl_terms[1]]]
#   b2 <- draws[[sl_terms[2]]]
#   
#   # Euclidean magnitude of the two linear components
#   sl_mag <- sqrt(b1^2 + b2^2)
#   
#   spatial_row <- tibble(
#     variable     = "s(lon, lat)",          # label used in the plot
#     Estimate     = 0,                      # centered at 0 (no sign)
#     Q2.5         = quantile(sl_mag, 0.025),
#     Q97.5        = quantile(sl_mag, 0.975)
#   )
#   
#   # Remove the raw slonlat_1/2 rows from the fixed table
#   posterior_summary_clean <- posterior_summary %>%
#     filter(!str_detect(variable, "^slonlat_\\d+$")) %>%
#     # Add the spatial magnitude row; it will plot as a line from 0 to Q97.5
#     bind_rows(spatial_row) %>%
#     # For ordering, use absolute magnitude: for spatial, use Q97.5; for others, abs(Estimate)
#     mutate(abs_estimate = if_else(variable == "s(lon, lat)", 
#                                   as.numeric(Q97.5), 
#                                   abs(Estimate))) %>%
#     arrange(desc(abs_estimate))
# } else {
#   # If slonlat terms not found, just drop any that might appear and continue
#   posterior_summary_clean <- posterior_summary %>%
#     filter(!str_detect(variable, "^slonlat_\\d+$")) %>%
#     mutate(abs_estimate = abs(Estimate)) %>%
#     arrange(desc(abs_estimate))
#   warning("Did not find both 'b_slonlat_1' and 'b_slonlat_2' in draws; spatial magnitude row skipped.")
# }

### 3) Save table
write_csv(posterior_summary, "outputs/11_niche_model/variable_importance.csv")


# 4) Plot -------------
nice_col  <- viridis(6, option = "D")[4]   # pleasant mid-tone
line_col  <- "grey60"
grid_col  <- "grey90"


df_var_imp <- lookup_table %>% 
  add_row(var="is_island",label="Island") %>% 
  add_row(var="is_forest",label="Forest") %>% 
  add_row(var="sforest_age_fill:is_forest_1",label="Forest age (yrs)") %>% 
  right_join(posterior_summary, by=c("var"="variable")) %>% 
  drop_na(label)

var_imp_plot <-  
  ggplot(df_var_imp,
    aes(x = reorder(label, Estimate),
        y = Estimate)) +
  
  geom_hline(yintercept = 0, col=nice_col,alpha=0.3)+
 # For most vars: draw symmetric CI; for spatial row: show 0 -> Q97.5 as a one-sided range.
  geom_linerange(
    aes(ymin = if_else(var == "s(lon, lat)", 0, Q2.5),
        ymax = if_else(var == "s(lon, lat)", Q97.5, Q97.5)),
    size = 0.9, color = line_col
  ) +
  geom_point(
    #data = subset(posterior_summary, var != "s(lon, lat)"),
    size = 3, color = nice_col
  ) +
  # #Put a distinct marker for the spatial magnitude at its upper CI (visual cue of magnitude)
  # geom_point(
  #   data = subset(posterior_summary, var == "s(lon, lat)"),
  #   aes(y = Q97.5),
  #   size = 3, shape = 18, color = nice_col
  # ) +
  coord_flip() +
  labs(title = "Variable Importance (Posterior Estimates)",
       x = "Predictor",
       y = "Effect Size (log-odds or magnitude)") +
  theme_minimal(base_size = 13) +
  my_theme14

var_imp_plot
ggsave("outputs/11_niche_model/fit_s_re_variable_importance_plot.png",
       var_imp_plot, width = 7, height = 5, dpi = 300)

# Manuscript summary:
# "We modeled accidental epiphyte occurrence using a Bayesian GLMM with a logit link,
# including 12 environmental predictors, latitude, island status, and landcover as fixed effects,
# and site as a random intercept. Sampling effort was incorporated as an offset for field data.
# Posterior estimates and LOO cross-validation indicated that temperature range, NPP,
# and forest age were among the strongest predictors. Partial dependence plots revealed
# positive associations with forest age and LAI, and negative associations with wind speed."


