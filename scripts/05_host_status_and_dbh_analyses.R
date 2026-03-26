# 04_host_status_and_dbh_analyses.R 
# Host size (DBH) effects on accidental epiphyte abundance and richness
# ---- 0) Setup & parameters ---------------------------------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(forcats)
  library(scales)
  library(patchwork)
  library(stringr)
  library(data.table)
  library(brms)
  library(emmeans)
  library(iNEXT)
  library(janitor)
  library(purrr)
  library(ggbreak)
})

IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
OUT_DIR_DBH   <- "outputs/05_host_status_and_dbh_analyses"
dir.create(OUT_DIR_DBH, showWarnings = FALSE, recursive = TRUE)

m_target <- 10L  # FIXED rarefaction target 

my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top",
    axis.line = element_blank()
  )

# ---- 1) Load & clean data ----------------------------------------------------
my_sites <- fread(IN_MY_SITES) %>% as_tibble() %>% clean_names()
inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% clean_names()

# ---- 2) Build site × DBH table (Script 1 original) --------------------------
df_dbh_raw <- my_sites %>%
  mutate(
    lon             = as.numeric(site_lon),
    lat             = as.numeric(site_lat),
    total_trees     = as.numeric(number_trees),
    obs_count_dbh   = suppressWarnings(as.numeric(n_epis_per_site_cat)),
    spp_count_dbh   = suppressWarnings(as.numeric(n_sp_per_site_cat)),
    obs_count       = suppressWarnings(as.numeric(n_epis_per_site)),
    spp_count       = suppressWarnings(as.numeric(n_sp_per_site)),
    dbh_cat         = tree_cat,
    setting         = setting
  ) %>%
  filter(is.finite(lat), is.finite(lon), dbh_cat != "0 - 10") %>%
  group_by(site, dbh_cat) %>%
  summarise(
    total_trees   = sum(total_trees,   na.rm = TRUE),
    obs_count_dbh = sum(obs_count_dbh, na.rm = TRUE),
    spp_count_dbh = sum(spp_count_dbh, na.rm = TRUE),
    obs_count     = first(obs_count,   na.rm = TRUE),
    spp_count     = first(spp_count,   na.rm = TRUE),
    lon           = first(lon), 
    lat           = first(lat),
    setting       = first(setting),
    .groups = "drop"
  ) %>%
  mutate(
    setting = fct_recode(setting,
                         "forest" = "forest",
                         "non-forest" = "non-forest") |> fct_drop(),
    dbh_cat = factor(dbh_cat, levels = c("10 - 30","30 - 60","60 - 90","90+"))
  )


# ---- 3) iNEXT BLOCK  -------------------------------------
# Standardise richness to 10 trees

# 3A. T_tbl 
T_tbl <- df_dbh_raw %>%
  transmute(site, dbh_cat, setting, T = as.integer(round(total_trees))) %>%
  filter(T > 0)

# 3B. incidence-frequency table, simplified
#     Count unique species occurrences per site × DBH
incidence_freq_tbl <- df_dbh_raw %>%
  filter(spp_count_dbh > 0) %>%        # only strata with species
  select(site, dbh_cat, setting, spp_count_dbh) %>%
  rowwise() %>%
  mutate(freq = list(rep(1L, spp_count_dbh))) %>%
  ungroup()

# 3C. Build inc_keys (size-based rarefaction, always m = 10)
inc_keys <- T_tbl %>%
  left_join(incidence_freq_tbl %>% select(site, dbh_cat, setting, freq),
            by = c("site","dbh_cat","setting")) %>%
  mutate(
    freq = replace_na(freq, list(integer(0))),
    level_size = m_target
  )

# 3D. iNEXT estimate (size-based)
spp_est_df <- inc_keys %>%
  rowwise() %>%
  mutate(
    tmp = list({
      f <- freq
      if (length(f) == 0) f <- integer(0)
      v <- as.numeric(c(T, f))
      out <- tryCatch(
        iNEXT::estimateD(list(v),
                         q = 0, datatype = "incidence_freq",
                         base = "size", level = level_size),
        error = function(e) NULL
      )
      if (is.null(out)) tibble(qD = 0) else tibble(qD = as.numeric(out$qD)[1])
    })
  ) %>%
  unnest(tmp) %>%
  ungroup() %>%
  transmute(site, dbh_cat, spp_exp_m = qD)

# 3E. Add back ALL strata (including true zeros)
df_dbh <- df_dbh_raw %>%
  left_join(spp_est_df, by = c("site","dbh_cat")) %>%
  mutate(
    spp_exp_m = replace_na(spp_exp_m, 0),
    obs_per10 = if_else(total_trees > 0,
                        m_target * obs_count_dbh / total_trees, 0),
    spp_per10 = if_else(total_trees > 0,
                        m_target * spp_count_dbh / total_trees, 0)
  )


# ---- 4) Bayesian models  -----------------------------------------

## 4.1) Abundance----------
df_counts <- df_dbh %>% filter(total_trees > 0)

fit_obs_zinb <- brm(
  formula = obs_count_dbh ~ setting * dbh_cat + (1 | site) + offset(log(total_trees)),
  data    = df_counts,
  family  = zero_inflated_negbinomial(),
  prior   = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd")
  ),
  chains  = 4, cores = 4, iter = 4000, warmup = 1000, seed = 123
)

## 4.2) Richness----------
df_rich <- df_dbh %>%
  transmute(site, setting, dbh_cat, spp_exp_m)

fit_spp_hg <- brm(
  formula = spp_exp_m ~ setting * dbh_cat + (1 | site),
  data    = df_rich,
  family  = hurdle_gamma(link = "log"),
  prior   = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd"),
    set_prior("exponential(1)", class = "shape")
  ),
  chains  = 4, cores = 4, iter = 4000, warmup = 1000, seed = 124
)

## 4.3) Host prob -----------------------------

# Build tree-level data: explode each site × DBH stratum into trees
fread("data/processed/inat.merged.csv") %>%  as_tibble()->inat_merged

df_host<-
inat_merged %>% 
  distinct(site,TreeID,dbh_cat=TreeCat) %>% 
  group_by(site,dbh_cat) %>% 
  summarise(host_count=n()) %>% 
  right_join(df_dbh_raw %>% select(site,dbh_cat,total_trees,setting)) %>% 
  mutate(host_count=if_else(is.na(host_count),0,host_count))

# Fit a binomial GLMM: probability that a tree is a host
fit_host_binom <- brm(
  host_count | trials(total_trees) ~ setting * dbh_cat + (1 | site),
  data = df_host,
  family = binomial(link = "logit"),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 1000, seed = 125
)


# ---- 6) Marginal means -------------------------------------------------------
## 6.2) Abundance------
emm_abund <- emmeans(fit_obs_zinb, ~ setting * dbh_cat,
                     type = "response", at = list(total_trees = 10))

abund_df <- emm_abund %>%
  as_tibble() %>%
  mutate(mean = prob,
         lwr = lower.HPD,
         upr = upper.HPD) %>%
  transmute(setting, dbh_cat, mean, lwr, upr) %>%
  mutate(setting = factor(setting, levels = c("forest","non-forest")),
         dbh_cat = factor(dbh_cat, levels = c("10 - 30","30 - 60","60 - 90","90+")))

## 6.2) Richness ------------
emm_rich <- emmeans(fit_spp_hg, ~ setting * dbh_cat, type = "response")

rich_df <- emm_rich %>%
  as_tibble() %>%
  mutate(mean = response,
         lwr = lower.HPD,
         upr = upper.HPD) %>%
  transmute(setting, dbh_cat, mean, lwr, upr) %>%
  mutate(setting = factor(setting, levels = c("forest","non-forest")),
         dbh_cat = factor(dbh_cat, levels = c("10 - 30","30 - 60","60 - 90","90+")))

## 6.3) Host prob--------
# Marginal predicted probabilities
emm_host <- emmeans(fit_host_binom, ~ setting * dbh_cat, type = "response")

host_df <- emm_host %>%
  as_tibble() %>%
  transmute(
    setting, dbh_cat,
    mean_prob = prob*100,
    lwr = lower.HPD*100,
    upr = upper.HPD*100
  ) %>%
  mutate(
    setting = factor(setting, levels = c("forest","non-forest")),
    dbh_cat = factor(dbh_cat, levels = c("10 - 30","30 - 60","60 - 90","90+"))
  )



# ---- 7) Figure ---------------------------------------------------
pos_dodge <- position_dodge2(width = 0.5)

## 7.1) Abundance-------
p_abund_means <- ggplot(abund_df, aes(x = setting, y = mean, colour = dbh_cat)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = pos_dodge, size = 0.9, alpha = 0.8) +
  scale_color_viridis_d("DBH category (cm)", option = "D", end = 0.9) +
  labs(x = NULL, y = paste0("Expected abundance per ", m_target, " trees")) +
  scale_y_break(breaks = c(7,20), space = 0, scales = c(1,5)) +
  my_theme14


## 7.2) Richness-------
p_rich_means <- ggplot(rich_df, aes(x = setting, y = mean, colour = dbh_cat)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = pos_dodge, size = 0.9, alpha = 0.8) +
  scale_color_viridis_d("DBH category (cm)", option = "D", end = 0.9, guide= "none") +
  labs(x = NULL, y = paste0("Expected richness per ", m_target, " trees")) +
  scale_y_break(breaks = c(7,10), space = 0, scales = c(1,5)) +
  my_theme14

## 7.3) host probability-------
p_host_prob <- ggplot(host_df, aes(x = setting, y = mean_prob, colour = dbh_cat)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = pos_dodge, size = 0.9, alpha = 0.8) +
  scale_color_viridis_d("DBH category (cm)", option = "D", end = 0.9, guide = "none") +
  labs(x = "Setting", y = "Host probability (%)") +
  scale_y_break(breaks = c(7,10), space = 0, scales = c(1,5)) +
  my_theme14

p_host_prob+plot_annotation(tag_levels = "A")
ggsave(file.path(OUT_DIR_DBH, "host_model_means_DBH.svg"),
       width = 4, height = 3, dpi = 500)
p_abund_means+plot_annotation(tag_levels = "A")
ggsave(file.path(OUT_DIR_DBH, "abund_model_means_DBH.svg"),
       width = 4, height =3.5, dpi = 500)
p_rich_means+plot_annotation(tag_levels = "A")
ggsave(file.path(OUT_DIR_DBH, "rich_model_means_DBH.svg"),
       width = 4, height = 3, dpi = 500)


p_abund_means+p_rich_means+p_host_prob+
  plot_annotation(tag_levels = "A")

# ---- 8) Export summaries -----------------------------------------------------
fx_obs  <- as_draws_df(fit_obs_zinb) %>% select(starts_with("b_"))
fx_rich <- as_draws_df(fit_spp_hg)   %>% select(starts_with("b_"))
fx_host <- as_draws_df(fit_host_binom)   %>% select(starts_with("b_"))

fx_obs_sum <- fx_obs %>%
  pivot_longer(everything(), names_to = "bterm", values_to = "draw") %>%
  mutate(term = str_remove(bterm, "^b_")) %>%
  group_by(term) %>%
  summarise(
    mean   = mean(draw),
    ci_low = quantile(draw, 0.025),
    ci_high= quantile(draw, 0.975),
    rr_mean= exp(mean),
    rr_ci_l = exp(ci_low),
    rr_ci_h = exp(ci_high),
    .groups = "drop"
  )

fx_rich_sum <- fx_rich %>%
  pivot_longer(everything(), names_to = "bterm", values_to = "draw") %>%
  mutate(term = str_remove(bterm, "^b_")) %>%
  group_by(term) %>%
  summarise(
    mean   = mean(draw),
    ci_low = quantile(draw, 0.025),
    ci_high= quantile(draw, 0.975),
    rr_mean= exp(mean),
    rr_ci_l = exp(ci_low),
    rr_ci_h = exp(ci_high),
    .groups = "drop"
  )


fx_host_sum <- fx_host %>%
  pivot_longer(everything(), names_to = "bterm", values_to = "draw") %>%
  mutate(term = str_remove(bterm, "^b_")) %>%
  group_by(term) %>%
  summarise(
    mean   = mean(draw),
    ci_low = quantile(draw, 0.025),
    ci_high= quantile(draw, 0.975),
    rr_mean= exp(mean),
    rr_ci_l = exp(ci_low),
    rr_ci_h = exp(ci_high),
    .groups = "drop"
  )

fwrite(fx_obs_sum,  file.path(OUT_DIR_DBH, "abundance_ZINB_fixed_effects_exp.csv"))
fwrite(fx_rich_sum, file.path(OUT_DIR_DBH, "richness_HG_fixed_effects_exp.csv"))
fwrite(fx_host_sum, file.path(OUT_DIR_DBH, "host_binomial.csv"))







# 9)Overall probability of being a host----------
# 9.1) Build prediction grid + weights
# weights = observed number of trees in each setting × dbh_cat
newdata_setting <- df_host %>%
  filter(total_trees > 0) %>%
  group_by(setting, dbh_cat) %>%
  summarise(total_trees = sum(total_trees), .groups = "drop") %>%
  group_by(setting) %>%
  mutate(weight = total_trees / sum(total_trees)) %>%
  ungroup()


# 9.2) Posterior predicted probabilities
#    re_formula = NA  -> population-level mean
#    (i.e. averaging over site random effects)
post_p <- posterior_linpred(
  fit_host_binom,
  newdata = newdata_setting,
  re_formula = NA,
  transform = TRUE
)

# post_p is: draws × rows(newdata_setting)
dim(post_p)


# 9.3) Weighted average probability by setting
#    within each posterior draw


idx_by_setting <- split(seq_len(nrow(newdata_setting)), newdata_setting$setting)

post_setting_mean <- lapply(names(idx_by_setting), function(s) {
  idx <- idx_by_setting[[s]]
  w   <- newdata_setting$weight[idx]
  
  tibble(
    setting = s,
    draw    = seq_len(nrow(post_p)),
    prob    = as.vector(post_p[, idx, drop = FALSE] %*% w)
  )
}) %>%
  bind_rows()

library(tidybayes)
setting_summary <- post_setting_mean %>%
  group_by(setting) %>%
  median_hdi(prob, .width = 0.95)

setting_summary
