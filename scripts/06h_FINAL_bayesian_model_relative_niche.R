# 30f_bayesian_model_use_availability.R
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(brms)
  library(posterior)
  library(ggplot2)
  library(viridis)
  library(data.table)
})

IN_DAT <- "outputs/06_bayesian_model_relative_niche/model_data/use_availability_rel01.csv"
OUT_DIR <- "outputs/06_bayesian_model_relative_niche/final_bayes"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- fread(IN_DAT) %>%
  filter(
    is.finite(moist_rel01),
    is.finite(maxT_rel01),
    is.finite(minT_rel01),
    is.finite(wind_rel01)
  ) %>%
  mutate(species_label = factor(species_label))


set.seed(123)

#0) prep----------
# Make sure species_label is a plain character vector
dat2 <- dat %>%
  mutate(species_label = as.character(species_label),
         species_label = str_replace(species_label,"_"," "))

# Count presences per species
pres_counts <- dat2 %>%
  filter(y == 1) %>%
  count(species_label, name = "n_pres")

# Count available background points per species
bg_counts <- dat2 %>%
  filter(y == 0) %>%
  count(species_label, name = "n_bg")

# Define target background size = 10 x presences (minimum 10)
bg_targets <- full_join(pres_counts, bg_counts, by = "species_label") %>%
  mutate(
    n_pres    = coalesce(n_pres, 0L),
    n_bg      = coalesce(n_bg, 0L),
    target_bg = pmax(10L, 10L * n_pres)
  ) %>%
  select(species_label, target_bg)

# Sample background points species-by-species using base R sampling
bg_small <- dat2 %>%
  filter(y == 0) %>%
  inner_join(bg_targets, by = "species_label") %>%
  group_split(species_label) %>%
  purrr::map_dfr(function(df_sp) {
    n_target <- df_sp$target_bg[1]
    n_avail  <- nrow(df_sp)
    
    idx <- sample.int(
      n = n_avail,
      size = n_target,
      replace = n_avail < n_target
    )
    
    df_sp[idx, , drop = FALSE]
  }) %>%
  select(-target_bg)

# Combine presences + sampled background 
dat_small <- bind_rows(
  dat2 %>% filter(y == 1),
  bg_small)


# Optional sanity check
dat_small %>%
  count(species_label, y) %>%
  arrange(species_label, y)



#plot
# compute extent from the data
x_rng <- range(dat_small$lon, na.rm = TRUE)
y_rng <- range(dat_small$lat, na.rm = TRUE)

ggplot(dat_small, aes(x = lon, y = lat, colour = species_label,shape=factor(y))) +
  geom_point(alpha = 0.7, size = 1.8) +
  coord_cartesian(
    xlim = x_rng,
    ylim = y_rng,
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  viridis::scale_colour_viridis(discrete=T)

#1) Model-------------
##priors

priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(lkj(2), class = "cor")
)

#Bayesian hierarchical logistic regression
fit <- brm(
  y ~ 0 + species_label +
    species_label:moist_rel01 +
    species_label:maxT_rel01 +
    species_label:minT_rel01 +
    species_label:wind_rel01,
  data   = dat_small,
  family = bernoulli(link = "logit"),
  prior  = c(
    prior(normal(0, 2), class = "b")
  ),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
summary(fit)

saveRDS(fit, file.path(OUT_DIR, "fit_use_availability.rds"))

# quick summaries
sink(file.path(OUT_DIR, "summary.txt"))
print(summary(fit))
sink()

message("DONE: model saved to ", file.path(OUT_DIR, "fit_use_availability_4chelsa.rds"))


#Posterior summaries-------------
fit<-readRDS(file.path(OUT_DIR, "fit_use_availability.rds"))
summary(fit)

dispersal <-
  fread("outputs/06_bayesian_model_relative_niche/environmental/dispersal_agents.csv") %>% 
  distinct(species_label=taxon_name,dispersal=dispersal_agent_summary)



# extract fixed effects
coefs <- as.data.frame(fixef(fit)) %>%
  rownames_to_column("term")

# make lookup table so species names stay nice
species_lookup <- tibble(
  species_label = unique(as.character(dat_small$species_label))
) %>%
  mutate(species_key = gsub("\\s+", "", species_label)) %>% 
  left_join(dispersal)

# keep only slope terms (exclude intercept-only species terms)
slopes_df <- coefs %>%
  filter(str_detect(term, ":")) %>%
  mutate(
    species_key = str_extract(term, "(?<=species_label)[^:]+"),
    predictor = str_extract(term, "(?<=:).+$")
  ) %>%
  left_join(species_lookup, by = "species_key") %>%
  mutate(
    predictor = recode(
      predictor,
      moist_rel01 = "Aridity index",
      maxT_rel01  = "Mean temperature of<br>warmest quarter (°C)",
      minT_rel01  = "Minimum temperature<br>of coldest month (°C)",
      wind_rel01  = "Mean near-surface<br>wind speed (ms<sup>-1</sup>)"),
    species_label = factor(species_label)
  ) %>% 
  mutate(
    predictor = factor(
      predictor,
      levels = c(
        "Mean temperature of<br>warmest quarter (°C)",
        "Minimum temperature<br>of coldest month (°C)",
        "Aridity index",
        "Mean near-surface<br>wind speed (ms<sup>-1</sup>)"
      )
    )
  )

#Plot--------------

#order species within each predictor by estimate
slopes_df <- slopes_df %>%
  group_by(predictor) %>%
  arrange(Estimate, .by_group = TRUE) %>%
  mutate(species_label = factor(species_label, levels = unique(species_label))) %>%
  ungroup()

# Standard ggplot theme 
my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = ggtext::element_markdown(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top",
    axis.text = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown()
  )



# make a colour-group variable that only keeps dispersal in the wind panel
slopes_df2 <- slopes_df %>%
  mutate(
    predictor_chr = trimws(as.character(predictor)),
    colour_group = if_else(
      str_detect(predictor_chr,"wind"),
      as.character(dispersal),
      "nonwind"
    )
  )

# colours for the real dispersal classes
disp_levels <- unique(as.character(slopes_df2$dispersal))
disp_cols <- setNames(viridis(length(disp_levels)), disp_levels)

# add a neutral colour for all non-wind panels
all_cols <- c(nonwind = "grey40", disp_cols)

ggplot(slopes_df2, aes(x = Estimate, y = species_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(
    aes(xmin = Q2.5, xmax = Q97.5, colour = colour_group),
    size = 1.5,alpha=0.7
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('", x, "')"))
  ) +
  facet_wrap(~ predictor, ncol = 4) +
  labs(
    x = "Effect on accidental epiphyte occurrence (log-odds)",
    y = NULL,
    colour = NULL
  ) +
  scale_color_manual("Dispersal mode",
    values = all_cols,
    breaks = names(disp_cols)   # only show dispersal classes in legend
  ) +
  my_theme14

ggsave(file.path(OUT_DIR,"model_results.jpeg"),width=11,height=5,dpi=500)


