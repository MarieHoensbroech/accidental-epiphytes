# 10b_dispersal_analyses.R
# Assess what dispersal syndromes accidental epiphytes use
# 1) using traditional stats
# 2) using Bayesian
# 3) Plotting


suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(viridis)
  library(forcats)
  library(scales)
  library(patchwork)
  library(stringr)
  library(data.table)
  library(cowplot)
  library(grid)
  library(ggimage)
  library(tidybayes)
  library(brms)
  library(rstatix)
  library(janitor)
  library(purrr)
  library(data.table)
  library(emmeans)
  library(broom.mixed)
})

#  Parameters -------------------------------
IN_SPECIES_MODES   <- "outputs/10_dispersal/species_dispersal_modes.csv"
IN_INAT_OBS_DISP   <- "outputs/10_dispersal/inat_obs_with_dispersal.csv"
IN_INAT_FIELD_DISP <- "outputs/10_dispersal/inat_field_with_dispersal.csv"

OUT_DIR_DISP  <- "outputs/10_dispersal"
dir.create(OUT_DIR_DISP, showWarnings = FALSE, recursive = TRUE)
update_geom_defaults("point", list(alpha = 0.8))
set.seed(2025)

# Theme -------------------------------
my_theme14 <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top"
  )

# Helpers -------------------------------
Mode_ <- function(x) {
  x <- x[is.finite(match(x, x)) & !is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}
collapse_values <- function(x, k = 3L) {
  ux <- unique(na.omit(x))
  if (length(ux) <= k) paste(ux, collapse = " | ") else Mode_(x)
}

# Load data -------------------------------
species_modes_final <- fread(IN_SPECIES_MODES) %>% clean_names() %>% as_tibble() %>% distinct()
inat_obs_with_disp  <- fread(IN_INAT_OBS_DISP) %>% clean_names() %>% as_tibble()%>% distinct()
inat_field_with_disp<- fread(IN_INAT_FIELD_DISP) %>% clean_names() %>% as_tibble()%>% distinct()

# Harmonise labels -------------------------------
# keep "mixed" as a class in descriptive plots, 
# but possibly EXCLUDE from the GLMM (or model with its own coefficient).
valid_agents <- c("animal","wind","water","ballistic", "unassisted","human","mixed","unknown","none")

# Combine datasets -------------------------------
df_dispersal <- bind_rows(
  inat_field_with_disp %>% 
    mutate(dataset = "field") %>% 
    select(-dbh_cat,-obs_count_dbh,-spp_count_dbh) %>% 
    unique(),
  inat_obs_with_disp   %>% mutate(dataset = "inat")
) %>%
  # Ensure consistent columns exist
  mutate(
    dispersal_agent = coalesce(dispersal_agent_summary, "unknown"),
    dispersal_agent = ifelse(is.na(dispersal_agent), "unknown", dispersal_agent),
    dispersal_agent = factor(dispersal_agent, levels = valid_agents),
    genus_name      = word(genus_name, 1),
    family          = coalesce(family, NA_character_)
  ) 


inat_field <- inat_field_with_disp %>% 
  mutate(
    dispersal_agent = coalesce(dispersal_agent_summary, "unknown"),
    dispersal_agent = ifelse(is.na(dispersal_agent), "unknown", dispersal_agent),
    dispersal_agent = factor(dispersal_agent, levels = valid_agents),
    genus_name      = word(genus_name, 1)
  ) 
# Summaries -------------------------------
sum_by_mode <- df_dispersal %>%
  mutate(dispersal_agent = fct_na_value_to_level(dispersal_agent, level = "unknown")) %>%
  group_by(dataset, dispersal_agent) %>%
  summarise(
    n_species = n_distinct(taxon_name),
    n_obs     = sum(suppressWarnings(as.numeric(epi_count)), na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(dataset, desc(n_species), desc(n_obs))

fwrite(sum_by_mode, file.path(OUT_DIR_DISP, "10b_summary_by_agent_dataset.csv"))
print(sum_by_mode, n = 100)

# __________________________----------------
# Plots

## bar proportion by souce-------------------------------
p_comp_species <-
  ggplot(sum_by_mode, aes(x = fct_reorder(dispersal_agent, n_species), y = n_species, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75)) +
  coord_flip() +
  scale_fill_viridis_d(begin=0.4,end = 0.8, option = "D") +
  labs(x = NULL, y = "Number of species", fill = "Dataset") +
  my_theme14

p_comp_obs <-
  ggplot(sum_by_mode, aes(x = fct_reorder(dispersal_agent, n_obs), y = n_obs, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75)) +
  coord_flip() +
  scale_fill_viridis_d(begin=0.4,end = 0.8, option = "D") +
  labs(x = NULL, y = "Number of epiphyte observations", fill = "Dataset") +
  my_theme14

g_comp <- p_comp_obs + p_comp_species + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR_DISP, "10b_dispersal~source_bars.png"), g_comp, width = 12, height = 6, dpi = 300)






## bar proportion by setting -------------------------------
sum_by_setting <- df_dispersal %>%
  drop_na(setting) %>% filter(!setting%in%c("Optional","")) %>% 
  mutate(dispersal_agent = fct_na_value_to_level(dispersal_agent, level = "unknown")) %>%
  group_by(dataset, dispersal_agent,setting) %>%
  summarise(
    n_species = n_distinct(taxon_name),
    n_obs     = sum(suppressWarnings(as.numeric(epi_count)), na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(dataset, desc(n_species), desc(n_obs))

fwrite(sum_by_setting, file.path(OUT_DIR_DISP, "10b_summary_by_agent+setting_dataset.csv"))

print(sum_by_mode, n = 100)

p_disp_species <-
  ggplot(sum_by_setting, aes(x = dispersal_agent, y = n_species, fill = setting)) +
  geom_col(position=position_dodge()) +
  coord_flip() +
  facet_wrap(~dataset,scales="free_x")+
  scale_fill_viridis_d(begin=0.4,end = 0.8, option = "D") +
  labs(x = "", y = "Number of species", fill = "Setting") +
  my_theme14

p_disp_obs <-
  ggplot(sum_by_setting, aes(x = dispersal_agent, y = n_obs, fill = setting)) +
  geom_col(position=position_dodge()) +
  coord_flip() +
  facet_wrap(~dataset,scales="free_x")+
  scale_fill_viridis_d(begin=0.4,end = 0.8, option = "D") +
  labs(x = "Dataset", y = "Number of epiphyte observations", fill = "Setting") +
  my_theme14

g_disp <- p_disp_obs + p_disp_species + plot_layout(ncol = 2)
g_disp
ggsave(file.path(OUT_DIR_DISP, "10b_dispersal~setting_bars.png"), g_disp, width = 12, height = 6, dpi = 300)



## box proportion field data -----------
inat_field <- inat_field %>%
  arrange(site) %>%
  filter(!is.na(setting), !setting %in% c("Optional", "")) %>%
  mutate(dispersal_agent = na_if(str_trim(dispersal_agent), "")) %>%
  mutate(
    dispersal_agent2 = case_when(
      is.na(taxon_name) ~ "none",
      TRUE                     ~ dispersal_agent
    ),
    
    dispersal_agent2 = fct_relevel(factor(dispersal_agent2), "none", "unknown")
  ) 

df_box <- inat_field %>%
  reframe(site, taxon_name,genus_name,dispersal_agent2, setting, dbh_cat,
          spp_count_dbh,spp_count,obs_count_dbh,obs_count) %>% 
  group_by(site,dbh_cat) 

levels_order <- c("none", "unknown", "ballistic","unassisted", "water", "wind", "animal")

df_box <- inat_field %>%
  reframe(
    site, taxon_name, dispersal_agent2, setting, 
    spp_count, obs_count
  ) %>%
  mutate(
    dispersal_agent2 = fct_relevel(factor(dispersal_agent2), levels_order),
    setting          = as.factor(setting)
  ) %>%
  drop_na(taxon_name) %>% 
  group_by(site, #dbh_cat, 
           setting, dispersal_agent2) %>%
  mutate(obs_disp = n(),
         sp_disp= n_distinct(taxon_name)) %>% 
  mutate(
    prop_obs    = ifelse(obs_disp > 0, obs_disp / obs_count, NA_real_),
    prop_spp    = ifelse(sp_disp > 0, sp_disp / spp_count, NA_real_)
  ) %>%
  ungroup()


p_prop_obs<-
ggplot(df_box, aes(x = dispersal_agent2, y = prop_obs, 
                   colour = setting)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth=1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Dispersal syndrome",
    y = "Proportion of observations per site"
  ) +
  my_theme14+
  scale_colour_viridis(discrete=T,begin = 0.4,end=0.8)+
  theme(legend.position = "top")


p_prop_sp<-
ggplot(df_box, aes(x = fct_reorder(dispersal_agent2, prop_spp), y = prop_spp, 
                   colour = setting)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth=1, 
               position=position_dodge2(preserve="single", padding=0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Dispersal agent",
    y = "Proportion of species per site"
  ) +
  my_theme14+
  scale_colour_viridis("Setting",discrete=T,begin = 0.4,end=0.8)+
  theme(legend.position = "top")

p_prop_sp

ggsave(file.path(OUT_DIR_DISP, "10b_prop~setting_box.png"), 
       p_prop_sp, width = 8, height = 6, dpi = 300)



# Traditional stats -------------------------------
# 1) Chi-square test: are agent distributions different between datasets?
## YES.......
tbl_agents <- df_dispersal %>%
  filter(!is.na(dispersal_agent), dispersal_agent != "unknown") %>%
  count(dataset, dispersal_agent, name = "n") %>%
  pivot_wider(names_from = dispersal_agent, values_from = n, values_fill = 0) %>%
  column_to_rownames("dataset") %>%
  as.matrix()

{
  chi_agents <- chisq.test(tbl_agents)
  cramer_v   <- rstatix::cramer_v(as.table(tbl_agents))
  capture.output({
    print(chi_agents)
    print(cramer_v)
  }, file = file.path(OUT_DIR_DISP, "10b_chisq_dataset_vs_agent.txt"))
}

# 2) Species-level proportions by agent (non-parametric)
sp_counts <- df_dispersal %>%
  group_by(dataset, dispersal_agent, taxon_name, genus_name) %>%
  summarise(epi_n = sum(suppressWarnings(as.numeric(epi_count)), na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(dispersal_agent), dispersal_agent != "unknown")

# Kruskal–Wallis within each dataset
kw_results <- df_box %>%
  group_by(dispersal_agent2) %>%
  filter(n_distinct(setting)>1) %>%
  kruskal_test(prop_spp ~ setting) %>% 
  
  mutate(
    p_value_label = case_when(
      is.na(p)             ~ NA_character_,
      p < 1e-4             ~ "< 0.0001",
      TRUE                       ~ sprintf("%.4f", p)
    )
  )


# Dunn post-hoc with BH correction
dunn_results <- df_box %>%
  group_by(dispersal_agent2) %>%
  filter(n_distinct(setting)>1) %>%
  dunn_test(prop_spp ~ setting, p.adjust.method = "bonferroni")

fwrite(kw_results,   file.path(OUT_DIR_DISP, "10b_kw_species_prop_by_setting.csv"))
fwrite(dunn_results, file.path(OUT_DIR_DISP, "10b_dunn_species_prop_by_setting.csv"))

#_____________________----------------
# Bayesian model proportion-binomial -------------------------------
# Negative Binomial GLMM: species-level epiphyte counts ~ agent (+ dataset) + (1|genus)
# Notes:
# - Overdispersion likely -> use negbinomial
# - We exclude "unknown"; keep or drop "mixed" (here we keep it as a class)
binom_obs <-df_box%>%
  filter(dispersal_agent2 != "unknown") %>%
  mutate(
    dispersal_agent = factor(dispersal_agent2,
                             levels=c("animal","wind","mixed","water","ballistic","unassisted")),
    setting = factor(setting)
  ) %>% 
  group_by(site, setting, dispersal_agent) %>%
  mutate(
    success = n_distinct(taxon_name,na.rm = T)) %>% # obs in this syndrome
  group_by(site,setting) %>% 
  reframe(success,dispersal_agent,                                              
    trials  = n_distinct(taxon_name, na.rm = TRUE)  # total obs at site
  ) %>%
  distinct() %>% 
  # guard against any accidental trials < success
  mutate(trials = pmax(trials, success))

# One joint model across syndromes with interaction:
m_obs <- brm(
  bf(success | trials(trials) ~ 0 + dispersal_agent * setting + (1 | site)),
  data    = binom_obs,
  family  = binomial(),
  prior   = c(
    set_prior("normal(0, 1.5)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd")  # random effect sd
    # no Intercept prior, because there is no Intercept
  ),
  chains  = 4, iter = 4000, seed = 2025,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

saveRDS(m_obs, file.path(OUT_DIR_DISP, "10b_brms_nbglmm_m1_sp_prop_by_setting.rds"))

print(summary(m_obs))



coef_tidy <- tidy(m_obs, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
fwrite(coef_tidy, file.path(OUT_DIR_DISP, "10b_brms_tidy_fixed_effects.csv"))


## Hypothesis tetsing-------
# Test Rural vs Urban difference per syndrome (example for 'wind')
hypothesis(m_obs, "dispersal_agentwind + dispersal_agentwind:settingUrban = 0")
#Urban vs Rural for animal
hypothesis(m_obs, "settingUrban + dispersal_agentanimal = 0")
hypothesis(m_obs, "dispersal_agentanimal = 0")
# For ANIMAL
hypothesis(m_obs, "settingUrban + dispersal_agentanimal = 0")

# TidyBayes: construct expected counts (on response scale) assuming average random effect ~ 0
emm <- emmeans(m_obs, ~ setting | dispersal_agent, type = "response") 
# 'type = "response"' returns posterior summaries on the probability scale

# Pairwise Urban - Rural contrasts within each syndrome (on the response scale):
contr <- contrast(emm, method = "revpairwise", by = "dispersal_agent", adjust = "none")
summary(contr, infer = c(TRUE, TRUE))  # add intervals
fwrite(as.data.frame(emm), file.path(OUT_DIR_DISP, "10b_brms_emmeans.csv"))



## Probability ---------
newdat <- expand.grid(
  dispersal_agent = levels(binom_obs$dispersal_agent),
  setting          = levels(binom_obs$setting),
  site             = NA  # population-level (marginal over site RE)
)


# trials is required for binomial link to compute mu, but on response scale it’s fine to set any positive number
newdat$trials <- 1
newdat$success <- 0

fits <- fitted(m_obs, newdata = newdat, re_formula = NA, summary = TRUE, scale = "response")
out  <- cbind(newdat[, c("dispersal_agent","setting")],
              as.data.frame(fits))  # columns: Estimate, Q2.5, Q97.5, etc.

### Plot------
p_prob<-
ggplot(out, aes(dispersal_agent, Estimate, color = setting)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.4)) +
  #facet_wrap(~ dispersal_agent, scales = "free_y") +
  labs(y = "Estimated proportion", x = NULL, 
       title = "Posterior fitted proportions by setting and syndrome") +
  scale_color_viridis(discrete=T, begin=0.4,end=0.8)+
  my_theme14

p_prob
ggsave(file.path(OUT_DIR_DISP, "10b_dispersal~source_bars.png"), 
       p_prob,width = 12, height = 6, dpi = 300)



## SUMMARY OF BAYESIAN BINOMIAL GLMM RESULTS---------------

# • Model: success | trials(trials) ~ dispersal_agent * setting + (1 | site)
# • Animal and wind dispersal remain dominant overall but effects are modest:
#     - Animal: Estimate ≈ 0.60 (95% CI: –0.09 to 1.47)
#     - Wind:   Estimate ≈ 0.41 (95% CI: –0.19 to 1.15)
# • Urban main effect: –0.64 (95% CI: –1.73 to 0.21) → no clear global increase.
# • Wind × Urban interaction: –1.03 (95% CI: –2.01 to –0.07) → strong decline in urban sites.
# • Odds ratios Urban vs Rural:
#     - Wind: OR ≈ 0.19 (95% HPD: 0.04–0.42)
#     - Animal: OR ≈ 0.55 (95% HPD: 0.11–1.11)
# • Rare syndromes (water, ballistic, unassisted) trend lower in Urban but with wide uncertainty.
# • Interpretation: Urban environments do NOT favor all syndromes; wind dispersal declines,
#   animal dispersal remains stable, others inconclusive.

#___________________________-------
# Bayesian: Composition model (multinomial) --------------------------------------------
# Goal: model the *composition* of dispersal syndromes jointly so class probs sum to 1.
# Data: species counts per site × setting × dispersal_agent (excluding 'unknown' & 'human').

## Data for multinomial (species-level counts) ----
classes <- c("animal","wind","mixed","water","ballistic","unassisted")

# Site-level data (richness per site with tree-based effort) ------------------
tree_effort <- inat_field %>%
  filter(!is.na(setting), !setting %in% c("Optional","")) %>%
  distinct(site,setting,total_trees) %>% 
  group_by(site, setting) %>%
  summarise(
    n_trees = sum(total_trees),
    .groups = "drop"
  )


rich_by_agent <- df_box %>%
  filter(!is.na(taxon_name), !dispersal_agent2 %in% c("unknown")) %>%
  mutate(
    dispersal_agent = fct_relevel(
      factor(dispersal_agent2),
      c("animal","wind","mixed","water","ballistic","unassisted")
    ),
    setting = droplevels(factor(setting, levels = c("Rural","Urban")))
  ) %>%
  group_by(site, setting, dispersal_agent) %>%
  summarise(sp_disp = n_distinct(taxon_name), .groups = "drop") %>%
  left_join(tree_effort, by = c("site","setting")) %>%
  mutate(
    n_trees    = coalesce(n_trees, 0L),
    n_trees    = ifelse(n_trees < 1, 1L, n_trees),  # guardrail
    log_ntrees = log(n_trees)
  )

fwrite(rich_by_agent, file.path(OUT_DIR_DISP, "10b_richness_by_syndrome_per_site.csv"))


# 1) SETUP & DATA PREP -------------------------------------------------
# Goal: Provide ONE robust solution that fixes the rank deficiency (ncol=12 vs rank=10)
# by using a cell-means parameterization, fits an NB-GLMM with a tree offset,
# tests ONE planned hypothesis, and then runs Bayesian post-hoc comparisons.


set.seed(2025)



rich_by_agent <- rich_by_agent %>%
  mutate(
    dispersal_agent = factor(dispersal_agent),
    setting         = factor(setting),
    site            = factor(site)
  )


# 2) MODEL: NB-GLMM WITH CELL-MEANS (IDENTIFIABLE) ---------------------
# Use cell-means only: 0 + dispersal_agent:setting
# This eliminates the rank deficiency caused by "0 + A * B".
# Priors are set tighter for stability (especially for shape), and we use a high adapt_delta.

m_nb <- brm(
  formula = sp_disp ~ 0 + dispersal_agent:setting + offset(log_ntrees) + (1 | site),
  data    = rich_by_agent,
  family  = negbinomial(),  # log link by default
  prior   = c(
    set_prior("normal(0, 1)",       class = "b"),
    set_prior("student_t(3, 0, 1)", class = "sd"),      # slightly tighter than 2.5
    set_prior("normal(0, 0.5)",     class = "shape")    # on log(shape): stabilizes overdispersion
  ),
  chains  = 4, iter = 3000, seed = 2025,
  control = list(adapt_delta = 0.99, max_treedepth = 14)
)

{
print(m_nb)
saveRDS(m_nb, file.path(OUT_DIR_DISP, "10b_brms_nbglmm_m1_sp_prop_by_setting.rds"))


# 3) ONE PLANNED HYPOTHESIS (AVERAGE SETTING EFFECT ON PER-TREE RATE) --
# Hypothesis: "Average per-tree richness differs between settings (marginalized over dispersal agents)."
# With the cell-means model, each coefficient is a log-rate-per-tree for an agent×setting cell.
# We compute the marginal (over agents) per setting on the link scale, then contrast settings.
# This is a direct posterior linear-contrast based test (fully Bayesian).

# Extract posterior draws of cell-means (on link scale) -----------------
draws <- as_draws_df(m_nb)  # includes .draw, .chain, .iteration and all 'b_' columns

# Identify the b_ columns for the interaction terms b_dispersal_agent{A}:setting{S}
coef_cols <- grep("^b_dispersal_agent[^:]+:setting[^:]+$", colnames(draws), value = TRUE)

# Parse agent and setting names out of term labels
coef_map <- tibble(term = coef_cols) %>%
  mutate(
    agent  = str_match(term, "^b_dispersal_agent([^:]+):setting([^:]+)$")[, 2],
    set    = str_match(term, "^b_dispersal_agent([^:]+):setting([^:]+)$")[, 3]
  )

# Long posterior table for the cell means (eta = log rate per tree)
draws_long <- draws %>%
  select(.draw, all_of(coef_cols)) %>%
  pivot_longer(-.draw, names_to = "term", values_to = "eta") %>%
  left_join(coef_map, by = "term")

# Compute marginal (over agents) per setting on the link scale
eta_by_setting <- draws_long %>%
  group_by(.draw, set) %>%
  summarise(eta = mean(eta), .groups = "drop")

# Build all pairwise setting contrasts on the link scale (differences in log-rate)
setting_levels <- sort(unique(eta_by_setting$set))
contrast_grid <- t(combn(setting_levels, 2)) %>%
  as.data.frame() %>%
  rename(set1 = V1, set2 = V2)

# For each pair (set1, set2), compute posterior diff = eta(set1) - eta(set2)
setting_contrasts <- contrast_grid %>%
  mutate(
    res = map2(set1, set2, ~{
      left  <- eta_by_setting %>% filter(set == .x) %>% select(.draw, eta) %>% rename(eta1 = eta)
      right <- eta_by_setting %>% filter(set == .y) %>% select(.draw, eta) %>% rename(eta2 = eta)
      out <- left %>% inner_join(right, by = ".draw") %>%
        mutate(diff_link = eta1 - eta2,
               ratio_resp = exp(diff_link))  # multiplicative effect on per-tree richness
      out
    })
  )

# Summarise ONE primary contrast: the first available pair of settings (edit if you have a specific pair)
primary_test <- setting_contrasts$res[[1]]
primary_pair <- setting_contrasts[1, c("set1", "set2")]

hyp_summary <- primary_test %>%
  summarise(
    mean_diff_link = mean(diff_link),
    l95_diff_link  = quantile(diff_link, 0.025),
    u95_diff_link  = quantile(diff_link, 0.975),
    P_greater_0    = mean(diff_link > 0),      # posterior probability set1 > set2 on per-tree rate
    mean_ratio     = mean(ratio_resp),
    l95_ratio      = quantile(ratio_resp, 0.025),
    u95_ratio      = quantile(ratio_resp, 0.975)
  ) %>%
  mutate(
    contrast = paste0("Setting ", primary_pair$set1, "  vs  ", primary_pair$set2,
                      " (marginal over dispersal_agent)")
  )

cat("\n--- PLANNED HYPOTHESIS: Average per-tree richness differs between settings ---\n")
print(hyp_summary)
cat("Interpretation: On the log scale, 'mean_diff_link' is the average difference.\n",
    "On the response scale, 'mean_ratio' is the multiplicative per-tree effect.\n",
    "P_greater_0 is the posterior probability that the first setting has higher per-tree richness.\n\n")


# 4) BAYESIAN POST-HOC TESTS VIA EMMEANS --------------------------------
# We perform pairwise contrasts (posterior) at offset = 0 (i.e., per-tree rate),
# so the marginal means are directly interpretable as per-tree expected richness.

# Build a reference grid with log_ntrees fixed to 0 (per-tree exposure)
rg <- ref_grid(m_nb, at = list(log_ntrees = 0))

# 4a) Pairwise comparisons of dispersal agents WITHIN each setting (per tree)
emm_agent_by_setting <- emmeans(rg, ~ dispersal_agent | setting, type = "response")
posthoc_agents_within_setting <- pairs(emm_agent_by_setting, adjust = "holm")

cat("\n--- POST-HOC: Pairwise differences among dispersal agents within each setting (per-tree) ---\n")
print(summary(posthoc_agents_within_setting))

# 4b) Pairwise comparisons of settings WITHIN each dispersal agent (per tree)
emm_setting_by_agent <- emmeans(rg, ~ setting | dispersal_agent, type = "response")
posthoc_settings_within_agent <- pairs(emm_setting_by_agent, adjust = "holm")

cat("\n--- POST-HOC: Pairwise differences among settings within each dispersal agent (per-tree) ---\n")
print(summary(posthoc_settings_within_agent))


# 5)--- Fitted per-tree richness (NB model) --------------------------------
# Adjusting your previous "probability" code to the new NB GLMM (m_nb) that models counts
# with an offset. We obtain posterior fitted *per-tree* expected richness by setting
# log_ntrees = 0 and marginalizing over site RE (re_formula = NA).

# Build prediction grid: all dispersal_agent × setting combinations, per tree
newdat <- expand.grid(
  dispersal_agent = levels(rich_by_agent$dispersal_agent),
  setting         = levels(rich_by_agent$setting)
)

# Per-tree exposure (offset = log(1) = 0) and population-level (no site RE)
newdat$log_ntrees <- 0
newdat$site       <- NA

# Get fitted values on the response scale (expected richness per tree)
fits <- fitted(
  m_nb,
  newdata     = newdat,
  re_formula  = NA,        # marginal over site random effects
  summary     = TRUE,
  scale       = "response" # expected count (per tree because log_ntrees=0)
)

out <- cbind(newdat[, c("dispersal_agent","setting")],
             as.data.frame(fits))  # columns: Estimate, Q2.5, Q97.5, Est.Error

# --- Plot ----------------------------------------------------------------
p_rate <-
  ggplot(out, aes(dispersal_agent, Estimate, color = setting)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                width = 0.2,
                position = position_dodge(width = 0.4)) +
  labs(
    y = "Estimated epiphyte richness per tree",
    x = NULL
    ,
    title = "Posterior fitted per-tree richness by setting and dispersal syndrome"
  ) +
  scale_color_viridis(discrete = TRUE, begin = 0.4, end = 0.8) +
  my_theme14

p_rate

# Save
ggsave(
  file.path(OUT_DIR_DISP, "10b_per_tree_richness_by_setting_and_syndrome.png"),
  p_rate, width = 8, height = 6, dpi = 300)
  



# 5) OPTIONAL: QUICK DIAGNOSTIC SUMMARIES (PRINT ONLY) ------------------
# (No figures to keep this a single continuous script)
summ <- summary(m_nb)
cat("\n--- MODEL DIAGNOSTICS (quick) ---\n")
cat("Max Rhat: ", max(summ$fixed$Rhat, na.rm = TRUE), "\n")
cat("Min Bulk-ESS (fixed): ", min(summ$fixed$Bulk_ESS, na.rm = TRUE), "\n")
cat("Min Tail-ESS (fixed): ", min(summ$fixed$Tail_ESS, na.rm = TRUE), "\n")

}
