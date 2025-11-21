# 09_hosttree_analyses.R
#
# Assess which trees accidental epiphytes occur in
# 1) using traditional stats
# 2) using Bayesian
# 3) Plotting


suppressPackageStartupMessages({
  library(tidyverse)     
  library(data.table)
  library(ggpubr)
  library(viridis)
  library(scales)
  library(maps)
  library(patchwork)
  library(cowplot)
  library(grid)          # base; optional to attach explicitly
  library(ggimage)
  library(rnaturalearth)
  library(brms)
  library(tidybayes)
  library(emmeans)
  library(multcompView)
  library(rstatix)
  library(janitor)
  library(tidytext)
  library(stringr)
  library(circlize)
})

# Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"

OUT_DIR_HOST <- "outputs/09_hosttree"
dir.create(OUT_DIR_HOST, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 12))
update_geom_defaults("point", list(alpha = 0.8))


clean_to_genus <- function(x) {
  # normalize first
  y <- x %>%
    str_to_lower() %>%
    str_replace_all("\\bcf\\b", "")  
  
  # map common names/typos to Latin genus (English + German hints)
  y <- case_when(
    
    is.na(y) ~ "unknown",
    str_detect(y, "\\bdead\\b|\\btot\\b")                                ~ "robinia",
    str_detect(y, "\\bpine\\b")                                          ~ "pinus",
    str_detect(y, "\\bplum\\b|\\bcherry\\b")                             ~ "prunus",
    str_detect(y, "\\bahorn\\b")                                         ~ "acer",
    str_detect(y, "\\bsvensk\\b")                                   ~ "sorbus",
    str_detect(y, "\\brubinia(s)?\\b|\\brobinie(n|s)?\\b|\\brubinie(n|s)?\\b") ~ "robinia",
    str_detect(y, "\\bvogelbeere\\b|\\beuropean\\s+ash\\b")              ~ "sorbus",
    str_detect(y, "\\bphoenix\\b|\\bcanariensis\\b")                     ~ "phoenix",
    str_detect(y, "\\bpalm(s)?\\b")                                      ~ "arecaceae",
    str_detect(y, "\\bpoplar(s)?\\b|\\bpappel(n|s)?\\b")                 ~ "populus",
    str_detect(y, "\\bbeech(es)?\\b|\\bbuche(n|s)?\\b")                  ~ "fagus",
    str_detect(y, "\\boak(s)?\\b|\\beiche(n|s)?\\b")                     ~ "quercus",
    str_detect(y, "\\bwillow(s)?\\b|\\bweide(n|s)?\\b")                  ~ "salix",
    str_detect(y, "\\bbirch(es)?\\b|\\bbirke(n)?\\b")                    ~ "betula",
    str_detect(y, "\\btila\\b|\\btil(l)?ia\\b|\\blinde(n)?\\b")          ~ "tilia",
    str_detect(y, "\\bsambuc(a|us)(es)?\\b|\\bholunder\\b")              ~ "sambucus",
    
    str_detect(y, "\\btree\\s*fern\\b") ~ "treefern",
    str_detect(y, "\\bseq[uio]{3}a\\b") ~ "sequoia",
    str_detect(
      y,
      "(^\\s*$|^\\s*\\?\\s*$|\\b(unknown|un|unclear|unsure|not\\s+sure|see\\s+obs|see|photo|picture|n/?a|na|cf)\\b)"
    ) ~ "unknown",
    
    
    TRUE ~ y
  )
  
  y %>%
    str_replace_all("[^\\p{L} ]+", " ") %>%  # keep only letters & spaces
    str_squish() %>%
    str_to_title() %>%
    word(1)                                   # keep the Genus only
}




my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% unique() %>% clean_names()
inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% unique() %>%  clean_names() %>% 
  reframe(site,id,taxon_name,tree_id,
          tree_genus=clean_to_genus(tree_sp),
          tree_sp, growing_site,
          epi_count = as.integer(str_extract(epi_count, "\\d+")),
          dbh_cat = tree_cat, taxon_rank
  ) %>% distinct()

inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  drop_na(latitude,longitude) %>%  clean_names() %>% 
  reframe(id=as.character(id),taxon_name,
          tree_sp,growing_site,
          tree_genus=clean_to_genus(tree_sp),
          epi_count = as.integer(str_extract(epi_count, "\\d+")),
          taxon_rank,  user_login,obs_date
  ) %>% distinct()

inat_obs %>% 
  distinct(id,tree_genus) %>% 
  fwrite("data/processed/tree_genus_clean.csv")
#inat_obs %>% distinct(tree_genus,tree_sp) %>% view()
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
# 0) Prep data -----------
df_host <- my_sites %>%
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
    spp_per10 = 10 * spp_count_dbh   / total_trees,
    dominanttree1 = word(dominanttree1, 1),
    dominanttree2 = word(dominanttree2, 1),
    dominanttree3 = word(dominanttree3, 1)
  ) %>%
  filter(is.finite(lat), is.finite(lon), !(lat == 0 & lon == 0),
         dbh_cat!="0 - 10") %>%
  distinct()


inat_field   <- inat_field %>%
  right_join(df_host) %>% 
  arrange(site)



# 1) KW-test and Dunn's post-hoc -------
## 1a) DBH~Observation count -----------
df_host %>%  
  kruskal_test(obs_per10~dominanttree1)
# dominant tree species has significant influence on observation counts in forests, not in willows

df_host %>% 
  dunn_test(obs_per10~dominanttree1) %>% arrange(p.adj)
# Fagus, Acer, Fraxinus, Pinus, Castanea vs SALIX


## 1b) DBH~species count
df_host %>% 
  kruskal_test(spp_per10~dominanttree1)
# dominant tree species has significant influence on species counts in forests, not in willows

df_host %>% 
  dunn_test(spp_per10~dominanttree1) %>% arrange(p.adj)
# Fagus, Acer, Fraxinus, Pinus, Castanea vs SALIX


# . ------------
# 2) Bayesian ------------
df_bayes <- df_host %>%
  filter(dbh_cat!="0 - 10",
         total_trees > 0) %>% 
  reframe(site,
          obs_count,
          spp_count,
          
          dbh_cat,
          dominanttree1,
          dominanttree2,
          dominanttree3,
          
          # Exposure offset on the log scale (adjusts for how many trees I sampled at a site)
          log_exposure = log(total_trees)
  )

priors_nb <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  set_prior("exponential(1)", class = "shape")  # negbin overdispersion
)


## 2.1) Observation counts -------------
fit_obs_host <- brm(
  formula = obs_count ~ dominanttree1 * dominanttree2 * dominanttree3 + (1 | site) + offset(log_exposure),
  data    = df_bayes,
  family  = negbinomial(),
  prior   = priors_nb,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  seed    = 42
)

# Posterior predictive checks
ppc_obs <- pp_check(fit_obs_host)     # visualize 
ppc_obs

summ_obs <- summary(fit_obs_host)     # summary
summ_obs

as_draws_df(fit_obs_host) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "bterm", values_to = "draw") %>%
  mutate(term = str_remove(bterm, "^b_")) %>%
  group_by(term) %>%
  summarise(
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
  arrange(term) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  fwrite(file.path(OUT_DIR_HOST,"obs.count~host_fixed_effects_summary.csv"))

# Observed zero rate
obs_z <- df_bayes %>%
  summarise(observed_zero_rate = mean(obs_count == 0, na.rm = TRUE))

# Predicted zero rate (overall)
pred_draws <- posterior_predict(fit_obs_host, ndraws = 1000)
pred_z <- tibble(predicted_zero_rate = mean(rowMeans(pred_draws == 0)))

# Bind and write
bind_cols(obs_z, pred_z) %>%
  fwrite(file.path(OUT_DIR_HOST,"obs.count~dbh+setting_zero_rates.csv"))

## 2.2) Species counts -----------
fit_spp_host <- brm(
  formula = spp_count ~ dominanttree1 * dominanttree2 * dominanttree3 + (1 | site) + offset(log_exposure),
  data    = df_bayes,
  family  = zero_inflated_negbinomial(),
  prior   = priors_nb,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  seed    = 43
)

ppc_spp <- pp_check(fit_spp_host)   
ppc_spp
spp_summary <- summary(fit_spp_host)
spp_summary

as_draws_df(fit_spp_host) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "bterm", values_to = "draw") %>%
  mutate(term = str_remove(bterm, "^b_")) %>%
  group_by(term) %>%
  summarise(
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
  arrange(term) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  fwrite(file.path(OUT_DIR_HOST,"spp.count~host_fixed_effects_summary.csv"))

# Observed zero rate
obs_z <- df_bayes %>%
  summarise(observed_zero_rate = mean(spp_count == 0, na.rm = TRUE))

# Predicted zero rate (overall)
pred_draws <- posterior_predict(fit_spp_host, ndraws = 1000)
pred_z <- tibble(predicted_zero_rate = mean(rowMeans(pred_draws == 0)))

# Bind and write
bind_cols(obs_z, pred_z) %>%
  fwrite(file.path(OUT_DIR_HOST,"spp.count~host_zero_rates.csv"))









# . ----------
# 3) Plots ------
## Prep ----------

## 3a) By dominant tree species-----------
df_long <- df_host %>%
  reframe(site,obs_per10,spp_per10,dominanttree1,dominanttree2) %>% 
  pivot_longer(obs_per10:spp_per10, names_to = "count_type", values_to = "count") %>%
  unique() %>% 
  mutate(
    count_type = if_else(str_detect(count_type,"spp"), "Species", "Individuals"),
    ratio = count
  ) %>%
  filter(is.finite(ratio)) %>%
  group_by(count_type) %>%
  mutate(
    # Reorder dominanttree1 *within* each count_type by the mean ratio
    dominanttree1 = fct_reorder(dominanttree1, ratio, .fun = median, .na_rm = TRUE)
  ) %>%
  ungroup()

# Keep legend globally ordered (by global mean) so it doesn't jump around
levels_global <- df_long %>%
  summarise(mu = median(ratio, na.rm = TRUE), .by = dominanttree1) %>%
  arrange(mu) %>%
  pull(dominanttree1) %>%
  as.character()


p_host_obs<-
  ggplot(
    df_long,
    aes(
      x = reorder(count_type, ratio, FUN = median, na.rm = TRUE),
      y = ratio,
      fill = dominanttree1
    )
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 5), # linear-ish below ~10, log-like above
    breaks = pretty_breaks(5)
  ) +
  geom_boxplot(position = position_dodge2(preserve = "single"), 
               width = 0.7, alpha = 0.9, color = "grey20", linewidth = 1) +
  scale_fill_viridis_d("Dominant tree", option = "D", end = 0.9, breaks=levels_global) +
  labs(x = "", y = "Observation count / ten trees") +
  my_theme14 +
  theme(legend.position = "top") 

p_host_obs
# 
# Values are plotted on a pseudo-logarithmic scale (linear below 1, log-like above) 
#to compress large values while preserving detail for small values. 
#This transformation avoids discontinuities and makes patterns across the full range easier to compare.
# 

ggsave(file.path(OUT_DIR_HOST, "count_by_dominant_tree.jpg"),
       p_host_obs, dpi = 500, width = 7, height = 5)


## 3b) By actual host tree-----------
# Helper: title case consistently
to_title <- function(x) str_to_sentence(x, locale = "en")

# Build per-site counts by host genus 
df_genus_counts <- inat_field %>%
  distinct(id, .keep_all = T) %>% 
  
  # CLEAN PUNCTUATION & WHITESPACE 
  mutate(tree_sp = str_replace_all(tree_sp, "[[:punct:]]+", " ")) %>%
  mutate(tree_sp = str_squish(tree_sp)) %>%
  mutate(tree_sp = str_to_sentence(tree_sp)) %>%
  group_by(site) %>% 
  reframe(
    site, tree_id,
    tree_sp = tolower(tree_sp),
    host_genus = to_title(word(tree_sp, 1)),
    id, epi_count
  ) %>%
  filter(!is.na(host_genus), host_genus != "") %>%
  group_by(site,tree_id, host_genus) %>%
  summarise(individuals = sum(epi_count, na.rm = TRUE), 
            .groups = "drop") 

# Lump rare genera 
df_genus_counts_lvl <- df_genus_counts %>% 
  group_by(host_genus) %>% 
  summarise(n=n_distinct(tree_id)) %>% 
  mutate(host_genus_2 = if_else(n <= 2,
                                "Other hosts",
                                host_genus) %>% 
           forcats::fct_inorder() %>%
           forcats::fct_relevel("Other hosts", after = Inf)) %>% 
  right_join(df_genus_counts) %>% 
  group_by(host_genus_2) %>% 
  reframe(n=n_distinct(tree_id),individuals)


# Global order of host_genus by *median* individuals across sites 
levels_global_host <- df_genus_counts_lvl %>%
  ungroup() %>% 
  summarise(med = median(individuals, na.rm = TRUE), .by = host_genus_2) %>%
  arrange(med) %>%
  pull(host_genus_2) %>%
  as.character()

df_genus_counts_lvl <- df_genus_counts_lvl %>% 
  mutate(host_genus_2 = factor(host_genus_2,levels=levels_global_host))



# plot
p_host_by_genus <-
  ggplot(
    df_genus_counts_lvl,
    aes(x = host_genus_2, y = individuals, fill = host_genus_2)
  ) +
  geom_boxplot(width = 0.7, alpha = 0.9, color = "grey20", linewidth = 1) +
  # Compress the tail but keep low counts readable
  scale_y_continuous(
    trans  = pseudo_log_trans(sigma = 5),
    breaks = pretty_breaks(6),
    expand = c(0.05,0)
  ) +
  scale_fill_viridis_d("Host genus", option = "D", end = 0.9, 
                       limits = levels_global_host, drop = FALSE, 
                       guide="none") +
  labs(x = "Host genus", y = "Number of epiphyte individuals") +
  geom_text(aes(x=host_genus_2,y=-1, label=n))+
  my_theme14 +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_host_by_genus

# Save
ggsave(file.path(OUT_DIR_HOST, "individuals_per_site_by_host_genus.jpg"),
       p_host_by_genus, dpi = 500, width = 7.5, height = 5)






## 3c) combine ----------
# Site -> dominant tree map 
df_dom <- df_host %>%
  distinct(site, dominanttree1, dominanttree2,dominanttree3
  )


df_counts <- inat_field %>%
  reframe(
    site,
    host_genus = clean_to_genus(tree_sp),
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    dominanttree1, dominanttree2,dominanttree3
  ) %>%
  filter(!is.na(site), !is.na(host_genus), host_genus != "", !is.na(individuals)) %>%
  group_by(site, host_genus,dominanttree1,dominanttree2,dominanttree3) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop")

# Add zeros for missing site × host_genus combos
all_sites  <- df_counts %>% distinct(site)
all_hosts  <- df_counts %>% distinct(host_genus)
full_grid  <- tidyr::crossing(all_sites, all_hosts)

df_counts_full <- full_grid %>%
  left_join(df_counts, by = c("site", "host_genus")) %>%
  mutate(individuals = replace_na(individuals, 0))

# Within-site totals & shares
df_comp <- df_counts_full %>%
  group_by(site) %>%
  mutate(site_total = sum(individuals, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(site_total > 0) %>%                    # sites with no epiphytes drop
  mutate(share = individuals / site_total)


# Long format for dominant ranks
df_comp_long <- df_comp %>%
  pivot_longer(
    c(dominanttree1, dominanttree2, dominanttree3),
    names_to  = "dom_rank",
    values_to = "dominant_tree"
  ) %>%
  filter(!is.na(dominant_tree), dominant_tree != "") %>%
  mutate(
    dom_rank = recode(dom_rank,
                      dominanttree1 = "1st dominant",
                      dominanttree2 = "2nd dominant",
                      dominanttree3 = "3rd dominant"
    )
  )

# Compute flags and medians
agg <- df_comp_long %>%
  summarise(
    ever      = any(individuals > 0, na.rm = TRUE),  # did this combo ever occur?
    med_share = median(share, na.rm = TRUE),
    n_sites   = n_distinct(site),
    .by = c(dom_rank, dominant_tree, host_genus)
  ) %>%
  mutate(
    med_share_plot = ifelse(ever, med_share, NA_real_)  # NA -> will render white
  )

# Order
host_levels <- agg %>%
  summarise(m = median(med_share, na.rm = TRUE), .by = host_genus) %>%
  arrange(desc(replace_na(m, -Inf))) %>%
  pull(host_genus)

agg <- agg %>%
  mutate(
    host_genus    = factor(host_genus, levels = host_levels),
    dominant_tree = forcats::fct_reorder(dominant_tree, n_sites, .desc = TRUE)
  )

if (!"ever" %in% names(agg)) {
  warning("No 'ever' column detected; deriving from med_share > 0 (conservative).")
  agg <- agg %>%
    mutate(ever = !is.na(med_share) & med_share > 0)
}

# X-axis: global host order by frequency + intensity
w_freq <- 0.7
w_int  <- 0.3

host_order_tbl <- agg %>%
  group_by(host_genus) %>%
  summarise(
    freq_host  = sum(ever, na.rm = TRUE),
    w_int_host = sum(med_share * n_sites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    freq_s  = rescale(freq_host,  to = c(0, 1)),
    inten_s = rescale(w_int_host, to = c(0, 1)),
    score_host = w_freq * freq_s + w_int * inten_s
  ) %>%
  arrange(desc(score_host))

host_levels <- host_order_tbl$host_genus

# Y-axis: per-facet (dom_rank) row order by frequency + intensity
row_order_tbl <- agg %>%
  group_by(dom_rank, dominant_tree) %>%
  summarise(
    freq_row  = sum(ever, na.rm = TRUE),
    w_int_row = sum(med_share * n_sites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(dom_rank) %>%
  mutate(
    freq_s  = rescale(freq_row,  to = c(0, 1)),
    inten_s = rescale(w_int_row, to = c(0, 1)),
    score_row = w_freq * freq_s + w_int * inten_s
  ) %>%
  arrange(dom_rank, desc(score_row))

# Join row scores back for plotting
agg2 <- agg %>%
  left_join(row_order_tbl %>% select(dom_rank, dominant_tree, score_row), by = c("dom_rank","dominant_tree")) %>%
  mutate(
    # Order hosts 
    host_genus = factor(host_genus, levels = host_levels),
    
    # order rows within facet using tidytextreorder_within()
    # (needs a numeric; larger = placed higher by default)
    row_score = score_row,
    
    # Optional: white-out never occurred cells
    med_share_plot = ifelse(ever, med_share, NA_real_)
  )

## 3d) heatmap  ------
p_heat <- ggplot(
  agg2,
  aes(
    x   = host_genus,
    y   = reorder_within(dominant_tree, row_score, within = dom_rank),
    fill = med_share_plot
  )
) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  scale_fill_viridis_c(
    "Median share\nof individuals",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    na.value = "white"  # white-out never-occurred cells
  ) +
  labs(x = "Host genus", y = "Dominant tree at site") +
  facet_grid(dom_rank ~ ., scales = "free_y", space = "free_y") +
  my_theme14 +
  scale_y_reordered() +  # converts the reorder_within() scale back to clean labels
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_heat
ggsave(file.path(OUT_DIR_HOST, "host_share_by_dominant_heatmap_ranked.jpg"),
       p_heat, dpi = 500, width = 9, height = 7)



# 4) Chord diagram host~epi----------
## 4a) Field data only ----------
inat_field %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    site, tree_id,
    host_genus=tree_genus,
    epi_genus = clean_to_genus(taxon_name),
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  filter(!is.na(site), !is.na(host_genus), host_genus != "", !is.na(individuals)) %>%
  group_by(site, host_genus,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, epi_genus) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(epi_genus,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(host_genus) %>% 
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(epi_genus = sub("^E: ", "", to)) %>%
  arrange(epi_genus) %>% 
  #arrange(match(epi_genus, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                      hosts_present)
grid.col  <- c(cols_host, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                              pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
png(file.path(OUT_DIR_HOST, "FIELD_host_epi_chord.png"), width = 5000, height = 5000, res=600)

circos.clear()
circos.par(start.degree = 90,
           track.margin = c(0.01, 0.01),
           gap.after = c(rep(2, length(hosts_present) - 1), 10,
                         rep(2, length(epis_present) - 1), 10))

# Pre-allocate an outer track for custom italic labels
chordDiagram(
  x = links_plot,
  grid.col = grid.col,
  col = links_plot$col,
  transparency = 0,              
  directional = 0,               
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Add italic sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  si   <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  theta_mid <- mean(xlim)
  label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
  circos.text(
    x = theta_mid, y = 0.1, labels = label,
    facing = "clockwise", niceFacing = TRUE,
    adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
  )
}, bg.border = NA)

dev.off()
}




## 4b) My iNat data only ----------

inat_obs %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section"),
         user_login=="marie-ho",
         obs_date>=as.Date("2025-04-01")) %>% 
  reframe(
    host_genus=tree_genus,
    epi_genus = clean_to_genus(taxon_name),
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  group_by(host_genus,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, epi_genus) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(epi_genus,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(host_genus) %>% 
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(epi_genus = sub("^E: ", "", to)) %>%
  arrange(epi_genus) %>% 
  #arrange(match(epi_genus, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                        pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
  png(file.path(OUT_DIR_HOST, "MYiNAT_host_epi_chord.png"), width = 5000, height = 5000, res=600)
  
  circos.clear()
  circos.par(start.degree = 90,
             track.margin = c(0.01, 0.01),
             gap.after = c(rep(2, length(hosts_present) - 1), 10,
                           rep(2, length(epis_present) - 1), 10))
  
  # Pre-allocate an outer track for custom italic labels
  chordDiagram(
    x = links_plot,
    grid.col = grid.col,
    col = links_plot$col,
    transparency = 0,              
    directional = 0,               
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # Add italic sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    theta_mid <- mean(xlim)
    label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
    circos.text(
      x = theta_mid, y = 0.1, labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
    )
  }, bg.border = NA)
  
  dev.off()
}





## 4c) iNat Europe ----------
fread("outputs/04_environmental/EUROPE_field+inat+bg_points.csv") %>% 
  mutate(id=as.character(id)) %>% 
  right_join(inat_obs) %>% 
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    host_genus=tree_genus,
    epi_genus = clean_to_genus(taxon_name),
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(host_genus,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, epi_genus) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(epi_genus,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>5) # keep only links with more than 2 occurrences


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  arrange(host_genus) %>% 
  #arrange(match(host_genus, host_order_all)) %>%
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(epi_genus = sub("^E: ", "", to)) %>%
  arrange(epi_genus) %>% 
  #arrange(match(epi_genus, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                      hosts_present)
grid.col  <- c(cols_host, cols_epi)

# # links coloured by epiphyte, with alpha by strength
# links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)],
#                               pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))

#same alpha always
links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)], 0.7)

# Plot 

png(file.path(OUT_DIR_HOST, "iNat_EUROPE_host_epi_chord.png"), width = 5000, height = 5000, res=600)

circos.clear()
circos.par(start.degree = 90,
           track.margin = c(0.01, 0.01),
           gap.after = c(rep(2, length(hosts_present) - 1), 10,
                         rep(2, length(epis_present) - 1), 10))

# Pre-allocate an outer track for custom italic labels
chordDiagram(
  x = links_plot,
  grid.col = grid.col,
  col = links_plot$col,
  transparency = 0,              
  directional = 0,               
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Add italic sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  si   <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  theta_mid <- mean(xlim)
  label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
  circos.text(
    x = theta_mid, y = 0.1, labels = label,
    facing = "clockwise", niceFacing = TRUE,
    adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
  )
}, bg.border = NA)

dev.off()
#



# 5) Chord diagram host~growing site----------
## 5a) Field data only ----------
inat_field %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    site, tree_id,
    host_genus=tree_genus,
    epi_genus = clean_to_genus(taxon_name),
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  filter(!is.na(site), !is.na(host_genus), host_genus != "", !is.na(individuals)) %>%
  group_by(host_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

gsite_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = gsite_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(host_genus) %>% 
  pull(from)

gsite_present <- links_plot %>%
  distinct(to) %>%
  mutate(growing_site = sub("^E: ", "", to)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = gsite_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(gsite_present)), gsite_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                        pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
  png(file.path(OUT_DIR_HOST, "FIELD_host_growing_site_chord.png"), 
      width = 6000, height = 5300, res=600)
  
  circos.clear()
  circos.par(start.degree = 90,
             track.margin = c(0.01, 0.01),
             gap.after = c(rep(2, length(hosts_present) - 1), 10,
                           rep(2, length(gsite_present) - 1), 10))
  
  # Pre-allocate an outer track for custom italic labels
  chordDiagram(
    x = links_plot,
    grid.col = grid.col,
    col = links_plot$col,
    transparency = 0,              
    directional = 0,               
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # Add italic sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    theta_mid <- mean(xlim)
    label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
    circos.text(
      x = theta_mid, y = 0.1, labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
    )
  }, bg.border = NA)
  
  dev.off()
}



## 5b) My iNat data only ----------

inat_obs %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section"),
         user_login=="marie-ho",
         obs_date>=as.Date("2025-04-01")) %>% 
  reframe(
    host_genus=tree_genus,
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  group_by(host_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(host_genus) %>% 
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(growing_site = sub("^E: ", "", to)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                        pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
  png(file.path(OUT_DIR_HOST, "MYiNAT_host_growing_site_chord.png"), width = 6000, height = 5300, res=600)
  
  circos.clear()
  circos.par(start.degree = 90,
             track.margin = c(0.01, 0.01),
             gap.after = c(rep(2, length(hosts_present) - 1), 10,
                           rep(2, length(epis_present) - 1), 10))
  
  # Pre-allocate an outer track for custom italic labels
  chordDiagram(
    x = links_plot,
    grid.col = grid.col,
    col = links_plot$col,
    transparency = 0,              
    directional = 0,               
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # Add italic sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    theta_mid <- mean(xlim)
    label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
    circos.text(
      x = theta_mid, y = 0.1, labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
    )
  }, bg.border = NA)
  
  dev.off()
}





## 5c) iNat Europe ----------
fread("outputs/04_environmental/EUROPE_field+inat+bg_points.csv") %>% 
  mutate(id=as.character(id)) %>% 
  right_join(inat_obs) %>% 
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    host_genus=tree_genus,
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(host_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>5) # keep only links with more than 2 occurrences


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  arrange(host_genus) %>% 
  #arrange(match(host_genus, host_order_all)) %>%
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(growing_site = sub("^E: ", "", to)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# # links coloured by epiphyte, with alpha by strength
# links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)],
#                               pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))

#same alpha always
links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)], 0.7)

# Plot 

png(file.path(OUT_DIR_HOST, "iNat_EUROPE_host_growing_site_chord.png"), width = 6000, height = 5300, res=600)

circos.clear()
circos.par(start.degree = 90,
           track.margin = c(0.01, 0.01),
           gap.after = c(rep(2, length(hosts_present) - 1), 10,
                         rep(2, length(epis_present) - 1), 10))

# Pre-allocate an outer track for custom italic labels
chordDiagram(
  x = links_plot,
  grid.col = grid.col,
  col = links_plot$col,
  transparency = 0,              
  directional = 0,               
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Add italic sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  si   <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  theta_mid <- mean(xlim)
  label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
  circos.text(
    x = theta_mid, y = 0.1, labels = label,
    facing = "clockwise", niceFacing = TRUE,
    adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
  )
}, bg.border = NA)

dev.off()
#



# 6) Chord diagram host~growing site----------
## 6a) Field data only ----------
inat_field %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    site, tree_id,
    host_genus=tree_genus,
    epi_genus = clean_to_genus(taxon_name),
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  filter(!is.na(site), !is.na(host_genus), host_genus != "", !is.na(individuals)) %>%
  group_by(epi_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(epi_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
epi_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

gsite_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    to = paste0("H: ", factor(epi_genus, levels = epi_order_all)),
    from   = paste0("E: ", factor(growing_site,  levels = gsite_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(epi_genus = sub("^H: ", "", to)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(epi_genus) %>% 
  pull(to)

gsites_present <- links_plot %>%
  distinct(from) %>%
  mutate(growing_site = sub("^E: ", "", from)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(from)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    to = factor(to, levels = epis_present),
    from   = factor(from,   levels = gsites_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_gsites  <- setNames(viridis(length(gsites_present), option = "D"),
                       gsites_present)
grid.col  <- c(cols_gsites, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                        pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
  png(file.path(OUT_DIR_HOST, "FIELD_epis_growing_site_chord.png"), 
      width = 6000, height = 5300, res=600)
  
  circos.clear()
  circos.par(start.degree = 90,
             track.margin = c(0.01, 0.01),
             gap.after = c(rep(2, length(epis_present) - 1), 10,
                           rep(2, length(gsites_present) - 1), 10))
  
  # Pre-allocate an outer track for custom italic labels
  chordDiagram(
    x = links_plot,
    grid.col = grid.col,
    col = links_plot$col,
    transparency = 0,              
    directional = 0,               
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # Add italic sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    theta_mid <- mean(xlim)
    label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
    circos.text(
      x = theta_mid, y = 0.1, labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
    )
  }, bg.border = NA)
  
  dev.off()
}



## 5b) My iNat data only ----------

inat_obs %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section"),
         user_login=="marie-ho",
         obs_date>=as.Date("2025-04-01")) %>% 
  reframe(
    host_genus=tree_genus,
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+"))
  ) %>%
  group_by(host_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(host_genus) %>% 
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(growing_site = sub("^E: ", "", to)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# links coloured by epiphyte, with alpha by strength
links_plot$col <- alpha(grid.col[as.character(links_plot$from)],
                        pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))



# Plot 
{
  png(file.path(OUT_DIR_HOST, "MYiNAT_host_growing_site_chord.png"), width = 6000, height = 5300, res=600)
  
  circos.clear()
  circos.par(start.degree = 90,
             track.margin = c(0.01, 0.01),
             gap.after = c(rep(2, length(hosts_present) - 1), 10,
                           rep(2, length(epis_present) - 1), 10))
  
  # Pre-allocate an outer track for custom italic labels
  chordDiagram(
    x = links_plot,
    grid.col = grid.col,
    col = links_plot$col,
    transparency = 0,              
    directional = 0,               
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # Add italic sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    theta_mid <- mean(xlim)
    label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
    circos.text(
      x = theta_mid, y = 0.1, labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
    )
  }, bg.border = NA)
  
  dev.off()
}





## 5c) iNat Europe ----------
fread("outputs/04_environmental/EUROPE_field+inat+bg_points.csv") %>% 
  mutate(id=as.character(id)) %>% 
  right_join(inat_obs) %>% 
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    host_genus=tree_genus,
    growing_site,
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(host_genus,growing_site) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(host_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(host_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(host_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(host_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>5) # keep only links with more than 2 occurrences


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(host_genus = sub("^H: ", "", from)) %>%
  arrange(host_genus) %>% 
  #arrange(match(host_genus, host_order_all)) %>%
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(growing_site = sub("^E: ", "", to)) %>%
  arrange(growing_site) %>% 
  #arrange(match(growing_site, epi_order_all)) %>%
  pull(to)

# # Sort by remaining abundance 
# host_weight <- links_plot %>% group_by(from) %>% summarise(w = sum(value), .groups = "drop")
# epi_weight  <- links_plot %>% group_by(to)   %>% summarise(w = sum(value), .groups = "drop")
# 
# hosts_present <- host_weight %>% right_join(tibble(from = hosts_present), by = "from") %>%
#   arrange(desc(w), match(from, hosts_present)) %>% pull(from)
# 
# epis_present  <- epi_weight %>% right_join(tibble(to = epis_present), by = "to") %>%
#   arrange(desc(w), match(to, epis_present)) %>% pull(to)

# 5) Apply factors to enforce the filtered order for plotting
links_plot <- links_plot %>%
  mutate(
    from = factor(from, levels = hosts_present),
    to   = factor(to,   levels = epis_present)
  )


# sectors
cols_epi <- setNames(rep("#BDBDBD", length(epis_present)), epis_present)  # light grey
cols_host  <- setNames(viridis(length(hosts_present), option = "D"),
                       hosts_present)
grid.col  <- c(cols_host, cols_epi)

# # links coloured by epiphyte, with alpha by strength
# links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)],
#                               pmin(0.9, 0.25 + 0.75 * (links_plot$value / max(links_plot$value))))

#same alpha always
links_plot$col <- scales::alpha(grid.col[as.character(links_plot$from)], 0.7)

# Plot 

png(file.path(OUT_DIR_HOST, "iNat_EUROPE_host_growing_site_chord.png"), width = 6000, height = 5300, res=600)

circos.clear()
circos.par(start.degree = 90,
           track.margin = c(0.01, 0.01),
           gap.after = c(rep(2, length(hosts_present) - 1), 10,
                         rep(2, length(epis_present) - 1), 10))

# Pre-allocate an outer track for custom italic labels
chordDiagram(
  x = links_plot,
  grid.col = grid.col,
  col = links_plot$col,
  transparency = 0,              
  directional = 0,               
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Add italic sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  si   <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  theta_mid <- mean(xlim)
  label <- gsub("^(H: |E: )", "", si)     # drop the H:/E: prefixes
  circos.text(
    x = theta_mid, y = 0.1, labels = label,
    facing = "clockwise", niceFacing = TRUE,
    adj = c(0, 0.1), cex = 0.8, font = 3  # font=3 → italic
  )
}, bg.border = NA)

dev.off()
#



#-------

