# 04b_latitude_analyses_expected_richness_abundance_maps.R
#
# iNEXT-based standardisation and curves
# - Incidence-based richness (trees as sampling units): expected richness at 100 trees
# - Individuals per 100 trees (rate)
# - Curves: incidence (q=0) and abundance (q=0/1/2), sample-size & coverage
# - Map: colour = richness at (up to) 100 trees; size = individuals / 100 trees; shape = setting

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(forcats)
  library(scales)
  library(maps)
  library(patchwork)
  library(stringr)
  library(cowplot)
  library(grid)
  library(sf)
  library(purrr)
})

# ---- Install/load iNEXT ------------------------
if (!requireNamespace("iNEXT", quietly = TRUE)) {
  install.packages("iNEXT") # requires internet
}
suppressPackageStartupMessages(library(iNEXT)) # ggiNEXT(), iNEXT(), estimateD()

# Parameters --------------------------------------------------------------
IN_MY_SITES     <- "data/processed/my.sites.csv"
IN_INAT_FIELD   <- "data/processed/inat.merged.csv"
OUT_DIR_LATITUDE <- "outputs/04_latitude_analyses"
dir.create(OUT_DIR_LATITUDE, showWarnings = FALSE, recursive = TRUE)

# Load data ---------------------------------------------------------------
my_sites <- fread(IN_MY_SITES) %>%
  as_tibble() %>%
  reframe(
    site,
    dbh_cat = TreeCat,
    setting = setting,                 
    obs_count = n_epis_per_site,       
    spp_count = n_sp_per_site,         
    obs_count_dbh = n_epis_per_site_cat,
    spp_count_dbh = n_sp_per_site_cat,
    total_trees = NumberTrees,         
    lat = siteLat,
    lon = siteLon
  ) %>%
  filter(abs(lat) > 0, abs(lon) > 0)

inat_field <- fread(IN_INAT_FIELD) %>%
  as_tibble() %>%
  reframe(
    site,
    dbh_cat = TreeCat,
    setting = Setting,
    obs_count = n_epis_per_site,
    spp_count = n_sp_per_site,
    obs_count_dbh = n_epis_per_site_cat,
    spp_count_dbh = n_sp_per_site_cat,
    total_trees = NumberTrees,
    lat = siteLat,
    lon = siteLon,
    id, taxon_id, taxon.name, taxon.rank,
    quality_grade,
    count = EpiCount_noNA,             
    tree_id = TreeID
  ) %>%
  filter(abs(lat) > 0, abs(lon) > 0) %>%
  unique()

# Standard ggplot theme ---------------------------------------------------
my_theme14 <-
  theme_classic(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(0.8, "lines"),
    legend.position  = "top",
    axis.line        = element_blank(),
    axis.ticks       = element_blank()
  )

# . ----------------------------------------------------------------------
# 0) Field work maps scaffolding (sites) ---------------------------------
world_df <- map_data("world")

# sum across DBH groups for site totals
sites_map <- my_sites %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(is.finite(lon), is.finite(lat)) %>%
  group_by(site) %>%
  summarise(
    lon         = first(lon),
    lat         = first(lat),
    n_sp        = first(spp_count, na.rm = TRUE),   
    n_epis      = first(obs_count, na.rm = TRUE),   
    total_trees = sum(total_trees, na.rm = TRUE), 
    setting     = first(na.omit(setting)),
    .groups = "drop"
  ) %>%
  filter(setting != "" & !is.na(setting)) %>%
  unique()

# . ----------------------------------------------------------------------
# 1) Build iNEXT inputs (incidence-frequency + abundance) -----------------
# Incidence-frequency: vector per site = c(T, y1, ..., yS)
incidence_freq_tbl <- inat_field %>%
  distinct(site, tree_id, taxon_id) %>%
  count(site, taxon_id, name = "y") # y = # of trees with species i

T_by_site <- sites_map %>% select(site, T = total_trees)

# Sites with at least one species incidence
sites_with_species <- unique(incidence_freq_tbl$site)

# iNEXT input for sites with species and T >= 2
inc_freq_list <- incidence_freq_tbl %>%
  group_by(site) %>%
  summarise(freq = list(y), .groups = "drop") %>%
  right_join(T_by_site, by = "site") %>%
  filter(site %in% sites_with_species, is.finite(T), T >= 2) %>%
  mutate(vec = map2(T, freq, ~ c(.x, if (is.null(.y)) numeric(0) else .y))) %>%
  { stats::setNames(.$vec, .$site) }

# Abundance: individuals per species per site (for q=0/1/2 curves)
abund_list <- inat_field %>%
  mutate(count_use = ifelse(is.na(count) | count <= 0, 1L, as.integer(round(count)))) %>%
  group_by(site, taxon_id) %>%
  summarise(n_ind = sum(count_use), .groups = "drop") %>%
  split(.$site) %>%
  lapply(function(df) { v <- df$n_ind; names(v) <- df$taxon_id; v })

# . ----------------------------------------------------------------------
# 2) Expected richness at 100 trees (no bootstrapping; avoid extrapolation when S <= 1)
m_target <- 100L

# Extract expected richness at m trees (incidence rarefaction) for m <= T
expected_richness_incidence <- function(T, y_vec, m) {
  if (length(y_vec) == 0L || m < 1L || T < 1L) return(0)
  m <- min(m, T)
  # E[S_m] = sum_i [ 1 - C(T - y_i, m) / C(T, m) ]; use lchoose for numeric stability
  sum(1 - exp(lchoose(T - y_vec, m) - lchoose(T, m)))
}

safe_estimate_incidence <- function(v, m_target) {
  # v is c(T, y1, ..., yS)
  T <- as.integer(v[1])
  y <- as.integer(v[-1])
  S <- length(y)
  
  # Rule: if S <= 1 → do NOT extrapolate; else allow up to 2T
  m_eff <- if (S <= 1L) min(m_target, T) else min(m_target, 2L * T)
  m_eff <- max(1L, m_eff)
  
  if (m_eff <= T) {
    # Exact rarefaction (no iNEXT call, no bootstrap)
    tibble::tibble(spp_exp = expected_richness_incidence(T, y, m_eff),
                   level_used = m_eff)
  } else {
    # Extrapolation via iNEXT but with NO bootstrap
    out <- iNEXT::estimateD(list(v), q = 0, datatype = "incidence_freq",
                            base = "size", level = as.integer(m_eff),
                            nboot = 0)
    tibble::tibble(spp_exp = out$qD, level_used = out$t)
  }
}

if (length(inc_freq_list)) {
  est_list <- purrr::imap(inc_freq_list, function(v, s) {
    res <- safe_estimate_incidence(v, m_target)
    tibble::tibble(site = s, spp_exp100 = res$spp_exp, level_used = res$level_used)
  })
  spp100_df <- dplyr::bind_rows(est_list)
} else {
  spp100_df <- tibble::tibble(site = character(), spp_exp100 = numeric(), level_used = integer())
}
# . ----------------------------------------------------------------------
# 3) Join estimates back to ALL sites (keep zero sites) -------------------
sites_map <- sites_map %>%
  left_join(spp100_df, by = "site") %>%
  mutate(
    spp_exp100    = ifelse(is.na(spp_exp100), 0, spp_exp100),
    level_used    = ifelse(is.na(level_used), 0L, level_used),
    indiv_per10   = 10  * n_epis / total_trees,
    indiv_per100  = 100 * n_epis / total_trees,
    reached_100   = level_used >= m_target  # TRUE if richness truly at 100 trees
  )

# . ----------------------------------------------------------------------
# 4) Map: colour = expected richness at (up to) 100 trees; size = indiv/100; shape = setting
sp_max <- max(sites_map$spp_exp100, na.rm = TRUE)
sp_breaks <- c(0, 1, 2, 5, 10, 30, 50)
sp_breaks <- sp_breaks[sp_breaks <= sp_max]
if (length(sp_breaks) == 0) sp_breaks <- c(0)

obs_max <- max(sites_map$indiv_per100, na.rm = TRUE)
obs_breaks <- c(0, 10, 25, 50, 100, 200, 500)
obs_breaks <- obs_breaks[obs_breaks <= obs_max]
if (length(obs_breaks) == 0) obs_breaks <- c(0)


p_site_inext100 <- ggplot() +
  geom_polygon(data = world_df, aes(long, lat, group = group),
               fill = "grey93", color = "grey80", linewidth = 0.2) +
  geom_point(
    data = sites_map,
    aes(lon, lat,
        color = spp_exp100,
        size  = indiv_per100,
        shape = setting),
    alpha = 0.8, stroke = 0.2
  ) +
  scale_color_viridis(
    begin = 0.20, end = 1,
    trans = scales::log1p_trans(),
    breaks = sp_breaks, labels = sp_breaks,
    name = "Richness\n(expected in \n100 trees)"
  ) +
  scale_size(
    range = c(2, 12),
    breaks = obs_breaks, labels = obs_breaks,
    name = "Individuals /\n100 trees",
    guide = "legend"
  ) +
  scale_shape(name = "Setting") +
  coord_quickmap(
    xlim = range(sites_map$lon, na.rm = TRUE) + c(-0.5, 0.5),
    ylim = range(sites_map$lat, na.rm = TRUE) + c(-0.5, 0.5)
  ) +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  guides(
    colour = guide_colourbar(
      direction      = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth       = unit(80, "pt"),
      barheight      = unit(10, "pt"),
      order          = 1
    ),
    size = guide_legend(
      title.position = "top",
      nrow           = 3, byrow = TRUE,
      order          = 2
    ),
    shape = guide_legend(
      title.position = "top",
      override.aes   = list(size = 4),
      order          = 3
    )
  ) +
  my_theme14 +
  theme(
    legend.position  = "right",
    legend.box       = "vertical",
    legend.direction = "vertical"
  )

# Draw & save map ---------------------------------------------------------
p_site_inext100

ggsave(file.path(OUT_DIR_LATITUDE, "Fig.1_sites_map_iNEXT_100trees.svg"),
       p_site_inext100, width = 8, height = 6, dpi = 500)

# . ----------------------------------------------------------------------
# 5) iNEXT curves ---------------------------------------------------------
# 5A) Incidence-based richness curves (q = 0)
# if (length(inc_freq_list)) {
#   inext_inc <- iNEXT(inc_freq_list, q = 0, datatype = "incidence_freq")  # default SE bands
#   
#   p_inc_size <- ggiNEXT(inext_inc, type = 1,
#                         aes(shape=setting)) +
#     scale_color_viridis_d(name = "Site") +
#     labs(x = "Number of trees (sampling units)", y = "Expected richness (q = 0)") +
#     my_theme14 +
#     theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical") +
#     geom_vline(xintercept = 100, linetype = 2, colour = "grey40")
#   
#   p_inc_cov <- ggiNEXT(inext_inc, type = 3) +
#     scale_color_viridis_d(name = "Site") +
#     labs(x = "Sample coverage", y = "Expected richness (q = 0)") +
#     my_theme14 +
#     theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical")
#   
#   ggsave(file.path(OUT_DIR_LATITUDE, "iNEXT_incidence_size_curves.png"),
#          p_inc_size, width = 10, height = 7, dpi = 300)
#   ggsave(file.path(OUT_DIR_LATITUDE, "iNEXT_incidence_coverage_curves.png"),
#          p_inc_cov, width = 10, height = 7, dpi = 300)
# }
# 
# # 5B) Abundance-based diversity curves (q = 0/1/2)
# if (length(abund_list)) {
#   inext_abund <- iNEXT(abund_list, q = c(0,1,2), datatype = "abundance")  # default SE bands
#   
#   p_abund_size <- ggiNEXT(inext_abund, type = 1) +
#     scale_color_viridis_d(name = "Site") +
#     labs(x = "Individuals (sample size)", y = "Diversity (Hill numbers, q = 0/1/2)") +
#     my_theme14 +
#     theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical")
#   
#   p_abund_cov <- ggiNEXT(inext_abund, type = 3) +
#     scale_color_viridis_d(name = "Site") +
#     labs(x = "Sample coverage", y = "Diversity (Hill numbers, q = 0/1/2)") +
#     my_theme14 +
#     theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical")
# 
# }

# . ----------------------------------------------------------------------
# 6) Export site summary --------------------------------------------------
sites_out <- sites_map %>%
  select(site, lon, lat, setting, total_trees, n_epis, n_sp,
         spp_exp100, level_used, reached_100, indiv_per10, indiv_per100)


fwrite(sites_out, file.path(OUT_DIR_LATITUDE, "sites_inext_summary_100trees.csv"))
