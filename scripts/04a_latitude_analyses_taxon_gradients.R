# 04_latitude_analyses_taxon_gradients.R
#

# Create maps of field sites
# Latitudinal patterns of species richness and observation counts of accidental epiphytes
# Determine dominant species and see how their distribution varies/shifts
#
# 0) Field work maps 
# 0a) "raw" numbers
# 0b) Effort-adjusted (per 10 trees)
#
# Analyses
# 1) Field work data only
# 2) All my 2025 season iNat observations
# 3) All iNat accidental epiphyte observations
#

rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)     # for stat_cor (optional)
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
  library(brms)
  library(sf)
  library(ggtext)
})

# Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"
IN_OCC_COVS   <- "outputs/04_environmental/env_vars_sel.csv"

OUT_DIR_LATITUDE  <- "outputs/04_latitude_analyses"
dir.create(OUT_DIR_LATITUDE,  showWarnings = FALSE, recursive = TRUE)

# Load data --------
my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% 
  reframe(site, dbh_cat=TreeCat,
          setting,
         obs_count=n_epis_per_site,
         spp_count=n_sp_per_site,
         obs_count_dbh=n_epis_per_site_cat,
         spp_count_dbh=n_sp_per_site_cat,
         total_trees=NumberTrees,
         lat=siteLat,
         lon=siteLon) %>% 
  filter(abs(lat)>0,abs(lon)>0)

# inat_field %>% distinct(taxon.rank)

inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% 
  filter(taxon.rank%in%c("species","genus","complex"),
         #!is.na(community_taxon_id) | quality_grade=="research"
         ) %>% #glimpse()
  reframe(site, dbh_cat=TreeCat,
          setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          tree_id = TreeID) %>% 
  filter(abs(lat)>0,abs(lon)>0) %>% 
  unique()


my_inat   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  filter(user_login=="marie-ho") %>% 
  reframe(lat=latitude,
          lon=longitude,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          family_name,
          community_taxon_id) %>% 
  unique()

inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5,
          lat=latitude,
          lon=longitude,
          setting=Setting) %>% 
  drop_na(lat,lon)

field_lat_mean <- mean(unique(my_sites$lat), na.rm = TRUE)


# Standard ggplot theme 
my_theme14 <- 
  theme_classic(base_size = 14) +
  
  theme(
    # Add a box around the panel(s)
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # Facet strip styling (kept from your snippet)
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(0.8, "lines"),
    
    # Legend
    legend.position  = "top",
    
    # Remove ALL axis lines (and ticks, if desired)
    axis.line        = element_blank(),
    axis.ticks       = element_blank()
  )

# . --------
#
# 1) iNat observations ~ latitude --------------
# Site-level metrics inside each 1° bin
inat_by_site <- inat_obs %>%
  filter(setting%in%c("Urban","Rural")) %>%
  filter(lat<55&lat>30, lon>-10&lon<6) %>% 
  mutate(
    lat_mid5 = floor(lat / 1) * 1 + 0.5,
    lon_mid5 = floor(lon / 1) * 1 + 0.5
  ) %>%
  group_by(lat_mid5) %>%
  reframe(lat,lon,setting,
          n_obs_site = n_distinct(id),
          n_count=sum(EpiCount_noNA),
          n_tax_site = n_distinct(taxon_id[taxon.rank == "species"], na.rm = TRUE)/n_obs_site,
          n_tax_site_all = n_distinct(taxon_id[taxon.rank == "species"], na.rm = TRUE)
   
  ) %>% unique()

# Summary (mean ± SE) per bin
inat_lat_se <- inat_by_site %>%
  group_by(lat_mid5) %>%
  reframe(
    n_sit   = n(),                              # number of sites in bin
    mean_tax = mean(n_tax_site),
    se_tax   = sd(n_tax_site) / sqrt(n_sit),
    mean_obs = mean(n_obs_site),
    se_obs   = sd(n_obs_site) / sqrt(n_sit)
  )


# PLOT ------
## 1A) Scatter: Species richness vs latitude ----------------------------
p_lat_species <- 
  ggplot(inat_by_site, aes(y = lat, x = n_tax_site,col=setting)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 6, stroke = 0.4,alpha=0.7) +
   labs(
    x = "Species / observation",
    y = "Latitude (°)"
  ) +
  scale_colour_viridis(discrete=T,option="G")+
  my_theme14+
  theme(legend.position = "top")



# ## 1B) Scatter: Observation count vs latitude ----------------------------
# 
# p_lat_obs <- 
#   ggplot(inat_by_site, aes(y = lat, x = n_obs_site, col=setting)) +
#   geom_smooth(method = "lm", se = FALSE) +
#   geom_point(size = 6, stroke = 0.4, alpha=0.7) +
#   labs(
#     x = "Observed individuals",
#     y = "Latitude (°)"
#   ) +
#   scale_colour_viridis(discrete=T,option="G")+
#   my_theme14
# 
# p_lat_combined <-
#   p_lat_obs + p_lat_species + theme(legend.position = "none")+
#   patchwork::plot_layout(axes = "collect")
# p_lat_combined
# 
# # ggsave(file.path(OUT_DIR_LATITUDE, "species_vs_latitude_merged.jpg"),
# #        p_lat_combined, width = 8, height = 6, dpi = 300)
# 

# 2) Species composition ~ latitude----------

## 2.1) All my iNat observations ---------
inat_field %>% 
  select(-setting) %>% 
  full_join(my_sites %>% select(site,setting) %>% distinct()) %>% 
  distinct() ->inat_field_join

my_inat_by_lat <-my_inat %>% 
  mutate(setting=if_else(setting=="Rural","forest","non-forest")) %>% 
  distinct() %>% 
  filter(!id%in%c(inat_field_join$id)) %>% 
  bind_rows(inat_field_join) %>% 
  mutate(
    lat_mid5 = floor(lat / 1) * 1 + 0.5,
    lon_mid5 = floor(lon / 1) * 1 + 0.5
  ) %>%
  filter(lat<55&lat>30, lon>-10&lon<6) %>% 
  group_by(lat_mid5,lon_mid5, setting) %>%
  mutate(
    n_tax_site_research = n_distinct(taxon_id[quality_grade=="research"], na.rm = TRUE),
    n_tax_site = n_distinct(taxon_id),
    n_obs_site = n_distinct(id),
    n_count = sum(count)
  )

my_inat_by_lat %>% 
  ungroup() %>% 
  select(taxon.rank) %>% 
  unique()

# Get dominant species

genus_species_n <- my_inat_by_lat %>%
  ungroup() %>% 
  mutate(taxon.name=if_else(family_name=="Poaceae","Poaceae",taxon.name)) %>% 
  filter(taxon.rank == "species"| family_name=="Poaceae") %>%
  mutate(genus = word(taxon.name, 1)) %>%
  distinct(genus, taxon_id) %>%
  count(genus, name = "n_species")


dominant_species_df <-
  my_inat_by_lat %>% 
  mutate(
    lat_mid5 = floor(lat / 1) * 1 + 0.5,
    lon_mid5 = floor(lon / 2) * 2 + 0.5) %>% 
  mutate(taxon.rank=if_else(family_name=="Poaceae","complex",taxon.rank),
         taxon.name=if_else(family_name=="Poaceae","Poaceae",taxon.name)) %>% 
  mutate(taxon.name=word(taxon.name,1)) %>% 
  group_by(lat_mid5,lon_mid5, setting, taxon.name) %>% 

  filter(taxon.rank%in%c("species","genus","complex","subgenus","subspecies"),
         #!is.na(community_taxon_id)|quality_grade=="research"
         ) %>% 
  reframe(sum=sum(count),taxon.name) %>% unique() %>% 
  group_by(lon_mid5,lat_mid5,setting) %>% 
  distinct() %>%
  mutate(rank = dense_rank(desc(sum))) %>% 
  group_by(lat_mid5,lon_mid5,setting) %>% 
  arrange(rank) %>% 
  slice_head(n=10) %>% 
  arrange(desc(lat_mid5), taxon.name) %>%
  ungroup() %>% 
  mutate(taxon.name = factor(taxon.name, levels = unique(taxon.name),ordered=T)) %>% 
  
  left_join(genus_species_n, by = c("taxon.name" = "genus")) %>%
  mutate(
      taxon.label = paste0(taxon.name, "<sup>",n_species,"</sup>"),
      taxon.label = factor(
        taxon.label,
        levels = unique(taxon.label),
        ordered = TRUE
      )
  )



dominant_species_df %>% 
  ungroup() %>% 
  select(taxon.label) %>% unique()

# ### 2.1a) MAP ------------
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# world_df <- map_data("world")
# 
# ### Create base map ------------
# base_map <- 
#   ggplot()+
#   geom_polygon(
#     data = world_df, aes(long, lat, group = group),
#     fill = "grey93", color = "grey80", linewidth = 0.2
#   ) +
#   coord_quickmap(
#     xlim = range(my_inat_by_lat$lon, na.rm = TRUE) + c(-0.5, 0.5),
#     ylim = range(my_inat_by_lat$lat, na.rm = TRUE) + c(-0.5, 0.5)
#   )  + 
#   labs(x="Longitude (°)",y="Latitude (°)",
#        caption = "''R'' = Rural, ''U'' = Urban")+
#   my_theme14+
#   
#   theme(
#     plot.background = element_blank(),  
#     panel.background = element_blank()
#   )
# 


### Create legend for all plots -------
# Global, fixed ordering for taxa (already set on dominant_species_df)
taxa_levels <- levels(unique(dominant_species_df$taxon.label))

# Create a named palette with viridis
taxa_pal <- setNames(
  viridis::viridis(length(taxa_levels), option = "D"),
  taxa_levels
)
# 
# 
# p_legend <- ggplot(dominant_species_df, aes(x = setting, y = sum, fill = taxon.name)) +
#   geom_col(position = "fill") +
#   
#   # Use the same global palette and order everywhere
#   scale_fill_manual("Genus name",
#                     values = taxa_pal,
#                     limits = taxa_levels,
#                     guide  = guide_legend(ncol = 3)
#   ) +
#   theme_minimal()+
#   theme(legend.position = "right",
#         legend.text = element_text(face = "italic"),
#         legend.margin        = margin(0, 0, 0, 0),
#         legend.box.margin    = margin(0, 0, 0, 0),
#         legend.box.spacing   = unit(0, "pt"),
#         legend.spacing       = unit(0, "pt"),
#         legend.spacing.y     = unit(0, "pt"),
#         plot.margin          = margin(0, 0, 0, 0)
#   )
# 
# 
# legend <- get_legend(p_legend)
# legend_plot <- ggdraw(legend)
# 
# 
# ### Create individual barplots ----------
# 
# make_bar <- function(df) {
#   df %>% 
#     mutate(setting=if_else(setting=="forest","F","N")) %>% 
#     
#     ggplot(aes(x = setting, y = sum, fill = taxon.name)) +
#     geom_col(position = "fill") +
#     
#     scale_fill_manual(
#       values = taxa_pal,
#       limits = taxa_levels,
#       drop = FALSE,
#       guide = "none"  
#     ) +
#     
#     theme_minimal(base_size = 8) +
#     labs(x = "", y = "") +
#     theme(
#       panel.background = element_rect(fill = NA, color = NA),  
#       plot.background = element_rect(fill = NA, color = NA),    
#       legend.position = "none",
#       axis.text.y = element_blank(),
#       axis.ticks = element_blank(),
#       panel.grid = element_blank()
#     )
# }
# 
# # Split data
# plot_data <- dominant_species_df %>%
#   group_by(lat_mid5, lon_mid5) %>%
#   group_split()
# 
# # Make the inset plots (one per group)
# plots <- purrr::map(plot_data, make_bar)
# 
# # Extract coords in the same order as 'plots'
# coords <- purrr::map_dfr(
#   plot_data,
#   ~ tibble(lon_mid5 = unique(.x$lon_mid5),
#            lat_mid5 = unique(.x$lat_mid5))
# )
# 
# # Build the tibble with a list-column of grobs
# inset_df <- coords %>%
#   mutate(subview = purrr::map(plots, cowplot::as_grob))
# 
# ### Combine & save --------
# map_no_legend<-
#   base_map +
#   geom_subview(
#     data = inset_df,
#     aes(x = lon_mid5, y = lat_mid5, subview = subview),
#     width  = rel(1),
#     height = rel(1)   # center on (x, y)
#   )+
#   theme(
#     plot.margin = margin(0, 0, 0, 0),
#     panel.spacing = unit(0, "pt")
#     
#   )
# 
# final <- plot_grid(
#   map_no_legend,
#   legend_plot,
#   nrow = 1,
#   rel_widths = c(1, 0.28),   # adjust the legend column width
#   align = "h"
# )
# 
# final
# 
# ggsave(file.path(OUT_DIR_LATITUDE, "myiNAT_species_gradient_map.jpeg"),
#        final, width = 12, height = 8, dpi = 600)
# 
# ggsave(file.path(OUT_DIR_LATITUDE, "myiNat_species_gradient_map.svg"),
#        final, width = 12, height = 8, dpi = 600)


# . --------------


### 2.1b) BAR PLOT########
df_segments <- dominant_species_df %>%
  filter(setting!="") %>% 
  group_by(lat_mid5, setting, taxon.label) %>%
  summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
  group_by(lat_mid5, setting) %>%
  mutate(
    total = sum(sum, na.rm = TRUE),
    prop  = if_else(total > 0, sum / total, 0)
  ) %>%
  ungroup() %>%
  # Fix factor levels so colors remain consistent
  mutate(taxon.label = factor(taxon.label, levels = taxa_levels, ordered = TRUE)) %>%
  # Order segments by taxon for reproducibility (or pick another order)
  arrange(setting, lat_mid5, taxon.label) %>%
  group_by(setting, lat_mid5) %>%
  # Build cumulative horizontal allocation [0..1]
  mutate(
    xmin = cumsum(dplyr::lag(prop, default = 0)),
    xmax = cumsum(prop)
  ) %>%
  ungroup() %>% 
  unique()


# Row thickness in latitude units (lat bins are 1°, so 0.8 leaves small gaps)
lat_bin_height <- 1

p_lat_segments <- ggplot(df_segments) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = lat_mid5 - lat_bin_height/2,
      ymax = lat_mid5 + lat_bin_height/2,
      fill = taxon.label
    ),
    color = NA
  ) +
  facet_grid(. ~ setting) +
  scale_x_continuous(
    name   = "Relative dominance",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous("Latitude",expand=F) +  
  scale_fill_manual("Taxon name",
                    values = taxa_pal,
                    limits = taxa_levels,
                    drop   = FALSE,
                    guide  = guide_legend(ncol = 3) 
  )+
  my_theme14+
  theme(
    legend.position     = "right",
    legend.text         = element_markdown(face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(5, 5, 5, 5)
  )+
  geom_text(aes(x=0.05,y=lat_mid5,label=total),col="white",size=3)

p_lat_segments


ggsave(file.path(OUT_DIR_LATITUDE, "Fig.2_MyiNat_species_latitude_bar.png"),
      p_lat_segments, width = 12, height = 6, dpi = 600)

ggsave(file.path(OUT_DIR_LATITUDE, "Fig.2_MyiNat_species_latitude_bar.svg"),
       p_lat_segments, width = 12, height = 6, dpi = 600)



# 
# 
# # .#######
# ## 2.2) FIELD SITES ----------
# dominant_species_field <-
#   inat_field %>% 
#   select(-setting) %>% 
#   full_join(my_sites %>% select(site,setting) %>% unique()) %>% 
#   drop_na(taxon.name) %>% 
#   unique() %>% 
#   mutate(
#     lat_mid5 = floor(lat / 1) * 1 + 0.5,
#     lon_mid5 = floor(lon/ 1) * 1 + 0.5) %>% 
#   mutate(genus_name=if_else(str_detect(taxon.name,"Polypodium")|
#                               str_detect(taxon.name,"Polypodiales"),
#                             "Polypodium sp.",stringr::word(taxon.name, 1))) %>% 
#   group_by(lat_mid5,lon_mid5, setting,taxon.name) %>% 
#   #filter(taxon.rank%in%c("species","genus","complex")) %>% 
#   reframe(sum=sum(count),genus_name) %>% unique() %>% 
#   group_by(lon_mid5,lat_mid5,setting) %>% 
#   distinct() %>%
#   mutate(rank = dense_rank(desc(sum))) %>% 
#   group_by(lat_mid5,lon_mid5,setting) %>% 
#   arrange(rank) %>% 
#   slice_head(n=10) %>% 
#   arrange(desc(lat_mid5), genus_name) %>%
#   ungroup() %>% 
#   mutate(genus_name = factor(genus_name, levels = unique(genus_name),ordered=T)) 
# 
# 
# dominant_species_field %>% 
#   ungroup() %>% 
#   select(genus_name) %>% unique()
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# ### 2.2a) MAP ---------
# ### Create base map ------------
# base_map <- 
#   ggplot()+
#   geom_polygon(
#     data = world_df, aes(long, lat, group = group),
#     fill = "grey93", color = "grey80", linewidth = 0.2
#   ) +
#   coord_quickmap(
#     xlim = range(my_inat_by_lat$lon, na.rm = TRUE) + c(-0.5, 0.5),
#     ylim = range(my_inat_by_lat$lat, na.rm = TRUE) + c(-0.5, 0.5)
#   )  + 
#   labs(x="Longitude (°)",y="Latitude (°)",
#        caption = "''F'' = Forest, ''N'' = Non-forest")+
#   my_theme14+
#   
#   theme(
#     plot.background = element_blank(),  
#     panel.background = element_blank()
#   )
# 
# 
# ### Create legend for all plots ---------
# # Global, fixed ordering for taxa (already set on dominant_species_df)
# taxa_levels_field <- levels(unique(dominant_species_field$genus_name))
# 
# # Create a named palette with viridis
# taxa_pal_field <- setNames(
#   viridis::viridis(length(taxa_levels_field), option = "D"),
#   taxa_levels_field
# )
# 
# 
# p_legend_field <- ggplot(dominant_species_field, aes(x = setting, y = sum, fill = genus_name)) +
#   geom_col(position = "fill") +
#   
#   # Use the same global palette and order everywhere
#   scale_fill_manual("Taxon name (by order of appearance from N to S)",
#                     values = taxa_pal_field,
#                     limits = taxa_levels_field,
#                     guide  = guide_legend(ncol = 2)
#   ) +
#   theme_minimal()+
#   theme(legend.position = "right",
#         legend.text = element_text(face = "italic"),
#         legend.margin        = margin(0, 0, 0, 0),
#         legend.box.margin    = margin(0, 0, 0, 0),
#         legend.box.spacing   = unit(0, "pt"),
#         legend.spacing       = unit(0, "pt"),
#         legend.spacing.y     = unit(0, "pt"),
#         plot.margin          = margin(0, 0, 0, 0)
#   )
# 
# 
# legend_field <- get_legend(p_legend_field)
# legend_field_plot <- ggdraw(legend_field)
# 
# 
# ###Create individual barplots ####
# # Function to create a bar plot for one location
# make_bar <- function(df) {
#   df %>% 
#     mutate(setting=if_else(setting=="Forest","F","N")) %>% 
#     
#     ggplot(aes(x = setting, y = sum, fill = genus_name)) +
#     geom_col(position = "fill") +
#     
#     scale_fill_manual(
#       values = taxa_pal_field,
#       limits = taxa_levels_field,
#       drop = FALSE,
#       guide = "none"  
#     ) +
#     
#     theme_minimal(base_size = 8) +
#     labs(x = "", y = "") +
#     theme(
#       panel.background = element_rect(fill = NA, color = NA),  
#       plot.background = element_rect(fill = NA, color = NA),    
#       legend.position = "none",
#       axis.text.y = element_blank(),
#       axis.ticks = element_blank(),
#       panel.grid = element_blank()
#     )
# }
# 
# # 1) Split data
# plot_data_field <- dominant_species_field %>%
#   group_by(lat_mid5, lon_mid5) %>%
#   group_split()
# 
# # 2) Make the inset plots (one per group)
# plots_field <- purrr::map(plot_data_field, make_bar)
# 
# # 3) Extract coords in the same order as 'plots'
# coords <- purrr::map_dfr(
#   plot_data_field,
#   ~ tibble(lon_mid5 = unique(.x$lon_mid5),
#            lat_mid5 = unique(.x$lat_mid5))
# )
# 
# # 4) Build the tibble with a list-column of grobs
# inset_df_field <- coords %>%
#   mutate(subview = purrr::map(plots_field, cowplot::as_grob))
# 
# ### Combine & save ---------
# map_no_legend_field<-
#   base_map +
#   geom_subview(
#     data = inset_df_field,
#     aes(x = lon_mid5, y = lat_mid5, subview = subview),
#     width  = rel(1),
#     height = rel(1)   # center on (x, y)
#   )+
#   theme(
#     plot.margin = margin(0, 0, 0, 0),
#     panel.spacing = unit(0, "pt")
#     
#   )
# 
# final_field <- plot_grid(
#   map_no_legend_field,
#   legend_field_plot,
#   nrow = 1,
#   rel_widths = c(1, 0.28),   # adjust the legend column width
#   align = "h"
# )
# 
# final_field
# # 
# # ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_genus_gradient_map.jpeg"),
# #        final_field, width = 12, height = 8, dpi = 600)
# # 
# # ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_genus_gradient_map.svg"),
# #        final_field, width = 12, height = 8, dpi = 600)
# # 
# 
# 
# # . -----------
# ## 2.2b) BAR PLOT ----------
# df_segments_field <- dominant_species_field %>%
#   group_by(lat_mid5, setting, genus_name) %>%
#   summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
#   group_by(lat_mid5, setting) %>%
#   mutate(
#     total = sum(sum, na.rm = TRUE),
#     prop  = if_else(total > 0, sum / total, 0)
#   ) %>%
#   ungroup() %>%
#   # Fix factor levels so colors remain consistent
#   mutate(taxon.name = factor(genus_name, levels = taxa_levels_field, ordered = TRUE)) %>%
#   # Order segments by taxon for reproducibility (or pick another order)
#   arrange(setting, lat_mid5, genus_name) %>%
#   group_by(setting, lat_mid5) %>%
#   # Build cumulative horizontal allocation [0..1]
#   mutate(
#     xmin = cumsum(dplyr::lag(prop, default = 0)),
#     xmax = cumsum(prop)
#   ) %>%
#   ungroup() %>% 
#   unique()
# 
# 
# # Row thickness in latitude units (lat bins are 2°, so 1.8 leaves small gaps)
# lat_bin_height <- 1
# 
# p_lat_segments_field <- ggplot(df_segments_field) +
#   geom_rect(
#     aes(
#       xmin = xmin,
#       xmax = xmax,
#       ymin = lat_mid5 - lat_bin_height/2,
#       ymax = lat_mid5 + lat_bin_height/2,
#       fill = genus_name
#     ),
#     color = NA
#   ) +
#   facet_grid(. ~ setting) +
#   scale_x_continuous(
#     name   = "Relative dominance",
#     limits = c(0, 1),
#     breaks = seq(0, 1, by = 0.25),
#     labels = label_percent(accuracy = 1),
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous("Latitude",expand=F) +  
#   scale_fill_manual("Genus",
#                     values = taxa_pal_field,
#                     limits = taxa_levels_field,
#                     drop   = FALSE,
#                     guide  = guide_legend(ncol = 2) 
#   )+
#   my_theme14+
#   theme(
#     legend.position     = "right",
#     legend.text         = element_text(face = "italic"),
#     panel.grid.major.x  = element_blank(),
#     panel.grid.minor    = element_blank(),
#     plot.margin         = margin(5, 5, 5, 5)
#   )
# 
# p_lat_segments_field
# 
# 
# # 
# # ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_genus_gradient_bar.jpeg"),
# #        p_lat_segments_field, width = 12, height = 6, dpi = 600)
# # 
# # ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_genus_gradient_bar.svg"),
# #        p_lat_segments_field, width = 12, height = 6, dpi = 600)
# # 
# # 
# 
# 
# # . ------------------
# 3) All iNat --------------
##  3.1) Richness/observations vs latitude ------------------------
inat_sf <- inat_obs %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Get countries and continents
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3, continent, name)

inat_pts <- inat_sf %>%
  filter(!is.na(longitude), !is.na(latitude), abs(longitude) <= 180, abs(latitude) <= 90) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Project points to meters, buffer, project back
# because points on small islands are otherwise part of the "continents"
pts_buf <- inat_pts %>%
  st_transform(3857) %>%      # meters
  st_buffer(dist = 100000) %>%  # m
  st_transform(4326)          # back to lon/lat

inat_joined_buf <- pts_buf %>%
  st_join(world, join = st_intersects, left = TRUE) %>%
  st_drop_geometry() %>%
  mutate(
    hemisphere = if_else(latitude >= 0, "North", "South"),
    continent = if_else(name=="Saint Helena", "Africa",continent),
    tax_world = case_when(
      continent %in% c("North America", "South America") ~ "New World",
      continent %in% c("Europe", "Africa", "Asia", "Oceania") ~ "Old World",
      TRUE ~ NA_character_
    )
  )

# inat_joined_buf %>%
#   filter(is.na(tax_world)) %>% View()


p_lat_idx <- inat_joined_buf %>%
  drop_na(tax_world) %>%
  filter(is.finite(latitude), is.finite(longitude),
         !is.na(taxon_id), !is.na(id)) %>%
  unique() %>%
  mutate( lat_mid5 = floor(latitude / 5) * 5 + 0.5) %>%
  group_by(lat_mid5, hemisphere, tax_world) %>%
  reframe(
    n_tax = n_distinct(taxon_id),
    n_obs = n_distinct(id),
    index = if_else(n_obs > 0, n_tax / n_obs, NA_real_)
  ) %>%
  unique() %>%
  ggplot(aes(x = index, y = lat_mid5)) +
  geom_point(aes(size = n_obs, alpha = n_tax), color = "forestgreen") +
  geom_smooth(method = "lm", orientation = "y", color = "black", se = TRUE) +
  facet_grid(rows = vars(hemisphere), cols = vars(tax_world), scales = "free_y") +
  scale_alpha_continuous("No. of species",range=c(0.3,1))+
  scale_size_continuous("No. of observations ", range = c(3,10),
                        breaks=c(1,10,100,1000))+
  labs(
    x = "Observed species per observation (n_tax / n_obs)",
    y = "Latitude (5° bin midpoint)")+

  stat_cor(
    method = "spearman",
    label.x.npc = "right",
    label.y.npc = "top",
    label.sep = ";",
    r.accuracy = 0.01,
    p.accuracy = 0.001,
    hjust=0.5,

    geom = "label",
    fill = "white",
    alpha = 0.7,
    label.padding = unit(2, "pt"),
    size = 3.5

  )+

  my_theme14+
  guides(
    alpha = guide_legend(override.aes = list(size = 5)))

p_lat_idx

ggsave(file.path(OUT_DIR_LATITUDE, "ALL_iNat_species_per_obs_vs_latitude.jpg"),
       p_lat_idx, width = 9, height = 8, dpi = 300)






# . -----------------------
# . ----------------------


# # Simple Bayesian analyses ----------------
# 
# species_by_tree <- inat_field %>%
#   distinct(site, dbh_cat, setting, lat, lon, taxon.name, tree_id)
# 
# # Count distinct trees with each species
# counts_species <- species_by_tree %>%
#   count(site, dbh_cat, setting, lat, lon, taxon.name, name = "count_trees_with_species")
# 
# # Join with total_trees per stratum
# df <- counts_species %>%
#     left_join(
#     inat_field %>%
#       select(site, dbh_cat, setting, lat, lon, total_trees) %>%
#       unique(),
#     by = c("site","dbh_cat","setting","lat","lon")
#   ) %>%
#   filter(total_trees > 0) %>%
#   mutate(
#     dbh_cat = factor(dbh_cat, levels = c("10 - 30","30 - 60","60 - 90","90+"), ordered = T),
#     taxon.name = as.factor(taxon.name),
#     setting = relevel(factor(setting), ref = "Rural"),
#     lat_cdeg = lat - field_lat_mean
# 
#   )# %>% #once off because of error in raw data
#   #mutate(total_trees=if_else(site=="S31" & count_trees_with_species == 3 & total_trees == 2 ,3,total_trees))
# 
# fit_additive <- brm(
#   count_trees_with_species | trials(total_trees) ~ lat + (1|setting),
#   data = df,
#   family = binomial(),
#   prior = c(
#     prior(normal(0, 1.5), class = "b"),
#     prior(normal(0, 1.5), class = "Intercept"),
#     prior(exponential(1), class = "sd")
#   ),
#   chains = 4, cores = 4, iter = 4000, seed = 2025
# )
# 
# summary(fit_additive)
# 
# 
# ##################
# fit_comp <- brm(
#   formula = count_trees_with_species | trials(total_trees) ~
#               lat_cdeg * lat + (1|setting),                    # site-level heterogeneity
#   data = df,
#   family = binomial(),
#   prior = c(
#     prior(normal(0, 1.5), class = "b"),
#     prior(normal(0, 1.5), class = "Intercept"),
#     prior(exponential(1), class = "sd")),
#   chains = 4, cores = 4, iter = 4000, seed = 2025
# )
# 
# summary(fit_comp)
# 
# #Posterior probability of direction
# hypothesis(fit_comp, "lat_cdeg < 0")
# hypothesis(fit_comp, "settingUrban > 0")
# hypothesis(fit_comp, "lat_cdeg:settingUrban > 0")
# 
# 
# 
