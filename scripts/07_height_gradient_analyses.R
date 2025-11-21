# 07_height_gradient_analyses.R
# 
# !!!!!!!!WORK IN PROGRESS!!!!!!!!
#
# Height distribution patterns of richness/observation ratios
# 1) For field data only
# 2) For my iNat data collected in 2025
# 3) For all accidental epiphyte observations
# 4) Does height change with temperature gradients?
# I use CHELSA bio1 (mean annual) and bio6 (lowest T in coldest quarter)
rm(list=ls())

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
  library(maps)
  library(patchwork)
  library(stringr)
  library(data.table)
  library(cowplot)
  library(grid)
  library(ggimage)
  library(rnaturalearth)
})

# Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"
IN_ENV        <- "outputs/04_environmental/env_vars_sel.csv"
IN_ENV_LOOKUP <- "outputs/04_environmental/var_lookup.csv"

#Temperature 
OUT_DIR_HEIGHT  <- "outputs/07_height_gradient"
dir.create(OUT_DIR_HEIGHT,  showWarnings = FALSE, recursive = TRUE)


#Load data ---------
my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% #glimpse()
  reframe(site, dbh_cat=TreeCat,setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon) %>% 
  filter(setting%in%c("Urban","Rural"))


occ_covs <- fread(IN_ENV) %>% distinct() %>% as_tibble() %>% 
  mutate(id=as.character(id))

var_lookup <- fread(IN_ENV_LOOKUP) %>% as_tibble()

inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% 
  reframe(site, dbh_cat=TreeCat,setting,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=siteLat,
          lon=siteLon,
          id=as.character(id), taxon_id,taxon.name,taxon.rank,
          quality_grade,
          count=EpiCount_noNA,
          epi_height=EpiHeight_num,
          genus_name) %>% 
  unique() %>% 
  filter(setting%in%c("Urban","Rural")) %>% 
  left_join(occ_covs)

inat_field<-inat_field %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name)) 

#Once off (maybe) because of missing obs field on iNat -fixed 05.11.
# inat_field %>%
#   mutate(setting=if_else(setting=="","Rural",setting),
#          dbh_cat=if_else(dbh_cat=="","90+",dbh_cat))->inat_field


all_inat   <- fread(IN_INAT_OBS) %>% as_tibble() %>% #glimpse()
  reframe(lat=latitude,
          lon=longitude,
          id=as.character(id), 
          taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          taxon.rank,
          epi_height = EpiHeight_num,
          genus_name,
          user_login, obs_date) %>% 
  unique() %>% 
  left_join(occ_covs) %>% 
  mutate(epi_height=if_else(epi_height>=15,epi_height/100,epi_height)) %>% 
  filter(setting%in%c("Urban","Rural"))

my_inat <- all_inat %>% filter(user_login=="marie-ho") %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name))  %>% 
  filter(setting%in%c("Urban","Rural"))


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

my_theme14_inside <- 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top",
    axis.ticks.length = unit(-0.2, "cm"),
    axis.text.y = element_text(margin = margin(l = 20, r = -32)),
    axis.text.x = element_text(margin = margin(t = -18, b = 30))
  )


# 1) Field data: height correlation --------------
## Boxplot: species~height ---------
p_box_height<-
  inat_field %>% 
  filter(taxon.rank%in%c("genus","complex","species","tribe","hybrid")) %>% 
  group_by(taxon.name) %>% 
  filter(n_distinct(id)>3) %>% 
  ggplot()+
  geom_boxplot(aes(x=fct_reorder(taxon.name, epi_height, .fun = median, .na_rm = T),
                   y=epi_height,
                   col=setting
  ),
  outlier.alpha = 0.7,
  linewidth=1, position = position_dodge2(preserve = "single",padding=0.2))+
  labs(x="",y="Growing height (m)")+
  scale_color_viridis("Setting",begin=0.4,end=0.8,discrete = T)+
  my_theme14+
  theme(axis.text.x = element_text(angle=45,  hjust = 1, face = "italic"))

ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_height_gradient_box.jpeg"),
       p_box_height,
       width = 7, height = 5, dpi = 600)
ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_height_gradient_box.svg"),
       p_box_height,
       width = 7, height = 5, dpi = 600)




## Height ~ latitude correlation -------

inat_field %>% 
  drop_na(epi_height) %>% 
  filter(taxon.rank%in%c("genus","species","complex","hybrid","tribe")) %>% 
  #filter(type=="Forest") %>% 
  #filter(genus_name!="") %>% 
  select(id,taxon.name,count) %>% 
  unique() %>% 
  group_by(taxon.name) %>% 
  reframe(sum=sum(count)) %>% 
  arrange(desc(sum)) %>% 
  unique() %>% 
  slice_head(n = 12) %>% 
  select(taxon.name) %>% deframe() -> top9_field


inat_field %>% 
  filter(taxon.name%in%top9_field) %>% 
  select(genus_name,id,epi_height,lat,lon,taxon.name,setting)->field_lat_height_df


### Plot: Most common species ---------
p_top9_lat<-
  field_lat_height_df %>% 
  ggplot(aes(x=lat, y=epi_height,col=setting))+
  scale_color_viridis("Setting",begin=0.4,end=0.8,discrete = T)+
  facet_wrap(~taxon.name) +
  geom_smooth(method="lm")+
  geom_point()+
  labs(x="Latitude (°)",y="Growing height (m)")+
  my_theme14 +
  theme(strip.text = element_text(face = "italic"))

ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_top9_lat_gradient_box.jpeg"),
       p_top9_lat,width=10,height=7,dpi=500)


### Polypodium --------
field_lat_height_df %>% 
  filter(str_detect(taxon.name,"Polypodium")) %>% 
  #filter(taxon.name!="Polypodium") %>% 
  ggplot(aes(x=lat, y=epi_height,col=taxon.name))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_viridis(discrete=T) +
  facet_wrap(~taxon.name,ncol=2)+
  my_theme14+
  stat_cor(
    method = "spearman",        
    label.x.npc = "right",        
    label.y.npc = "top",
    label.sep = ";",             
    r.accuracy = 0.01,            
    p.accuracy = 0.001,         
    hjust=1,
    show.legend = F) -> poly_lat_height_p

poly_lat_height_p

## Barplot: Species-Height-distribution ----------
df_height_field <- field_lat_height_df %>%
  mutate(epi_height_mid = floor(epi_height / 1) * 1 + 0.5) %>% 
  group_by(epi_height_mid, setting, taxon.name) %>%
  mutate(sum=1) %>%  #every row is one individual
  summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
  group_by(epi_height_mid, setting) %>%
  mutate(
    total = sum(sum, na.rm = TRUE),
    prop  = if_else(total > 0, sum / total, 0)
  ) %>%
  ungroup() %>%
  # Fix factor levels so colors remain consistent
  #mutate(taxon.name = factor(taxon.name, levels = taxa_levels_field, ordered = TRUE)) %>%
  # Order segments by taxon for reproducibility (or pick another order)
  arrange(setting, epi_height_mid, taxon.name) %>%
  group_by(setting, epi_height_mid) %>%
  # Build cumulative horizontal allocation [0..1]
  mutate(
    xmin = cumsum(dplyr::lag(prop, default = 0)),
    xmax = cumsum(prop)
  ) %>%
  ungroup() %>% 
  unique()


# Row thickness in height meters
height_bin_height <- 1

p_height_segments_field <- 
  ggplot(df_height_field) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = epi_height_mid - height_bin_height/2,
      ymax = epi_height_mid + height_bin_height/2,
      fill = taxon.name
    ),
    color = NA
  ) +
  facet_grid(~ setting) +
  scale_x_continuous(
    name   = "Relative dominance",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous("Height (m)") +  
  # scale_fill_manual("Taxon name",
  #                   values = taxa_pal_field,
  #                   limits = taxa_levels_field,
  #                   drop   = FALSE,
  #                   guide  = guide_legend(ncol = 2) 
  # )+
  scale_fill_viridis("Taxon name", begin=1,end=0,discrete=T)+
  my_theme14+
  theme(
    legend.position     = "right",
    legend.text         = element_text(face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(5, 5, 5, 5)
  )

p_height_segments_field



ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_species_gradient_bar.jpeg"),
       p_height_segments_field, width = 12, height = 8, dpi = 600)

ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_species_gradient_bar.svg"),
       p_height_segments_field, width = 12, height = 8, dpi = 600)




# ___________ -----------------



# 2) My iNat data -------------

## boxplot: species~height ---------------
my_inat %>% 
  filter(obs_date>=lubridate::dmy("01/04/2025")) %>% 
  filter(taxon.rank%in%c("genus","complex","species","tribe","section","hybrid"),
         setting != "Optional",
         genus_name!="") %>% 
  group_by(genus_name,setting) %>% 
  filter(n_distinct(id)>3) %>% 
  drop_na(epi_height) %>% 
  ggplot()+
  geom_boxplot(aes(x=fct_reorder(genus_name, epi_height, .fun = max, .na_rm = T) ,
                   y=epi_height,
                   col=setting),
               linewidth=1, width=1,
               position=position_dodge2(preserve = "single", padding=0.2)
  )+
  # facet_wrap(~setting, scales="free_x",
  #            space="free_x")+
  scale_color_viridis("Setting",discrete=T, begin=0.3,end=0.7)+
  labs(x="",y="Growing height (m)")+
  my_theme14+
  theme(axis.text.x = element_text(angle=45,  hjust = 1, face = "italic"))

ggsave(file.path(OUT_DIR_HEIGHT, "MYINAT_height_gradient_box.jpeg"),
       width = 13, height = 6, dpi = 600)

## Most common species ----------
my_inat %>% 
  drop_na(epi_height) %>% 
  filter(taxon.rank%in%c("genus","species","complex","hybrid","tribe")) %>% 
  #filter(type=="Forest") %>% 
  #filter(genus_name!="") %>% 
  select(id,taxon.name,count) %>% 
  unique() %>% 
  group_by(taxon.name) %>% 
  reframe(sum=sum(count)) %>% 
  arrange(desc(sum)) %>% 
  unique() %>% 
  slice_head(n = 12) %>% 
  select(taxon.name) %>% deframe() -> top9_myinat


my_inat %>% 
  filter(taxon.name%in%top9_myinat) %>% 
  select(genus_name,id,epi_height,lat,lon,taxon.name,setting)->myinat_lat_height_df



## Barplot: Species-Height-distribution ----------
myinat_height_field <- myinat_lat_height_df %>%
  filter(!setting%in%c("","Optional")) %>% 
  mutate(epi_height_mid = floor(epi_height / 1) * 1 + 0.5) %>% 
  group_by(epi_height_mid, setting, taxon.name) %>%
  mutate(sum=1) %>%  #every row is one individual
  summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
  group_by(epi_height_mid, setting) %>%
  mutate(
    total = sum(sum, na.rm = TRUE),
    prop  = if_else(total > 0, sum / total, 0)
  ) %>%
  ungroup() %>%
  # Fix factor levels so colors remain consistent
  #mutate(taxon.name = factor(taxon.name, levels = taxa_levels_field, ordered = TRUE)) %>%
  # Order segments by taxon for reproducibility (or pick another order)
  arrange(setting, epi_height_mid, taxon.name) %>%
  group_by(setting, epi_height_mid) %>%
  # Build cumulative horizontal allocation [0..1]
  mutate(
    xmin = cumsum(dplyr::lag(prop, default = 0)),
    xmax = cumsum(prop)
  ) %>%
  ungroup() %>% 
  unique()


# Row thickness in height meters
height_bin_height <- 1

p_height_segments_myinat <- 
  ggplot(myinat_height_field) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = epi_height_mid - height_bin_height/2,
      ymax = epi_height_mid + height_bin_height/2,
      fill = taxon.name
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
  scale_y_continuous("Height (m)") +  
  # scale_fill_manual("Taxon name",
  #                   values = taxa_pal_field,
  #                   limits = taxa_levels_field,
  #                   drop   = FALSE,
  #                   guide  = guide_legend(ncol = 2) 
  # )+
  scale_fill_viridis("Taxon name", begin=1,end=0,discrete=T)+
  my_theme14+
  theme(
    legend.position     = "right",
    legend.text         = element_text(face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(5, 5, 5, 5)
  )

p_height_segments_myinat
ggsave(file.path(OUT_DIR_HEIGHT, "MYINAT_species_gradient_bar.jpeg"),
       width = 13, height = 6, dpi = 600)

# 3) All iNat -------------

all_inat %>% 
  drop_na(epi_height) %>% 
  filter(genus_name!="") %>% 
  select(id,genus_name,taxon.name,count) %>% 
  unique() %>% 
  group_by(genus_name) %>% 
  reframe(sum=sum(count)) %>% 
  arrange(desc(sum)) %>% 
  unique() %>% 
  slice_head(n = 12) %>% 
  select(genus_name) %>% deframe() -> top9_inat


all_inat %>% 
  drop_na(epi_height) %>% 
  filter(genus_name%in%top9_inat) %>% 
  #select(taxon.name,id,epi_height,lat,lon,taxon.name) %>% 
  mutate(epi_height=if_else(epi_height>15, epi_height/100, epi_height))->inat_lat_height_df


all_inat %>%
  filter(str_detect(taxon.name,"Polypodium")) %>% 
  mutate(epi_height=if_else(epi_height>15, epi_height/100, epi_height)) %>% 
  group_by(taxon.name) %>% 
  filter(n_distinct(id)>20) %>% 
  ggplot(aes(y=lat, x=epi_height, col=taxon.name))+
  facet_wrap(~taxon.name,
             #scales="free_y"
  ) +
  geom_point()+
  geom_smooth(method="lm", col="grey")+
  my_theme14+
  scale_color_viridis(discrete=T)+
  stat_cor(
    method = "spearman",        
    label.x.npc = "right",        
    label.y.npc = "top",
    label.sep = ";",             
    r.accuracy = 0.01,            
    p.accuracy = 0.001,         
    hjust=1,
    show.legend = F)



## HEIGHT VS ENVIRONMENT-----------
# 1) Long format + labels
inat_env_height_long <- 
  inat_lat_height_df %>% 
  tidyr::pivot_longer(CHELSA_bio7:forest_age,
                      names_to = "var", values_to = "value") %>% 
  left_join(var_lookup)

# 2) Get the plotting order of variables based on 'label' (optional)
vars_tbl <- inat_env_height_long %>%
  dplyr::distinct(var, label) %>%
  dplyr::arrange(label)

# 3) Open a multi-page PDF and print one page per environmental variable
pdf(file.path(OUT_DIR_HEIGHT, "INAT_TOP9_growing_height~environmental.pdf"),
    width = 10, height = 8)

for (i in seq_len(nrow(vars_tbl))) {
  v   <- vars_tbl$var[i]
  lab <- vars_tbl$label[i]
  
  dat <- inat_env_height_long %>% dplyr::filter(var == v)
  
  p <- ggplot(dat, aes(x = value, y = epi_height, col = genus_name)) +
    geom_point(alpha=0.4,size=2) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_viridis("", discrete = TRUE, guide = "none", begin = 1, end = 0) +
    my_theme14 +
    labs(
      x = lab,
      y = "Growing height (m)"
    ) +
    facet_wrap(~ genus_name) +
    # Spearman’s rho displayed inside each facet
    ggpubr::stat_cor(
      aes(x = value, y = epi_height),
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ";",
      r.accuracy = 0.01,
      p.accuracy = 0.001,
      vjust = 1,
      show.legend = FALSE
    ) +
    # Make facet titles italic (genus names)
    theme(strip.text = element_text(face = "italic"))
  
  print(p)
}

dev.off()



## Most common species ----------
all_inat %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name)) %>% 
  drop_na(epi_height) %>% 
  filter(taxon.rank%in%c("genus","species","complex","hybrid")) %>% 
  #filter(type=="Forest") %>% 
  #filter(genus_name!="") %>% 
  select(id,taxon.name,count) %>% 
  unique() %>% 
  group_by(taxon.name) %>% 
  reframe(sum=sum(count)) %>% 
  arrange(desc(sum)) %>% 
  unique() %>% 
  slice_head(n = 15) %>% 
  select(taxon.name) %>% deframe() -> top9_allinat


all_inat %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name)) %>% 
  filter(taxon.name%in%top9_allinat) %>% 
  select(genus_name,id,epi_height,lat,lon,taxon.name,setting)->allinat_lat_height_df



## Barplot: Species-Height-distribution ----------
allinat_segments_height_field <- allinat_lat_height_df %>%
  filter(!setting%in%c("","Optional")) %>% 
  mutate(epi_height_mid = floor(epi_height / 1) * 1 + 0.5) %>% 
  group_by(epi_height_mid, setting, taxon.name) %>%
  mutate(sum=1) %>%  #every row is one individual
  summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
  group_by(epi_height_mid, setting) %>%
  mutate(
    total = sum(sum, na.rm = TRUE),
    prop  = if_else(total > 0, sum / total, 0)
  ) %>%
  ungroup() %>%
  # Fix factor levels so colors remain consistent
  #mutate(taxon.name = factor(taxon.name, levels = taxa_levels_field, ordered = TRUE)) %>%
  # Order segments by taxon for reproducibility (or pick another order)
  arrange(setting, epi_height_mid, taxon.name) %>%
  group_by(setting, epi_height_mid) %>%
  # Build cumulative horizontal allocation [0..1]
  mutate(
    xmin = cumsum(dplyr::lag(prop, default = 0)),
    xmax = cumsum(prop)
  ) %>%
  ungroup() %>% 
  unique()


# Row thickness in height meters
height_bin_height <- 1


p_height_segments_allinat <- 
  ggplot(allinat_segments_height_field) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = epi_height_mid - height_bin_height/2,
      ymax = epi_height_mid + height_bin_height/2,
      fill = taxon.name
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
  scale_y_continuous("Height (m)") +  
  # scale_fill_manual("Taxon name",
  #                   values = taxa_pal_field,
  #                   limits = taxa_levels_field,
  #                   drop   = FALSE,
  #                   guide  = guide_legend(ncol = 2) 
  # )+
  scale_fill_viridis("Taxon name", begin=1,end=0,discrete=T)+
  my_theme14+
  theme(
    legend.position     = "right",
    legend.text         = element_text(face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(5, 5, 5, 5)
  )

p_height_segments_allinat


ggsave(file.path(OUT_DIR_HEIGHT, "ALLINAT_species_gradient_bar.jpeg"),
       p_height_segments_allinat,
       width = 13, height = 6, dpi = 600)

