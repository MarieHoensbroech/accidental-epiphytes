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
})

# Parameters -----------------------------------
IN_MY_SITES   <- "data/processed/my.sites.csv"
IN_INAT_FIELD <- "data/processed/inat.merged.csv"
IN_INAT_OBS   <- "data/processed/inat_observations.csv"
IN_ENV        <- "data/processed/10_occ_env.csv"
  
#Temperature 
OUT_DIR_HEIGHT  <- "outputs/figures/07_height_gradient"
dir.create(OUT_DIR_HEIGHT,  showWarnings = FALSE, recursive = TRUE)


#Load data ---------
my_sites   <- fread(IN_MY_SITES) %>% as_tibble() %>% #glimpse()
  reframe(site=Site, dbh_cat=TreeCat,type=Type,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=SiteLat,
          lon=SiteLon)


occ_covs <- fread(IN_ENV) %>% 
  filter(str_detect(site,"N")) %>% 
  select(lat,lon,CHELSA_bio1, CHELSA_bio6) %>% 
  unique()


inat_field <- fread(IN_INAT_FIELD) %>% as_tibble() %>% #glimpse()
  reframe(site=Site, dbh_cat=TreeCat,type=Type,
          obs_count=n_epis_per_site,
          spp_count=n_sp_per_site,
          obs_count_dbh=n_epis_per_site_cat,
          spp_count_dbh=n_sp_per_site_cat,
          total_trees=NumberTrees,
          lat=SiteLat,
          lon=SiteLon,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          epi_height=EpiHeight_num,
          genus_name) %>% 
  unique() %>% 
  left_join(occ_covs)


all_inat   <- fread(IN_INAT_OBS) %>% as_tibble() %>% #glimpse()
  reframe(lat=latitude,
          lon=longitude,
          id, taxon_id,taxon.name,taxon.rank,
          quality_grade,
          setting=Setting,
          count=EpiCount_noNA,
          taxon.rank,
          epi_height = EpiHeight_num,
          genus_name, user_login) %>% 
  unique() %>% 
  left_join(occ_covs) %>% 
  mutate(epi_height=if_else(epi_height>=15,epi_height/100,epi_height))

my_inat <- all_inat %>% filter(user_login=="marie-ho") 


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




# 1) Field data: height correlation --------------
## Boxplot: species~height ---------
inat_field %>% 
  filter(taxon.rank%in%c("genus","complex","species","tribe")) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name)) %>% 
  group_by(taxon.name) %>% 
  filter(n_distinct(id)>3) %>% 
  ggplot()+
  geom_boxplot(aes(x=fct_reorder(taxon.name, epi_height, .fun = mean, .na_rm = T) ,y=epi_height),
               linewidth=1)+
  facet_wrap(~type, scales="free_x")+
  labs(x="",y="Growing height (m)")+
  my_theme14+
  theme(axis.text.x = element_text(angle=45,  hjust = 1, face = "italic"))

ggsave(file.path(OUT_DIR_HEIGHT, "FIELD_height_gradient_box.jpeg"),
       width = 8, height = 6, dpi = 600)




## Height ~ latitude correlation -------

inat_field %>% 
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
  slice_head(n = 12) %>% 
  select(taxon.name) %>% deframe() -> top9_field


inat_field %>% 
  filter(taxon.name%in%top9_field) %>% 
  select(genus_name,id,epi_height,lat,lon,taxon.name,type)->field_lat_height_df


### Most common species ---------
field_lat_height_df %>% 
  ggplot(aes(y=lat, x=epi_height))+
  facet_wrap(~taxon.name) +
  geom_point()+
  geom_smooth(method="lm",colour="grey")+
  my_theme14


### Polypodium --------
field_lat_height_df %>% 
  filter(str_detect(taxon.name,"Polypodium")) %>% 
  #filter(taxon.name!="Polypodium") %>% 
  ggplot(aes(y=lat, x=epi_height,col=taxon.name))+
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
 


## Barplot: Species-Height-distribution ----------
df_height_field <- field_lat_height_df %>%
  mutate(epi_height_mid = floor(epi_height / 1) * 1 + 0.5) %>% 
  group_by(epi_height_mid, type, taxon.name) %>%
  mutate(sum=1) %>%  #every row is one individual
  summarise(sum = sum(sum, na.rm = TRUE), .groups = "drop") %>%
  group_by(epi_height_mid, type) %>%
  mutate(
    total = sum(sum, na.rm = TRUE),
    prop  = if_else(total > 0, sum / total, 0)
  ) %>%
  ungroup() %>%
  # Fix factor levels so colors remain consistent
  #mutate(taxon.name = factor(taxon.name, levels = taxa_levels_field, ordered = TRUE)) %>%
  # Order segments by taxon for reproducibility (or pick another order)
  arrange(type, epi_height_mid, taxon.name) %>%
  group_by(type, epi_height_mid) %>%
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
  facet_grid(. ~ type) +
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



ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_species_gradient_bar.jpeg"),
       p_lat_segments_field, width = 12, height = 8, dpi = 600)

ggsave(file.path(OUT_DIR_LATITUDE, "FIELD_species_gradient_bar.svg"),
       p_lat_segments_field, width = 12, height = 8, dpi = 600)














# 2) My iNat data -------------

## boxplot: species~height ---------------
my_inat %>% 
  #select(taxon.rank) %>% unique()
  filter(taxon.rank%in%c("genus","complex","species","tribe","section"),
         setting != "Optional",
         genus_name!="") %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Polypodium")|
                              str_detect(taxon.name,"Polypodiales"),
                            "Polypodium sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Acer"),
                            "Acer sp.",taxon.name)) %>% 
  mutate(taxon.name=if_else(str_detect(taxon.name,"Quercus"),
                            "Quercus sp.",taxon.name)) %>% 
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
  scale_color_viridis(discrete=T, begin=0.3,end=0.7)+
  labs(x="",y="Growing height (m)")+
  my_theme14+
  theme(axis.text.x = element_text(angle=45,  hjust = 1, face = "italic"))

ggsave(file.path(OUT_DIR_HEIGHT, "ALL_height_gradient_box.jpeg"),
       width = 13, height = 6, dpi = 600)

## Most common species ----------
my_inat %>% 
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
  slice_head(n = 12) %>% 
  select(taxon.name) %>% deframe() -> top9_field


my_inat %>% 
  filter(taxon.name%in%top9_field) %>% 
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



inat_lat_height_df %>% 
  ggplot(aes(x=CHELSA_bio1, y=epi_height,col=genus_name))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_viridis(discrete=T) +
  my_theme14+
  facet_wrap(~genus_name)+
  stat_cor(
    method = "spearman",        
    label.x.npc = "left",        
    label.y.npc = "top",
    label.sep = ";",             
    r.accuracy = 0.01,            
    p.accuracy = 0.001,         
    vjust=1,
    show.legend = F)

## Most common species ----------
all_inat %>% 
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
  select(taxon.name) %>% deframe() -> top9_field


all_inat %>% 
  filter(taxon.name%in%top9_field) %>% 
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


# . ------------
# . -----------


# 4) Temperature -------------

all_inat %>% 
  drop_na(epi_height) %>% 
  filter(genus_name!="") %>% 
  select(id,genus_name,taxon.name,count) %>% 
  unique() %>% 
  group_by(taxon.name) %>% 
  reframe(sum=sum(count)) %>% 
  arrange(desc(sum)) %>% 
  unique() %>% 
  slice_head(n = 9) %>% 
  select(taxon.name) %>% deframe() -> top9_inat


all_inat %>% 
  drop_na(epi_height) %>% 
  filter(taxon.name%in%top9_inat) %>% 
  select(taxon.name,id,epi_height,lat,lon,taxon.name) %>% 
  mutate(epi_height=if_else(epi_height>15, epi_height/100, epi_height))->inat_lat_height_df


all_inat %>%
  filter(str_detect(taxon.name,"Polypodium")) %>% 
  mutate(epi_height=if_else(epi_height>15, epi_height/100, epi_height)) %>% 
  group_by(taxon.name) %>% 
  filter(n_distinct(id)>20) %>% 
  ggplot(aes(x=CHELSA_bio6, y=epi_height, col=taxon.name))+
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



lat_height_df %>% 
  filter(str_detect(taxon.name,"Polypodium")) %>% 
  #filter(taxon.name!="Polypodium") %>% 
  ggplot(aes(y=lat, x=epi_height,col=taxon.name))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_viridis(discrete=T) +
  my_theme14+
  stat_cor(
    method = "spearman",        
    label.x.npc = "right",        
    label.y.npc = "top",
    label.sep = ";",             
    r.accuracy = 0.01,            
    p.accuracy = 0.001,         
    hjust=1,
    show.legend = F)
