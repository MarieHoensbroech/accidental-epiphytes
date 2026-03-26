# 07_APPENDIX_epi_host_growing_site_associations.R
#
# Assess which trees accidental epiphytes occur in
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

OUT_DIR_HOST <- "outputs/07_host_growing_sites"
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
          dbh_cat = tree_cat, taxon_rank,
          moss=if_else(moss=="yes","mossy","")
  ) %>% distinct()

inat_obs   <- fread(IN_INAT_OBS) %>% as_tibble() %>% 
  drop_na(latitude,longitude) %>%  clean_names() %>% 
  reframe(id=as.character(id),taxon_name,
          tree_sp,growing_site,
          tree_genus=clean_to_genus(tree_sp),
          epi_count = as.integer(str_extract(epi_count, "\\d+")),
          taxon_rank,  user_login,obs_date,
          moss=if_else(moss=="yes","mossy","")
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

# Plots ------
# 1) Chord diagram host~epi----------
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
inat_obs %>% 
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



# 5) Chord diagram epi~growing site----------
## 5b) Field data only ----------

inat_field %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    epi_genus = clean_to_genus(taxon_name),
    growing_site=if_else(growing_site%in%c(""," ","  "),"unknown",growing_site),
    growing_site=paste(moss, growing_site,sep=" "),
    growing_site=str_trim(growing_site),
    
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(growing_site,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(epi_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(epi_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(epi_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(epi_genus) %>% 
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
  png(file.path(OUT_DIR_HOST, "FIELD_epi_growing_site_chord.png"), width = 6000, height = 5300, res=600)
  
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




## 5b) My iNat data only ----------

inat_obs %>%
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section"),
         user_login=="marie-ho",
         obs_date>=as.Date("2025-04-01")) %>% 
  reframe(
    epi_genus = clean_to_genus(taxon_name),
    growing_site=if_else(growing_site%in%c(""," ","  "),"unknown",growing_site),
    growing_site=paste(moss, growing_site,sep=" "),
    growing_site=str_trim(growing_site),
    
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(growing_site,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(epi_genus, growing_site) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

epi_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(epi_genus, levels = host_order_all)),
    to   = paste0("E: ", factor(growing_site,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>1) # keep only links with more than 1 occurrence


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(epi_genus = sub("^H: ", "", from)) %>%
  #arrange(match(host_genus, host_order_all)) %>%
  arrange(epi_genus) %>% 
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
  png(file.path(OUT_DIR_HOST, "MYiNAT_epi_growing_site_chord.png"), width = 6000, height = 5300, res=600)
  
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
inat_obs %>% 
  filter(taxon_rank%in%c("genus","species","variety","subspecies","complex","hybrid","form",
                         "subgenus","section")) %>% 
  reframe(
    epi_genus = clean_to_genus(taxon_name),
    growing_site=if_else(growing_site%in%c(""," ","  "),"unknown",growing_site),
    growing_site=paste(moss, growing_site,sep=" "),
    growing_site=str_trim(growing_site),
    
    individuals = as.integer(str_extract(epi_count, "\\d+")),
    individuals = if_else(is.na(individuals)|individuals=="",1,individuals)
  ) %>%
  filter(!is.na(individuals)) %>%
  group_by(growing_site,epi_genus) %>%
  summarise(individuals = sum(individuals, na.rm = TRUE), .groups = "drop") ->chord_df



# Aggregate individuals per host × epiphyte pair
links <- chord_df %>%
  group_by(growing_site,epi_genus) %>%
  summarise(n_indiv = sum(individuals, na.rm = TRUE), .groups = "drop") %>%
  filter(n_indiv > 0)

# Order host/epi by total abundance to bring dominant taxa forward
host_order_all <- links %>%
  group_by(growing_site) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(growing_site)

epi_order_all <- links %>%
  group_by(epi_genus) %>%
  summarise(total = sum(n_indiv), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(epi_genus)

# Make the two partitions explicit so it's strictly bipartite
links_plot <- links %>%
  mutate(
    from = paste0("H: ", factor(growing_site, levels = host_order_all)),
    to   = paste0("E: ", factor(epi_genus,  levels = epi_order_all)),
    value = n_indiv
  ) %>%
  select(from, to, value) %>% 
  filter(value>5) # keep only links with more than 2 occurrences


# Recompute the sector orders only for remaining sectors,
#    but keep the relative order from the global ordering
hosts_present <- links_plot %>%
  distinct(from) %>%
  mutate(growing_site = sub("^H:", "", from)) %>%
  arrange(growing_site) %>% 
  #arrange(match(host_genus, host_order_all)) %>%
  pull(from)

epis_present <- links_plot %>%
  distinct(to) %>%
  mutate(epi_genus = sub("^E: ", "", to)) %>%
  arrange(epi_genus) %>% 
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

png(file.path(OUT_DIR_HOST, "iNat_EUROPE_epi_growing_site_chord.png"), width = 6000, height = 5300, res=600)

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




