# 08_pca.R 
#
# 1-2) clean, summarise, and join field+inat+bg data with select environment variables
# 3-5) prep, run, plot PCA 
# 6) post hoc permutation tests


suppressPackageStartupMessages({
  library(tidyverse)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(scales)
  library(vegan)
  library(sf)
  library(rnaturalearth)
  library(units)
  library(rnaturalearthdata)
  library(tibble)
  library(patchwork) 
  library(data.table)
  library(forcats)
  library(vegan)
  library(ggnewscale) 
  library(viridisLite)
  
})

set.seed(2025)


# 0) LOAD DATA ---------------
my_sites<-fread("data/processed/my.sites.csv") %>% as_tibble()

env_vars_sel<-fread("outputs/04_environmental/env_vars_sel.csv") %>% as_tibble() %>% 
  mutate(taxon.name=if_else(taxon.name%in%c("",NA), "no_taxon",taxon.name)) %>% 
  mutate(LAI=if_else(is.na(LAI),0,LAI),
         forest_age=if_else(is.na(forest_age),0,forest_age))

var_lookup <- fread("outputs/04_environmental/var_lookup.csv")



# 2) Clean names, scale variables----------
sapply(env_vars_sel, function(x) sum(is.na(x))) #check how many NAs per column
# Drop rows with NA only for PCA vars, not grouping cols
dat_pca <- 
  env_vars_sel %>%
  select(where(~ sum(is.na(.x)) < 1000)) %>% 
  drop_na() %>% 
  mutate(source=case_when(str_detect(site,"S")~"Field",
                          str_detect(site,"B")~"iNat",
                          .default="iNat")
  )

sapply(dat_pca, function(x) sum(is.na(x))) #check how many NAs per column
# Standardize PCA vars (constant-safe)
scale_const_safe <- function(v){
  v <- as.numeric(v)
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) v - mean(v, na.rm = TRUE) else as.numeric(scale(v))
}


pca_vars <- c(dat_pca %>% select(CHELSA_bio7:CHELSA_bio15) %>% 
                names(), "lat_cdeg")


# __________------- 
#3) Run PCA----------------
pca_df <- dat_pca %>%
  mutate(across(any_of(setdiff(pca_vars, "is_island")), scale_const_safe),
         is_island = is_island)

X <- pca_df %>% select(any_of(pca_vars)) %>% as.data.frame()

pca <- prcomp(X, center = TRUE, scale. = TRUE)

# 4) Plot prep ---------------
# Scores & explained variance
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x[, 1:2]) %>%
  setNames(c("PC1","PC2")) %>%
  bind_cols(
    pca_df %>% select(source, presence, is_island, taxon.name, quality_grade)
  )

# Env variable set + pretty labels
meta_cols <- c("lat","lon","is_island","source","presence","taxon.name")
env_vars_clean <- pca_df %>%
  select(where(is.numeric), -any_of(c("PC1","PC2"))) %>%
  select(-any_of(intersect(names(.), meta_cols))) %>%
  names()

# Adjust labels
env_labels <-var_lookup %>% deframe()

# ENVFIT on PC1–PC2; filter R² ≥ 20% and p ≤ 0.05
ord_mat <- as.matrix(scores[, c("PC1","PC2")])
env_dat <- pca_df[, env_vars_clean, drop = FALSE]

ef_env <- vegan::envfit(ord_mat, env_dat, permutations = 999)

env_arrows <- as.data.frame(ef_env$vectors$arrows) %>%
  tibble::rownames_to_column("variable") %>%
  mutate(r2 = ef_env$vectors$r,
         p  = ef_env$vectors$pvals)

r2_thresh <- 0 #0.150
p_thresh  <- 0.05

arrow_df <- env_arrows %>%
  filter(r2 >= r2_thresh, p <= p_thresh) %>%
  reframe(
    variable,
    PC1    = PC1 * 3,   # scale arrows for visibility
    PC2    = PC2 * 3,
    r2, p,
    label  = dplyr::recode(variable, !!!env_labels, .default = (variable))
  )


# Data for the iNat & Field plots (points + 95% CI ellipses)

scores2 <- scores %>%
  mutate(pres_f  = factor(presence, c(0,1), c("absent","present")))

inat_pts <- scores2 %>% filter(source == "iNat")
field_pts <- scores2 %>% filter(source == "Field")


# taxon ellipses (research-grade presences only), Top 10 by frequency

TOP_N_TAXA <- 10

taxon_base <- scores %>%
  filter(presence == 1,
         quality_grade == "research",
         !is.na(taxon.name))

tax_counts <- taxon_base %>% count(taxon.name, sort = TRUE, name = "n")
top_taxa   <- tax_counts %>% slice_head(n = TOP_N_TAXA) %>% pull(taxon.name)

taxon_ell <- taxon_base %>%
  filter(taxon.name %in% top_taxa) %>%
  group_by(taxon.name) %>%
  filter(n() >= 5) %>%   # ensure enough points per ellipse
  ungroup() %>%
  mutate(
    taxon_col = factor(
      taxon.name,
      levels = tax_counts$taxon.name[tax_counts$taxon.name %in% top_taxa]  
    )
  ) %>%
  droplevels()

taxon_pts <- taxon_ell

tax_levels <- levels(taxon_ell$taxon_col)
pal_taxon  <- setNames(viridisLite::viridis(length(tax_levels), option = "D"), tax_levels)

x_lab <- paste0("PC1 (", percent(var_exp[1]), ")")
y_lab <- paste0("PC2 (", percent(var_exp[2]), ")")


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

#5) PLOT-------------
{
  inat_pts<-inat_pts %>% 
    mutate(is_island=if_else(is_island==1,"island","mainland"))
  field_pts<-field_pts %>% 
    mutate(is_island=if_else(is_island==1,"island","mainland"))
  
  ## 5A) iNat: points (alpha by presence), ellipses (95% CI by presence × is_island)--------
  p_inat <- ggplot(inat_pts, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey90") +
    geom_vline(xintercept = 0, color = "grey90") +
    geom_point(aes(color = pres_f, shape = is_island, alpha = pres_f), size = 1.6) +
    stat_ellipse(aes(color = pres_f, group = interaction(pres_f, is_island),
                     linetype=is_island),
                 type = "norm", level = 0.95, linewidth = 1.15) +
    scale_color_viridis_d(name = "Presence", begin = 0.5, end = 1) +
    scale_linetype("Setting")+
    scale_alpha_manual(values = c("absent" = 0.10, "present" = 0.70)) +
    guides(alpha = "none", color = "none", shape = "none", linetype = "none") +  # hide legends for this subplot
    labs(title = "iNat data", x = x_lab, y = y_lab) +
    my_theme14
  
  ## 5B) Field: points (fixed alpha), ellipses------------
  p_field <- ggplot(field_pts, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey90") +
    geom_vline(xintercept = 0, color = "grey90") +
    geom_point(aes(color = pres_f, shape = is_island), size = 1.6, alpha = 0.7) +
    stat_ellipse(aes(color = pres_f, group = interaction(pres_f, is_island),
                     linetype=is_island),
                 type = "norm", level = 0.95, linewidth = 1.15) +
    scale_color_viridis_d(name = "Presence", begin = 0.5, end = 1) +
    scale_linetype("Setting")+
    scale_shape("Setting")+
    labs(title = "Field data", x = x_lab, y = y_lab) +
    my_theme14
  
  ## 5C) Env arrows --------
  p_env <- ggplot(arrow_df, aes(0, 0)) +
    geom_hline(yintercept = 0, color = "grey90") +
    geom_vline(xintercept = 0, color = "grey90") +
    geom_segment(aes(xend = PC1, yend = PC2), color = "grey30",
                 arrow = arrow(length = grid::unit(0.18, "cm")),
                 alpha=0.8) +
    ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = label),
                             color = "grey20", size = 3.1) +
    labs(title = "Environmental correlations", x = x_lab, y = y_lab) +
    theme(legend.position = "none")+
    my_theme14
  
  ## 5D) Taxon points----------
  p_taxon <- ggplot(taxon_pts, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey90") +
    geom_vline(xintercept = 0, color = "grey90") +
    # points under ellipses
    geom_point(aes(color = taxon_col), size = 1.4, alpha = 0.35) +
    # ellipses (same colour outline, slightly thicker)
    stat_ellipse(aes(color = taxon_col, group = taxon_col),
                 type = "norm", level = 0.95, linewidth = 1.15) +
    scale_color_manual(values = pal_taxon, name = "Taxon",
                       limits = tax_levels, breaks = tax_levels) +
    # Put the legend INSIDE this panel (top-left; tweak to taste)
    my_theme14 +
    theme(
      legend.position       = c(0.7, 0.98),
      legend.justification  = c(0, 1),
      legend.background     = element_rect(fill = scales::alpha("white", 0.75), color = NA),
      legend.text = element_text(face="italic")
    ) +
    labs(title = "Dominant taxon ellipses", x = x_lab, y = y_lab)
  
  
  ## 5E) Global axis limits ----------
  # Collect all PC1/PC2 used across panels: scores + arrow endpoints
  {
    pc_all <- dplyr::bind_rows(
      dplyr::select(scores, PC1, PC2),
      dplyr::select(arrow_df, PC1, PC2)
    )
    
    # xth and yth percentiles for scores
    qx <- quantile(scores$PC1, probs = c(0.001, 0.999), na.rm = TRUE)
    qy <- quantile(scores$PC2, probs = c(0.001, 0.999), na.rm = TRUE)
    
    # Combine with arrow ranges so arrows are not clipped
    xr <- range(c(qx, arrow_df$PC1), na.rm = TRUE)
    yr <- range(c(qy, arrow_df$PC2), na.rm = TRUE)
    
    # Small padding
    pad_x <- diff(xr) *0.3
    pad_y <- diff(yr) *0.3
    
    xlim_all <- c(xr[1] - pad_x, xr[2] + pad_x)
    ylim_all <- c(yr[1] - pad_y, yr[2] + pad_y)
  }
}

## 5F) combine and save----------
p_four <- (p_inat + p_field) / (p_env + p_taxon) &
  # Consistent axes across all facets
  coord_cartesian(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

p_four

ggsave("outputs/08_pca/p_biplot.jpeg",
       p_four, width = 14, height = 10, dpi = 450)
ggsave("outputs/08_pca/p_biplot.svg",
       p_four, width = 14, height = 10, dpi = 450)



#6) RDA: presence~is_island+source ------------------

X_std <- as.matrix(X)  # already standardized
meta  <- pca_df %>% select(source, presence, is_island) %>% 
  mutate(across(where(is.factor), droplevels))

# Multivariate linear model in RDA form
mod_rda <- rda(X_std ~ presence * is_island + source, data = meta)

## POST HOC---------------
# Overall permutation test
perm_overall <- anova.cca(mod_rda, permutations = 999)   
# Marginal (type III–like) tests for each term
perm_margin  <- anova.cca(mod_rda, permutations = 999, by = "margin")
perm_terms <- anova.cca(mod_rda, permutations = 999, by = "terms")

perm_overall
perm_margin
perm_terms

capture.output(perm_overall, file = "outputs/08_pca/permanova_overall.txt")
capture.output(perm_margin,  file = "outputs/08_pca/permanova_margin.txt")
capture.output(perm_terms,  file = "outputs/08_pca/permanova_terms.txt")


## Logistic model on iNat PC scores (presence ~ PC1 * PC2 + setting)-------------------------------------------------------------
inat_scores <- scores %>% #filter(source == "iNat") %>% 
  mutate(presence = factor(presence, c(0,1), c("absence","presence"))) %>% 
  bind_cols(X %>% select(-is_island)) 

# Sanitize column names
names(inat_scores) <- make.names(names(inat_scores))

inat_scores %>% 
  select(CHELSA_bio7:ncol(.)) %>% names() ->glm_vars

# Create formula 
glm_formula <- as.formula(
  paste("I(presence == 'presence') ~ PC1 * PC2 +", paste(glm_vars, collapse = " + "))
)


glm <- glm(glm_formula, data = inat_scores, family = binomial())

summary_glm <- summary(glm)


capture.output(summary_glm, file = "outputs/08_pca/glm_presence_on_pcs.txt")


cor(inat_scores$lat_cdeg, inat_scores$PC1)
cor(inat_scores$lat_cdeg, inat_scores$PC2)

