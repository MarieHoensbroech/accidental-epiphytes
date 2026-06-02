# 07_APPENDIX_SPECIES_LIST.R
#
# Species list (family, genus, species) + site matrix
# formatted for Word appendix

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(stringr)
  library(officer)
  library(flextable)
})

# Parameters -----------------------------------
IN_INAT_FIELD <- "data/processed/inat.merged.csv"

OUT_DIR_LIST <- "outputs/08_species_list"
dir.create(OUT_DIR_LIST, showWarnings = FALSE, recursive = TRUE)

OUT_DOC <- file.path(OUT_DIR_LIST, "species_list.docx")

# Load data ------------------------------------

inat_field <- fread(IN_INAT_FIELD) %>%
  as_tibble() %>%
  clean_names() %>%
  distinct() %>%
  transmute(
    site,
    taxon_name,
    family = family_name,
    epi_count = as.integer(str_extract(epi_count, "\\d+"))
  )

# Prepare taxonomy ------------------------------

species_df <- inat_field %>%
  separate(taxon_name,
           into = c("genus", "species"),
           sep = " ",
           fill = "right") %>%
  mutate(
    species = replace_na(species, "sp."),
    present = "X"
  )

# Split into identified vs unidentified --------

identified <- species_df %>%
  filter(!is.na(family) & family != "")

unidentified <- species_df %>%
  filter(is.na(family) | family == "")

# Build matrix function -------------------------

make_matrix <- function(df) {
  df %>%
    distinct(family, genus, species, site, present) %>%
    pivot_wider(
      names_from = site,
      values_from = present,
      values_fill = ""
    )
}

identified_mat   <- make_matrix(identified)
unidentified_mat <- make_matrix(unidentified)

# Sort: family, genus, species, with "sp." last
identified_mat <- identified_mat %>%
  arrange(
    family,
    genus,
    species == "sp.",   # TRUE goes last
    species
  )

# Count sites
identified_mat <- identified_mat %>%
  mutate(
    n_sites = rowSums(across(where(is.character), ~ . == "X"))
  )

# Word-style grouping (blank repeated family/genus)
identified_mat <- identified_mat %>%
  group_by(family) %>%
  mutate(family = if_else(row_number() == 1, family, "")) %>%
  group_by(genus, .add = TRUE) %>%
  mutate(genus = if_else(row_number() == 1, genus, "")) %>%
  ungroup()

# Format unidentified block ---------------------

if (nrow(unidentified_mat) > 0) {
  
  unidentified_mat <- unidentified_mat %>%
    mutate(
      family = "Unidentified",
      genus = replace_na(genus, ""),
      species = replace_na(species, "")
    ) %>%
    arrange(genus, species) %>%
    group_by(family) %>%
    mutate(family = if_else(row_number() == 1, family, "")) %>%
    group_by(genus, .add = TRUE) %>%
    mutate(genus = if_else(row_number() == 1, genus, "")) %>%
    ungroup()
  
  unidentified_mat <- unidentified_mat %>%
    mutate(
      n_sites = rowSums(across(where(is.character), ~ . == "X"))
    )
  
}

# Combine tables --------------------------------

final_table <- bind_rows(
  identified_mat,
  unidentified_mat
)

# Reorder columns
final_table <- final_table %>%
  rename(Family=family, Genus=genus, Species=species) %>% 
  reframe(Family, Genus, Species, n_sites)

# Create Word table -----------------------------

ft <- flextable(final_table)

ft <- ft %>%
  autofit() %>%
  theme_booktabs() %>%
  align(j = 1:3, align = "left") %>%
  align(j = 4:ncol(final_table), align = "center") %>%
  valign(valign = "center")

ft <- compose(compose(ft, j = "Genus", value = as_paragraph(as_i(Genus))), 
              j = "Species", value = as_paragraph(as_i(Species)))

ft
# Save to Word ----------------------------------

doc <- read_docx() %>%
  body_add_par("Appendix 2: Species list", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = OUT_DOC)

message("Word species list saved to: ", OUT_DOC)

