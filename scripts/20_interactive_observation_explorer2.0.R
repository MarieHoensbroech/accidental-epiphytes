# 20_interactive_observation_explorer.R
#
# My first ever shiny app 
# (AI taught me a lot here)

# One can browse accidental epiphyte observations,
# filter by quality grade / user / rank / family / taxon, and visualize:
#  - Observation points
#  - iNaturalist taxon ranges (exported in step 2)
#  - Observation-derived ranges (generated in Step 2)
#  
# See my field sites and our sailed track 
#
# Every click will result in the ranges showing (if layer activated),
# and a side panel that shows the observation details including photos.
#
# Below the main map appears a brief summary of total observations/count and 
# the %overlap of the iNat taxon range and the obs-derived range
#
# Check it out!
# https://mariehoensbroech.shinyapps.io/accidental_epiphytes_interactive_observation_explorer/

rm(list=ls())

# Libraries ---------------
suppressPackageStartupMessages({
  library(httr2)
  library(units)
  library(tidyverse)
  library(data.table)
  library(shiny)
  library(bslib)
  library(leaflet)
  library(sf)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(lubridate)
  library(htmltools)
  library(leaflet)
})
# --- Helpers------------
#GitHub directory listing with headers + token fallback
# Uses httr2 if available; otherwise safely returns NULL.
github_list_dir <- function(owner = "MarieHoensbroech",
                            repo  = "accidental-epiphytes",
                            path) {
  api_url <- sprintf(
    "https://api.github.com/repos/%s/%s/contents/%s",
    owner, repo, utils::URLencode(path, reserved = TRUE)
  )
  
  # Prefer httr2 for better headers & auth
  if (requireNamespace("httr2", quietly = TRUE)) {
    # Build request step-by-step (no braces inside |>):
    req <- httr2::request(api_url)
    req <- httr2::req_method(req, "GET")
    req <- httr2::req_user_agent(req, "accidental-epiphytes-shiny-app")
    req <- httr2::req_headers(req, Accept = "application/vnd.github.v3+json")
    
    # Optional bearer token (local or shinyapps.io env var)
    tok <- Sys.getenv("GITHUB_TOKEN", unset = "")
    if (nzchar(tok)) {
      req <- httr2::req_auth_bearer_token(req, tok)
    }
    
    # Perform safely
    resp <- try(httr2::req_perform(req), silent = TRUE)
    if (inherits(resp, "try-error")) return(NULL)
    if (httr2::resp_status(resp) != 200L) return(NULL)
    
    return(httr2::resp_body_json(resp, simplifyVector = TRUE))
  }
  
  # Fallback: plain JSON fetch (may hit 403/429; then return NULL)
  out <- try(jsonlite::fromJSON(api_url), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

# regex & selection helper (was collapsing with "\n") 
pick_first_download <- function(df, exts) {
  if (is.null(df) || NROW(df) == 0) return(NA_character_)
  if (!all(c("name","download_url","type") %in% names(df))) return(NA_character_)
  dd <- df[df$type == "file", , drop = FALSE]
  if (NROW(dd) == 0) return(NA_character_)
  dd <- dd[order(tolower(dd$name)), , drop = FALSE]
  # Use OR '|' between extensions (previously "\n" prevented any match)
  rx <- paste0("\\.(", paste(exts, collapse="|"), ")$", collapse = "")
  hit <- dd[grepl(rx, dd$name, ignore.case = TRUE), , drop = FALSE]
  if (NROW(hit) == 0) return(NA_character_)
  hit$download_url[1]
}


# Parameters -----------------------------------
POINTS_LYR  <- "observations"
INAT_CSV    <- "data/processed/inat_observations.csv"  # optional (photos)
LEGEND_COLS <- c(research = "#1a9641", casual = "#fdae61", needs_id = "#d7191c") 
arrow_col <- "black"
arrow_col_url <- gsub("#", "%23", arrow_col, fixed = TRUE)

# Data source: my GitHub repo under /appdata
GITHUB_MODE     <- TRUE
GITHUB_RAW_BASE <- "https://raw.githubusercontent.com/MarieHoensbroech/accidental-epiphytes/main/appdata"

# Relative filenames under that folder
INAT_REMOTE_REL <- "inat_observations.csv"
MERGED_REMOTE   <- "inat.merged.csv"
SITES_REMOTE    <- "my.sites.csv"

# _______________-----------------
# Load and prep data ------------------------------------
## Ranges-------------
range_lookup <- fread(file.path(GITHUB_RAW_BASE,"range_lookup.csv")) %>% as_tibble() %>% 
  mutate(taxon_id=as.character(taxon_id))

## Field sites---------
my_sites <- fread(file.path(GITHUB_RAW_BASE,"my.sites.csv")) %>% as_tibble() %>% distinct()

## Field observations ----------
inat_merged <- fread(file.path(GITHUB_RAW_BASE,"inat.merged.csv")) %>% as_tibble() %>% 
  reframe(TreeID,id,EpiID,site,id=as.character(id)) %>% 
  distinct() %>% 
  right_join(my_sites %>% group_by(site) %>% 
               reframe(siteLat,siteLon,TotalTrees=sum(NumberTrees)) %>% 
               distinct()) %>% 
  right_join(my_sites %>% select(site,Name) %>% unique())

## All iNat observations -----------
obs_sf <- fread(file.path(GITHUB_RAW_BASE,"inat_observations.csv")) %>% as_tibble() %>% 
  reframe(latitude, longitude, id, taxon_id, uri, taxon.name, taxon.common, family_name,
          photo_urls_concat,
          quality_grade, user_login, taxon.rank, description,
          EpiSize, n_obs, sum_individuals,
          EpiCount, GrowingSite,
          EpiHeight=EpiHeight_num, Moss, TreeSp, DBH, TreeHeight, Setting) %>% 
  mutate(id = as.character(id),
         taxon_id = as.character(taxon_id),
         photo_urls_list = str_split(photo_urls_concat, "\\|")) %>% 
  distinct() %>% 
  full_join(inat_merged) %>% 
  mutate(source=if_else(is.na(TreeID),"iNat data","Field data")) %>% 
  distinct(id,EpiHeight,.keep_all = T) %>% 
  left_join(range_lookup) %>% 
  drop_na(latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  distinct()

## Quality grade colour coding -----------
obs_shiny <- obs_sf %>% 
  unique() %>% 
  mutate(grade_key = case_when(
    quality_grade %in% names(LEGEND_COLS) ~ quality_grade,
    TRUE ~ "needs_id"
  )) %>%
  mutate(col = unname(LEGEND_COLS[grade_key])) %>% 
  distinct(id,.keep_all = T)


## Journey tracks and field site polygons-----------------
zip_tracks <- "https://raw.githubusercontent.com/MarieHoensbroech/accidental-epiphytes/main/appdata/fieldsitesGIS/tracks/tracks.zip"
zip_sites <- "https://raw.githubusercontent.com/MarieHoensbroech/accidental-epiphytes/main/appdata/fieldsitesGIS/polygons/fieldsites.zip"

tracks <- st_read(paste0("/vsizip/vsicurl/", zip_tracks), 
                  layer = "tracks3110hondarribia", quiet = FALSE) 

sites  <- st_read(paste0("/vsizip/vsicurl/", zip_sites),  
                  layer = "fieldsites", quiet = FALSE)

tracks_wgs84 <- if (!is.null(tracks)) tryCatch(st_transform(tracks, 4326), error = function(e) tracks) else NULL
sites_wgs84  <- if (!is.null(sites))  tryCatch(st_transform(sites,  4326), error = function(e) sites)  else NULL

if (!is.null(sites_wgs84)) {
  
  sites_wgs84$Name_safe <- ifelse(is.na(sites_wgs84$Name) | !nzchar(sites_wgs84$Name),
                                  paste0("Site_", seq_len(nrow(sites_wgs84))),
                                  sites_wgs84$Name)
  
  # areas (ha) with numeric fallback; also keep mean for imputation
  sites_m <- tryCatch(st_transform(sites_wgs84, 3857), error = function(e) sites_wgs84)
  sites_wgs84$area_ha_calc <- suppressWarnings(as.numeric(st_area(sites_m)) / 10000)
  sites_wgs84$area_ha_mean <- suppressWarnings(mean(sites_wgs84$area_ha_calc[is.finite(sites_wgs84$area_ha_calc)], na.rm = TRUE))
}

sites_wgs84<- sites_wgs84 %>% 
  mutate(
    Name_safe = if_else(Name_safe=="Site5_LeidenCanalWllows","Site4_LeidenCanalWillows",Name_safe),
    Name_safe = if_else(Name_safe=="Site4_LeidenDunes","Site5_LeidenDunes",Name_safe)
  )



tracks_wgs84 <- tracks_wgs84 %>%
  mutate(
    dt_txt = str_extract(Name, "\\b\\d{1,2} [A-Z]{3} \\d{4} \\d{2}:\\d{2}\\b"),
    dt = dmy_hm(dt_txt, locale = "C"),
    date_str = format(dt, "%d %b %Y"),
    km_val = as.numeric(distance_k),
    len_m = st_length(geometry),
    len_km = as.numeric(set_units(len_m, "km")),
    Name_safe = sprintf("%s; %s km",
                        date_str,
                        format(round(km_val, 2), nsmall = 1))
  ) %>%
  select(-dt_txt, -dt, -date_str, -km_val)

tracks_wgs84 <-tracks_wgs84 %>% 
  distinct(Name_safe,.keep_all = T) 

sites_wgs84_clean <- sites_wgs84 %>%
  st_make_valid() %>%
  filter(!st_is_empty(geometry))

# Centroids for star markers so sites remain visible when zoomed out
sites_centroids <- if (!is.null(sites_wgs84)) st_centroid(sites_wgs84) else NULL

sites_centroids_clean <- sites_centroids %>%
  st_make_valid() %>%
  filter(!st_is_empty(geometry))

tracks_wgs84_clean <- tracks_wgs84 %>%
  st_make_valid() %>%
  filter(!st_is_empty(geometry))

tracks_wgs84_clean <- tracks_wgs84_clean %>%
  mutate(
    dt_txt = str_extract(Name, "\\b\\d{1,2} [A-Z]{3} \\d{4} \\d{2}:\\d{2}\\b"),
    dt = dmy_hm(dt_txt, locale = "C"),
    date_str = ifelse(!is.na(dt), format(dt, "%d %b %Y"), NA_character_),
    km_val = as.numeric(distance_k),
    date_str = replace_na(date_str, "Unknown date"),
    km_val = replace_na(km_val, 0),
    Name_safe = sprintf("%s; %s km",
                        date_str,
                        format(round(km_val, 1), nsmall = 1))
  ) %>%
  select(-dt_txt, -dt, -date_str, -km_val)

sites_wgs84_clean <- sites_wgs84_clean %>%
  mutate(
    Name_safe = as.character(replace_na(Name_safe, "Unnamed site")),
    layer_id  = paste0("site_poly_", make.names(Name_safe, unique = FALSE))
  )

sites_centroids_clean <- sites_centroids_clean %>%
  mutate(
    Name_safe = as.character(replace_na(Name_safe, "Unnamed site")),
    layer_id  = paste0("site_marker_", make.names(Name_safe, unique = FALSE))
  )

tracks_wgs84_clean <- tracks_wgs84_clean %>%
  mutate(
    Name_safe = as.character(replace_na(Name_safe, "Unknown track")),
    layer_id  = paste0("track_", make.names(Name_safe, unique = FALSE))
  )


#_____________________________------------------------
## "Last updated" text  -----------------------
# Get "Last updated" for a GitHub raw file:
url <- file.path(GITHUB_RAW_BASE, INAT_REMOTE_REL)
resp <- request(url) |> req_method("HEAD") |> req_user_agent("shiny-app") |> req_perform()
hdrs <- resp_headers(resp)
lm <- hdrs[["date"]] %||% hdrs[["Date"]]


# Extract and format the date
d <- format(as.POSIXct(lm, format = "%a, %d %b %Y %H:%M:%S", tz = "GMT"), "%d %b %Y")

last_updated_text <- paste0("Last updated on ",d,"by Marie Hoensbroech.")

#____________________----------------
# SITE SUMMARY POP-UP (MP4 ONLY) -------------------------------------------
site_summary_html <- function(site_key) {
  # -- Helpers (scoped) --------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a
  
  build_video_tag_mp4 <- function(poster, mp4_url) {
    # Returns a <video> tag with a single MP4 <source>. Poster is optional but recommended.
    poster_attr <- htmltools::htmlEscape(poster %||% "")
    mp4_attr    <- htmltools::htmlEscape(mp4_url)
    paste0(
      "<video controls preload='metadata' crossorigin='anonymous' playsinline ",
      "style='width:100%; max-height:240px; margin:6px 0;' poster='", poster_attr, "'>",
      "<source src='", mp4_attr, "' type='video/mp4'/>",
      "Your browser does not support the video tag.",
      "</video>"
    )
  }
  
  # -- 1) Tables & counts ------------------------------------------------------
  ms <- my_sites |> dplyr::filter(Name == site_key)
  
  total_trees <- if (nrow(ms) > 0) sum(ms$NumberTrees, na.rm = TRUE) else NA_integer_
  
  dbh_tbl <- if (nrow(ms) > 0) {
    ms |>
      dplyr::group_by(TreeCat) |>
      dplyr::summarise(Trees = sum(NumberTrees, na.rm = TRUE), .groups = "drop")
  } else {
    tibble::tibble(TreeCat = character(), Trees = integer())
  }
  
  site_codes <- unique(ms$site)
  obs_here <- tryCatch(
    obs_shiny |>
      sf::st_drop_geometry() |>
      dplyr::filter(source == "Field data", site %in% site_codes, !is.na(EpiID)),
    error = function(e) tibble::tibble()
  )
  n_epis  <- if (nrow(obs_here) > 0) dplyr::n_distinct(obs_here$EpiID)  else 0
  n_hosts <- if (nrow(obs_here) > 0) dplyr::n_distinct(obs_here$TreeID) else 0
  
  host_species_vec <- if (nrow(obs_here) > 0) {
    obs_here |>
      dplyr::distinct(TreeSp) |>
      dplyr::filter(!is.na(TreeSp), nzchar(TreeSp)) |>
      dplyr::pull(TreeSp) |>
      unique() |>
      sort()
  } else character(0)
  host_species_str <- if (length(host_species_vec) > 0)
    paste(host_species_vec, collapse = ", ")
  else
    "none recorded"
  
  # -- 2) Match site polygon row ----------------------------------------------
  sites_tbl <- if (exists("sites_wgs84_clean", inherits = TRUE) && !is.null(get("sites_wgs84_clean"))) {
    get("sites_wgs84_clean")
  } else {
    get("sites_wgs84")
  }
  
  srow <- try(
    sites_tbl |>
      dplyr::filter(.data$Name == site_key | .data$Name_safe == site_key),
    silent = TRUE
  )
  
  if (inherits(srow, "try-error") || nrow(srow) < 1) {
    return(paste0(
      "<div style='font-size:12px;'>",
      "<b>", htmltools::htmlEscape(site_key), "</b><br/>",
      "<b>Total trees:</b> ",
      ifelse(is.na(total_trees), "(n/a)", format(total_trees, big.mark = ",")),
      "<br/><em>Site geometry not found.</em>",
      "</div>"
    ))
  }
  
  # Area (fallback to mean)
  area_val <- suppressWarnings(as.numeric(srow$area_ha_calc[1]))
  used_mean <- FALSE
  if (!is.finite(area_val) || is.na(area_val)) {
    area_val <- suppressWarnings(as.numeric(sites_tbl$area_ha_mean[1]))
    used_mean <- TRUE
  }
  
  # Dominant species list
  dom_cols <- intersect(paste0("dominanttree", 1:6), names(ms))
  dom_vals <- ms |>
    dplyr::select(dplyr::all_of(dom_cols)) |>
    tidyr::pivot_longer(dplyr::everything(), values_to = "sp") |>
    dplyr::mutate(sp = trimws(sp)) |>
    dplyr::filter(!is.na(sp), nzchar(sp), sp != '""') |>
    dplyr::count(sp, sort = TRUE)
  dom_list <- if (nrow(dom_vals) > 0)
    paste(dom_vals$sp[1:min(5, nrow(dom_vals))], collapse = ", ")
  else
    "(not recorded)"
  
  dbh_html <- if (nrow(dbh_tbl) > 0) {
    paste0(
      "<ul>",
      paste(
        sprintf(
          "<li>%s: %s</li>",
          htmltools::htmlEscape(dbh_tbl$TreeCat),
          format(dbh_tbl$Trees, big.mark = ",")
        ),
        collapse = ""
      ),
      "</ul>"
    )
  } else "<em>No DBH data for this site.</em>"
  
  # -- 3) MEDIA (MP4 only) -----------------------------------------------------
  folder_name <- dplyr::coalesce(
    if ("Name_clean" %in% names(srow)) srow$Name_clean[1] else NA_character_,
    if ("Name_safe"  %in% names(srow)) srow$Name_safe[1]  else NA_character_,
    site_key
  )
  folder_rel <- sprintf("fieldwork_media/%s", folder_name)
  
  contents <- github_list_dir(path = folder_rel)   # may be NULL on 403/429
  # MP4 ONLY:
  vid_mp4 <- pick_first_download(contents, exts = c("mp4"))
  poster  <- pick_first_download(contents, exts = c("jpg", "jpeg", "png"))
  
  raw_base <- sprintf(
    "https://raw.githubusercontent.com/MarieHoensbroech/accidental-epiphytes/main/fieldwork_media/%s/",
    utils::URLencode(folder_name, reserved = TRUE)
  )
  # If no poster is listed by API, we still try conventional cover.jpg (ok if 404; it’s just an <img> attr)
  if (is.na(poster)) poster <- paste0(raw_base, "cover.jpg")
  
  # Build the video HTML: only if we actually found an MP4 in the folder
  if (!is.na(vid_mp4)) {
    video_html <- build_video_tag_mp4(poster, vid_mp4)
  } else {
    # No MP4 available → show poster + small note (no MOV embed)
    poster_attr <- htmltools::htmlEscape(poster %||% "")
    video_html <- paste0(
      "<img src='", poster_attr, "' alt='Site cover' ",
      "style='width:100%; max-height:240px; object-fit:cover; margin:6px 0;'/>",
      "<div style='font-size:11px; color:#666;'>No MP4 available in this folder yet.</div>"
    )
  }
  
  # GitHub media folder link (proper <a>)
  gh_folder_link <- sprintf(
    "https://github.com/MarieHoensbroech/accidental-epiphytes/tree/main/fieldwork_media/%s",
    utils::URLencode(folder_name, reserved = TRUE)
  )
  gh_folder_anchor <- paste0(
    "<a href='", gh_folder_link, "' target='_blank' rel='noopener'>",
    "Open media folder (photos &amp; videos)</a>"
  )
  
  # -- 4) Final HTML -----------------------------------------------------------
  paste0(
    "<div style='font-size:12px;'>",
    "<b>", htmltools::htmlEscape(site_key), "</b><br/>",
    "<b>Total trees:</b> ",
    ifelse(is.na(total_trees), "(n/a)", format(total_trees, big.mark = ",")), "<br/>",
    "<b>Trees per DBH category (cm):</b> ", dbh_html,
    "<b>Epiphytes found:</b> ", format(n_epis, big.mark = ","), " on <b>", format(n_hosts, big.mark = ","), "</b> ",
    "host trees of species: ", htmltools::htmlEscape(host_species_str), "<br/>",
    "<b>Area:</b> ", round(area_val, 3), ifelse(used_mean, " (imputed mean)", ""), "<br/>",
    "<b>Dominant tree species:</b> ", htmltools::htmlEscape(dom_list),
    "<div style='margin-top:6px;'></div>",
    video_html,
    "<div style='margin-top:4px;'>", gh_folder_anchor, "</div>",
    "</div>"
  )
}




`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

#_________________________-----------------
#Custom marker----------------
pin_svg <- paste0(
  "data:image/svg+xml;utf8,",
  utils::URLencode(
    sprintf(
      '<svg xmlns="http://www.w3.org/2000/svg" width="28" height="38" viewBox="0 0 28 38">
         <!-- Pin background with thin beige outline -->
         <path d="M14 0
                  C6.3 0 0 6.3 0 14
                  c0 9.5 12.2 23.3 13.4 24.6
                  a1 1 0 0 0 1.2 0
                  C15.8 37.3 28 23.5 28 14
                  C28 6.3 21.7 0 14 0z"
               fill="#B2DFDB"
               stroke="%s"
               stroke-width="1.5"
               stroke-linejoin="round"/>
         <!-- Star (light beige) -->
         <polygon points="14,8 16.2,12.9 21.6,13.3 17.2,16.6 18.7,21.8 14,19 9.3,21.8 10.8,16.6 6.4,13.3 11.8,12.9"
                  fill="%s"/>
       </svg>',
      "#F5E6CC", "#F5E6CC"
    ),
    reserved = TRUE)
)

site_icon <- icons(
  iconUrl     = pin_svg,
  iconWidth   = 28, iconHeight = 38,
  iconAnchorX = 14, iconAnchorY = 38,   # pin tip anchors to point
  popupAnchorX = 0,  popupAnchorY = -38
)

# __________________________ ---------------------
# BUILD THE APP ------------------------------
# UI -----------------------
ui <- fluidPage(
  theme = bs_theme(version = 5),
  
  
  # Global font size reduction (~2pt)
  
  tags$head(
    tags$style(HTML("
    h6 {
      margin-bottom: 1px;
      margin-top: 1px;
    }

    :root { --app-font-reduce: 2pt; }
    body { font-size: calc(1rem - var(--app-font-reduce)); }

    /* Keep common UI bits in sync with base scaling */
    .form-label, .selectize-input, .selectize-dropdown, .btn,
    .card, .card-header, .leaflet-control, label, p, h6, h5 {
      font-size: calc(1rem - var(--app-font-reduce));
    }

    /* Ensure Leaflet tooltips are visible above Bootstrap elements */
    .leaflet-tooltip {
      opacity: 1 !important;
      display: block !important;
      z-index: 10000 !important;
    }
  "))
  ),
  
  
  
  tags$head(
    tags$meta(
      `http-equiv` = "Content-Security-Policy",
      content = paste(
        "default-src 'self' blob: data: https:;",
        # allow images and video posters from GitHub raw
        "img-src 'self' data: https: blob: https://raw.githubusercontent.com https://*.githubusercontent.com;",
        # allow media/video to load from these
        "media-src 'self' https://raw.githubusercontent.com https://*.githubusercontent.com data: blob:;",
        "style-src 'self' 'unsafe-inline' https: data:;",
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https:;",
        "font-src 'self' data: https:;",
        sep = " "
      )
    )
  ),
  
  
  titlePanel("Sailing for Science: mapping accidental epiphytes along the Eastern Atlantic temperate coast"),
  
  sidebarLayout(
    # Narrower filter column
    sidebarPanel(
      width = 2,
      actionButton("reset_filters", "Reset filters", icon = icon("redo")),
      selectizeInput("grade_filter",  "Select quality grade:", choices = NULL, selected = "All"),
      selectizeInput("user_filter",   "Select user:",         choices = NULL, selected = "All"),
      selectizeInput("rank_filter",   "Select taxon rank:",   choices = NULL, selected = "All"),
      selectizeInput("family_filter", "Select family:",       choices = NULL, selected = "All"),
      selectizeInput("taxon_filter",  "Select taxon:",        choices = NULL, selected = "All")
    ),
    
    mainPanel(
      width = 10,
      fluidRow(
        column(
          width = 8,
          leafletOutput("map", height = "600px")
        ),
        column(
          width = 4,
          bslib::card(
            bslib::card_header("Observation Details"),
            uiOutput("obs_details")
          )
        )
      ),
      uiOutput("selected_taxon_summary"),
      tags$p(last_updated_text, style = "font-size: 12px; color: gray; margin-top: 10px;")
    )
  )
)



# Server --------------------------------
server <- function(input, output, session) {
  
  # Welcome message -------------------
  observeEvent(TRUE, {
    welcome_img_file <- "welcome.jpg"  
    inat_link <- "https://www.inaturalist.org/projects/accidental-epiphytes"
    inat_ref <- paste0(
      "<a href='", inat_link, "' target='_blank' rel='noopener'>",
      "Accidental Epiphyte Project</a> on iNaturalist (at the time of publishing this app)."
    )
    
    # Image URLs
    img_src <- sprintf(
      "https://raw.githubusercontent.com/MarieHoensbroech/accidental-epiphytes/main/fieldwork_media/000_fieldwork_highlights/%s",
      utils::URLencode(welcome_img_file, reserved = TRUE)
    )
    
    showModal(modalDialog(
      title = "Welcome aboard!",
      easyClose = TRUE, footer = NULL, size = "l",
      HTML(paste0(
        # CSS tweak for title spacing
        "<style>.modal-title {margin-bottom:1px !important;}</style>",
        
        # Subtitle block
        "<div style='font-size:16px; font-weight:500; color:#444; margin-top:0; margin-bottom:8px;'>",
        "🗺️ This map documents a 2025 sailing expedition combined with ecological fieldwork across multiple sites between Lelystad, Netherlands, and Hondarribia, Basque Country. ",
        "The project blends science, art, and open data to explore accidental epiphyte distribution—currently in Europe—while promoting sustainable research practices through citizen science platforms like iNaturalist.<br/><br/>",
        "🌿 Follow our sailing route, field locations, and all accidental epiphyte observations submitted to our ",
        inat_ref,
        "</div>",
        
        # Image block with caption
        "<div style='margin:0 0 12px 0; overflow:hidden;'>",
        "<img src='", img_src, "' alt='Welcome image' ",
        "style='width:100%; height:auto; display:block; object-fit:cover; max-height:400px;'>",
        "<div style='font-size:12px; color:#555; text-align:center; margin-top:4px;'>",
        "S/V Salka Valka is a Norwegian-built Hero 107 CS Ketch.",
        "</div>",
        "</div>",
        
        # Body text
        "<div style='line-height:1.4'>",
        "<b>Layers & basemaps</b><br/>",
        "• Use the top-right layer control to toggle <em>All observations</em>, <em>Selected taxon</em>, <em>Ranges</em>, <em>Field sites</em>, and <em>Salka Valka's track</em>.<br/><br/>",
        
        "<b>Filtering observations</b><br/>",
        "• Left panel lets you filter by quality grade, user, rank, family, and taxon. Click <em>Reset filters</em> anytime.<br/><br/>",
        
        "<b>Exploring observations</b><br/>",
        "• Click a point to see details (with photos) of the observation on the right.<br/>",
        "• The <em>Ranges</em> overlay shows the iNat taxon range and the accidental-epiphyte range.<br/><br/>",
        
        "<b>Field sites</b><br/>",
        "• Hover over tracks to see when and how far we sailed.<br/>",
        "• Click a site polygon (zoom in!) to open a popup with totals, DBH structure, dominant species, and a link to photos & videos.<br/><br/>",
        
        "<b>Tips</b><br/>",
        "• Switch basemaps (Minimal, Topography, Satellite) from the layer control.<br/>",
        "• Popups can be closed by clicking the map background.",
        "</div>"
      ))
    ))
  }, once = TRUE)
  
  # Auto-reset filters on app start----------
  observe({
    isolate({
      session$sendCustomMessage("resetFilters", "All")
    })
  })
  
  observe({
    updateSelectizeInput(
      session, "taxon_filter",
      choices = NULL,
      selected = "All",
      server = TRUE
    )
  })
  observe({ updateSelectizeInput(session, "user_filter",   selected = "All") })
  observe({ updateSelectizeInput(session, "grade_filter",  selected = "All") })
  observe({ updateSelectizeInput(session, "family_filter", selected = "All") })
  
  # Reactive Values --------------
  selected_taxon    <- reactiveVal(NULL)
  clicked_obs       <- reactiveVal(NULL)
  clicked_obs_data  <- reactiveVal(NULL)
  matching_points   <- reactiveVal(NULL)
  
  # Filter choices ----------------------
  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "grade_filter",
                         choices = c("All", sort(unique(obs_shiny$quality_grade))),
                         selected = "All")
    updateSelectizeInput(session, "user_filter",
                         choices = c("All", sort(unique(obs_shiny$user_login))),
                         selected = "All")
    updateSelectizeInput(session, "rank_filter",
                         choices = c("All", sort(unique(obs_shiny$taxon.rank))),
                         selected = "All")
    updateSelectizeInput(session, "family_filter",
                         choices = c("All", sort(unique(obs_shiny$family_name))),
                         selected = "All")
    updateSelectizeInput(session, "taxon_filter",
                         choices = c("All", sort(unique(obs_shiny$taxon.name))),
                         selected = "All")
  })
  
  # Trigger filter reset on startup (once)
  session$onFlushed(function() {
    updateSelectizeInput(session, "grade_filter",
                         choices = c("All", sort(unique(obs_shiny$quality_grade))),
                         selected = "All")
    updateSelectizeInput(session, "user_filter",
                         choices = c("All", sort(unique(obs_shiny$user_login))),
                         selected = "All")
    updateSelectizeInput(session, "rank_filter",
                         choices = c("All", sort(unique(obs_shiny$taxon.rank))),
                         selected = "All")
    updateSelectizeInput(session, "family_filter",
                         choices = c("All", sort(unique(obs_shiny$family_name))),
                         selected = "All")
    updateSelectizeInput(session, "taxon_filter",
                         choices = c("All", sort(unique(obs_shiny$taxon.name))),
                         selected = "All")
  }, once = TRUE)
  
  # Filtering --------------
  filtered_obs <- reactive({
    data <- obs_shiny %>% ungroup()
    if (input$grade_filter  != "All") { data <- data %>% filter(quality_grade      == input$grade_filter)  }
    if (input$user_filter   != "All") { data <- data %>% filter(user_login         == input$user_filter)   }
    if (input$rank_filter   != "All") { data <- data %>% filter(taxon.rank         == input$rank_filter)   }
    if (input$family_filter != "All") { data <- data %>% filter(family_name  == input$family_filter) }
    if (input$taxon_filter  != "All") { data <- data %>% filter(taxon.name         == input$taxon_filter)  }
    data
  })
  
  # Output 1: Map -------------------
  output$map <- renderLeaflet({
    req(input$grade_filter, input$user_filter, input$rank_filter, input$family_filter, input$taxon_filter)
    data <- filtered_obs()

    if (nrow(data) == 0) {
      return(leaflet() %>% addTiles() %>% addControl("No observations to display", position = "topright"))
    }

    coords <- st_coordinates(data)

    map<-
      leaflet(data = data) %>%

      ### Basemaps --------------
      addProviderTiles(providers$CartoDB.Positron, group = "Minimal") %>%

      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%

      addProviderTiles(providers$OpenTopoMap, group = "Topography") %>%

      
      # Create custom panes with explicit z-index -------------
      addMapPane("rangesPane", zIndex = 400) %>%
      addMapPane("sitesPane", zIndex = 410) %>%
      addMapPane("tracksPane", zIndex = 420) %>%
      addMapPane("obsPane", zIndex = 430) %>%
      

      ### Add observations ------------------
    addCircleMarkers(
      lng = as.numeric(coords[,1]),
      lat = as.numeric(coords[,2]),
      layerId = ~as.character(id),
      radius = 5,
      color = ~col,
      fillOpacity = 0.6,
      group = "All observations",
      options = pathOptions(pane = "obsPane")
    ) %>%

      ### Add legend -------------

    addLegend("bottomright", colors = unique(obs_shiny$col),
              labels = unique(obs_shiny$grade_key), title = "Quality grade") %>%

      addLegend("bottomleft", colors = c("#90C987", "#F7F056"),
                labels = c("Taxon range", "Accidental epiphyte range"), title = "Range type") %>%

      ### Add field sites ------------
    addPolygons(
      data = sites_wgs84_clean,
      layerId = ~layer_id,
      color = "#3b6fb6", weight = 2, fill = TRUE, fillOpacity = 0.15, dashArray = "5,5",
      highlightOptions = highlightOptions(weight = 3, color = "#2c5282", bringToFront = TRUE),
      label = ~Name_safe,
      group = "Marie's field work: Field sites",
      options = pathOptions(pane = "sitesPane")
    ) %>%

      addMarkers(
        data = sites_centroids_clean,
        layerId = ~layer_id,
        icon = site_icon,
        label = ~Name_safe,
        group = "Marie's field work: Field sites",
        options = pathOptions(pane = "obsPane")
      ) %>%
      ### Add track ----------------
      addPolylines(
        data = tracks_wgs84_clean,
        layerId = ~layer_id,
        label = ~Name_safe,
        labelOptions = labelOptions(
          # sticky keeps the tooltip near the cursor
          sticky = TRUE, direction = "auto", opacity = 0.95,
          textsize = "12px", offset = c(0, -5)
        ),
        color = "cadetblue", weight = 3, opacity = 0.8,
        group = "Marie's field work: Salka Valka's track",
        options = pathOptions(pane = "tracksPane")
      ) %>%


      ### Add layers control --------------------
    addLayersControl(
      baseGroups = c("Minimal","Topography","Satellite"),
      overlayGroups = c("Ranges","All observations",
                        "Observations of the same taxon",
                        "Selected observation",
                        "Marie's field work: Field sites",
                        "Marie's field work: Salka Valka's track"),
      options = layersControlOptions(collapsed = FALSE))
  })
  
  
  
  # Output 2: Observation Details Pane ------------
  normalize_inat_size <- function(u) gsub("/(square|small|medium|large|original)/", "/medium/", u, perl = TRUE)
  
  output$obs_details <- renderUI({
    obs.clicked <- clicked_obs_data()
    if (is.null(obs.clicked) || !is.data.frame(obs.clicked) || nrow(obs.clicked) == 0) {
      return(h6("No observation selected."))
    }
    
    ## Photos  ---------
    # CSV uses 'photo_urls_concat' (pipe-separated). Split on "|" and trim.
    raw_urls <- unique(unlist(strsplit(obs.clicked$photo_urls_concat[1] %||% "", "\\|", fixed = FALSE)))
    raw_urls <- raw_urls[!is.na(raw_urls) & nzchar(raw_urls)]
    
    # Force 'medium' size on iNaturalist (fast & lightweight)
    photo_urls <- normalize_inat_size(raw_urls)
    
    n_photos <- length(photo_urls)
    carousel_id <- "photoCarousel"
    
    # Build carousel items without forced upscaling (no 'w-100')
    photo_items <- if (n_photos > 0) {
      tagList(lapply(seq_along(photo_urls), function(i) {
        tags$div(
          class = paste("carousel-item", if (i == 1) "active"),
          tags$img(
            src   = photo_urls[i],
            # Avoid blur: never upscale beyond container; keep aspect ratio; neutral bg
            style = paste(
              "display:block; margin:0 auto;",
              "max-height:320px; max-width:100%; height:auto;",
              "object-fit:contain; background:#f8f9fa;"
            )
          )
        )
      }))
    } else NULL
    
    # Scoped CSS to ensure arrows are visible and spacing is tidy
    ui_carousel <- if (n_photos > 0) {
      tagList(
        singleton(tags$head(
          tags$style(HTML(sprintf("
          #%s { --bs-carousel-control-color: #222; --bs-carousel-control-opacity: .9; }
          #%s .carousel-control-prev-icon, #%s .carousel-control-next-icon { filter: drop-shadow(0 0 1px rgba(0,0,0,.4)); }
          /* Reduce gap above/below the carousel */
          #%s { margin-top:6px; margin-bottom:8px; }
        ", carousel_id, carousel_id, carousel_id, carousel_id)))
        )),
        tags$div(
          id = carousel_id, class = "carousel slide", `data-bs-ride` = if (n_photos > 1) "carousel" else NULL,
          tags$div(class = "carousel-inner", photo_items),
          if (n_photos > 1) tags$button(
            class = "carousel-control-prev", type = "button",
            `data-bs-target` = paste0("#", carousel_id), `data-bs-slide` = "prev",
            tags$span(class = "carousel-control-prev-icon", `aria-hidden` = "true"),
            tags$span(class = "visually-hidden", "Previous")
          ),
          if (n_photos > 1) tags$button(
            class = "carousel-control-next", type = "button",
            `data-bs-target` = paste0("#", carousel_id), `data-bs-slide` = "next",
            tags$span(class = "carousel-control-next-icon", `aria-hidden` = "true"),
            tags$span(class = "visually-hidden", "Next")
          )
        ),
        if (n_photos > 1) tags$script(HTML(sprintf("
        (function(){
          var el = document.getElementById('%s');
          if (el && window.bootstrap && bootstrap.Carousel) {
            try { new bootstrap.Carousel(el, { interval: 3000, pause: 'hover', wrap: true }); } catch(e){}
          }
        })();", carousel_id)))
      )
    } else NULL
    
    # Details panel ----
    tagList(
      h6(
        tags$a(href = unique(obs.clicked$uri)[1], "View on iNaturalist", target = "_blank", rel = "noopener")
      ),
      ui_carousel,  # <- the carousel goes here
      h6(strong("Taxon name:"), tagList(em(unique(obs.clicked$taxon.name)[1]),
                                        paste0(" (", unique(obs.clicked$quality_grade)[1], ")"))),
      h6(strong("Common name:"), unique(obs.clicked$taxon.common)[1]),
      h6(strong("Identified to:"), paste(unique(obs.clicked$taxon.rank)[1], " level")),
      h6(strong("Family:"), unique(obs.clicked$family_name)[1]),
      h5(strong("Observation fields")),
      h6(strong("Size of epiphyte (m):"), unique(obs.clicked$EpiSize)[1]),
      h6(strong("Count of individuals:"), unique(obs.clicked$EpiCount)[1]),
      h6(strong("Growing site:"), unique(obs.clicked$GrowingSite)[1]),
      h6(strong("Height of epiphyte (m):"), unique(obs.clicked$EpiHeight)[1]),
      h6(strong("Moss:"), unique(obs.clicked$Moss)[1]),
      h6(strong("Host species:"), unique(obs.clicked$TreeSp)[1]),
      h6(strong("Host DBH (m):"), unique(obs.clicked$DBH)[1]),
      h6(strong("Host height (m):"), unique(obs.clicked$TreeHeight)[1]),
      h6(strong("Setting:"), unique(obs.clicked$Setting)[1]),
      h6(strong(tagList("Notes from observer (", tags$em(unique(obs.clicked$user_login)[1]), "): ")),
         unique(obs.clicked$description)[1])
    )
  })
  
  
  # Observe event 1: Click handling & map updates -------------
  observeEvent(input$map_marker_click, {
    click <- input$map_marker_click
    clicked_id <- as.character(click$id)
    clicked <- obs_shiny[obs_shiny$id == clicked_id, ]
    if (nrow(clicked) == 0) return()
    
    clicked_obs_data(clicked)
    
    # Defensive extraction of taxon_id
    clicked_taxon_id <- clicked$taxon_id[1]
    if (is.na(clicked_taxon_id)) {
      warning("clicked_taxon_id is NA")
      return()
      
    }
    
    selected_taxon(unique(clicked$taxon.name)[1])
    
    mp <- obs_shiny %>% filter(taxon_id == clicked_taxon_id)
    matching_points(mp)  
    
    # Clear & show
    leafletProxy("map") %>%
      clearGroup("Observations of the same taxon") %>%
      clearGroup("All observations") %>%
      clearGroup("Ranges") %>%
      clearGroup("Selected observation") %>%
      showGroup("Ranges")
    
    # Observe event 2: Range polygons ----
    {
      # Convert obs_range to valid MULTIPOLYGON
      obs_range_geom  <- clicked$obs_range_wkt  %>%
        st_as_sfc(., crs = 4326)
      
      inat_range_geom <- clicked$inat_range_wkt %>%
        st_as_sfc(., crs = 4326)
      
      if (length(obs_range_geom)  > 0 && !st_is_empty(obs_range_geom)) {
        leafletProxy("map") %>%
          addPolygons(data = obs_range_geom, 
                      color = "#F7F056", weight = 2, 
                      fillOpacity = 0.4, group = "Ranges",
                      options = pathOptions(pane = "rangesPane"))
      }
      if (length(inat_range_geom) > 0 && !st_is_empty(inat_range_geom)) {
        leafletProxy("map") %>%
          addPolygons(data = inat_range_geom, 
                      color = "#90C987", weight = 2, 
                      fillOpacity = 0.3, group = "Ranges",
                      options = pathOptions(pane = "rangesPane"))
      }
    }
    
    # Observe event 3: (Highlighted) markers ----
    ## All observations (dimmed)--------------
    leafletProxy("map") %>%
      addCircleMarkers(
        data = obs_shiny,
        lng = as.numeric(st_coordinates(obs_shiny)[,1]),
        lat = as.numeric(st_coordinates(obs_shiny)[,2]),
        layerId = ~as.character(id),
        radius = 5,
        color = ~col,
        fillOpacity = 0.1,
        group = "All observations"
      )
    
    ## Observations of the same taxon (highlighted)-------------
    mp_now <- matching_points()
    if (!is.null(mp_now) && nrow(mp_now) > 0) {
      leafletProxy("map") %>%
        addCircleMarkers(
          data = mp_now,
          lng = st_coordinates(mp_now)[,1],
          lat = st_coordinates(mp_now)[,2],
          layerId = ~as.character(id),
          radius = 7,
          color = "darkgoldenrod",
          fillOpacity = 0.9,
          group = "Observations of the same taxon"
        )
    }
    
    ## Selected observation---------------
    leafletProxy("map") %>%
      addCircleMarkers(
        data = clicked,
        lng = as.numeric(st_coordinates(clicked)[,1]),
        lat = as.numeric(st_coordinates(clicked)[,2]),
        radius = 10,
        color = "white",
        weight = 3,
        fillOpacity = 1,
        group = "Selected observation"
      )
  
      
  }))
  
  ### Observe event 4: Field Site popup --------------
  
  observeEvent(input$map_shape_click, {
    e <- input$map_shape_click
    # Abort early if missing click info or data
    req(e$id, e$lng, e$lat)
    req(exists("sites_wgs84_clean"), !is.null(sites_wgs84_clean))
    
    srow <- dplyr::filter(sites_wgs84_clean, layer_id == e$id)
    req(nrow(srow) >= 1)
    
    site_key <- dplyr::coalesce(srow$Name[1], srow$Name_safe[1])
    
    html <- site_summary_html(site_key)
    leafletProxy("map") |>
      clearPopups() |>
      
      clearPopups() %>%
      addPopups(lng = e$lng, lat = e$lat,
                popup = htmltools::HTML(html),
                options = popupOptions(maxWidth = 420))
    
    
    # Reset marker click after handling
    isolate({
      session$sendCustomMessage("resetMarkerClick", NULL)
    
  })
  
  ### Observe event 5: Zoom to site when clicking a centroid marker -----------------------------
  observeEvent(input$map_marker_click, {
    click <- input$map_marker_click
    req(click$id)
    
    if (!startsWith(click$id, "site_marker_")) return()
    
    # Find the corresponding centroid row (to get Name_safe)
    req(exists("sites_centroids_clean"), !is.null(sites_centroids_clean))
    cent_row <- try(dplyr::filter(sites_centroids_clean, layer_id == click$id), silent = TRUE)
    if (inherits(cent_row, "try-error") || nrow(cent_row) < 1) return()
    
    # Match the polygon by Name_safe and get its bounding box
    site_key  <- cent_row$Name_safe[1]
    req(exists("sites_wgs84_clean"), !is.null(sites_wgs84_clean))
    poly_row  <- try(dplyr::filter(sites_wgs84_clean, Name_safe == site_key), silent = TRUE)
    if (inherits(poly_row, "try-error") || nrow(poly_row) < 1) return()
    
    bb <- sf::st_bbox(poly_row$geometry)
    
    # Optional: add a tiny padding so the polygon isn’t flush with the edges
    pad <- 0.0001  # ~100 m at mid-lats; increase if you like
    lng1 <- as.numeric(bb["xmin"]) - pad
    lat1 <- as.numeric(bb["ymin"]) - pad
    lng2 <- as.numeric(bb["xmax"]) + pad
    lat2 <- as.numeric(bb["ymax"]) + pad
    
    leafletProxy("map") %>%
      fitBounds(lng1 = lng1, lat1 = lat1, lng2 = lng2, lat2 = lat2)
    
    # (Optional) also open the site popup immediately at the centroid:
    # html <- site_summary_html(site_key)
    # leafletProxy("map") %>%
    #   clearPopups() %>%
    #   addPopups(lng = click$lng, lat = click$lat,
    #             popup = htmltools::HTML(html),
    #             options = popupOptions(maxWidth = 420))
  })
  
  
  ## Observe event 6: Summary: observation and individual count -----------
  output$selected_taxon_summary <- renderUI({
    req(selected_taxon(), matching_points(), clicked_obs_data)
    data <- matching_points()
    obs.clicked <- clicked_obs_data()
    req(NROW(data) > 0)
    
    # Use provided columns if present; otherwise fallbacks
    n_obs   <- unique(data$n_obs)
    n_obs   <- n_obs[!is.na(n_obs)]
    n_obs   <- if (length(n_obs) == 1) n_obs else nrow(data)  # fallback to count of rows
    
    sum_individuals <- unique(data$sum_individuals)
    sum_individuals <- sum_individuals[!is.na(sum_individuals)]
    sum_individuals <- if (length(sum_individuals) == 1) sum_individuals else suppressWarnings(sum(as.numeric(data$sum_individuals), na.rm = TRUE))
    
    req(!is.na(n_obs), !is.na(sum_individuals))
    
    taxon_name <- selected_taxon()
    
    
    # Base message 
    message_base <- if ((n_obs == 1 && sum_individuals > 1) || (n_obs == 1 && sum_individuals == 0)) {
      paste0("<i>", taxon_name, "</i> was observed ", n_obs,
             " time with a total of at least ", sum_individuals, " individuals.")
    } else if (n_obs == 1 && sum_individuals == 1) {
      paste0("<i>", taxon_name, "</i> was observed ", n_obs,
             " time with a total of at least ", sum_individuals, " individual.")
    } else {
      paste0("<i>", taxon_name, "</i> was observed ", n_obs,
             " times with a total of at least ", sum_individuals, " individuals.")
    }
    
    # Decide the overlap sentence
    # Treat additional labels as "range missing" 
    is_range_missing <- isTRUE(obs.clicked$status %in% c("range_missing"))
    
    # Be tolerant to tiny floating errors around zero
    is_zero_pct <- !is.na(obs.clicked$pct_inat) && abs(obs.clicked$pct_inat) < 1e-12
    
    message_extra <- if (is_range_missing) {
      "Overlap could not be calculated as one range is missing."
    } else if (is.na(obs.clicked$pct_inat)) {
      "Overlap could not be calculated."
    } else if (is_zero_pct || obs.clicked$pct_inat == 0) {
      paste0("The ranges of <i>", taxon_name,
             "</i> as accidentally epiphytic or terrestrial lifeform do not currently overlap.")
    } else {
      paste0("As accidental epiphyte, <i>", taxon_name, "</i> has been spotted in ",
             round(obs.clicked$pct_inat, 2), "% of its total range.")
    }
    
    # Combine and return as HTML
    HTML(paste(message_base, message_extra))})
  
  
  ### Observe event 7: Deselect everything when clicking on empty map space---------------
    observeEvent(input$map_click, {
      # Ignore initial render
      if (is.null(input$map_click$lat)) return()
      
      # Clear everything
      clicked_obs(NULL)
      clicked_obs_data(NULL)
      selected_taxon(NULL)
      matching_points(NULL)
      
      leafletProxy("map") %>%
        clearPopups() %>%
        clearGroup("Selected observation") %>%
        clearGroup("Observations of the same taxon") %>%
        clearGroup("Ranges")
      
      output$obs_details <- renderUI({ h6("No observation selected.") })
      output$selected_taxon_summary <- renderUI(NULL)
    })
}
# ______________________________ ---------------
# RUN APP-------------
shinyApp(ui, server)


# ___________________ -------------------
# DEPLOY ONLINE ----
# Currently deployed at:
# https://mariehoensbroech.shinyapps.io/accidental_epiphytes_interactive_observation_explorer/

# library(rsconnect)
# 
# # Path to the directory that contains app script
# proj_root <- "C:/Users/Marie/Documents/LatGradientsAccEpis/github"
# app_dir   <- file.path(proj_root, "scripts", "20_shiny_interactive_explorer")
# 
# # Deploy the app directory (the bundle will include app.R + appdata/* we just copied)
# rsconnect::deployApp(
#   appDir  = app_dir,
#   appPrimaryDoc = "20_shiny_interactive_observation_explorer.R",
#   appName = "accidental_epiphytes_interactive_observation_explorer"
# )

