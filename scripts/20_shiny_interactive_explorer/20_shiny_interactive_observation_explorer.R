#!/usr/bin/env Rscript
# =============================================================================
# 05_interactive_observation_explorer.R
# -----------------------------------------------------------------------------
# Purpose:
#   Lightweight Shiny explorer to browse accidental epiphyte observations,
#   filter by quality grade / user / rank / family / taxon, and visualize:
#     - Observation points
#     - iNaturalist taxon ranges (geomodel export)
#     - Observation-derived ranges (from Step 2)
#
# Inputs (created by Steps 1–3; with fallbacks for older runs):
#   - data/spatial/obs_points.gpkg   (layer: "observations")
#   - data/spatial/obs_ranges.gpkg   (layer: "obs_range")   [taxon_id polygons]
#   - data/spatial/inat_ranges.gpkg  (layer: "inat_range")  [taxon_id polygons]
#   - data/processed/inat_observations.csv (for photos & details)  [optional]
#   - (fallback) obs_sf.rds  — legacy object if GPKG not present
# =============================================================================

suppressPackageStartupMessages({
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
})

# ------------------------------- Parameters -----------------------------------
POINTS_LYR  <- "observations"
INAT_CSV    <- "appdata/inat_observations.csv"  # optional (photos)
LEGEND_COLS <- c(research = "#1a9641", casual = "#fdae61", needs_id = "#d7191c") 
arrow_col <- "black"
arrow_col_url <- gsub("#", "%23", arrow_col, fixed = TRUE)


# ------------------------------- Load data ------------------------------------
range_lookup <- fread("appdata/range_lookup.csv") %>% as_tibble() %>% 
  mutate(taxon_id=as.character(taxon_id))

obs_sf <- fread("appdata/inat_observations.csv") %>% as_tibble() %>% 
  select(latitude, longitude, id, taxon_id, uri, taxon.name, taxon.common, family_name,
         photo_urls_concat,
         quality_grade, user_login, taxon.rank, description,
         EpiSize, n_obs, sum_individuals,
         EpiCount, GrowingSite,
         EpiHeight, Moss, TreeSp, DBH, TreeHeight, Setting) %>% 
  mutate(id = as.character(id),
         taxon_id = as.character(taxon_id),
         photo_urls_list = str_split(photo_urls_concat, "\\|")) %>% 
  left_join(range_lookup) %>% 
  drop_na(latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Map a color per point by quality_grade
obs_shiny <- obs_sf %>% 
  unique() %>% 
  mutate(grade_key = case_when(
  quality_grade %in% names(LEGEND_COLS) ~ quality_grade,
  TRUE ~ "needs_id"
)) %>%
  mutate(col = unname(LEGEND_COLS[grade_key]))



# Compute "Last updated" from the freshest input
last_updated_text <- paste0("Last updated on: ",
                            format(as_date(
                              file.info(INAT_CSV)$mtime), 
                              "%d %B %Y"))
##### UI #####
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

    /* Keep common UI bits in sync with base scaling 
    .form-label, .selectize-input, .selectize-dropdown, .btn,
    .card, .card-header, .leaflet-control, label, p, h6, h5 {
      font-size: calc(1rem - var(--app-font-reduce));
    }
  "))
  ),
  
  
  titlePanel("Accidental epiphytes"),
  
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

##### Server #####
server <- function(input, output, session) {
  
  # Auto-reset filters on app start
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
  
  ##### Reactive Values #####
  selected_taxon    <- reactiveVal(NULL)
  clicked_obs       <- reactiveVal(NULL)
  clicked_obs_data  <- reactiveVal(NULL)
  matching_points   <- reactiveVal(NULL)
  
  # Reset filter choices from data (use obs_shiny, not an undefined 'obs')
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
  
  ##### Filtering #####
  filtered_obs <- reactive({
    data <- obs_shiny %>% ungroup()
    if (input$grade_filter  != "All") { data <- data %>% filter(quality_grade      == input$grade_filter)  }
    if (input$user_filter   != "All") { data <- data %>% filter(user_login         == input$user_filter)   }
    if (input$rank_filter   != "All") { data <- data %>% filter(taxon.rank         == input$rank_filter)   }
    if (input$family_filter != "All") { data <- data %>% filter(family_name  == input$family_filter) }
    if (input$taxon_filter  != "All") { data <- data %>% filter(taxon.name         == input$taxon_filter)  }
    data
  })
  
  ##### Output 1: Base map #####
  output$map <- renderLeaflet({
    req(input$grade_filter, input$user_filter, input$rank_filter, input$family_filter, input$taxon_filter)
    data <- filtered_obs()
    if (nrow(data) == 0) {
      return(leaflet() %>% addTiles() %>% addControl("No observations to display", position = "topright"))
    }
    
    coords <- st_coordinates(data)
    leaflet(data = data) %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Minimal") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      addProviderTiles(providers$OpenTopoMap, group = "Topography") %>%
      addCircleMarkers(
        lng = as.numeric(coords[,1]),
        lat = as.numeric(coords[,2]),
        layerId = ~as.character(id),
        radius = 5,
        color = ~col,
        fillOpacity = 0.6,
        group = "All observations"
      ) %>%
      addLegend("bottomright", colors = unique(obs_shiny$col), labels = unique(obs_shiny$grade_key), title = "Quality grade") %>%
      addLegend("bottomleft", colors = c("#90C987", "#F7F056"), labels = c("Taxon range", "Accidental epiphyte range"), title = "Range type") %>%
      addLayersControl(
        baseGroups = c("Minimal", "Topography", "Satellite"),
        overlayGroups = c("Ranges", "All observations", "Observations of the same taxon", "Selected observation"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  # (Optional) to later add an input$basemap, toggle here.
  # observe({
  #   if (!is.null(input$basemap)) {
  #     leafletProxy("map") %>%
  #       hideGroup(c("Minimal", "Topography", "Satellite")) %>%
  #       showGroup(input$basemap)
  #   }
  # })
  
  ##### Output 2: Observation Details Pane (with carousel + arrow color CSS) #####
  output$obs_details <- renderUI({
    obs.clicked <- clicked_obs_data()
    if (is.null(obs.clicked) || !is.data.frame(obs.clicked) || nrow(obs.clicked) == 0) {
      return(h6("No observation selected."))
    }
    
    # ---- Photos ----
    photo_urls <- unique(unlist(obs.clicked$photo_urls_list))
    photo_urls <- photo_urls[!is.na(photo_urls) & nzchar(photo_urls)]
    n_photos   <- length(photo_urls)
    carousel_id <- "photoCarousel"
    
    # Build carousel items
    photo_items <- tagList(lapply(seq_along(photo_urls), function(i) {
      tags$div(
        class = paste("carousel-item", if (i == 1) "active"),
        tags$img(
          src = photo_urls[i],
          class = "d-block w-100",
          style = "max-height:300px; object-fit:contain; background:#f8f9fa;"
        )
      )
    }))
    
    # --- Scoped CSS + Carousel UI (your snippet) ---
    ui_carousel <- if (n_photos > 0) {
      tagList(
        # Scoped CSS (Bootstrap 5.2+/5.3)
        singleton(tags$head(
          tags$style(HTML(sprintf("
        /* Scope the variables to this carousel ID 
        #%s {
          --bs-carousel-control-color: %s;      /* arrow color 
          --bs-carousel-control-opacity: .95;   /* visible by default 
          --bs-carousel-control-hover-opacity: 1;
          --bs-carousel-control-icon-width: 1.75rem; /* optional: slightly bigger icons 
        }

        /* Optional: subtle affordance on hover 
        #%s .carousel-control-prev:hover .carousel-control-prev-icon,
        #%s .carousel-control-next:hover .carousel-control-next-icon {
          filter: brightness(0.95) contrast(1.05);
        }

        /* ---- Fallback for older Bootstrap builds ----
           If the variables above don't take, uncomment the two rules below. 

        /*
        #%s .carousel-control-prev-icon {
          background-image: url(\"data:image/svg+xml,%%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16' fill='%s'%%3E%%3Cpath d='M11.354 1.646a.5.5 0 0 1 0 .708L6.707 7l4.647 4.646a.5.5 0 0 1-.708.708l-5-5a.5.5 0 0 1 0-.708l5-5a.5.5 0 0 1 .708 0z'%%3E%%3C/svg%%3E\");
        }
        #%s .carousel-control-next-icon {
          background-image: url(\"data:image/svg+xml,%%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16' fill='%s'%%3E%%3Cpath d='M4.646 1.646a.5.5 0 0 1 .708 0l5 5a.5.5 0 0 1 0 .708l-5 5a.5.5 0 0 1-.708-.708L9.293 7 4.646 2.354a.5.5 0 0 1 0-.708z'%%3E%%3C/svg%%3E\");
        }
        
      ",
                                  carousel_id, arrow_col,
                                  carousel_id, carousel_id,
                                  carousel_id, arrow_col_url,
                                  carousel_id, arrow_col_url
          )))
        )),
        
        # The carousel itself
        tags$div(
          id = carousel_id,
          class = "carousel slide",
          `data-bs-ride` = if (n_photos > 1) "carousel" else NULL,  # auto-start only if >1
          
          # Slides
          tags$div(class = "carousel-inner", photo_items),
          
          # Controls: render ONLY if more than one photo
          if (n_photos > 1) tags$button(
            class = "carousel-control-prev",
            type  = "button",
            `data-bs-target` = paste0("#", carousel_id),
            `data-bs-slide`  = "prev",
            tags$span(class = "carousel-control-prev-icon", `aria-hidden` = "true"),
            tags$span(class = "visually-hidden", "Previous")
          ),
          if (n_photos > 1) tags$button(
            class = "carousel-control-next",
            type  = "button",
            `data-bs-target` = paste0("#", carousel_id),
            `data-bs-slide`  = "next",
            tags$span(class = "carousel-control-next-icon", `aria-hidden` = "true"),
            tags$span(class = "visually-hidden", "Next")
          )
        )
      )
    } else {
      NULL
    }
    
    
    # ---- textbox and carousel UI ----
    tagList(
      h6(
        tags$a(
          href   = unique(obs.clicked$uri)[1],
          "View on iNaturalist",
          target = "_blank"
        )
      ),
      
      # Carousel (only if there is at least one photo)
      ui_carousel,
      
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
      h6(
        strong(tagList("Notes from observer (", tags$em(unique(obs.clicked$user_login)[1]), "): ")),
        unique(obs.clicked$description)[1]
      ),
      
      # Initialize the carousel in BS5 (only if it can actually slide)
      if (n_photos > 1) {
        tags$script(HTML(sprintf("
          (function(){
            var el = document.getElementById('%s');
            if (el && window.bootstrap && bootstrap.Carousel) {
              try { new bootstrap.Carousel(el, { interval: 3000, pause: 'hover', wrap: true }); } catch(e){}
            }
          })();
        ", carousel_id)))
      }
    )
  })
  
  ##### Click handling & map updates #####
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
    matching_points(mp)  # <-- update the reactiveVal (this was the missing link)
    
    # Clear & show
    leafletProxy("map") %>%
      clearGroup("Observations of the same taxon") %>%
      clearGroup("All observations") %>%
      clearGroup("Ranges") %>%
      clearGroup("Selected observation") %>%
      showGroup("Ranges")
    
    ## Range polygons ----
    {
      # Convert obs_range to valid MULTIPOLYGON
      obs_range_geom  <- clicked$obs_range_wkt  %>%
        st_as_sfc(., crs = 4326)
      
      inat_range_geom <- clicked$inat_range_wkt %>%
        st_as_sfc(., crs = 4326)
      
      if (length(obs_range_geom)  > 0 && !st_is_empty(obs_range_geom)) {
        leafletProxy("map") %>%
          addPolygons(data = obs_range_geom, color = "#F7F056", weight = 2, fillOpacity = 0.5, group = "Ranges")
      }
      if (length(inat_range_geom) > 0 && !st_is_empty(inat_range_geom)) {
        leafletProxy("map") %>%
          addPolygons(data = inat_range_geom, color = "#90C987", weight = 2, fillOpacity = 0.3, group = "Ranges")
      }
    }
    
    #### Output 4: markers ----
    # All observations (dimmed)
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
    
    # Observations of the same taxon (highlighted)
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
    
    # Selected observation
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
  })
  
  ##### Output 5: sum observations and individuals #####
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
    
    
    # Base message (your existing grammar)
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
    # Treat additional labels as "range missing" if you use them
    is_range_missing <- isTRUE(obs.clicked$status %in% c("range_missing"))
    
    # (Optional) be tolerant to tiny floating errors around zero
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
}

##### Run the Shiny App #####
shinyApp(ui, server)





