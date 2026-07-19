# Original exploratory plotting notebook; contains optional packages and alternatives.
################################
#––– 0. Tweak as needed:
csv_folder <- "/Users/dr.amirsaeed/Downloads/NTNU Research/CO2 Research/CO2 storage BO with real data/Sleipner 2019/Sleipner plume boundaries/data/CSV_exports"
# pattern to catch csv files
file_pattern <- "_parsed\\.csv$"        

#––– 1. Load libraries
library(sf)       # for spatial data types & plotting via ggplot2
library(ggplot2)  # for plotting

#––– 2. List and read all parsed CSVs into a list
parsed_files <- list.files(csv_folder, pattern = file_pattern, full.names = TRUE)
if (length(parsed_files)==0) stop("No parsed CSVs found in ", csv_folder)

# function to turn one CSV into an sf MULTIPOLYGON
csv_to_sf <- function(fpath) {
  df <- read.csv(fpath, stringsAsFactors = FALSE)
  
  # split by segment ID
  segs <- split(df, df$SEG_ID)
  
  # for each segment, build one coordinate ring (ensure it's closed)
  rings <- lapply(segs, function(d) {
    coords <- as.matrix(d[, c("X_EASTING","Y_NORTHING")])
    if (!all(coords[1,] == coords[nrow(coords),])) {
      coords <- rbind(coords, coords[1,])
    }
    coords
  })
  
  # each layer may have multiple rings → one MULTIPOLYGON
  # we wrap the list of rings in one extra list:
  mp  <- st_multipolygon(list(rings))
  # build a one-row sf object
  st_sf(
    Layer    = unique(df$Layer),
    geometry = st_sfc(mp),
    crs      = 25832   # or whatever CRS your easting/northing use
  )
}

#––– 3. Apply to all files, then combine into one sf
sf_list <- lapply(parsed_files, csv_to_sf)
plume_sf <- do.call(rbind, sf_list)

#––– 4. Plot with ggplot2
ggplot(plume_sf) +
  geom_sf(aes(fill = Layer), alpha = 0.5, color = "black") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Sleipner CO2 Plume Boundaries by Layer (2019)",
    fill  = "Layer"
  )

library(sf)
library(ggplot2)

ggplot(plume_sf) +
  geom_sf(fill = NA, color = "steelblue") +
  facet_wrap(~ Layer, ncol = 3) +
  coord_sf() +
  theme_minimal() +
  labs(title = "Sleipner CO2 Plume by Layer (separate panels)")


library(tmap)

# static map
tm_shape(plume_sf) +
  tm_fill("Layer", alpha = 0.5, palette = "Set2") +
  tm_borders(col = "black") +
  tm_facets(by = "Layer") +
  tm_layout(title = "Sleipner Plume by Layer (tmap facets)")

library(sp)
library(raster)
library(gganimate)

# Create base plot
p <- ggplot(plume_sf) +
  geom_sf(aes(fill = Layer), alpha = 0.7) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "CO2 Plume Evolution: {closest_state}", 
       subtitle = "Sleipner Storage Site")

# Animate through layers
anim <- p +
  transition_states(Layer, transition_length = 1, state_length = 2) +
  enter_fade() +
  exit_shrink()

# Save animation
animate(anim, fps = 30, duration = 20, renderer = gifski_renderer())


library(leaflet)
library(viridisLite)
library(viridis)

# Convert to geographic CRS
plume_wgs84 <- st_transform(plume_sf, 4326)

# Create color palette
pal <- colorFactor(viridis_pal(option = "C")(9), domain = plume_wgs84$Layer)

# Create interactive map
leaflet(plume_wgs84) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons(
    fillColor = ~pal(Layer),
    weight = 1,
    opacity = 1,
    color = "white",
    fillOpacity = 0.7,
    highlight = highlightOptions(weight = 3, color = "yellow"),
    popup = ~Layer,
    group = "Layers"
  ) %>%
  addLayersControl(
    overlayGroups = "Layers",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(pal = pal, values = ~Layer, title = "Plume Layers")

ggplot(plume_sf) +
  geom_sf(aes(fill = Layer), show.legend = FALSE) +
  facet_wrap(~Layer, ncol = 3) +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Sleipner CO2 Plume Layers")


library(sf)
library(lwgeom)
library(dplyr)
library(purrr)
library(leaflet)
library(leaflet.extras)
library(viridis)
library(osmdata)   # for opq() et al.

# — Repair & reproject plume —
plume_wgs84 <- plume_sf %>%
  st_make_valid() %>%
  st_transform(4326)

# — Define simple block-rectangles —
bbox_156 <- matrix(c(
  1.80, 58.30,
  2.00, 58.30,
  2.00, 58.40,
  1.80, 58.40,
  1.80, 58.30
), ncol = 2, byrow = TRUE)

bbox_159 <- matrix(c(
  1.90, 58.35,
  2.10, 58.35,
  2.10, 58.45,
  1.90, 58.45,
  1.90, 58.35
), ncol = 2, byrow = TRUE)

blocks_sf <- st_sf(
  Block    = c("15/6", "15/9"),
  geometry = st_sfc(
    st_polygon(list(bbox_156)),
    st_polygon(list(bbox_159))
  ),
  crs      = 4326
)

# — Color palettes —
pal_plume  <- colorFactor(viridis(9, option = "C"), domain = plume_wgs84$Layer)
pal_blocks <- colorFactor(c("blue","red"), domain = blocks_sf$Block)

# — Centroids (no warning) —
centroids_df <- st_centroid(plume_wgs84, byid = TRUE) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(
    Layer      = plume_wgs84$Layer,
    col        = pal_plume(Layer),
    label_html = sprintf(
      "<span style='color:%s; font-weight:bold;'>%s</span>",
      col, Layer
    )
  )

# — Fetch nearby places from OSM —
bb <- st_bbox(blocks_sf)
bbox_expanded <- c(
  xmin = bb["xmin"] - 0.5,
  ymin = bb["ymin"] - 0.5,
  xmax = bb["xmax"] + 0.5,
  ymax = bb["ymax"] + 0.5
)

osm_places <- opq(bbox = bbox_expanded) %>%
  add_osm_feature(key = "place",
                  value = c("city","town","village","hamlet")) %>%
  osmdata_sf()

places_pts <- data.frame()
if (!is.null(osm_places$osm_points)) {
  places_pts <- osm_places$osm_points %>%
    mutate(
      name  = if ("name"  %in% names(.)) as.character(name)  else NA_character_,
      place = if ("place" %in% names(.)) as.character(place) else NA_character_
    ) %>%
    filter(!is.na(name)) %>%
    mutate(
      lng = st_coordinates(geometry)[,1],
      lat = st_coordinates(geometry)[,2]
    ) %>%
    as_tibble() %>%
    dplyr::select(name, place, lng, lat)
}

# — Prepare map center as a list —
sleipner_center <- list(lng = 1.917, lat = 58.367)

# — Build Leaflet map —
map <- leaflet() %>%
  addTiles() %>%
  addPolygons(
    data        = plume_wgs84,
    fillColor   = ~pal_plume(Layer),
    fillOpacity = 0.7,
    color       = "#444444",
    weight      = 1
  ) %>%
  addLabelOnlyMarkers(
    data         = centroids_df,
    lng          = ~X, lat = ~Y,
    label        = ~lapply(label_html, htmltools::HTML),
    labelOptions = labelOptions(noHide = TRUE,
                                direction = "center",
                                textOnly  = TRUE)
  ) %>%
  addPolygons(
    data        = blocks_sf,
    fillColor   = ~pal_blocks(Block),
    fillOpacity = 0.2,
    color       = ~pal_blocks(Block),
    weight      = 2
  ) %>%
  { if (nrow(places_pts) > 0) {
    addMarkers(
      .,
      data  = places_pts,
      lng   = ~lng,
      lat   = ~lat,
      label = ~name,
      popup = ~paste0("<b>", name, "</b><br>Type: ", place)
    )
  } else . } %>%
  addLegend(
    position = "bottomright",
    pal      = pal_plume,
    values   = plume_wgs84$Layer,
    title    = "CO₂ Plume Layers"
  ) %>%
  setView(
    lng  = sleipner_center$lng,
    lat  = sleipner_center$lat,
    zoom = 10
  )

map

