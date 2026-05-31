# Map of field sites — Maxillariinae epiphytes, NW Ecuador
# ============================================================

# install.packages(c("sf", "ggplot2", "dplyr",
#                    "rnaturalearth", "rnaturalearthdata", "ggspatial"))

library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# ── 1. Load GeoJSON files ─────────────────────────────────────────────────────
# The files have slightly different property columns ('stroke' is absent in some).
# dplyr::bind_rows() silently drops the sf class when columns mismatch,
# so we keep only the columns we need before combining.

# Set this to the folder that contains your GeoJSON files
data_dir <- paste0(BASE_DIR, "/geojson_to_csv/raw")

geojson_files <- c(
  Maquipucuna   = "Maquipucuna.geojson",
  Mashpi        = "Mashpi.geojson",
  MindoTarabita = "MindoTarabita.geojson",   # merged with TarabitaMindo below
  TarabitaMindo = "TarabitaMindo.geojson",
  MiradorMindo  = "MiradorMindo.geojson",
  Yanayacu      = "Yanayacu.geojson"
)

load_geojson <- function(filename, site_name, dir = data_dir) {
  sf::st_read(file.path(dir, filename), quiet = TRUE) %>%
    # Keep only what we need — avoids column-mismatch when rbinding
    transmute(site = site_name,
              description = if ("description" %in% names(.)) description else NA_character_)
}

site_list <- mapply(load_geojson, geojson_files, names(geojson_files),
                    SIMPLIFY = FALSE)

# sf-safe row-bind: do.call(rbind) preserves the sf class
all_features <- do.call(rbind, site_list)

# Merge TarabitaMindo into MindoTarabita
all_features <- all_features %>%
  mutate(site = if_else(site == "TarabitaMindo", "MindoTarabita", site))

# ── 2. Split by geometry type ─────────────────────────────────────────────────
obs_points <- all_features %>%
  filter(sf::st_geometry_type(geometry) == "POINT")

transects <- all_features %>%
  filter(sf::st_geometry_type(geometry) == "LINESTRING")

# ── 3. Site label positions ───────────────────────────────────────────────────
# Labels sit at transect centroids where transects exist; fall back to
# obs_points centroid for Yanayacu (no transect recorded in the GeoJSON).
transect_labels <- transects %>%
  group_by(site) %>%
  summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
  mutate(geometry = sf::st_centroid(geometry))

fallback_labels <- obs_points %>%
  filter(!site %in% transect_labels$site) %>%
  group_by(site) %>%
  summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
  mutate(geometry = sf::st_centroid(geometry))

site_labels <- rbind(transect_labels, fallback_labels)

# Per-site label offsets (degrees) — keeps long names off the dots
# MiradorMindo / MindoTarabita are <0.07° apart so need opposite directions
label_nudges <- data.frame(
  site    = c("Mashpi", "Maquipucuna", "MiradorMindo", "MindoTarabita", "Yanayacu"),
  nudge_x = c( 0.00,     0.08,         -0.14,           0.15,            0.10),
  nudge_y = c( -0.04,     0.07,          0.07,           -0.05,           0.06)
)
label_pos <- site_labels |>
  dplyr::mutate(X = sf::st_coordinates(geometry)[, 1],
                Y = sf::st_coordinates(geometry)[, 2]) |>
  sf::st_drop_geometry() |>
  dplyr::left_join(label_nudges, by = "site") |>
  dplyr::mutate(X = X + nudge_x, Y = Y + nudge_y)

# ── 4. Basemap ────────────────────────────────────────────────────────────────
ecuador    <- ne_states(country = "Ecuador", returnclass = "sf")
nw_provs <- c("Pichincha", "Esmeraldas", "Imbabura",
              "Santo Domingo de los Tsáchilas", "Cotopaxi",
              "Napo", "Sucumbios", "Tungurahua")
map_area <- ecuador %>% filter(name %in% nw_provs)

bbox <- sf::st_bbox(obs_points)
xlim <- c(bbox["xmin"] - 0.30, bbox["xmax"] + 0.55)  # extra east for Yanayacu
ylim <- c(bbox["ymin"] - 0.40, bbox["ymax"] + 0.30)

# ── 5. Colours ────────────────────────────────────────────────────────────────
site_colours <- setNames(
  scico::scico(5, palette = "lipari", begin = 0.10, end = 0.88),
  c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu")
)

# ── 6. Plot ───────────────────────────────────────────────────────────────────
p <- ggplot() +
  
  geom_sf(data = map_area,
          fill = "#f5f0e8", colour = "grey55", linewidth = 0.35) +
  
  geom_sf(data = obs_points,
          aes(colour = site),
          size = 1.8, alpha = 0.55, shape = 16) +
  
  geom_sf(data = transects,
          aes(colour = site),
          linewidth = 1.3, alpha = 0.85) +
  
  geom_text(data = label_pos,
            aes(x = X, y = Y, label = site, colour = site),
            size = 3, fontface = "bold",
            show.legend = FALSE) +
  
  scale_colour_manual(values = site_colours, name = "Site") +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.2, "cm"), width = unit(1.2, "cm")
  ) +
  
  labs(
    title    = "Field sites",
    subtitle = "NW Ecuador · Chocó Andino ",
    x = "Longitude", y = "Latitude"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey40")
  )

print(p)

ggsave("site_map.png", plot = p,
       width = 10, height = 8, dpi = 300, bg = "white")

message("Saved: site_map.png")