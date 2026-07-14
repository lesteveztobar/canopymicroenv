# plot_functions.R
# Reusable plotting functions for the whole canopymicroenv pipeline. Each
# function loads its own inputs from data/processed (or geojson_to_csv) and
# saves its output(s) into OUTPUT_DIR — nothing here plots-and-forgets.
# If an input file doesn't exist yet, the function skips with a message
# instead of erroring, so plot_all.R can be re-run at any point in the
# pipeline and it just draws whatever is available so far.
#
# Refactored from the standalone plotmap.R / plot_temperatures.R /
# plot_bestfit_3d.R, plus new export wrappers around the plot_abundance() /
# plot_3d_abundance() / plot_experiment() functions already defined in
# get_colonization.R.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(abind)
library(scatterplot3d)

source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

# ── Field site map ─────────────────────────────────────────────────────────────
# NW Ecuador Maxillariinae field sites, from GeoJSON transects/points.
# sf/dplyr/rnaturalearth(data)/ggspatial/scico are loaded here rather than at
# file level: this is the only function in the file that needs them (all are
# geospatial/mapping-specific, several with heavy system-library dependencies
# -- GDAL/PROJ/GEOS/UDUNITS -- that every other function in this file has no
# reason to require just to be sourced).
plot_site_map <- function(geojson_dir = file.path(BASE_DIR, "geojson_to_csv", "raw"),
                           out_dir = OUTPUT_DIR) {
  library(sf)
  library(dplyr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
  library(scico)

  geojson_files <- c(
    Maquipucuna   = "Maquipucuna.geojson",
    Mashpi        = "Mashpi.geojson",
    MindoTarabita = "MindoTarabita.geojson",   # merged with TarabitaMindo below
    TarabitaMindo = "TarabitaMindo.geojson",
    MiradorMindo  = "MiradorMindo.geojson",
    Yanayacu      = "Yanayacu.geojson"
  )

  load_geojson <- function(filename, site_name, dir = geojson_dir) {
    sf::st_read(file.path(dir, filename), quiet = TRUE) %>%
      transmute(site = site_name,
                description = if ("description" %in% names(.)) description else NA_character_)
  }

  site_list    <- mapply(load_geojson, geojson_files, names(geojson_files), SIMPLIFY = FALSE)
  all_features <- do.call(rbind, site_list)
  all_features <- all_features %>%
    mutate(site = if_else(site == "TarabitaMindo", "MindoTarabita", site))

  obs_points <- all_features %>% filter(sf::st_geometry_type(geometry) == "POINT")
  transects  <- all_features %>% filter(sf::st_geometry_type(geometry) == "LINESTRING")

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

  # Per-site label offsets (degrees) — keeps long names off the dots.
  label_nudges <- data.frame(
    site    = c("Mashpi", "Maquipucuna", "MiradorMindo", "MindoTarabita", "Yanayacu"),
    nudge_x = c( 0.00,     0.08,         -0.14,           0.15,            0.10),
    nudge_y = c(-0.04,     0.07,          0.07,           -0.05,           0.06)
  )
  label_pos <- site_labels |>
    dplyr::mutate(X = sf::st_coordinates(geometry)[, 1],
                  Y = sf::st_coordinates(geometry)[, 2]) |>
    sf::st_drop_geometry() |>
    dplyr::left_join(label_nudges, by = "site") |>
    dplyr::mutate(X = X + nudge_x, Y = Y + nudge_y)

  ecuador  <- ne_states(country = "Ecuador", returnclass = "sf")
  nw_provs <- c("Pichincha", "Esmeraldas", "Imbabura",
                "Santo Domingo de los Tsáchilas", "Cotopaxi",
                "Napo", "Sucumbios", "Tungurahua")
  map_area <- ecuador %>% filter(name %in% nw_provs)

  bbox <- sf::st_bbox(obs_points)
  xlim <- c(bbox["xmin"] - 0.30, bbox["xmax"] + 0.55)  # extra east for Yanayacu
  ylim <- c(bbox["ymin"] - 0.40, bbox["ymax"] + 0.30)

  site_colours <- setNames(
    scico::scico(5, palette = "lipari", begin = 0.10, end = 0.88),
    c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu")
  )

  p <- ggplot() +
    geom_sf(data = map_area, fill = "#f5f0e8", colour = "grey55", linewidth = 0.35) +
    geom_sf(data = obs_points, aes(colour = site), size = 1.8, alpha = 0.55, shape = 16) +
    geom_sf(data = transects, aes(colour = site), linewidth = 1.3, alpha = 0.85) +
    geom_text(data = label_pos, aes(x = X, y = Y, label = site, colour = site),
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = site_colours, name = "Site") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr", style = north_arrow_fancy_orienteering(),
      height = unit(1.2, "cm"), width = unit(1.2, "cm")
    ) +
    labs(title = "Field sites", subtitle = "NW Ecuador · Chocó Andino ",
         x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey40")
    )

  out_path <- file.path(out_dir, "site_map.png")
  ggsave(out_path, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Microclimate temperature profile ────────────────────────────────────────────
# 3D temperature volume (interactive, plotly) plus a side-view heatmap +
# boxplot with Kruskal-Wallis / Dunn post-hoc height comparison, for one
# site's microenv RDS.
.build_temp_volume <- function(env_heights, day_type = "tmax", var = "Tz") {
  layers <- lapply(env_heights, function(h) {
    r <- h[[day_type]][[var]]
    if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
  })
  abind::abind(layers, along = 3)
}

plot_temperature_profile <- function(site_name, out_dir = OUTPUT_DIR,
                                      processed_dir = PROCESSED_DIR) {
  microenv_path <- file.path(processed_dir, sprintf("microenv_%s.rds", site_name))
  if (!file.exists(microenv_path)) {
    message("Skipping ", site_name, " temperature profile -- no ", microenv_path)
    return(invisible(NULL))
  }
  env       <- readRDS(microenv_path)
  heights_m <- microenv_heights(env)
  # load_height()/microenv_heights() (get_colonization.R) handle both the new
  # manifest-only format (.heights + .height_dir, per-height RDS on scratch)
  # and the old format where each height is embedded as its own list element.
  env_heights <- lapply(heights_m, function(h) load_height(env, h))
  vol_tmax    <- .build_temp_volume(env_heights, "tmax", "Tz")

  nrow_r  <- dim(vol_tmax)[1]
  ncol_r  <- dim(vol_tmax)[2]
  nheight <- dim(vol_tmax)[3]

  if (!is.null(env$.spatial)) {
    e        <- env$.spatial$ext
    x_coords <- seq(e$xmin, e$xmax, length.out = ncol_r)
    y_coords <- seq(e$ymax, e$ymin, length.out = nrow_r)  # north→south
    x_lab <- "Longitude"; y_lab <- "Latitude"
  } else {
    x_coords <- seq_len(ncol_r)
    y_coords <- seq_len(nrow_r)
    x_lab <- "col (W→E)"; y_lab <- "row (S→N)"
  }

  # ── 3D volume (plotly) ────────────────────────────────────────────────────
  grid <- expand.grid(xi = x_coords, yi = y_coords, zi = seq_len(nheight))
  grid$z  <- heights_m[grid$zi]
  grid$Tz <- vol_tmax[cbind(nrow_r + 1 - grid$yi, grid$xi, grid$zi)]
  grid    <- grid[!is.na(grid$Tz), ]

  fig <- plot_ly(
    data = grid, x = ~xi, y = ~yi, z = ~z, value = ~Tz, type = "volume",
    isomin = quantile(grid$Tz, 0.05), isomax = quantile(grid$Tz, 0.95),
    opacity = 0.15, surface = list(count = 10),
    colorscale = list(
      c(0,    "#313695"), c(0.25, "#74add1"), c(0.5, "#ffffbf"),
      c(0.75, "#f46d43"), c(1,    "#a50026")
    ),
    colorbar = list(title = "T (°C)")
  ) |>
    layout(
      title = sprintf("Air temperature 3D volume — %s (warmest representative day)", site_name),
      scene = list(
        xaxis = list(title = x_lab), yaxis = list(title = y_lab),
        zaxis = list(title = "Height (m)")
      )
    )

  volume_path <- file.path(out_dir, sprintf("temp_volume_%s.html", site_name))
  # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
  # writes a small "<name>_files/" dependency folder alongside the HTML
  # instead -- keep the two together when copying/viewing elsewhere.
  htmlwidgets::saveWidget(fig, volume_path, selfcontained = FALSE)
  message("Saved: ", volume_path)

  # ── Side-view heatmap + per-height boxplot ────────────────────────────────
  side_df <- do.call(rbind, lapply(seq_along(env_heights), function(i) {
    r <- env_heights[[i]]$tmax$Tz
    m <- if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
    col_means <- colMeans(m, na.rm = TRUE)
    data.frame(col = seq_along(col_means), height = heights_m[i], Tz = col_means)
  }))
  side_df <- side_df[!is.na(side_df$Tz), ]

  # All pixel values per height — more power for the test than col means alone
  all_px <- do.call(rbind, lapply(seq_along(env_heights), function(i) {
    r <- env_heights[[i]]$tmax$Tz
    m <- if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
    data.frame(height = factor(heights_m[i]), Tz = as.vector(m))
  }))
  all_px <- all_px[!is.na(all_px$Tz), ]

  # Kruskal-Wallis: are any height levels different?
  kw <- kruskal.test(Tz ~ height, data = all_px)
  kw_label <- sprintf("Kruskal-Wallis: χ²(%.0f) = %.2f, p %s",
                      kw$parameter, kw$statistic,
                      ifelse(kw$p.value < 0.001, "< 0.001",
                             sprintf("= %.3f", kw$p.value)))
  message(kw_label)

  # Dunn post-hoc pairwise comparisons (Holm correction)
  if (!requireNamespace("dunn.test", quietly = TRUE)) install.packages("dunn.test")
  dunn_res <- dunn.test::dunn.test(all_px$Tz, all_px$height,
                                    method = "holm", kw = FALSE, label = FALSE)
  dunn_df <- data.frame(
    comparison = dunn_res$comparisons,
    p_adj      = dunn_res$P.adjusted,
    sig        = ifelse(dunn_res$P.adjusted < 0.001, "***",
                 ifelse(dunn_res$P.adjusted < 0.01,  "**",
                 ifelse(dunn_res$P.adjusted < 0.05,  "*", "ns")))
  )
  message("Pairwise Dunn tests (Holm-adjusted):")
  print(dunn_df[order(dunn_df$p_adj), ], row.names = FALSE)

  # Compact letter display for the boxplot annotation
  if (!requireNamespace("multcompView", quietly = TRUE)) install.packages("multcompView")
  p_mat  <- setNames(dunn_res$P.adjusted, dunn_res$comparisons)
  cld    <- multcompView::multcompLetters(p_mat, threshold = 0.05)$Letters
  cld_df <- data.frame(height = as.numeric(names(cld)), letter = cld)

  p_heat <- ggplot(side_df, aes(x = col, y = height, fill = Tz)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"),
      name = "T (°C)", na.value = "grey90"
    ) +
    labs(x = "Horizontal position (col, W→E)", y = "Height (m)",
         title = sprintf("Vertical cross-section — %s (warmest representative day)", site_name),
         subtitle = kw_label) +
    theme_minimal(base_size = 11)

  p_box <- ggplot(all_px, aes(x = Tz, y = height, group = height)) +
    geom_boxplot(aes(fill = after_stat(middle)), width = 0.06, outlier.size = 0.5) +
    geom_text(data = cld_df,
              aes(x = max(all_px$Tz, na.rm = TRUE), y = height, label = letter),
              hjust = -0.2, size = 3.5, inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"), guide = "none"
    ) +
    labs(x = "T (°C)", y = NULL,
         caption = "Letters = Dunn post-hoc (Holm); shared letter → no significant difference") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank())

  p_combined  <- p_heat + p_box + plot_layout(widths = c(3, 1))
  profile_path <- file.path(out_dir, sprintf("temp_profile_%s.png", site_name))
  ggsave(profile_path, plot = p_combined, width = 11, height = 5, dpi = 300, bg = "white")
  message("Saved: ", profile_path)

  invisible(list(volume = fig, profile = p_combined))
}

# ── Best-fit vs. realistic 3D comparison (simple standalone model) ─────────────
plot_bestfit_3d_comparison <- function(out_dir = OUTPUT_DIR) {
  if (!exists("best_run", inherits = TRUE)) {
    source("scripts/simple_colonization.R")
    best_params <- params
    best_params$establishment_prob <- 0.10
    best_params$repro_rate         <- 500
    best_params$survival_A         <- 0.95
    best_params$maxDisp            <- 4
    best_params$n_founders         <- 100
    best_params$carCap             <- 40
    best_params$timesteps          <- 30
    message("Re-running best fit...")
    best_run <- run_simple_colonization(best_params, seed = 42)
  }

  p      <- best_run$params
  T_last <- p$timesteps

  abundA_last <- best_run$abundA[, , , T_last]
  zone_arr    <- best_run$zone

  idx <- which(abundA_last > 0, arr.ind = TRUE)
  df  <- data.frame(
    x = idx[, 1], y = idx[, 2], z = idx[, 3],
    n = abundA_last[idx], zone = zone_arr[idx]
  )
  df$x_m <- (df$x - 0.5) * p$resolution
  df$y_m <- (df$y - 0.5) * p$resolution
  df$z_m <- (df$z - 0.5) * (p$max_height / p$zDim)

  zone_cols <- c(
    "1" = "#8B4513", "2" = "#A0522D", "3" = "#6B8E23",
    "4" = "#228B22", "5" = "#32CD32"
  )
  df$col     <- zone_cols[as.character(df$zone)]
  df$pt_size <- pmin(df$n / max(df$n) * 6 + 1, 7)

  # Realistic-params run for comparison (above extinction threshold, not saturating)
  source("scripts/simple_colonization.R")
  base_p <- params
  base_p$establishment_prob <- 0.05
  base_p$repro_rate         <- 150
  base_p$carCap             <- 8
  base_p$n_founders         <- 30
  base_p$timesteps          <- 30
  message("Running baseline simulation...")
  invisible(capture.output(
    base_out <- run_simple_colonization(base_p, seed = 7),
    type = "output"
  ))

  make_df <- function(run_out, timestep) {
    abA  <- run_out$abundA[, , , timestep]
    zArr <- run_out$zone
    idx  <- which(abA > 0, arr.ind = TRUE)
    if (nrow(idx) == 0) return(NULL)
    data.frame(
      x_m = (idx[, 1] - 0.5) * run_out$params$resolution,
      y_m = (idx[, 2] - 0.5) * run_out$params$resolution,
      z_m = (idx[, 3] - 0.5) * (run_out$params$max_height / run_out$params$zDim),
      n = abA[idx], zone = zArr[idx]
    )
  }

  df_base <- make_df(base_out, base_p$timesteps)
  df_best <- make_df(best_run, T_last)

  set.seed(42)
  df_best_plot <- df_best[sample(nrow(df_best), min(nrow(df_best), 5000)), ]
  df_base_plot <- df_base   # usually sparse enough to keep all

  add_cols <- function(df) { df$col <- zone_cols[as.character(df$zone)]; df }
  df_base_plot <- add_cols(df_base_plot)
  df_best_plot <- add_cols(df_best_plot)

  out_path <- file.path(out_dir, "bestfit_3d.png")
  png(out_path, width = 3200, height = 1400, res = 180, type = "cairo")
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(10, 10, 3))

  plot_panel <- function(df_p, title_str, cex_sym = 0.55) {
    par(mar = c(2, 2, 3, 1))
    scatterplot3d(
      x = df_p$x_m, y = df_p$y_m, z = df_p$z_m,
      color = df_p$col, pch = 16, cex.symbols = cex_sym,
      xlab = "East–West (m)", ylab = "South–North (m)", zlab = "Height (m)",
      main = title_str, angle = 35, scale.y = 0.6, grid = TRUE, box = FALSE,
      col.axis = "grey40", col.grid = "grey88", col.lab = "grey20",
      cex.axis = 0.8, cex.lab = 0.9
    )
  }

  plot_panel(df_base_plot,
    sprintf("Realistic params  (p_e=0.05, λ=150, carCap=8)\nyear %d — %d adults",
            base_p$timesteps, sum(df_base$n)))

  plot_panel(df_best_plot,
    sprintf("Best-fit params  (p_e=0.10, λ=500, carCap=40)\nyear %d — %d adults  [5k sample]",
            T_last, sum(df_best$n)))

  par(mar = c(2, 0, 3, 1))
  plot.new()
  legend(
    "center",
    legend = c("Zone 1 – trunk base", "Zone 2 – lower trunk",
               "Zone 3 – upper trunk", "Zone 4 – inner crown",
               "Zone 5 – outer crown"),
    col = unname(zone_cols), pch = 16, pt.cex = 1.5, bty = "n", cex = 0.95,
    title = "Johansson zone", title.col = "grey20", title.font = 2
  )

  message("Saved: ", out_path)
  invisible(out_path)
}

# ── Colonization sensitivity-experiment sweeps ──────────────────────────────────
# One RDS per site x experiment tag, produced by run_colonization_onesite.R
# via batch_exp.sh. This mapping must stay in sync with the EXP/PARAMS
# pairing in scripts/batch_exp.sh.
EXP_PARAM_MAP <- c(
  pollination_success       = "p_poll",
  adult_survival_intercept  = "beta0A",
  germination_probability   = "p_germ",
  reproduction_cost         = "cost_repro",
  climate_sensitivity_rh    = "beta_rh",
  precipitation_sensitivity = "beta_precip",
  niche_tolerance           = "niche_pad",
  founder_number            = "n_founders"
)

plot_colonization_experiment <- function(site_name, exp_tag, out_dir = OUTPUT_DIR,
                                          processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (!is.data.frame(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- single-run result, use plot_default_colonization_run() instead")
    return(invisible(NULL))
  }
  if (!"param_value" %in% names(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- factorial result (no param_value column), use plot_factorial_experiment() instead")
    return(invisible(NULL))
  }

  param_name <- EXP_PARAM_MAP[[exp_tag]] %||% exp_tag
  p <- plot_experiment(result, param_name,
                       title = sprintf("Effect of %s — %s", param_name, site_name))

  out_path <- file.path(out_dir, sprintf("experiment_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 12, height = 5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Factorial experiment (e.g. n_founders x p_poll x p_germ x p_s1) ────────────
# Extinction-rate heatmap across two swept parameters, faceted by however many
# more are present — generalizes to any factorial from run_factorial_experiment(),
# not just a specific 3- or 4-parameter design. Swept columns are detected as
# whatever's left after removing the fixed output columns, with n_founders
# (if present) placed on the primary x-axis since it's usually the parameter
# of most direct interest. Mirrors the simple model's Fig. 5 (report.pdf
# sec. 3.3), generalized past a single facet variable via facet_grid.
plot_factorial_experiment <- function(site_name, exp_tag = "reproduction_factorial",
                                      out_dir = OUTPUT_DIR, processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (!is.data.frame(result) || "param_value" %in% names(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- not a factorial result, use plot_colonization_experiment() instead")
    return(invisible(NULL))
  }
  fixed_cols <- c("rep", "t", "totalS", "totalJ", "totalA", "total", "extinct")
  swept <- setdiff(names(result), fixed_cols)
  if ("n_founders" %in% swept) swept <- c("n_founders", setdiff(swept, "n_founders"))
  if (length(swept) < 2) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- fewer than 2 swept columns found (", paste(swept, collapse = ", "), ")")
    return(invisible(NULL))
  }

  t_max <- max(result$t)
  final <- result[result$t == t_max, ]
  form  <- as.formula(paste("cbind(extinct, total) ~", paste(swept, collapse = " + ")))
  combo_summary <- aggregate(form, data = final, FUN = mean)

  x_var <- swept[1]; y_var <- swept[2]
  facet_vars <- swept[-(1:2)]

  p <- ggplot(combo_summary,
             aes(x = factor(.data[[x_var]]), y = factor(.data[[y_var]]), fill = extinct)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#08519c", "#f1a340", "#a50026"),
                         limits = c(0, 1), name = "Extinction\nrate") +
    labs(x = x_var, y = y_var,
         title = sprintf("Factorial — %s", site_name),
         subtitle = if (length(facet_vars) > 0)
           sprintf("Extinction rate at year %d, faceted by %s",
                   t_max, paste(facet_vars, collapse = " x "))
         else
           sprintf("Extinction rate at year %d", t_max)) +
    theme_minimal(base_size = 11)

  if (length(facet_vars) == 1) {
    p <- p + facet_wrap(as.formula(paste("~", facet_vars[1])), labeller = label_both)
  } else if (length(facet_vars) >= 2) {
    p <- p + facet_grid(as.formula(paste(facet_vars[1], "~", facet_vars[2])), labeller = label_both)
    if (length(facet_vars) > 2)
      message("Note: faceting by ", facet_vars[1], " and ", facet_vars[2],
              " only; ", paste(facet_vars[-(1:2)], collapse = ", "),
              " collapsed via aggregation.")
  }

  out_path <- file.path(out_dir, sprintf("factorial_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Single default colonization run (no swept parameter) ───────────────────────
plot_default_colonization_run <- function(site_name, exp_tag = "default",
                                           out_dir = OUTPUT_DIR,
                                           processed_dir = PROCESSED_DIR,
                                           animate = TRUE) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (is.data.frame(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- sweep result, use plot_colonization_experiment() instead")
    return(invisible(NULL))
  }
  # run_replicated() output: list(runs = <one runcolonization() per
  # replicate>, summary = <tidy data frame>). Falls back to treating `result`
  # itself as a single run for any older RDS saved before that change.
  runs <- if (is.list(result) && !is.null(result$runs)) result$runs else list(result)
  runs <- Filter(Negate(is.null), runs)
  if (length(runs) == 0) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no successful replicates")
    return(invisible(NULL))
  }
  message(length(runs), " replicate(s) found for ", site_name, " / ", exp_tag)

  saved <- lapply(seq_along(runs), function(i) {
    suffix <- if (length(runs) > 1) sprintf("_rep%d", i) else ""
    abundance_path <- file.path(out_dir, sprintf("abundance_%s_%s%s.png", site_name, exp_tag, suffix))
    png(abundance_path, width = 1800, height = 900, res = 150, type = "cairo")
    plot_abundance(runs[[i]])
    dev.off()
    message("Saved: ", abundance_path)
    abundance_path
  })

  # 3D: static snapshot (final year) + animated (all years) for the first
  # replicate only, to avoid generating one large HTML per replicate by default.
  fig_3d <- plot_3d_abundance(runs[[1]])
  volume_path <- file.path(out_dir, sprintf("abundance_3d_%s_%s.html", site_name, exp_tag))
  if (!is.null(fig_3d)) {
    # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
    # writes a small "<name>_files/" dependency folder alongside the HTML
    # instead -- keep the two together when copying/viewing elsewhere.
    htmlwidgets::saveWidget(fig_3d, volume_path, selfcontained = FALSE)
    message("Saved: ", volume_path)
  }

  anim_path <- NULL
  if (animate) {
    anim_path <- file.path(out_dir, sprintf("abundance_3d_animated_%s_%s.html", site_name, exp_tag))
    plot_3d_abundance_animated(runs[[1]], out_path = anim_path)
  }

  invisible(list(abundance_png = saved, abundance_3d = fig_3d, abundance_3d_animated = anim_path))
}
