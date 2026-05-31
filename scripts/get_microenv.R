# get_microenv.R
# Microenvironment input into model
# Lizeth Estévez Tobar — University of Bonn, 2026
# ── Niche characterisation ────────────────────────────────────────────────────
# Prepares the observation dataframe for niche extraction:
# - snaps each observed height to the nearest modelled height step per site
# - builds hKey for lookup into the models list (e.g. "Mashpi_h5.0")
# - assigns each observation to its nearest ERA5 grid cell using per-site
#   cell coordinates — valid cells and their locations vary by site depending
#   on the raster extent used during data acquisition
# ERA5 cells are ~0.25° (~28km) — coarse macroclimate input, but microclimf
# downscales to fine-scale microclimate using DTM + vegetation structure
prepare_observations <- function(csv_path, models) {
  # get modelled heights — strip site prefix e.g. "Mashpi_h5.0" -> 5.0
  modelled_heights <- as.numeric(sub(".*_h", "", names(models)))

  # build per-site cell coordinate lookup from valid cells in models
  # each site has different valid cell indices and coordinates
  site_cell_coords <- list()
  for (key in names(models)) {
    site_name <- sub("_h.*", "", key) # e.g. "Mashpi_h5.0" -> "Mashpi"
    if (!site_name %in% names(site_cell_coords)) {
      m <- models[[key]]
      valid <- which(sapply(m, function(x) inherits(x, "micropoint")))
      site_cell_coords[[site_name]] <- do.call(rbind, lapply(valid, function(j) {
        data.frame(cell = j, lat = m[[j]]$lat, lon = m[[j]]$long)
      }))
    }
  }
  # build valid cell indices per model key — needed by extract_niches()
  valid_per_model <- lapply(models, function(m) {
    which(sapply(m, function(x) inherits(x, "micropoint")))
  })
  # snap observed height to nearest modelled height step
  snap_to_nearest <- function(h, available) {
    available[which.min(abs(available - h))]
  }

  # find nearest ERA5 cell for a given observation using site-specific coords
  nearest_cell <- function(obs_lat, obs_lon, cells) {
    dists <- sqrt((cells$lat - obs_lat)^2 + (cells$lon - obs_lon)^2)
    cells$cell[which.min(dists)]
  }

  obs <- read_csv(csv_path, col_types = COMBINED_COL_TYPES) |>
    filter(!is.na(Height_m), !is.na(Genus)) |>
    mutate(
      hSnapped = sapply(Height_m, snap_to_nearest, available = modelled_heights),
      hKey = sprintf("%s_h%.1f", Area_or_Site, hSnapped),
      datetime = as.POSIXct(datetime, tz = "UTC"),
      nearest_cell = mapply(function(lat, lon, site) {
        nearest_cell(lat, lon, site_cell_coords[[site]])
      }, lat, lon, Area_or_Site)
    ) |>
    select(
      Genus, species, FinalID, Height_m, hSnapped, hKey,
      CanopyHeight_m, lat, lon, Area_or_Site, datetime, nearest_cell
    )

  log_msg(sprintf(
    "Observations prepared: %d rows, %d unique heights, %d sites",
    nrow(obs), length(unique(obs$hKey)), length(unique(obs$Area_or_Site))
  ))
  return(list(obs = obs, valid_per_model = valid_per_model))
}

# Extracts microclimate niche summary for each observation:
# for each individual, retrieves the modelled hourly weather at its
# height and nearest ERA5 cell, splits into day (swdown > 0) and night,
# and summarises mean + sd for temp, relhum, swdown, lwdown, windspeed, precip.
# valid_per_model: named list of valid cell indices per model key, from
#   prepare_observations() — needed to match each observation to the correct
#   ERA5 cell for its site
# Day/night split captures ecologically distinct conditions:
# daytime = photosynthesis, transpiration, convective cooling;
# nighttime = dew formation, radiative cooling, cold stress
extract_niches <- function(obs, models, valid_per_model) {
  niches <- data.frame()

  for (i in 1:nrow(obs)) {
    h <- obs$hKey[i]
    c <- obs$nearest_cell[i]

    # retrieve hourly weather for this observation's height and ERA5 cell
    # use the correct valid cells for this specific model key
    valid_cells_for_h <- valid_per_model[[h]]
    c_actual <- valid_cells_for_h[which.min(abs(valid_cells_for_h - c))]
    weather <- models[[h]][[c_actual]]$weather
    weather$is_day <- weather$swdown > 0 # day = any incoming solar radiation
    day <- weather[weather$is_day, ]
    night <- weather[!weather$is_day, ]

    # build one-row summary of day and night microclimate conditions
    # as.numeric() used defensively in case of unexpected type coercion
    avgs <- data.frame(
      # daytime variables (including radiation)
      day_temp_mean = mean(as.numeric(day$temp)),
      day_temp_sd = sd(as.numeric(day$temp)),
      day_relhum_mean = mean(as.numeric(day$relhum)),
      day_relhum_sd = sd(as.numeric(day$relhum)),
      day_swdown_mean = mean(as.numeric(day$swdown)),
      day_swdown_sd = sd(as.numeric(day$swdown)),
      day_lwdown_mean = mean(as.numeric(day$lwdown)),
      day_lwdown_sd = sd(as.numeric(day$lwdown)),
      day_windspeed_mean = mean(as.numeric(day$windspeed)),
      day_windspeed_sd = sd(as.numeric(day$windspeed)),
      day_precip_mean = mean(as.numeric(day$precip)),
      day_precip_sd = sd(as.numeric(day$precip)),
      # nighttime variables (no swdown — always zero at night)
      night_temp_mean = mean(as.numeric(night$temp)),
      night_temp_sd = sd(as.numeric(night$temp)),
      night_relhum_mean = mean(as.numeric(night$relhum)),
      night_relhum_sd = sd(as.numeric(night$relhum)),
      night_lwdown_mean = mean(as.numeric(night$lwdown)),
      night_lwdown_sd = sd(as.numeric(night$lwdown)),
      night_windspeed_mean = mean(as.numeric(night$windspeed)),
      night_windspeed_sd = sd(as.numeric(night$windspeed)),
      night_precip_mean = mean(as.numeric(night$precip)),
      night_precip_sd = sd(as.numeric(night$precip)),
      # observation metadata for joining back to species/site information
      Genus = obs$Genus[i],
      species = obs$species[i],
      FinalID = obs$FinalID[i],
      Height_m = obs$Height_m[i],
      hSnapped = obs$hSnapped[i],
      Area_or_Site = obs$Area_or_Site[i],
      lat = obs$lat[i],
      lon = obs$lon[i]
    )
    niches <- rbind(niches, avgs)
  }

  log_msg(sprintf(
    "Niche extraction complete: %d observations, %d unique species",
    nrow(niches), length(unique(niches$FinalID))
  ))
  return(niches)
}
# Retrieves hourly weather data for a given height and site from the models list.
# Used as the primary microclimate lookup inside all submodel functions.
# Returns the full weather dataframe (672 rows × 10 columns) for the height step
# nearest to the requested height.

# TODO: match to nearest ERA5 cell for [x,y] position rather than always
#       using first valid cell — adequate for prototype
get_clim <- function(site_name, height, models, valid_per_model) {
  h_key  <- sprintf("%s_h%.1f", site_name, height)
  cell_c <- valid_per_model[[h_key]][1]
  if (is.null(cell_c) || is.na(cell_c)) return(NULL)
  models[[h_key]][[cell_c]]$weather
}

# Extracts a monthly slice of hourly weather data.
# Divides the 672-hour weather dataframe into 12 equal chunks of 56 hours.
# TODO: replace with real calendar month slicing once annual ERA5 data available
get_clim_month <- function(clim, month) {
  if (is.null(clim) || nrow(clim) == 0) {
    return(NULL)
  }
  hours_per_month <- floor(nrow(clim) / 12)
  month_hours <- ((month - 1) * hours_per_month + 1):(month * hours_per_month)
  month_hours <- month_hours[month_hours <= nrow(clim)]
  clim[month_hours, ]
}
# Downloads ETH Global Canopy Height 2020 (Lang et al.) for the site extent
# via Google Earth Engine and resamples to simulation resolution.
# Returns a [xDim × yDim] numeric matrix of canopy heights in metres.
# Source: ETH GlobalCanopyHeight_2020_10m_v1 (users/nlang/...)
# References: Lang et al. 2023
get_canopy_grid <- function(site, xDim, yDim, resolution, out_dir,
                            google_drive_folder = "rgee_backup") {
  cache_file <- file.path(out_dir, paste0(site$Site, "_canopy_grid.rds"))
  if (file.exists(cache_file)) {
    message("Canopy grid cache found, loading...")
    return(readRDS(cache_file))
  }

  message("Downloading canopy height raster from GEE...")
  e <- c(site$lon_min, site$lat_min, site$lon_max, site$lat_max)
  aoi <- ee$Geometry$Rectangle(e)
  img <- ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")$
    select("b1")$clip(aoi)

  # Check Drive first — skip export if already there
  googledrive::drive_auth(email = "lizethestevezt@gmail.com")
  folder <- googledrive::drive_find(pattern = google_drive_folder, type = "folder", n_max = 1)
  drive_files <- googledrive::drive_ls(folder)
  drive_file <- drive_files[grepl("canopy_height", drive_files$name), ]

  if (nrow(drive_file) == 0) {
    message("Exporting canopy height to Drive...")
    task <- ee$batch$Export$image$toDrive(
      image          = img,
      description    = "canopy_height_export",
      folder         = google_drive_folder,
      fileNamePrefix = paste0(site$Site, "_canopy_height"),
      region         = aoi$bounds()$getInfo()$coordinates[[1]],
      scale          = resolution,
      crs            = "EPSG:4326"
    )
    task$start()
    rgee::ee_monitoring(task, max_attempts = 200, quiet = FALSE)
    drive_files <- googledrive::drive_ls(folder)
    drive_file <- drive_files[grepl("canopy_height", drive_files$name), ]
  } else {
    message("Canopy height found on Drive, downloading...")
  }

  tmp_path <- tempfile(fileext = ".tif")
  googledrive::drive_download(file = drive_file[1, ], path = tmp_path, overwrite = TRUE)
  canopy_rast <- terra::rast(tmp_path)

  target_rast <- terra::rast(
    nrows = yDim, ncols = xDim,
    xmin = site$lon_min, xmax = site$lon_max,
    ymin = site$lat_min, ymax = site$lat_max,
    crs = "EPSG:4326"
  )
  canopy_rast <- terra::resample(canopy_rast, target_rast, method = "bilinear")

  canopy_mat <- t(as.matrix(canopy_rast, wide = TRUE))
  canopy_mat[is.na(canopy_mat)] <- median(canopy_mat, na.rm = TRUE)

  saveRDS(canopy_mat, cache_file)
  message(
    "Canopy grid saved: ", xDim, " × ", yDim, " cells, range ",
    round(min(canopy_mat)), "–", round(max(canopy_mat)), "m"
  )
  return(canopy_mat)
}
