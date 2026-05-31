# get_climate_inputs.R
# Data acquisition pipeline for canopymicroenv
# Lizeth Estévez Tobar — University of Bonn, 2026

# ── Site preparation ──────────────────────────────────────────────────────────

# Reads the combined CSV and derives a single-row site summary with:
# - bounding box (lat/lon min/max) for ERA5 and raster construction
# - time window (tme_start/tme_end) for ERA5 request
# - observed height range (hObs_min/hObs_max) for model height sequence
make_site <- function(csv_path) {
  message("Reading combined CSV...")
  df <- read_csv(csv_path,
                 na        = c("", "NA", "N/A"),
                 col_types = COMBINED_COL_TYPES) |>
    dplyr::filter(!is.na(Source), !is.na(Area_or_Site))
  
  time_windows <- df |>
    dplyr::filter(!is.na(datetime)) |>
    dplyr::mutate(datetime = as.POSIXlt(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) |>
    dplyr::summarise(
      tme_start = min(datetime),
      tme_end   = max(datetime),
      hObs_min  = min(Height_m, na.rm = TRUE),
      hObs_max  = max(Height_m, na.rm = TRUE),
      .groups   = "drop"
    )
  
  coord_windows <- df |>
    dplyr::filter(!is.na(lat), !is.na(lon)) |>
    dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon)) |>
    dplyr::summarise(
      lat_min = min(lat), lat_max = max(lat),
      lon_min = min(lon), lon_max = max(lon),
      .groups = "drop"
    )
  
  site <- dplyr::bind_cols(data.frame(Site = "AllSites"), coord_windows, time_windows)
  message("Site created: AllSites | lon [", round(site$lon_min, 3), ", ", round(site$lon_max, 3),
          "] lat [", round(site$lat_min, 3), ", ", round(site$lat_max, 3), "]")
  message("Time window: ", format(site$tme_start), " to ", format(site$tme_end))
  return(site)
}

# Reads the combined CSV and derives a per-site summary dataframe with one row
# per field site (Area_or_Site), each with:
# - bounding box padded by pad° for ERA5 cell coverage
# - time window from actual observation datetimes at that site
# - observed height range for the model height sequence
# Use this for the per-site loop in getmicroenv.R; use make_site() when running
# a single AllSites bounding box instead.
make_sites <- function(csv_path, pad = 0.01) {
  message("Reading combined CSV...")
  df <- read_csv(csv_path,
                 na        = c("", "NA", "N/A"),
                 col_types = COMBINED_COL_TYPES) |>
    dplyr::filter(!is.na(Source), !is.na(Area_or_Site))
  
  df |>
    dplyr::filter(!is.na(lat), !is.na(lon), !is.na(datetime), !is.na(Height_m)) |>
    dplyr::mutate(
      lat      = as.numeric(lat),
      lon      = as.numeric(lon),
      datetime = as.POSIXlt(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    ) |>
    dplyr::group_by(Area_or_Site) |>
    dplyr::summarise(
      lat_min   = min(lat)      - pad,
      lat_max   = max(lat)      + pad,
      lon_min   = min(lon)      - pad,
      lon_max   = max(lon)      + pad,
      tme_start = min(datetime),
      tme_end   = max(datetime),
      hObs_min  = min(Height_m, na.rm = TRUE),
      hObs_max  = max(Height_m, na.rm = TRUE),
      .groups   = "drop"
    ) |>
    dplyr::rename(Site = Area_or_Site)
}
# ── ERA5 data acquisition ─────────────────────────────────────────────────────

# Merges the three ERA5 stepType netCDF files (accum, avg, instant) that CDS
# delivers separately into one combined file. Renames radiation variables to
# match the names expected by microclimdata::era5_process()
merge_era5_steptype_files <- function(pathin, pathout) {
  library(ncdf4)
  
  nc_files <- list.files(pathin, pattern = "stepType.*\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) stop("No stepType .nc files found in ", pathin)
  message("Merging ", length(nc_files), " ERA5 stepType files...")
  
  rename_map <- c(avg_snlwrf = "msnlwrf", avg_sdlwrf = "msdwlwrf")
  
  datasets  <- lapply(nc_files, nc_open)
  src       <- datasets[[1]]
  out_dims  <- lapply(src$dim, function(d) ncdim_def(d$name, d$units, d$vals, unlim = d$unlim))
  names(out_dims) <- names(src$dim)
  
  seen_vars <- names(src$dim)
  out_vars  <- list()
  var_data  <- list()
  
  for (ds in datasets) {
    for (vname in names(ds$var)) {
      if (vname %in% seen_vars) next
      seen_vars <- c(seen_vars, vname)
      var       <- ds$var[[vname]]
      out_name  <- ifelse(vname %in% names(rename_map), rename_map[vname], vname)
      var_dims  <- lapply(var$dim, function(d) out_dims[[d$name]])
      out_vars[[out_name]] <- ncvar_def(name = out_name, units = var$units,
                                        dim = var_dims, missval = var$missval)
      var_data[[out_name]] <- ncvar_get(ds, vname)
      message("  ", vname, " -> ", out_name)
    }
  }
  
  nc_out <- nc_create(pathout, vars = out_vars)
  for (vname in names(var_data)) ncvar_put(nc_out, vname, var_data[[vname]])
  nc_close(nc_out)
  lapply(datasets, nc_close)
  message("Merged file written to ", pathout)
  return(pathout)
}

# ERA5 land-sea mask (lsm) sometimes has near-land cells with values just below 1
# (e.g. 0.96) which microclimdata treats as ocean and excludes.
# This fix rounds near-land cells up to 1 so they are included in processing.
fix_lsm <- function(nc_path) {
  nc  <- ncdf4::nc_open(nc_path, write = TRUE)
  lsm <- ncdf4::ncvar_get(nc, "lsm")
  n   <- sum(lsm < 1 & lsm >= 0.95)
  lsm[lsm >= 0.95] <- 1
  ncdf4::ncvar_put(nc, "lsm", lsm)
  ncdf4::nc_close(nc)
  message("LSM fix: ", n, " near-land cells set to 1")
}

# Downloads ERA5 hourly climate data for the site bounding box and time window,
# merges stepType files, fixes LSM, and processes to grid format for microclimf.
# Skips download if merged file already exists (overwrite = FALSE).
get_weather <- function(site, credentials, r, tme, dir, overwrite = FALSE) {
  message("")
  merged_file <- file.path(dir, paste0(site$Site, ".nc"))
  
  if (!file.exists(merged_file) || overwrite) {
    expected_files <- c(
      "data_stream-oper_stepType-accum.nc",
      "data_stream-oper_stepType-avg.nc",
      "data_stream-oper_stepType-instant.nc"
    )
    existing <- file.exists(file.path(dir, expected_files))
    
    if (!all(existing) || overwrite) {
      message("Downloading ERA5 data for site ", site$Site, "...")
      req <- mcera5::build_era5_request(
        xmin         = site$lon_min, xmax = site$lon_max,
        ymin         = site$lat_min, ymax = site$lat_max,
        start_time   = site$tme_start, end_time = site$tme_end,
        by_month     = TRUE, outfile_name = site$Site
      )
      ecmwfr::wf_request(
        request  = req[[1]],
        user     = credentials$username[credentials$Site == "CDS"],
        transfer = TRUE, path = paste0(dir, "/"), retry = 120, verbose = TRUE
      )
      for (z in list.files(dir, pattern = "\\.zip$", full.names = TRUE)) {
        unzip(z, exdir = dir)
        unlink(z)
      }
    } else {
      message("ERA5 stepType files already exist, skipping download.")
    }
    
    message("Merging ERA5 stepType files...")
    merge_era5_steptype_files(pathin = dir, pathout = merged_file)
    message("Fixing land-sea mask...")
    fix_lsm(merged_file)
    files <- list.files(dir, pattern = "^data_stream", full.names = TRUE)
    if (length(files) > 0) { unlink(files); message("Deleted ", length(files), " stepType files.") }
    
  } else {
    message("ERA5 merged file already exists for ", site$Site, ", skipping download.")
    fix_lsm(merged_file)
  }
  
  message("Processing ERA5 data to grid format...")
  weatherdata <- microclimdata::era5_process(
    tme = tme, req = NA, pathin = paste0(dir, "/"), r = r, out = "grid"
  )
  message("Weather data ready.")
  return(weatherdata)
}


# ── Terrain and landcover ─────────────────────────────────────────────────────

# Downloads a digital elevation model for the site extent via elevatr.
# Reprojects to EPSG:32717 (UTM 17S) for metric calculations, then back to
# EPSG:4326 for use with microclimf. Caches as dtm.tif to avoid re-downloading.
get_dtm <- function(r, dir, mask = FALSE) {
  message("")
  cache_file <- file.path(dir, "dtm.tif")
  
  if (file.exists(cache_file)) {
    message("DTM cache found, loading...")
    return(terra::rast(cache_file))
  }
  
  message("Downloading digital elevation model...")
  terra::values(r) <- 1
  raster_utm <- terra::project(r, "EPSG:32717")
  terra::values(raster_utm) <- 1
  dtm <- microclimdata::dem_download(r = raster_utm, msk = mask)
  
  terra::writeRaster(dtm, cache_file)
  message("DTM downloaded and cached: ", nrow(dtm), " x ", ncol(dtm),
          " pixels at ", round(terra::res(dtm)[1], 1), "m resolution")
  return(dtm)
}

# Downloads ESA WorldCover 10m landcover via Google Earth Engine, exports to
# Google Drive, and downloads to disk. Checks Drive first to avoid re-exporting.
# Requires rgee initialisation before calling.
get_landcover <- function(site, r, out_dir, type = "ESA",
                          overwrite = FALSE, google_drive_folder = "rgee_backup") {
  message("")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_path    <- file.path(out_dir, paste0(site$Site, "_landcover_", type, ".tif"))
  drive_prefix <- paste0(site$Site, "_ESA_WorldCover")
  
  if (file.exists(save_path) && !overwrite) {
    message("Landcover already exists on disk, loading...")
    return(terra::rast(save_path))
  }
  
  message("Checking Google Drive for landcover...")
  googledrive::drive_auth(email = "lizethestevezt@gmail.com")
  folder      <- googledrive::drive_find(pattern = google_drive_folder, type = "folder", n_max = 1)
  drive_files <- googledrive::drive_ls(folder)
  drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]
  
  if (nrow(drive_file) == 0 || overwrite) {
    message("Exporting landcover from GEE to Drive (this takes ~15 min)...")
    e          <- terra::ext(r)
    epsg_code  <- paste0("EPSG:", terra::crs(r, describe = TRUE)$code)
    aoi        <- ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
    aoi_coords <- aoi$bounds()$getInfo()$coordinates[[1]]
    img        <- ee$ImageCollection("ESA/WorldCover/v100")$first()
    task <- ee$batch$Export$image$toDrive(
      image = img, description = paste0(site$Site, "_landcover_export"),
      folder = google_drive_folder, fileNamePrefix = drive_prefix,
      region = aoi_coords, scale = 10, crs = epsg_code
    )
    task$start()
    rgee::ee_monitoring(task, max_attempts = 200, quiet = FALSE)
    drive_file <- googledrive::drive_ls(folder) |> 
      dplyr::filter(grepl(drive_prefix, name))
  } else {
    message("Landcover found on Google Drive, downloading...")
  }
  
  googledrive::drive_download(file = drive_file[1, ], path = save_path, overwrite = TRUE)
  message("Landcover saved to ", save_path)
  return(terra::rast(save_path))
}


# ── Vegetation and soil parameters ────────────────────────────────────────────

# Reclassifies ESA landcover to microclimdata habitat types, then builds
# the full vegetation parameter object (canopy structure, LAI, reflectance)
# needed by runpointmodela(). Aggregates to coarser resolution for speed.
get_vegetation <- function(r, lcover, tme, lat, lon) {
  message("")
  message("Reclassifying ESA landcover to habitat types...")
  rcl <- matrix(c(
    10,  2,  20,  6,  30, 10,  40, 13,  50, 14,
    60, 16,  70, 16,  80, 12,  90, 12,  95,  2,  100, 10
  ), ncol = 2, byrow = TRUE)
  
  hbts        <- terra::classify(lcover, rcl, others = 2)
  hbts        <- terra::unwrap(hbts)
  hbts_coarse <- terra::aggregate(hbts, fact = 10, fun = "modal")
  message("Habitat raster: ", nrow(hbts_coarse), " x ", ncol(hbts_coarse),
          " pixels | unique types: ", length(unique(terra::values(hbts_coarse), na.rm = TRUE)))
  
  message("Running vegpfromhab...")
  vegetation <- vegpfromhab(habitats = hbts_coarse, tme = as.POSIXlt(tme, tz = "UTC"),
                            lat = lat, long = lon)
  message("Vegetation parameters ready.")
  return(vegetation)
}

# Downloads MODIS LAI (500m) for the site extent and time period via NASA
# Earthdata. Mosaics tiles, crops to template extent, and reprojects.
# For higher resolution (10m HRVPP) set reso = 10 — requires WEkEO credentials.
get_lai <- function(r, tme, reso = 500, pathout, credentials, template) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  if (!reso %in% c(10, 500)) stop("reso must be one of 10 or 500")
  
  year <- tme$year[1] + 1900
  
  if (reso == 10) {
    # HRVPP via WEkEO — placeholder, see legacy scripts for implementation
  } else {
    library(luna)
    if (tme[length(tme)] < as.POSIXlt("2000-02-18", tz = "UTC"))
      stop("No data available prior to 2000-02-18")
    
    existing_hdf <- list.files(pathout, pattern = "\\.hdf$", full.names = TRUE)
    
    if (length(existing_hdf) == 0) {
      message("Downloading MODIS LAI...")
      e  <- terra::ext(r)
      r2 <- terra::project(terra::rast(e, crs = terra::crs(r)), "EPSG:4326")
      e2 <- terra::ext(r2)
      st <- substr(as.character(tme[1]), 1, 10)
      ed <- substr(as.character(tme[length(tme)]), 1, 10)
      mf <- luna::getNASA("MOD15A2H", st, ed, aoi = e2, version = "061", download = FALSE)
      if (length(mf) == 0) stop("No data for specified location or time period")
      luna::getNASA("MOD15A2H", st, ed, aoi = e2, version = "061",
                    download = TRUE, path = pathout,
                    username = credentials$username[credentials$Site == "NASA"],
                    password = credentials$password[credentials$Site == "NASA"],
                    server   = "LPDAAC_ECS")
    } else {
      message("MODIS LAI files already exist (", length(existing_hdf), " files), skipping download.")
    }
  }
  
  lai_files <- list.files(pathout, pattern = "\\.hdf$", full.names = TRUE)
  message("Mosaicing ", length(lai_files), " LAI tiles...")
  lai <- terra::rast(lai_files[1])
  if (length(lai_files) > 1) {
    for (f in lai_files[-1]) {
      lai_new <- terra::rast(f)
      lai     <- terra::mosaic(lai, lai_new)
    }
  }
  message("Cropping and reprojecting LAI to template extent...")
  template_in_lai_crs <- terra::project(template, terra::crs(lai))
  lai <- terra::crop(lai, terra::ext(template_in_lai_crs))
  lai <- terra::project(lai, terra::crs(template))
  message("LAI ready: ", nlyr(lai), " layers")
  return(lai)
}

# Downloads and processes MODIS BRDF/albedo for the site.
# Fills NA values with the spatial mean and resamples to template resolution.
# Caches processed result as albedo_processed.rds to avoid re-downloading.
get_albedo <- function(template, tme, pathout, credentials) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  alb_cache <- file.path(pathout, "albedo_processed.rds")
  
  if (file.exists(alb_cache)) {
    message("Albedo cache found, loading...")
    return(readRDS(alb_cache))
  }
  
  message("Downloading MODIS albedo...")
  albedo_download(r = template, tme = tme, pathout = paste0(pathout, "/"), credentials = credentials)
  message("Processing albedo...")
  alb <- albedo_process(r = template, pathin = paste0(pathout, "/"))
  
  if (terra::crs(alb) != terra::crs(template)) {
    message("Reprojecting albedo to template CRS...")
    template_in_alb_crs <- terra::project(template, terra::crs(alb))
    alb <- terra::crop(alb, terra::ext(template_in_alb_crs))
    alb <- terra::project(alb, terra::crs(template))
  }
  
  message("Filling NAs and resampling albedo...")
  alb_val <- mean(terra::values(alb), na.rm = TRUE)
  alb[is.na(alb)] <- alb_val
  alb <- terra::resample(alb, template)
  
  saveRDS(alb, alb_cache)
  message("Albedo ready and cached.")
  return(alb)
}

# Builds the full soil/ground parameter object for microclimf:
# downloads SoilGrids physical properties, derives soil type,
# computes ground reflectance from LAI + albedo, and assembles
# a soilcharac object. Caches to groundparams.rds.
get_soil <- function(r, template, tme, credentials, landcover, dir, albedodir, laidir) {
  message("")
  cache_file <- file.path(dir, "groundparams.rds")
  if (file.exists(cache_file)) {
    message("Soil cache found, loading...")
    return(readRDS(cache_file))
  }
  
  message("[1/5] Downloading soil data from SoilGrids...")
  soil_r <- r
  terra::values(soil_r) <- 1
  soilproperties <- microclimdata::soildata_download(r = soil_r, pathdir = paste0(dir, "/"),
                                                     deletefiles = FALSE)
  message("Soil properties downloaded: ", paste(names(soilproperties), collapse = ", "))
  
  message("[2/5] Deriving soil type...")
  soiltype <- microclimdata:::soildata_gettype(soilproperties)
  message("Soil type range: ", min(terra::values(soiltype), na.rm = TRUE),
          " - ", max(terra::values(soiltype), na.rm = TRUE))
  
  message("[3/5] Getting LAI...")
  lai <- get_lai(r = r, tme = tme, pathout = laidir,
                 credentials = credentials, template = template)
  
  message("[4/5] Getting leaf inclination coefficient and albedo...")
  x   <- x_calc(landcover = landcover, lctype = "ESA")
  alb <- get_albedo(template = template, tme = tme,
                    pathout = albedodir, credentials = credentials)
  
  message("[5/5] Computing ground reflectance...")
  # Reduce to single layer and aggregate to coarser resolution to save memory
  alb_coarse <- terra::aggregate(mean(alb, na.rm = TRUE), fact = 10, fun = "mean")
  lai_coarse <- terra::aggregate(mean(lai, na.rm = TRUE), fact = 10, fun = "mean")
  x_coarse   <- terra::aggregate(x,                       fact = 10, fun = "mean")
  
  alb_coarse <- terra::resample(alb_coarse, x_coarse)
  lai_coarse <- terra::resample(lai_coarse, x_coarse)
  
  groundr <- reflectance_calc(alb = alb_coarse, lai = lai_coarse,
                              x = x_coarse, plotprogress = FALSE)$gref
  
  message("Building soilcharac object...")
  soilc <- list(
    soiltype = terra::wrap(terra::resample(soiltype, template, method = "near")),
    groundr  = terra::wrap(terra::resample(groundr,  template, method = "bilinear"))
  )
  class(soilc) <- "soilcharac"
  
  saveRDS(soilc, cache_file)
  message("Soil parameters ready and cached.")
  return(soilc)
}