# get_climateinputs.R
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
# - bounding box padded by pad° to ensure ERA5 cells overlap the study area
# - time window from actual observation datetimes at that site
# - observed height range for the model height sequence
# Use this for the per-site loop in runmicroenv.R; use make_site() when running
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
# match the names expected by microclimdata::era5_process().
merge_era5_steptype_files <- function(pathin, pathout) {
  library(ncdf4)

  nc_files <- list.files(pathin, pattern = "stepType.*\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) stop("No stepType .nc files found in ", pathin)
  message("Merging ", length(nc_files), " ERA5 stepType files...")

  # CDS renamed two radiation variables between API versions
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
# merges stepType files, fixes LSM, and processes to a point climate data frame
# for runpointmodel(). Skips download if merged file already exists.
get_weather <- function(site, credentials, r, tme, dir, overwrite = FALSE, output = "point") {
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

  message("Processing ERA5 data to point climate data frame...")
  weatherdata <- microclimdata::era5_process(
    tme = tme, req = NA, pathin = paste0(dir, "/"), r = r, out = output
  )
  message("Weather data ready.")
  return(weatherdata)
}


# ── Terrain and landcover ─────────────────────────────────────────────────────

# Downloads a digital elevation model for the site extent via elevatr.
# Projects r to UTM first so the downloaded DEM has a sensible metric resolution,
# then reprojects the result back to WGS84 (EPSG:4326) to match all other inputs.
# Caches as dtm.tif to avoid re-downloading.
get_dtm <- function(r, dir, mask = FALSE) {
  message("")
  cache_file <- file.path(dir, "dtm.tif")

  if (file.exists(cache_file)) {
    message("DTM cache found, loading...")
    return(terra::rast(cache_file))
  }

  message("Downloading digital elevation model...")
  # dem_download resamples to the template raster, so a 2×2 template gives a
  # 2×2 DTM. Instead, build a ~90m resolution UTM template to get a usable DEM,
  # then reproject to WGS84. EPSG:32717 = UTM zone 17S covers western Ecuador.
  e_wgs <- terra::ext(r)
  r_utm <- terra::project(terra::rast(e_wgs, crs = "EPSG:4326"), "EPSG:32717")
  e_utm <- terra::ext(r_utm)
  r_utm90 <- terra::rast(
    xmin = e_utm$xmin, xmax = e_utm$xmax,
    ymin = e_utm$ymin, ymax = e_utm$ymax,
    res  = 90, crs = "EPSG:32717"
  )
  terra::values(r_utm90) <- 1
  dtm <- microclimdata::dem_download(r = r_utm90, msk = FALSE)
  dtm <- terra::project(dtm, "EPSG:4326")

  terra::writeRaster(dtm, cache_file)
  message("DTM downloaded and cached: ", nrow(dtm), " x ", ncol(dtm),
          " pixels at ", round(terra::res(dtm)[1] * 111320, 0), "m resolution")
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

# Downloads MODIS LAI (500m) for the site extent and time period via NASA
# Earthdata, then mosaics tiles using microclimdata::lai_mosaic().
# Skips download if HDF files already exist in pathout.
get_lai <- function(r, tme, pathout, credentials, reso = 500) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  if (!reso %in% c(10, 500)) stop("reso must be one of 10 or 500")
  if (tme[length(tme)] < as.POSIXlt("2000-02-18", tz = "UTC"))
    stop("No MODIS data available prior to 2000-02-18")

  lai_file     <- file.path(pathout, "lai_mosaic.tif")
  existing_hdf <- list.files(pathout, pattern = "\\.hdf$", full.names = TRUE)

  if (!file.exists(lai_file)) {
    if (length(existing_hdf) == 0) {
      message("Downloading MODIS LAI...")
      microclimdata::lai_download(
        r           = r,
        tme         = tme,
        reso        = reso,
        pathout     = paste0(pathout, "/"),
        credentials = data.frame(
          username = credentials$username[credentials$Site == "NASA"],
          password = credentials$password[credentials$Site == "NASA"]
        )
      )
    } else {
      message("MODIS LAI files already exist (", length(existing_hdf), " files), skipping download.")
    }
    message("Mosaicing LAI tiles...")
    laidata <- microclimdata::lai_mosaic(r = r, pathin = paste0(pathout, "/"), reso = reso)
    terra::writeRaster(laidata, lai_file)
  } else {
    message("LAI mosaic cache found, loading...")
    laidata <- terra::rast(lai_file)
  }

  message("LAI ready: ", nrow(laidata), "x", ncol(laidata), " @ ", nlyr(laidata), " layers")
  return(laidata)
}

# Downloads and processes MODIS BRDF/albedo (WSA shortwave, band 30) for the
# site. Caches the processed SpatRaster as albedo_processed.rds to avoid
# re-downloading on subsequent runs.
get_albedo <- function(r, tme, pathout, credentials) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  alb_cache <- file.path(pathout, "albedo_processed.rds")

  if (file.exists(alb_cache)) {
    message("Albedo cache found, loading...")
    return(readRDS(alb_cache))
  }

  message("Downloading MODIS albedo...")
  microclimdata::albedo_download(
    r           = r,
    tme         = tme,
    pathout     = paste0(pathout, "/"),
    credentials = credentials
  )
  message("Processing albedo...")
  albedodata <- microclimdata::albedo_process(r = r, pathin = paste0(pathout, "/"))

  saveRDS(albedodata, alb_cache)
  message("Albedo ready and cached.")
  return(albedodata)
}

# Computes ground and leaf reflectance from LAI, albedo, and the leaf inclination
# coefficient derived from landcover. Aggregates all inputs to a common coarse
# grid before calling reflectance_calc() to avoid geometry mismatch errors.
# Returns a list with $gref (ground reflectance) and $lref (leaf reflectance).
get_reflectance <- function(lai, alb, landcover, cachefile = NULL) {
  message("")

  if (!is.null(cachefile) && file.exists(cachefile)) {
    message("Reflectance cache found, loading...")
    cached <- readRDS(cachefile)
    return(list(gref = terra::unwrap(cached$gref), lref = terra::unwrap(cached$lref)))
  }

  # x_calc maps ESA landcover codes to per-pixel leaf inclination coefficients
  message("Computing leaf inclination coefficients...")
  x_lc <- microclimdata::x_calc(landcover = landcover, lctype = "ESA")

  # Collapse multi-layer LAI to a single mean layer, then align all three inputs
  # to x_lc's grid before passing to reflectance_calc
  message("Aggregating inputs to common grid for reflectance_calc...")
  lai_agg <- terra::resample(terra::app(lai, mean, na.rm = TRUE), x_lc)
  alb_agg <- terra::resample(alb, x_lc)

  message("Computing reflectance...")
  refldata <- microclimdata::reflectance_calc(
    lai          = lai_agg,
    alb          = alb_agg,
    x            = x_lc,
    plotprogress = FALSE
  )
  message("refldata$gref range: ", round(min(terra::values(refldata$gref), na.rm = TRUE), 3),
          " – ", round(max(terra::values(refldata$gref), na.rm = TRUE), 3))
  message("refldata$lref range: ", round(min(terra::values(refldata$lref), na.rm = TRUE), 3),
          " – ", round(max(terra::values(refldata$lref), na.rm = TRUE), 3))

  if (!is.null(cachefile)) {
    saveRDS(list(gref = terra::wrap(refldata$gref), lref = terra::wrap(refldata$lref)), cachefile)
    message("Reflectance cached.")
  }
  return(refldata)
}

# Builds the vegparams object for microclimf using microclimdata::create_veggrid().
# Downloads vegetation height from GEE if not already cached in dir.
# Resamples vhgt and lai to the landcover grid before creating vegp, as
# create_veggrid() requires all inputs to share the same geometry.
get_vegetation <- function(r, lcover, lai, refldata, dir, site_name) {
  message("")
  vhgt_file    <- file.path(dir, "vhgt.tif")
  drive_prefix <- paste0("canopy_height_", site_name)

  if (!file.exists(vhgt_file)) {
    googledrive::drive_auth(email = "lizethestevezt@gmail.com")
    folder      <- googledrive::drive_find(pattern = "rgee_backup", type = "folder", n_max = 1)
    drive_files <- googledrive::drive_ls(folder)
    drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]

    if (nrow(drive_file) == 0) {
      message("Vegetation height not found on Drive — exporting from GEE for ", site_name, "...")
      # patch vegheight_download to use a site-specific Drive filename
      .orig_vhgt <- get("vegheight_download", envir = getNamespace("microclimdata"))
      assignInNamespace("vegheight_download",
        function(r, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
          reticulate::use_python(paste0(pathtopython, "python"), required = TRUE)
          if (!is.na(projectname)) rgee::ee$Initialize(project = projectname)
          e  <- terra::ext(r)
          r2 <- terra::rast(e); terra::crs(r2) <- terra::crs(r)
          r2 <- terra::project(r2, "EPSG:4326"); e <- terra::ext(r2)
          proj_string <- terra::crs(r, describe = TRUE)
          epsg_code   <- paste0("EPSG:", proj_string$code)
          aoi         <- rgee::ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
          aoi_coords  <- aoi$bounds()$getInfo()$coordinates[[1]]
          canopy_height <- rgee::ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
          task <- rgee::ee$batch$Export$image$toDrive(
            image          = canopy_height,
            description    = paste0("canopy_height_", site_name),
            folder         = GoogleDrivefolder,
            fileNamePrefix = drive_prefix,
            region         = aoi_coords,
            scale          = 10,
            crs            = epsg_code
          )
          task$start()
          if (!silent) microclimdata:::.monitor_task(task$id)
        },
        ns = "microclimdata"
      )
      microclimdata::vegheight_download(
        r                 = r,
        GoogleDrivefolder = "rgee_backup",
        pathtopython      = "/Users/lizethestevezt/miniforge3/bin/",
        projectname       = "ee-lizethestevezt"
      )
      # re-check Drive after export
      drive_files <- googledrive::drive_ls(folder)
      drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]
      if (nrow(drive_file) == 0)
        stop("GEE export completed but ", drive_prefix, " not found on Drive")
    } else {
      message("Vegetation height found on Drive for ", site_name, ", downloading...")
    }

    googledrive::drive_download(file = drive_file[1, ], path = vhgt_file, overwrite = TRUE)
    message("Vegetation height downloaded and cached.")
  }
  vhgt <- terra::rast(vhgt_file)

  # All inputs to create_veggrid must share the same geometry — resample to lcover grid
  message("Resampling vhgt, lai, and lref to landcover grid...")
  vhgt_rs <- terra::resample(vhgt,         lcover,        method = "bilinear")
  lai_rs   <- terra::resample(lai,          lcover,        method = "bilinear")
  lref_rs  <- terra::resample(refldata$lref, lcover,       method = "bilinear")
  gref_rs  <- terra::resample(refldata$gref, lcover,       method = "bilinear")

  message("Creating vegparams...")
  vegetationdata <- microclimdata::create_veggrid(
    landcover = lcover,
    vhgt      = vhgt_rs,
    lai       = lai_rs,
    refldata  = list(gref = gref_rs, lref = lref_rs),
    lctype    = "ESA"
  )
  message("vegparams ready | pai layers: ", nlyr(terra::rast(vegetationdata$pai)))
  return(vegetationdata)
}

# Builds the soilcharac object for microclimf using microclimdata::create_soilgrid().
# Downloads SoilGrids physical properties if not already cached.
# refldata provides ground reflectance (gref); the soil type is derived internally
# by create_soilgrid() from the physical properties.
get_soil <- function(r, dir, landcover, refldata) {
  message("")
  soil_cache <- file.path(dir, "soilproperties.rds")

  if (!file.exists(soil_cache)) {
    message("Downloading SoilGrids data...")
    soil_r <- r
    terra::values(soil_r) <- 1
    soilprops <- microclimdata::soildata_download(
      r           = soil_r,
      pathdir     = paste0(dir, "/"),
      deletefiles = FALSE
    )
    saveRDS(soilprops, soil_cache)
  } else {
    message("Soil cache found, loading...")
    soilprops <- readRDS(soil_cache)
  }

  # gref and lref must match the landcover grid geometry for create_soilgrid
  message("Creating soilcharac...")
  gref_rs <- terra::resample(refldata$gref, landcover, method = "bilinear")
  lref_rs <- terra::resample(refldata$lref, landcover, method = "bilinear")
  soildata <- microclimdata::create_soilgrid(
    soildata  = soilprops,
    refldata  = list(gref = gref_rs, lref = lref_rs),
    landcover = landcover
  )
  # soildata_downscale assigns 0 to water/masked pixels; checkinputs requires 1–11
  st <- terra::rast(soildata$soiltype)
  st[st < 1 | st > 11] <- NA
  soildata$soiltype <- terra::wrap(st)

  message("soilcharac ready")
  return(soildata)
}
