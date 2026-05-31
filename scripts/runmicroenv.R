# runmicroenv.R
# Data acquisition and point model execution for canopy microenvironment characterisation
# Runs ERA5 download, terrain, landcover, vegetation and soil preparation,
# then executes runpointmodela() across the full observed height range.
# Output: named list of micropoint objects saved to data/processed/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Patches to external packages ─────────────────────────────────────────────
# These fix known bugs in dependencies without modifying package source.
# Must be applied before loading the packages they patch.

# ecmwfr::wf_request — newer versions added a `job_name` argument that the
# CDS API endpoint rejects. This patch wraps the original and drops job_name.
if (!exists("original_wf_request")) {
  original_wf_request <- ecmwfr::wf_request
  assignInNamespace("wf_request", 
                    function(request, user = "ecmwfr", transfer = TRUE, 
                             path = tempdir(), time_out = 7200, retry = 240, 
                             job_name, verbose = TRUE) {
                      original_wf_request(request = request, user = user, transfer = transfer, 
                                          path = path, time_out = time_out, retry = retry, 
                                          verbose  = verbose)}, ns = "ecmwfr")}

# microclimdata:::reflectance_calc — variable name bug: `wgt` used instead of
# `bwgt` in the original source, causing ground reflectance calculation to fail.
assignInNamespace("reflectance_calc",
  function(alb, lai, x, plotprogress = TRUE, maxiter = 50, tol = 0.001, bwgt = 0.5) {
    e1 <- terra::intersect(terra::ext(lai), terra::ext(alb))
    e   <- terra::intersect(e1, terra::ext(x))
    lai <- terra::crop(lai, e)
    alb <- terra::crop(alb, e)
    x   <- terra::crop(x, e)
    all_same <- terra::compareGeom(lai, alb, x)
    if (all_same) {
      tst  <- exp(-mean(as.vector(lai), na.rm = TRUE))
      lref <- (x * 0 + 0.5) * (1 - bwgt) + bwgt * alb  # fix: wgt -> bwgt
      gref <- x * 0 + 0.15
      mxdif <- tol * 10
      paim  <- as.matrix(lai, wide = TRUE)
      xm    <- as.matrix(x,   wide = TRUE)
      albm  <- as.matrix(alb, wide = TRUE)
      itr   <- 1
      while (mxdif > tol) {
        if (tst < 0.5) {
          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref2, wide = TRUE), paim, xm, albm), x)
          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
        } else {
          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref, wide = TRUE), paim, xm, albm), x)
          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
        }
        gref  <- bwgt * gref + (1 - bwgt) * gref2
        lref  <- bwgt * lref + (1 - bwgt) * lref2
        mxdif1 <- mean(abs(as.vector(gref) - as.vector(gref2)), na.rm = TRUE)
        mxdif2 <- mean(abs(as.vector(lref) - as.vector(lref2)), na.rm = TRUE)
        mxdif  <- max(mxdif1, mxdif2)
        itr <- itr + 1
        if (itr > maxiter) mxdif <- 0
      }
    } else {
      stop("Geometries of input rasters do not match")
    }
    return(list(gref = gref, lref = lref))
  },
  ns = "microclimdata"
)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)

# Python environment for rgee / Earth Engine
reticulate::use_python("/Users/lizethestevezt/miniforge3/bin/python", required = TRUE)

# ── Setup ─────────────────────────────────────────────────────────────────────
source("scripts/get_microenv.R")
source("scripts/paths.R")

mycredentials <- readRDS("/Users/lizethestevezt/canopymicroenv/credentials.rds")

# Build per-site table — one row per field site with bounding box,
# time window, and observed height range derived from the combined CSV.
sites <- make_sites("data/csv/combinedv3.csv", pad = 0.15)

# Important: Run these ONE BY ONE before running the loop
# they require authentication and a code to paste
ee$Authenticate()

# ── Logging ───────────────────────────────────────────────────────────────────
# Timestamped log file written to logs/ alongside console output
dir.create("/Users/lizethestevezt/canopymicroenv/logs", recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR, sprintf("runmicroenv_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("runmicroenv.R started")

# ── Site loop ─────────────────────────────────────────────────────────────────
for (site in split(sites, seq_len(nrow(sites)))) {
  log_msg(sprintf("── Site %s (%d/%d) ──", site$Site, which(sites$Site == site$Site), nrow(sites)))
  
  # Re-initialize Earth Engine for each site — session can expire between sites
  ee$Initialize(project = "ee-lizethestevezt")
  
  # Creating data storage directories 
  site_dir        <- file.path(RAW_DIR, site$Site)
  site_dtm_dir    <- file.path(site_dir, "dtm")
  site_soil_dir   <- file.path(site_dir, "soil")
  site_era5_dir   <- file.path(site_dir, "era5")
  site_alb_dir    <- file.path(site_dir, "albedo")
  site_lai_dir    <- file.path(site_dir, "lai")
  site_lcover_dir <- file.path(site_dir, "landcover")
  
  for (d in c(site_dir, site_dtm_dir, site_soil_dir, 
              site_era5_dir, site_alb_dir, site_lcover_dir, site_lai_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  # ── 0. Spatial and temporal objects ──────────────────────────────────────────
  # 2×2 raster covering site bounding box — used for landcover and vegetation
  raster <- terra::rast(
    nrows = 2, ncols = 2,
    xmin  = site$lon_min, xmax = site$lon_max,
    ymin  = site$lat_min, ymax = site$lat_max,
    crs   = "EPSG:4326"
  )
  raster_utm <- terra::project(raster, "EPSG:32717")
  
  # 50×50 raster for DTM — finer resolution needed for terrain modelling
  r_dtm <- terra::rast(
    nrows = 50, ncols = 50,
    xmin  = site$lon_min, xmax = site$lon_max,
    ymin  = site$lat_min, ymax = site$lat_max,
    crs   = "EPSG:4326"
  )
  
  # Buffered raster for ERA5 download — 0.5° padding ensures we capture the
  # full ERA5 cells overlapping the study area (ERA5 resolution ~0.25°)
  r_buffered <- terra::rast(
    nrows = 2, ncols = 2,
    xmin  = site$lon_min - 0.5, xmax = site$lon_max + 0.5,
    ymin  = site$lat_min - 0.5, ymax = site$lat_max + 0.5,
    crs   = "EPSG:4326"
  )
  
  # High-resolution template raster for soil parameter resampling
  TemplateRaster <- terra::rast(
    extent     = terra::ext(site$lon_min, site$lon_max, site$lat_min, site$lat_max),
    resolution = 0.0001,
    crs        = "EPSG:4326"
  )
  
  # Time objects — ERA5 is requested for full calendar months containing the
  # observation window (month-padded) to ensure complete diurnal cycles.
  # tmeHourly: hourly sequence for the full padded month (used for ERA5 + model)
  # tmeMonth:  daily sequence for the padded month (used for soil/albedo)
  tme_start <- as.POSIXlt(site$tme_start, tz = "UTC")
  tme_end   <- as.POSIXlt(site$tme_end,   tz = "UTC")
  
  tme <- as.POSIXlt(seq(from = tme_start, to = tme_end, by = "day"))
  
  tmeMonth <- as.POSIXlt(seq(
    from = as.POSIXlt(format(tme_start, "%Y-%m-01 00:00:00"), tz = "UTC"),
    to   = as.POSIXlt(format(seq(as.POSIXlt(format(tme_start, "%Y-%m-01"), tz = "UTC"),
                                 by = "month", length.out = 2)[2] - 3600,
                             "%Y-%m-%d %H:00:00"), tz = "UTC"),
    by = "day"
  ))
  
  tmeHourly <- as.POSIXlt(seq(
    from = as.POSIXlt(format(tme_start, "%Y-%m-01 00:00:00"), tz = "UTC"),
    to   = as.POSIXlt(format(seq(as.POSIXlt(format(tme_start, "%Y-%m-01"), tz = "UTC"),
                                 by = "month", length.out = 2)[2] - 3600,
                             "%Y-%m-%d %H:00:00"), tz = "UTC"),
    by = "hour"
  ))
  
  # ── 1. Data acquisition ───────────────────────────────────────────────────────
  
  # ERA5 hourly climate data for the buffered site extent and padded time window.
  # Skips download if merged file already exists.
  # TODO: once weatherdata structure is confirmed with str(weatherdata), implement
  # clipping to expedition window ±1 day (tme_start - 1day to tme_end + 1day)
  # to reduce model runtime without losing diurnal cycle completeness.
  log_msg("Acquiring weather data...")
  weatherdata <- get_weather(
    site        = site,
    credentials = mycredentials,
    r           = r_buffered,
    tme         = tmeHourly,
    dir         = site_era5_dir
  )
  
  # Digital elevation model — 50×50 grid, no coastal masking needed
  log_msg("Acquiring DTM...")
  dtmdata <- get_dtm(r = r_dtm, dir = site_dtm_dir, mask = FALSE)
  dtmdata <- terra::project(dtmdata, "EPSG:4326")
  
  # ESA WorldCover 10m landcover via Google Earth Engine
  # Required for both vegetation and soil parameter derivation
  log_msg("Acquiring landcover...")
  landcoverdata <- get_landcover(site = site, r = raster, out_dir = site_lcover_dir)
  
  # Vegetation parameters (canopy structure, LAI, reflectance) derived from
  # landcover reclassification via microclimdata::vegpfromhab()
  # PackedSpatRaster objects are unwrapped after loading (terra serialisation issue)
  log_msg("Deriving vegetation parameters...")
  vegetationdata <- get_vegetation(
    r      = raster,
    lcover = landcoverdata,
    tme    = tmeHourly,
    lat    = mean(c(site$lat_min, site$lat_max)),
    lon    = mean(c(site$lon_min, site$lon_max))
  )
  for (n in names(vegetationdata)) {
    if (class(vegetationdata[[n]])[1] == "PackedSpatRaster")
      vegetationdata[[n]] <- terra::unwrap(vegetationdata[[n]])
  }
  
  # Soil parameters: SoilGrids physical properties + MODIS LAI + albedo +
  # ground reflectance. Cached to groundparams.rds after first run.
  log_msg("Deriving soil parameters...")
  soildata <- get_soil(
    r         = raster,
    template  = TemplateRaster,
    tme       = tmeMonth,
    credentials = mycredentials,
    landcover = landcoverdata,
    dir       = site_soil_dir,
    albedodir = site_alb_dir,
    laidir    = site_lai_dir
  )
  for (n in names(soildata)) {
    if (class(soildata[[n]])[1] == "PackedSpatRaster")
      soildata[[n]] <- terra::unwrap(soildata[[n]])
  }
  
  # ── 2. Point model height loop ────────────────────────────────────────────────
  # Runs runpointmodela() once per height step across the full observed height
  # range of the dataset (site$hObs_min to site$hObs_max, step 0.2m).
  # Results stored as named list: models[["h3.4"]] = micropoint object at 3.4m.
  # runpointmodela() is used (not runpointmodel()) because climate forcing is
  # supplied as a gridded array (climarrayr), not a single-column dataframe.
  tme_model <- as.POSIXlt(tmeHourly, tz = "UTC")
  heights   <- seq(from = site$hObs_min, to = site$hObs_max, by = 0.1)
  models    <- list()
  
  log_msg(sprintf("Running point model: %.1f – %.1f m (%d steps)",
                  min(heights), max(heights), length(heights)))
  
  for (h in heights) {
    message("")
    log_msg(sprintf("runpointmodela @ %.1f m", h))
    message("")
    models[[sprintf("%s_h%.1f", site$Site, h)]] <- microclimf::runpointmodela(
      climarrayr = weatherdata,
      tme        = tme_model,
      reqhgt     = h,
      dtm        = dtmdata,
      vegp       = vegetationdata,
      soilc      = soildata
    )
    message("")
  }
  
  log_msg(sprintf(paste0("Saving %d models to pointmodel_", site$Site, ".rds"), length(models)))
  saveRDS(models, file.path(PROCESSED_DIR, paste0("pointmodel_", site$Site, ".rds")))
}
log_msg("runmicroenv.R complete")