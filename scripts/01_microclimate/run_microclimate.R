# run_microclimate.R
# Data acquisition and microclimate modelling for canopy microenvironment characterisation.
# For each site: downloads ERA5, DTM, landcover, LAI, albedo, and soil data;
# builds vegparams and soilcharac; runs the point model and then the grid model
# across the full observed height range.
# Output: per-site microenvironment lists saved to data/processed/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/patches.R")
# ── Libraries ─────────────────────────────────────────────────────────────────
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
library(luna)

# ── Patches to external packages ─────────────────────────────────────────────


# Python environment for rgee / Earth Engine
PYTHON_PATH <- Sys.getenv("CANOPY_PYTHON",
  unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python")
reticulate::use_python(PYTHON_PATH, required = TRUE)

# ── Setup ─────────────────────────────────────────────────────────────────────
source("scripts/01_microclimate/lib.R")
source("scripts/02_model/config/paths.R")
source("scripts/00_data_conversion/helper_functions.R")

mycredentials <- readRDS(file.path("/home/s38leste_hpc/canopymicroenv/credentials.rds"))

# Build per-site table — one row per field site with bounding box,
# time window, and observed height range derived from the combined CSV.
sites <- make_sites(OBSERVATIONS_CSV, pad = 0.15)

# Important: Run these ONE BY ONE before running the loop —
# they require authentication and a code to paste into the console.
ee$Authenticate()

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR, sprintf("runmicroenv_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("runmicroenv.R started")
models   <- list()
microenv <- list()

# ── Site loop ─────────────────────────────────────────────────────────────────
for (site in split(sites, seq_len(nrow(sites)))) {
  log_msg(sprintf("── Site %s (%d/%d) ──", site$Site, which(sites$Site == site$Site), nrow(sites)))

  # Re-initialize Earth Engine for each site — session can expire between sites
  ee$Initialize(project = "ee-lizethestevezt")

  # Data storage directories for this site
  site_dir      <- file.path(RAW_DIR, site$Site)
  site_dtm_dir  <- file.path(site_dir, "dtm")
  site_soil_dir <- file.path(site_dir, "soil")
  site_era5_dir <- file.path(site_dir, "era5")
  site_alb_dir  <- file.path(site_dir, "albedo")
  site_lai_dir  <- file.path(site_dir, "lai")
  site_lcover_dir <- file.path(site_dir, "landcover")

  for (d in c(site_dir, site_dtm_dir, site_soil_dir,
              site_era5_dir, site_alb_dir, site_lcover_dir, site_lai_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  site_env_path   <- file.path(PROCESSED_DIR, sprintf("microenv_%s.rds",   site$Site))
  site_model_path <- file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site$Site))

  if (file.exists(site_env_path)) {
    log_msg(sprintf("Microenv already exists for %s, skipping.", site$Site))
    next
  }

  if (file.exists(site_model_path)) {
    log_msg(sprintf("Point model already exists for %s, will load after data prep.", site$Site))
  }

  # ── 0. Spatial and temporal objects ──────────────────────────────────────────
  # Single 2×2 raster in WGS84 — the common spatial anchor passed to all
  # download functions. The packages crop and reproject internally.
  raster <- terra::rast(
    nrows = 2, ncols = 2,
    xmin  = site$lon_min, xmax = site$lon_max,
    ymin  = site$lat_min, ymax = site$lat_max,
    crs   = "EPSG:4326"
  )
  terra::values(raster) <- 1

  tme_start <- as.POSIXlt(site$tme_start, tz = "UTC")
  tme_end   <- as.POSIXlt(site$tme_end,   tz = "UTC")

  # Expand to full calendar month(s) to ensure complete diurnal cycles for ERA5.
  # era5_process() and runpointmodel() both need hourly sequences.
  tme <- as.POSIXlt(seq(
    from = as.POSIXlt(format(tme_start, "%Y-%m-01 00:00:00"), tz = "UTC"),
    to   = as.POSIXlt(format(
      seq(as.POSIXlt(format(tme_start, "%Y-%m-01"), tz = "UTC"),
          by = "month", length.out = 2)[2] - 3600,
      "%Y-%m-%d %H:00:00"), tz = "UTC"),
    by = "hour"
  ), tz = "UTC")

  # ── 1. Data acquisition ────────────────────────────────────────────────────

  log_msg("Acquiring weather data...")
  weatherdata <- get_weather(
    site        = site,
    credentials = mycredentials,
    r           = raster,
    tme         = tme,
    dir         = site_era5_dir,
    output      = "grid"
  )

  log_msg("Acquiring DTM...")
  dtmdata <- get_dtm(r = raster, dir = site_dtm_dir, mask = FALSE)

  log_msg("Acquiring landcover...")
  landcoverdata <- get_landcover(site = site, r = raster, out_dir = site_lcover_dir)

  log_msg("Acquiring LAI...")
  laidata <- get_lai(
    r           = raster,
    tme         = tme,
    pathout     = site_lai_dir,
    credentials = mycredentials
  )

  log_msg("Acquiring albedo...")
  albedodata <- get_albedo(
    r           = raster,
    tme         = tme,
    pathout     = site_alb_dir,
    credentials = mycredentials
  )

  # runmicro() requires dtm, vegp, and soilc to share an identical grid.
  # ESA WorldCover (10m) is far finer than the DTM (~90m), so resample
  # landcover to the DTM grid before building vegparams and soilcharac.
  # nearest-neighbour preserves discrete land-cover class values.
  landcover_dtm <- terra::resample(landcoverdata, dtmdata, method = "near")

  # Ground and leaf reflectance — derived from LAI, albedo, and landcover.
  # Must be computed before vegparams and soilcharac, both of which need refldata.
  log_msg("Computing reflectance...")
  refldata <- get_reflectance(
    lai       = laidata,
    alb       = albedodata,
    landcover = landcover_dtm,
    cachefile = file.path(site_dir, "reflectance.rds")
  )

  log_msg("Deriving vegetation parameters...")
  vegetationdata <- get_vegetation(
    r         = raster,
    lcover    = landcover_dtm,
    lai       = laidata,
    refldata  = refldata,
    dir       = site_dir,
    site_name = site$Site
  )

  log_msg("Deriving soil parameters...")
  soildata <- get_soil(
    r         = raster,
    dir       = site_soil_dir,
    landcover = landcover_dtm,
    refldata  = refldata
  )

  # ── 2. Point model ────────────────────────────────────────────────────────

  if (file.exists(site_model_path)) {
    log_msg(sprintf("Found existing point model for %s, loading from disk...", site$Site))
    models[[site$Site]] <- readRDS(site_model_path)
  } else {
    log_msg(sprintf("Running point model for site %s...", site$Site))
    # runpointmodela runs the point model for each ERA5 grid cell separately,
    # preserving spatial variation in macroclimate across the site.
    models[[site$Site]] <- microclimf::runpointmodela(
      climarrayr = weatherdata,
      tme        = tme,
      reqhgt     = 0.05,
      dtm        = dtmdata,
      vegp       = vegetationdata,
      soilc      = soildata
    )
    saveRDS(models[[site$Site]], site_model_path)
    log_msg(sprintf("Saved point model to %s", site_model_path))
  }

  # ── 3. Grid microclimate model (runmicro) ─────────────────────────────────

  # Subset the point model to the hottest and coldest hour of each month —
  # runmicro() runs the full spatial grid for those representative timesteps only,
  # which is far cheaper than running it for every hourly timestep.
  log_msg("Subsetting point model to monthly max/min temperature days...")
  micropoint_mx <- microclimf::subsetpointmodela(models[[site$Site]], tstep = "month", what = "tmax")
  micropoint_mn <- microclimf::subsetpointmodela(models[[site$Site]], tstep = "month", what = "tmin")

  heights        <- seq(0.1, site$hObs_max, by = 0.1)
  height_dir     <- file.path(PROCESSED_DIR, sprintf("microenv_%s_heights", site$Site))
  dir.create(height_dir, recursive = TRUE, showWarnings = FALSE)

  # dtmc must match the ERA5 grid exactly (same cells as micropointa).
  # weatherdata[[1]] is a PackedSpatRaster — unwrap one layer to get the template.
  # Keep as plain SpatRaster: method="R" passes dtmc to .cca() which calls dim()
  # on it directly; dim(PackedSpatRaster) returns NULL and breaks the array build.
  era5_template <- terra::rast(weatherdata[[1]])[[1]]
  dtmc <- terra::resample(dtmdata, era5_template, method = "bilinear")

  # Save each height to its own small file as it is computed so that:
  #   (a) memory never accumulates across heights, and
  #   (b) a killed run can resume from where it left off.
  # terra SpatRasters are not fork-safe — wrap before forking, unwrap inside.
  dtmdata_w  <- terra::wrap(dtmdata)
  dtmc_w     <- terra::wrap(dtmc)

  n_cores    <- max(1L, parallel::detectCores() - 1L)
  n_heights  <- length(heights)
  log_msg(sprintf("  Launching parallel height loop: %d heights on %d cores...",
                  n_heights, n_cores))

  # Workers can't call the parent log_msg — capture log_file in closure and
  # write directly. flock-style append is safe across forked processes on macOS.
  .log_file <- log_file
  .wlog <- function(msg) {
    stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "][worker] ", msg)
    cat(stamped, "\n", file = .log_file, append = TRUE)
  }

  parallel::mclapply(seq_along(heights), function(i) {
    h          <- heights[i]
    h_key      <- sprintf("h%.2f", h)
    h_rds_path <- file.path(height_dir, sprintf("%s.rds", h_key))

    if (file.exists(h_rds_path)) {
      .wlog(sprintf("[%d/%d] %.2f m — already done, skipping.", i, n_heights, h))
      return(invisible(NULL))
    }

    .wlog(sprintf("[%d/%d] %.2f m — starting runmicro...", i, n_heights, h))

    dtm_  <- terra::unwrap(dtmdata_w)
    dtmc_ <- terra::unwrap(dtmc_w)

    mout_mx <- microclimf::runmicro(
      micropoint = micropoint_mx,
      reqhgt     = h,
      vegp       = vegetationdata,
      soilc      = soildata,
      dtm        = dtm_,
      dtmc       = dtmc_,
      altcorrect = 1,
      method     = "R"
    )
    mout_mn <- microclimf::runmicro(
      micropoint = micropoint_mn,
      reqhgt     = h,
      vegp       = vegetationdata,
      soilc      = soildata,
      dtm        = dtm_,
      dtmc       = dtmc_,
      altcorrect = 1,
      method     = "R"
    )

    saveRDS(list(tmax = mout_mx, tmin = mout_mn), h_rds_path)
    .wlog(sprintf("[%d/%d] %.2f m — done. Tz max: [%.1f, %.1f] | min: [%.1f, %.1f]",
                  i, n_heights, h,
                  min(mout_mx$Tz, na.rm = TRUE), max(mout_mx$Tz, na.rm = TRUE),
                  min(mout_mn$Tz, na.rm = TRUE), max(mout_mn$Tz, na.rm = TRUE)))
    invisible(NULL)
  }, mc.cores = n_cores)

  # Report how many heights actually completed
  n_done <- sum(file.exists(file.path(height_dir, sprintf("h%.2f.rds", heights))))
  log_msg(sprintf("  Parallel height loop done: %d/%d heights complete.", n_done, n_heights))

  # All heights done — assemble the final list and write the combined RDS.
  log_msg("  All heights complete. Assembling combined microenv RDS...")
  site_height_envs <- list()
  for (h in heights) {
    h_key <- sprintf("h%.2f", h)
    site_height_envs[[h_key]] <- readRDS(file.path(height_dir, sprintf("%s.rds", h_key)))
  }

  # Attach spatial reference so downstream scripts can recover lon/lat axes.
  site_height_envs$.spatial <- list(
    ext = terra::ext(dtmdata),
    crs = terra::crs(dtmdata)
  )

  # Attach the ERA5 weather time series from the first valid point model cell.
  # This gives the colonization model access to precip, winddir, and all other
  # macro-climate variables that runmicro() doesn't output spatially.
  # Cell 1 is representative — macro variables (precip, wind direction) have
  # negligible within-site variation at ERA5 resolution (~9km grid).
  site_height_envs$.weather <- models[[site$Site]][[1]]$weather

  microenv[[site$Site]] <- site_height_envs

  saveRDS(site_height_envs, site_env_path)
  log_msg(sprintf("Saved combined microenvironment to %s", site_env_path))

  # Clean up per-height files now that the combined RDS is safely written.
  unlink(height_dir, recursive = TRUE)
  log_msg("  Per-height cache cleaned up.")
}

# ── 4. Save collated outputs ──────────────────────────────────────────────────
log_msg(sprintf("Saving %d site point models to master file...", length(models)))
saveRDS(models, file.path(PROCESSED_DIR, "pointmodel.rds"))

# Collect all per-site microenvs from disk (includes sites skipped this run).
all_site_envs <- lapply(
  setNames(sites$Site, sites$Site),
  function(s) {
    p <- file.path(PROCESSED_DIR, sprintf("microenv_%s.rds", s))
    if (file.exists(p)) readRDS(p) else NULL
  }
)
all_site_envs <- Filter(Negate(is.null), all_site_envs)
log_msg(sprintf("Saving %d site microenvironments to master file...", length(all_site_envs)))
saveRDS(all_site_envs, file.path(PROCESSED_DIR, "microenv.rds"))

log_msg("runmicroenv.R complete")
