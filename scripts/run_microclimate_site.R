# run_microclimate_site.R
# Runs the full microclimate pipeline for ONE site, specified as a command-line arg.
# Called by hpc_job.sh: Rscript run_microclimate_site.R Maquipucuna
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript run_microclimate_site.R <SiteName> [n_months] [height_step]")
TARGET_SITE <- args[1]
N_MONTHS    <- if (length(args) >= 2) as.integer(args[2]) else 12L
# Vertical spacing (m) between height tiers. Default 0.1 matches production
# sites. A coarser step (e.g. 0.5) is for the height-resolution efficiency
# diagnostic — it reuses the cached weather/DTM/point-model data (those don't
# depend on height spacing) and only redoes the per-height loop below, so
# testing a coarser resolution is proportionally cheaper, not free.
HEIGHT_STEP <- if (length(args) >= 3) as.numeric(args[3]) else 0.1

source("scripts/patches.R")
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
library(luna)

PYTHON_PATH <- Sys.getenv("CANOPY_PYTHON",
  unset = "/home/lestevez/miniforge3/bin/python")
reticulate::use_python(PYTHON_PATH, required = TRUE)

source("scripts/get_microenv.R")
source("scripts/get_climateinputs.R")
source("scripts/paths.R")
source("scripts/helper_functions.R")

mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))
cds_row <- mycredentials[mycredentials$Site == "CDS", ]
ecmwfr::wf_set_key(key = cds_row$password, user = cds_row$username)
sites         <- make_sites(file.path(CSV_DIR, "combinedv3.csv"), pad = 0.15)
site          <- sites[sites$Site == TARGET_SITE, ]
if (nrow(site) == 0) stop(sprintf("Site '%s' not found in sites table.", TARGET_SITE))
site <- split(site, seq_len(nrow(site)))[[1]]

ee$Initialize(project = "ee-lizethestevezt")

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
  sprintf("microclim_%s_%s.log", TARGET_SITE, format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg(sprintf("run_microclimate_site.R started for %s", TARGET_SITE))

# ── Directories ───────────────────────────────────────────────────────────────
site_dir        <- file.path(RAW_DIR,  site$Site)
site_dtm_dir    <- file.path(site_dir, "dtm")
site_soil_dir   <- file.path(site_dir, "soil")
site_era5_dir   <- file.path(site_dir, "era5")
site_alb_dir    <- file.path(site_dir, "albedo")
site_lai_dir    <- file.path(site_dir, "lai")
site_lcover_dir <- file.path(site_dir, "landcover")
for (d in c(site_dir, site_dtm_dir, site_soil_dir,
            site_era5_dir, site_alb_dir, site_lcover_dir, site_lai_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Resolution-suffixed file names only kick in for a non-default height step,
# so production (0.1m) runs keep their original, unsuffixed file names and
# coarser diagnostic runs never collide with or overwrite them.
res_suffix      <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""
site_env_path   <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site$Site, res_suffix))
site_model_path <- file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site$Site))

# Per-height temp files go to Lustre scratch (set via CANOPY_SCRATCH in hpc_job.sh).
# Falls back to PROCESSED_DIR for local/interactive runs.
scratch_base <- Sys.getenv("CANOPY_SCRATCH", unset = "")
scratch_base <- if (nchar(scratch_base) > 0 && dir.exists(scratch_base)) scratch_base else PROCESSED_DIR

if (file.exists(site_env_path)) {
  log_msg(sprintf("Microenv already exists for %s — nothing to do.", site$Site))
  quit(status = 0)
}

# ── Spatial / temporal setup ──────────────────────────────────────────────────
raster <- terra::rast(
  nrows = 2, ncols = 2,
  xmin  = site$lon_min, xmax = site$lon_max,
  ymin  = site$lat_min, ymax = site$lat_max,
  crs   = "EPSG:4326"
)
terra::values(raster) <- 1

tme_obs_end <- as.POSIXlt(site$tme_end, tz = "UTC")
tme_to <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "month", length.out = 2L)[2L] - 3600L,
  "%Y-%m-%d %H:00:00"), tz = "UTC")
tme_from <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "-1 month", length.out = N_MONTHS)[N_MONTHS],
  "%Y-%m-01 00:00:00"), tz = "UTC")
tme <- as.POSIXlt(seq(from = tme_from, to = tme_to, by = "hour"), tz = "UTC")
site$tme_start <- tme_from
site$tme_end   <- tme_to

# ── 1. Data acquisition ───────────────────────────────────────────────────────
log_msg("Acquiring weather data...")
weatherdata <- get_weather(site = site, credentials = mycredentials,
                           r = raster, tme = tme, dir = site_era5_dir, output = "grid")

log_msg("Acquiring DTM...")
dtmdata <- get_dtm(r = raster, dir = site_dtm_dir, mask = FALSE)

log_msg("Acquiring landcover...")
landcoverdata <- get_landcover(site = site, r = raster, out_dir = site_lcover_dir)

log_msg("Acquiring LAI...")
laidata <- get_lai(r = raster, tme = tme, pathout = site_lai_dir, credentials = mycredentials)

log_msg("Acquiring albedo...")
albedodata <- get_albedo(r = raster, tme = tme, pathout = site_alb_dir, credentials = mycredentials)

landcover_dtm <- terra::resample(landcoverdata, dtmdata, method = "near")

log_msg("Computing reflectance...")
refldata <- get_reflectance(lai = laidata, alb = albedodata, landcover = landcover_dtm,
                            cachefile = file.path(site_dir, "reflectance.rds"))

log_msg("Deriving vegetation parameters...")
vegetationdata <- get_vegetation(r = raster, lcover = landcover_dtm, lai = laidata,
                                 refldata = refldata, dir = site_dir, site_name = site$Site)

log_msg("Deriving soil parameters...")
soildata <- get_soil(r = raster, dir = site_soil_dir,
                     landcover = landcover_dtm, refldata = refldata)

# ── 2. Point model ────────────────────────────────────────────────────────────
if (file.exists(site_model_path)) {
  log_msg("Loading existing point model...")
  model <- readRDS(site_model_path)
} else {
  log_msg("Running point model...")
  model <- microclimf::runpointmodela(
    climarrayr = weatherdata, tme = tme, reqhgt = 0.05,
    dtm = dtmdata, vegp = vegetationdata, soilc = soildata)
  saveRDS(model, site_model_path)
  log_msg(sprintf("Saved point model to %s", site_model_path))
}

# ── 3. Grid model (parallel heights) ─────────────────────────────────────────

# runpointmodela returns NA for cells where the model fails (e.g. ocean/masked
# cells). runmicro iterates weather[[k]] up to length(dtmc) — the ERA5 cell
# count — so the model list must stay that length. Replace NA cells with a copy
# of the first valid cell; their spatial output is masked/unused anyway.
n_before  <- length(model)
valid_mp  <- Filter(function(x) inherits(x, "micropoint"), model)
if (length(valid_mp) == 0) stop("All grid cells failed in runpointmodela — check ERA5 coverage and LSM.")
model <- lapply(model, function(x) if (inherits(x, "micropoint")) x else valid_mp[[1]])
log_msg(sprintf("Valid grid cells: %d / %d", length(valid_mp), n_before))

log_msg("Subsetting point model to monthly max/min days...")
micropoint_mx <- microclimf::subsetpointmodela(model, tstep = "month", what = "tmax")
micropoint_mn <- microclimf::subsetpointmodela(model, tstep = "month", what = "tmin")

heights    <- seq(0.1, site$hObs_max, by = HEIGHT_STEP)
height_dir <- file.path(scratch_base, sprintf("microenv_%s%s_heights", site$Site, res_suffix))
dir.create(height_dir, recursive = TRUE, showWarnings = FALSE)

era5_template <- terra::rast(weatherdata[[1]])[[1]]
dtmc          <- terra::resample(dtmdata, era5_template, method = "bilinear")
dtmdata_w     <- terra::wrap(dtmdata)
dtmc_w        <- terra::wrap(dtmc)

n_cores   <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1L)))
n_heights <- length(heights)
log_msg(sprintf("Launching parallel height loop: %d heights on %d cores...", n_heights, n_cores))

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
    .wlog(sprintf("[%d/%d] %.2f m — skipped.", i, n_heights, h))
    return(invisible(NULL))
  }
  .wlog(sprintf("[%d/%d] %.2f m — starting...", i, n_heights, h))
  dtm_  <- terra::unwrap(dtmdata_w)
  dtmc_ <- terra::unwrap(dtmc_w)
  mout_mx <- microclimf::runmicro(micropoint = micropoint_mx, reqhgt = h,
    vegp = vegetationdata, soilc = soildata, dtm = dtm_, dtmc = dtmc_,
    altcorrect = 1, method = "R")
  mout_mn <- microclimf::runmicro(micropoint = micropoint_mn, reqhgt = h,
    vegp = vegetationdata, soilc = soildata, dtm = dtm_, dtmc = dtmc_,
    altcorrect = 1, method = "R")
  saveRDS(list(tmax = mout_mx, tmin = mout_mn), h_rds_path)
  .wlog(sprintf("[%d/%d] %.2f m — done.", i, n_heights, h))
  invisible(NULL)
}, mc.cores = n_cores)

n_done <- sum(file.exists(file.path(height_dir, sprintf("h%.2f.rds", heights))))
log_msg(sprintf("Height loop done: %d/%d complete.", n_done, n_heights))

# ── 4. Save manifest ─────────────────────────────────────────────────────────
# Each height file is ~6 GB in memory (10 spatial arrays × 374×372×288).
# Loading all heights at once would require hundreds of GB — instead save a
# small manifest that points to the per-height files in scratch. Downstream
# code reads individual heights on demand via load_height().
log_msg("Saving microenv manifest...")
manifest <- list(
  .heights    = heights,
  .height_dir = height_dir,
  .spatial    = list(ext = terra::ext(dtmdata), crs = terra::crs(dtmdata)),
  .weather    = model[[1]]$weather
)
saveRDS(manifest, site_env_path)
log_msg(sprintf("Manifest saved to %s  (height_dir=%s)", site_env_path, height_dir))
log_msg("Done. Height files remain in scratch — do not ws_release until downstream analysis is complete.")
