# allsites.R
# Per-site colonization runs — vertical niche partitioning of epiphytic Maxillariinae
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(rgee)
source("scripts/helper_functions.R")
source("scripts/get_climateinputs.R")
source("scripts/get_microenv.R")
source("scripts/get_colonization.R")
source("scripts/paths.R")

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
                      sprintf("allsites_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("allsites.R started")

# ── Load models and observations ──────────────────────────────────────────────
model_files <- list.files(PROCESSED_DIR, pattern = "^pointmodel_.*\\.rds$",
                          full.names = TRUE)
models <- list()
for (f in model_files) models <- c(models, readRDS(f))
log_msg(sprintf("Loaded %d height models across %d sites",
                length(models), length(model_files)))

prep   <- prepare_observations("data/csv/combinedv3.csv", models)
obs    <- prep$obs
niches <- extract_niches(obs, models, prep$valid_per_model)
saveRDS(niches, file.path(PROCESSED_DIR, "niches.rds"))

# ── Site list ─────────────────────────────────────────────────────────────────
sites_df <- make_sites("data/csv/combinedv3.csv", pad = 0.001)
site_names <- c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu")

res <- 5  # metres per cell

# ── GEE init ─────────────────────────────────────────────────────────────────
ee_Initialize(project = "ee-lizethestevezt",
              user    = "lizethestevezt@gmail.com", drive = TRUE)

# ── Parameters ────────────────────────────────────────────────────────────────
# canopy_z set per site in loop below
superforgivingparams <- list(
  # ── Growth (Mondragón et al. 2007) ──
  beta0GrowthS = 0.22,  beta0GrowthJ = 0.45,  beta0GrowthA = 1.40,
  beta1Growth  = 0.3,   # literature: 0.8
  sigma        = 0.8,   # literature: 0.5
  
  # ── Survival (Raventós et al. 2015) ──
  beta0Seedling = 2.0,  # literature: -1.5
  beta0Juvenile = 2.5,  # literature: -0.5
  beta0Adult    = 3.0,  # literature:  1.0
  beta1         = 0.105,
  
  # ── Reproduction ──
  p_flower = 0.80, p_poll = 0.50, p_germ = 0.30, p_s1 = 0.83,
  # repRate = 0.80 * 0.50 * 0.30 * 0.83 ≈ 0.10
  
  # ── Dispersal (Murren & Ellison 1998) ──
  Ut       = 1,
  lambda   = 1,
  canopy_z = 10  # placeholder — overwritten per site below
)

# ── Per-site loop ─────────────────────────────────────────────────────────────
results <- list()

for (sn in site_names) {
  log_msg(sprintf("══ Starting site: %s ══", sn))
  
  site      <- sites_df[sites_df$Site == sn, ]
  site_niches <- niches[niches$Area_or_Site == sn, ]
  
  if (nrow(site_niches) == 0) {
    log_msg(sprintf("No observations for %s — skipping", sn))
    next
  }
  
  # compute grid dimensions
  lat_range_m <- (max(site_niches$lat) - min(site_niches$lat)) * 111000
  lon_range_m <- (max(site_niches$lon) - min(site_niches$lon)) * 111000 *
    cos(mean(site_niches$lat) * pi / 180)
  xDim <- max(round(lon_range_m / res), 10) + 4
  yDim <- max(round(lat_range_m / res), 10) + 4
  
  log_msg(sprintf("Grid: %d x %d cells | %d observations | %d species",
                  xDim, yDim, nrow(site_niches),
                  length(unique(site_niches$FinalID))))
  
  # canopy grid
  out_dir <- file.path(RAW_DIR, sn)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  canopy_grid <- get_canopy_grid(site, xDim, yDim,
                                 resolution = res, out_dir = out_dir)
  med_val <- median(canopy_grid, na.rm = TRUE)
  if (is.na(med_val)) med_val <- 10
  canopy_grid[is.na(canopy_grid)] <- med_val
  log_msg(sprintf("Canopy grid: %.0f-%.0fm, NAs filled with %.0fm",
                  min(canopy_grid), max(canopy_grid), med_val))
  
  # update canopy_z for this site
  superforgivingparams$canopy_z <- mean(canopy_grid, na.rm = TRUE)
  
  # run
  results[[sn]] <- runcolonization(
    site            = site,
    niches          = site_niches,
    canopy_grid     = canopy_grid,
    models          = models,
    valid_per_model = prep$valid_per_model,
    timesteps       = 50,
    resolution      = res,
    carCap          = 1,
    maxDisp         = 5,
    spinup          = 15,
    stochastic      = FALSE,
    Visualize       = TRUE,
    sleeptime       = 0.1,
    visualize_dispersion = FALSE,
    parameters      = superforgivingparams,
    allsites        = TRUE
  )
  
  # save individual site result
  saveRDS(results[[sn]],
          file.path(PROCESSED_DIR, paste0("colonization_", sn, ".rds")))
  log_msg(sprintf("Site %s complete — saved.", sn))
}

# save combin results
saveRDS(results, file.path(PROCESSED_DIR, "colonization_allsites.rds"))
log_msg("All sites complete.")

# ── Post-hoc 3D inspection ────────────────────────────────────────────────────
# call for any site after the loop:
# plot_3d_abundance(results[["Mashpi"]])
# plot_3d_abundance(results[["Maquipucuna"]], t = 25)