# run_colonization_onesite.R
# Canopy colonization model — vertical niche partitioning of epiphytic Maxillariinae
# Runs one site × one parameter set. Designed to be called interactively or
# from a SLURM job array:
#   Rscript run_colonization_onesite.R <site> <params_file> <experiment_tag>
#
# Arguments (all optional — fall back to defaults if omitted):
#   site           Site name matching Area_or_Site in combinedv3.csv
#                  Default: "Maquipucuna"
#   params_file    Path to an RDS file containing the params list
#                  Default: uses the literature params defined below
#   experiment_tag Label appended to output file names for identification
#                  Default: "default"
#
# Output: data/processed/colonization_<site>_<tag>.rds
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(plotly)
library(ggplot2)
source("scripts/paths.R")
source("scripts/get_colonization.R")

# ── Command-line arguments (for cluster / batch runs) ────────────────────────
args         <- commandArgs(trailingOnly = TRUE)
site_name    <- if (length(args) >= 1) args[1] else "Maquipucuna"
params_file  <- if (length(args) >= 2) args[2] else NULL
exp_tag      <- if (length(args) >= 3) args[3] else "default"

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
  sprintf("colonization_%s_%s_%s.log", site_name, exp_tag,
          format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("run_colonization_onesite.R started")

# ── 1. Load microenvironment ──────────────────────────────────────────────────
microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s.rds", site_name))
microenv      <- readRDS(microenv_path)

# If this RDS was saved before .weather was added, patch it once here.
if (is.null(microenv$.weather)) {
  log_msg("Patching microenv with ERA5 weather from point model...")
  pm <- readRDS(file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site_name)))
  microenv$.weather <- pm[[1]]$weather
  saveRDS(microenv, microenv_path)
  rm(pm)
  log_msg("Patch saved.")
}
log_msg(sprintf("Heights available: %s",
  paste(names(microenv)[!names(microenv) %in% c(".spatial", ".weather")],
        collapse = ", ")))

# ── 2. Field observations ─────────────────────────────────────────────────────
niches <- read.csv("data/csv/combinedv3.csv")
niches <- niches[
  !is.na(niches$lat) & !is.na(niches$lon) &
  !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
log_msg(sprintf("Observations for %s: %d across %d species", site_name,
  sum(niches$Area_or_Site == site_name),
  length(unique(niches$FinalID[niches$Area_or_Site == site_name]))))

# ── 3. Canopy grid (fallback) ─────────────────────────────────────────────────
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name],
                    na.rm = TRUE)
log_msg(sprintf("Mean canopy height: %.1f m", mean_canopy))
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)  # used only if forestparams=NULL

# ── 4. Site object ────────────────────────────────────────────────────────────
site <- list(Site = site_name)

# ── 4b. Forest structure parameters ──────────────────────────────────────────
# Structural values from Myster (2017) primary Maquipucuna cloud forest plots:
#   mean dsh = 22.7 cm → trunk_r = 0.114 m
#   stem density = 272–324 trees/ha (≥10 cm dsh) → stems_per_ha = 298
# Crown geometry from pantropical allometry (Williams et al. 2019 review).
# Epiphyte footprint: ~4×5 cm pseudobulb cluster = 0.02 m² (field estimate).
forestparams <- list(
  stems_per_ha          = 298,    # Myster (2017) Table 5, mean of 4 primary MR plots
  mean_hgt              = 8.4,    # m — mean canopy height at this elevation
  sd_hgt                = 3.5,
  mean_crown_r          = 2.0,    # m — crown radius
  sd_crown_r            = 0.8,
  trunk_r               = 0.114,  # m — mean dsh 22.7 cm / 2 (Myster 2017)
  branch_density        = 3.0,    # m² branch surface per m² projected crown area
  epiphyte_footprint_m2 = 0.02    # m² bark area per Maxillariinae individual
)

# ── 5. Parameters ─────────────────────────────────────────────────────────────
# If a params file was passed, load it — otherwise use literature defaults.
# To run a sensitivity experiment:
#   params <- readRDS("params_experiment1.rds")  (or build programmatically)
#   saveRDS(params, "params_experiment1.rds")
#   Rscript run_colonization_onesite.R Maquipucuna params_experiment1.rds exp1
if (!is.null(params_file) && file.exists(params_file)) {
  params <- readRDS(params_file)
  log_msg(sprintf("Loaded params from %s", params_file))
} else {
  log_msg("Using default literature params.")
  params <- list(
    # ── s(z, e): survival ────────────────────────────────────────────────────
    beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
    beta1   =  0.10,
    z_S_min =  0.0,  z_S_max =  1.0,
    z_J_min =  1.0,  z_J_max =  7.0,
    z_A_min =  7.0,  z_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions ────────────────────────────
    psi0S        = -3.30,  psi0J        = -2.70,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_z_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
    p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ──────────────────────────────────────────────
    canopy_z = mean_canopy,  lambda = 1,  Ut = 1
  )
}

# ── 6. Sanity check ───────────────────────────────────────────────────────────
clim_test <- get_clim(as.numeric(sub("h", "",
  names(microenv)[!names(microenv) %in% c(".spatial", ".weather")][1])), microenv)
stopifnot(
  is.data.frame(clim_test),
  nrow(clim_test) == 48,
  all(c("temp", "relhum", "windspeed", "swdown", "precip", "winddir") %in% names(clim_test))
)
log_msg(sprintf("Climate check passed: %.1f°C mean temp, %.0f mm/yr precip",
  mean(clim_test$temp, na.rm = TRUE),
  mean(clim_test$precip, na.rm = TRUE) * 8760))

# ── 7. Run colonization model ─────────────────────────────────────────────────
log_msg(sprintf("Starting colonization run [%s | %s]...", site_name, exp_tag))

result <- runcolonization(
  site         = site,
  niches       = niches,
  canopy_grid  = canopy_grid,
  microenv     = microenv,
  timesteps    = 30,
  resolution   = 10,
  carCap       = 5,           # fallback scalar (used only if forestparams = NULL)
  maxDisp      = 10,
  spinup       = 5,
  Visualize    = !interactive(),
  sleeptime    = 0.2,
  parameters   = params,
  forestparams = forestparams
)

out_path <- file.path(PROCESSED_DIR,
  sprintf("colonization_%s_%s.rds", site_name, exp_tag))
saveRDS(result, out_path)
log_msg(sprintf("Done. Results saved to %s", out_path))

# ── 8. Plots (interactive only) ───────────────────────────────────────────────
if (interactive()) {
  plot_abundance(result)
  plot_3d_abundance(result)
}
