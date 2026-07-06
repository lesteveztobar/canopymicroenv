# run_colonization_allsites.R
# Runs the colonization model for all sites × one parameter set.
# Can be called interactively or submitted as a SLURM job array where each
# task handles one site:
#
#   sbatch --array=1-5 run_slurm.sh
#   (run_slurm.sh calls: Rscript run_colonization_allsites.R <params_file> <tag>)
#
# Or run all sites sequentially in one session:
#   Rscript run_colonization_allsites.R params_experiment1.rds exp1
#
# Arguments (optional):
#   params_file    Path to RDS with params list. Default: literature params.
#   experiment_tag Label appended to output filenames. Default: "default"
#   slurm_task_id  If set, runs only that site index (1-based). Used by SLURM.
#                  Reads from SLURM_ARRAY_TASK_ID env var automatically.
#
# Output: data/processed/colonization_<site>_<tag>.rds for each site
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

# ── Arguments ─────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
params_file <- if (length(args) >= 1) args[1] else NULL
exp_tag     <- if (length(args) >= 2) args[2] else "default"

# SLURM array task ID — if set, run only that site
slurm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))

# ── Sites ─────────────────────────────────────────────────────────────────────
all_sites <- c("Maquipucuna", "MiradorMindo", "Mashpi", "MindoTarabita", "Yanayacu")
sites_to_run <- if (!is.na(slurm_id)) all_sites[slurm_id] else all_sites

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
  sprintf("allsites_%s_%s.log", exp_tag, format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg(sprintf("run_colonization_allsites.R started | tag=%s | sites=%s",
                exp_tag, paste(sites_to_run, collapse = ", ")))

# ── Field observations (loaded once, filtered per site inside loop) ────────────
niches_all <- read.csv("data/csv/combinedv3.csv")
niches_all <- niches_all[
  !is.na(niches_all$lat) & !is.na(niches_all$lon) &
  !is.na(niches_all$Height_m) & !is.na(niches_all$FinalID), ]

# ── Parameters ────────────────────────────────────────────────────────────────
if (!is.null(params_file) && file.exists(params_file)) {
  base_params <- readRDS(params_file)
  log_msg(sprintf("Loaded params from %s", params_file))
} else {
  log_msg("Using default literature params.")
  base_params <- list(
    beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
    beta1   =  0.10,
    z_S_min =  0.0,  z_S_max =  1.0,
    z_J_min =  1.0,  z_J_max =  7.0,
    z_A_min =  7.0,  z_A_max = 20.0,
    psi0S        = -3.30,  psi0J        = -2.70,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_z_base =  0.80,
    cost_repro   =  0.50,
    p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
    lambda  = 1,     Ut      = 1
  )
}

# ── Site loop ─────────────────────────────────────────────────────────────────
results <- list()

for (site_name in sites_to_run) {
  log_msg(sprintf("── %s ──", site_name))

  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s.rds", site_name))
  if (!file.exists(microenv_path)) {
    log_msg(sprintf("  microenv not found for %s — skipping", site_name))
    next
  }

  microenv <- readRDS(microenv_path)
  if (is.null(microenv$.weather)) {
    pm_path <- file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site_name))
    if (file.exists(pm_path)) {
      pm <- readRDS(pm_path)
      microenv$.weather <- pm[[1]]$weather
      saveRDS(microenv, microenv_path)
      rm(pm)
      log_msg("  Patched .weather into microenv.")
    } else {
      log_msg("  No point model found for .weather patch — skipping")
      next
    }
  }

  niches      <- niches_all[niches_all$Area_or_Site == site_name, ]
  mean_canopy <- mean(niches$CanopyHeight_m, na.rm = TRUE)
  canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)  # fallback only

  # canopy_z is site-specific — update it in params
  params <- base_params
  params$canopy_z <- mean_canopy

  # Forest structure — site-specific parameters.
  # Sources:
  #   Maquipucuna (~1340m): Myster (2017) primary MR plots — 298 stems/ha, mean dsh 22.7cm
  #   Mashpi (~920m): our observed mean canopy height 17.5m; stem density from
  #     Homeier et al. (2010) for lower montane wet forest at ~1000m in Ecuador:
  #     ~350 stems/ha (≥10cm dbh), mean dbh ~24cm
  #   MiradorMindo (~1313m): elevation close to Maquipucuna; slightly adjusted
  #     toward lower-montane values (Homeier et al. 2010)
  #   MindoTarabita (~1450m, estimated): upper montane transition; Valencia et al.
  #     (1994) mid-elevation Ecuador — ~320 stems/ha, smaller mean dbh ~20cm
  #   Yanayacu (~1600m, estimated): upper cloud forest; Homeier et al. (2010)
  #     upper montane band — ~380 stems/ha, shorter trees, dbh ~18cm
  #
  # branch_density and epiphyte_footprint_m2 are kept constant across sites
  # (no site-specific literature available; used as calibration parameters).
  forest_by_site <- list(
    Maquipucuna = list(
      stems_per_ha = 298,   mean_hgt = 8.4,  sd_hgt = 3.5,
      mean_crown_r = 2.0,   sd_crown_r = 0.8,
      trunk_r = 0.114,                          # dsh 22.7cm / 2 / 100
      branch_density = 3.0, epiphyte_footprint_m2 = 0.02
    ),
    Mashpi = list(
      stems_per_ha = 350,   mean_hgt = 17.0, sd_hgt = 5.0,
      mean_crown_r = 2.5,   sd_crown_r = 1.0,
      trunk_r = 0.120,                          # dsh ~24cm
      branch_density = 3.0, epiphyte_footprint_m2 = 0.02
    ),
    MiradorMindo = list(
      stems_per_ha = 310,   mean_hgt = 9.5,  sd_hgt = 3.5,
      mean_crown_r = 2.1,   sd_crown_r = 0.8,
      trunk_r = 0.115,
      branch_density = 3.0, epiphyte_footprint_m2 = 0.02
    ),
    MindoTarabita = list(
      stems_per_ha = 320,   mean_hgt = 7.5,  sd_hgt = 3.0,
      mean_crown_r = 1.9,   sd_crown_r = 0.7,
      trunk_r = 0.100,                          # dsh ~20cm
      branch_density = 3.0, epiphyte_footprint_m2 = 0.02
    ),
    Yanayacu = list(
      stems_per_ha = 380,   mean_hgt = 6.0,  sd_hgt = 2.5,
      mean_crown_r = 1.7,   sd_crown_r = 0.6,
      trunk_r = 0.090,                          # dsh ~18cm
      branch_density = 3.0, epiphyte_footprint_m2 = 0.02
    )
  )
  forestparams <- forest_by_site[[site_name]]
  if (is.null(forestparams)) {
    log_msg(sprintf("  No site-specific forestparams for %s — using Maquipucuna defaults", site_name))
    forestparams <- forest_by_site[["Maquipucuna"]]
  }

  out_path <- file.path(PROCESSED_DIR,
    sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (file.exists(out_path)) {
    log_msg(sprintf("  Result already exists (%s) — skipping", out_path))
    next
  }

  result <- tryCatch(
    runcolonization(
      site         = list(Site = site_name),
      niches       = niches_all,
      canopy_grid  = canopy_grid,
      microenv     = microenv,
      timesteps    = 30,
      resolution   = 10,
      carCap       = 5,
      maxDisp      = 10,
      spinup       = 5,
      Visualize    = FALSE,
      parameters   = params,
      forestparams = forestparams
    ),
    error = function(e) {
      log_msg(sprintf("  ERROR: %s", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(result)) {
    saveRDS(result, out_path)
    log_msg(sprintf("  Saved to %s", out_path))
    results[[site_name]] <- result
  }
}

log_msg(sprintf("All done. Completed %d/%d sites.",
                length(results), length(sites_to_run)))
