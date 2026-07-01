# run_experiments_thesis.R
# One-at-a-time sensitivity experiments for the full thesis colonization model.
# Requires microenv_Maquipucuna.rds in data/processed/.
# Runs experiments in parallel across parameter values (one core per value×rep).
#
# Usage: Rscript scripts/run_experiments_thesis.R
#   or:  source("scripts/run_experiments_thesis.R")  (interactive)
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(parallel)
source("scripts/paths.R")
source("scripts/get_colonization.R")

N_CORES <- max(1L, detectCores() - 1L)
cat(sprintf("Using %d parallel cores\n", N_CORES))

# ── Load site data (shared across all runs) ───────────────────────────────────

microenv_path <- file.path(PROCESSED_DIR, "microenv_Maquipucuna.rds")
if (!file.exists(microenv_path))
  stop("microenv_Maquipucuna.rds not found — run run_microclimate.R first")

cat("Loading microenv...\n")
microenv <- readRDS(microenv_path)
if (is.null(microenv$.weather)) {
  pm <- readRDS(file.path(PROCESSED_DIR, "pointmodel_Maquipucuna.rds"))
  microenv$.weather <- pm[[1]]$weather; rm(pm)
}

niches      <- read.csv("data/csv/combinedv3.csv")
niches      <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                      !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == "Maquipucuna"],
                    na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site        <- list(Site = "Maquipucuna")

forestparams <- list(
  stems_per_ha = 298, mean_hgt = 8.4, sd_hgt = 3.5,
  mean_crown_r = 2.0, sd_crown_r = 0.8, trunk_r = 0.114,
  branch_density = 3.0, epiphyte_footprint_m2 = 0.02
)

# ── Baseline parameters ───────────────────────────────────────────────────────

base_params <- list(
  beta0S  = -0.24, beta0J  =  0.41, beta0A  =  1.73, beta1  =  0.10,
  z_S_min =  0.0,  z_S_max =  1.0,
  z_J_min =  1.0,  z_J_max =  7.0,
  z_A_min =  7.0,  z_A_max = 20.0,
  psi0S        = -3.30, psi0J        = -2.70,
  beta_precip  =  3e-4, beta_rh      =  0.010,
  sigma        =  0.10, delta_z_base =  0.80,
  cost_repro   =  0.50,
  p_poll  = 0.30, p_germ  = 0.001, p_s1 = 0.45,
  canopy_z = mean_canopy, lambda = 1, Ut = 1
)

# ── Logging (suppressed in parallel workers) ──────────────────────────────────

log_msg <- function(msg) message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

# ── Single-run wrapper ────────────────────────────────────────────────────────



