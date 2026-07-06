# height_resolution_experiment.R
# Outcome diagnostic (as opposed to resolution_diagnostics.R's timing-only
# diagnostic): runs the actual colonization model at each microenv height
# resolution (0.1 / 0.25 / 0.5 / 1.0m) with real timesteps/spinup and several
# replicates each, to check whether height resolution changes RESULTS
# (population trajectories, extinction rate) and not just runtime.
#
# Uses best_case.rds by default (generous params, persists reliably) rather
# than literature defaults, so trajectories stay non-zero and comparable
# across height steps -- isolates the resolution question from the separate
# "does the model persist under realistic params" question.
#
# Horizontal resolution is held fixed at 10 (matches production run_one()) --
# only height resolution varies here.
#
# Requires microenv_<site>.rds (0.1m) and microenv_<site>_h<step>.rds for
# coarser steps to already exist (height_res_array.sh).
#
# Usage: Rscript scripts/complex_model/height_resolution_experiment.R <site> [height_steps] [params_file]
#   e.g.: Rscript scripts/complex_model/height_resolution_experiment.R Maquipucuna "0.1,0.25,0.5,1.0" data/params/best_case.rds
# Output: data/processed/height_resolution_<site>.rds (tidy: height_step, n_heights, rep, t, totalS/J/A, total, extinct)
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(parallel)
source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

args         <- commandArgs(trailingOnly = TRUE)
site_name    <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
height_steps <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(strsplit(args[2], ",")[[1]]) else c(0.1, 0.25, 0.5, 1.0)
params_file  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else file.path(PARAMS_DIR, "best_case.rds")

TIMESTEPS  <- 30  # matches production run_colonization_onesite.R
SPINUP     <- 5
RESOLUTION <- 10  # fixed -- isolates height-resolution effect only

N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)

log_msg <- function(msg) message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

# ── Site setup (mirrors run_colonization_onesite.R) ─────────────────────────
niches <- read.csv("data/csv/combinedv3.csv")
niches <- niches[
  !is.na(niches$lat) & !is.na(niches$lon) &
  !is.na(niches$Height_m) & !is.na(niches$FinalID), ]

mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site <- list(Site = site_name)

forestparams <- list(
  stems_per_ha = 298, mean_hgt = 8.4, sd_hgt = 3.5,
  mean_crown_r = 2.0, sd_crown_r = 0.8, trunk_r = 0.114,
  branch_density = 3.0, epiphyte_footprint_m2 = 0.02
)

params <- readRDS(params_file)
params$canopy_z <- mean_canopy
N_REPS <- if (!is.null(params$n_reps)) params$n_reps else 5
log_msg(sprintf("Loaded params from %s (n_reps=%d)", params_file, N_REPS))

# ── Run one height step's worth of replicates ───────────────────────────────
run_height_step <- function(step) {
  suffix <- if (step != 0.1) sprintf("_h%.2f", step) else ""
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, suffix))
  if (!file.exists(microenv_path)) {
    log_msg(sprintf("Skipping height step %.2fm -- no %s (run height_res_array.sh first)",
                    step, microenv_path))
    return(NULL)
  }

  # run_one()/run_replicated() (get_colonization.R) read `microenv` as a
  # global via lexical scoping -- update it here so they pick up this step's.
  microenv <<- readRDS(microenv_path)
  n_heights <- length(microenv_heights(microenv))
  log_msg(sprintf("Height step %.2fm: %d tiers. Building climate cache...", step, n_heights))
  clim_cache <- build_clim_cache(microenv)

  log_msg(sprintf("Height step %.2fm: running %d replicate(s) x %d timesteps (spinup %d)...",
                  step, N_REPS, TIMESTEPS, SPINUP))
  out <- run_replicated(params, n_reps = N_REPS, timesteps = TIMESTEPS,
                        spinup = SPINUP, clim_cache = clim_cache)

  df <- out$summary
  if (is.null(df) || nrow(df) == 0) {
    log_msg(sprintf("Height step %.2fm: no replicate succeeded.", step))
    return(NULL)
  }
  df$height_step <- step
  df$n_heights   <- n_heights
  df
}

all_rows   <- lapply(height_steps, run_height_step)
summary_df <- do.call(rbind, Filter(Negate(is.null), all_rows))

if (is.null(summary_df) || nrow(summary_df) == 0) {
  stop("No height-step variants produced usable results for ", site_name)
}

out_path <- file.path(PROCESSED_DIR, sprintf("height_resolution_%s.rds", site_name))
saveRDS(summary_df, out_path)

# ── Descriptive comparison ──────────────────────────────────────────────────
final_df <- summary_df[summary_df$t == TIMESTEPS, ]

cat("\n== Final total abundance (t=", TIMESTEPS, ") by height step ==\n", sep = "")
print(aggregate(total ~ height_step, data = final_df,
                FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))))

# `extinct` is constant within a (height_step, rep) -- one row per replicate.
extinct_per_rep <- aggregate(extinct ~ height_step + rep, data = summary_df,
                             FUN = function(x) x[1])
cat("\n== Extinction rate by height step ==\n")
print(aggregate(extinct ~ height_step, data = extinct_per_rep, FUN = mean))

# Supplementary omnibus check -- indicative only given small n_reps.
if (length(unique(final_df$height_step)) > 1 && nrow(final_df) >= 4) {
  kt <- tryCatch(kruskal.test(total ~ factor(height_step), data = final_df),
                error = function(e) NULL)
  if (!is.null(kt)) {
    cat("\n== Kruskal-Wallis test: final total abundance ~ height_step ==\n")
    print(kt)
  }
}

cat(sprintf("\nSaved: %s\n", out_path))
