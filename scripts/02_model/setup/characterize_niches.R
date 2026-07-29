# characterize_niches.R
# Precomputes each species' realized climate niche (temp/RH/light) pooled
# across EVERY site where it was observed, not just the one site.rds file
# happens to be modeling — a species with 2-3 records at one site may have
# several more at others. Saves data/processed/species_niches.rds, loaded by
# init_colonization() (get_colonization.R) in preference to the this-site-
# only get_niche() fallback. Also saves data/processed/niche_background_density.rds
# (the pooled reference distribution every species is scored against), reused
# by check_niche_suitability.R / plot_niche_suitability() for diagnostics.
#
# Each observed individual is augmented with the 12 within-year monthly
# conditions at its recorded height, rather than collapsed to that height's
# single annual mean — a species realistically tolerates whatever seasonal
# range its occupied height experiences over a year, not just that height's
# average, and pooling the full monthly record turns each raw observation
# into 12 data points instead of 1 (helpful for the kernel density fit, not
# just the old box), giving a far better-supported niche estimate than the
# handful of field observations alone could. Each monthly value already
# averages a representative warm day + cold day (48 hourly values), so it
# reflects a month-scale typical condition rather than a single transient
# hour. swdown (light) additionally drops zero/night-time readings before
# averaging, consistent with how light is used elsewhere in the model (see
# run_pass2_establish()/survival_logit() in get_colonization.R) — a monthly
# mean that included every dark hour would just be diluted by a fixed
# day/night ratio rather than reflecting daytime irradiance.
#
# The background reference is every height tier's monthly climate at every
# site (not just observed heights) — "available but not necessarily
# occupied" conditions, pooled once and shared by every species so all
# niches are scored against the same yardstick.
#
# Rerun this whenever you add observations to OBSERVATIONS_CSV (paths.R;
# data/csv/combined_with_identification.csv by default, overridable via
# CANOPY_OBS_CSV) — it's the only step that needs to change; every
# colonization run downstream picks up the refined niches automatically
# without needing to recompute anything itself.
#
# Requires microenv_<site>[_h<step>].rds to already exist for every site with
# observations (run_microenv.sh) — reads each site's climate once.
#
# Usage: Rscript scripts/02_model/setup/characterize_niches.R [height_step]
#   height_step defaults to 0.25 (the production resolution) — pass 0.1 to
#   use the unsuffixed manifests instead, if you have those for every site.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
HEIGHT_STEP <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
manifest_suffix <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""

VARS <- NICHE_VARS  # c("temp", "relhum", "swdown") — see get_colonization.R

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
sites <- sort(unique(niches$Area_or_Site))

# Monthly mean of variable `v` in climate table `cl` for calendar month `m`,
# via .niche_var_mean() (get_colonization.R) — swdown drops zero/night-time
# readings first (see header note above); temp/relhum use every hourly
# reading in the month.
.monthly_var_mean <- function(cl, v, m) .niche_var_mean(cl[cl$month == m, , drop = FALSE], v)

# ── Pass 1: per-observation ("presence") climate, and every height tier's
# monthly climate at every site ("background", pooled below) ────────────────
obs_clim_rows <- list()
bg_clim_rows  <- list()

for (s in sites) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", s, manifest_suffix))
  if (!file.exists(microenv_path)) {
    message("Skipping ", s, " -- no microenv at ", microenv_path, " (run_microenv.sh first)")
    next
  }
  message("Reading climate for ", s, "...")
  microenv <- readRDS(microenv_path)
  heights  <- microenv_heights(microenv)
  cc       <- build_clim_cache(microenv)

  # Monthly-mean climate per height (12 rows x length(VARS)): each row
  # averages that month's representative warm + cold day (48 hourly values),
  # capturing within-year seasonal variation instead of one collapsed annual
  # mean.
  clim_monthly <- lapply(cc$clim_by_height, function(cl) {
    if (is.null(cl)) return(matrix(NA_real_, nrow = 12, ncol = length(VARS), dimnames = list(NULL, VARS)))
    sapply(VARS, function(v) vapply(1:12, function(m) .monthly_var_mean(cl, v, m), numeric(1)))
  })

  # Every height tier's monthly climate feeds the background pool, whether
  # or not any species was observed at that height — this is "what's
  # available", not "what's occupied".
  for (h_idx in seq_along(clim_monthly)) {
    for (m in 1:12) {
      row <- clim_monthly[[h_idx]][m, ]
      if (all(!is.na(row))) bg_clim_rows[[length(bg_clim_rows) + 1]] <- row
    }
  }

  obs_site <- niches[niches$Area_or_Site == s, ]
  for (i in seq_len(nrow(obs_site))) {
    h    <- obs_site$Height_m[i]
    sp   <- obs_site$FinalID[i]
    h_idx <- which.min(abs(heights - h))
    # Augment this observation with all 12 within-year monthly conditions at
    # its recorded height, instead of just that height's annual mean.
    for (m in 1:12) {
      obs_clim_rows[[length(obs_clim_rows) + 1]] <-
        c(species = sp, clim_monthly[[h_idx]][m, ])
    }
  }
  message(sprintf("  %s: %d observations, %d height tiers", s, nrow(obs_site), length(heights)))
}

if (length(obs_clim_rows) == 0) stop("No usable observations across any site with a microenv file.")
if (length(bg_clim_rows) == 0) stop("No usable background climate across any site with a microenv file.")

obs_df <- as.data.frame(do.call(rbind, obs_clim_rows), stringsAsFactors = FALSE)
for (v in VARS) obs_df[[v]] <- as.numeric(obs_df[[v]])

bg_df <- as.data.frame(do.call(rbind, bg_clim_rows), stringsAsFactors = FALSE)
message(sprintf("Background pool: %d voxel-months across %d site(s)", nrow(bg_df), length(sites)))

bg_density <- build_background_density(bg_df, VARS)
bg_out_path <- file.path(PROCESSED_DIR, "niche_background_density.rds")
saveRDS(bg_density, bg_out_path)
message("Saved background density to ", bg_out_path)

# ── Pass 2: pooled niche model per species, across every site it was
# observed at ────────────────────────────────────────────────────────────────
species_ids <- sort(unique(obs_df$species))
niche_cache <- lapply(species_ids, function(sp) {
  clim_vals <- as.matrix(obs_df[obs_df$species == sp, VARS])
  n_obs     <- nrow(clim_vals) / 12L  # each field observation contributes 12 monthly rows
  n_sites   <- length(unique(niches$Area_or_Site[niches$FinalID == sp]))
  message(sprintf("%-30s %3d observations (%4d monthly data points) across %d site(s)",
                  sp, n_obs, nrow(clim_vals), n_sites))
  niche_density_model(clim_vals, bg_density, VARS)
})
names(niche_cache) <- species_ids

out_path <- file.path(PROCESSED_DIR, "species_niches.rds")
saveRDS(niche_cache, out_path)
message(sprintf("\nSaved %d species niches to %s", length(niche_cache), out_path))
