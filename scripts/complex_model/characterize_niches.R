# characterize_niches.R
# Precomputes each species' realized climate niche (temp/RH) pooled across
# EVERY site where it was observed, not just the one site.rds file happens
# to be modeling — a species with 2-3 records at one site may have several
# more at others. Saves data/processed/species_niches.rds, loaded by
# init_colonization() (get_colonization.R) in preference to the this-site-
# only get_niche() fallback.
#
# Rerun this whenever you add observations to data/csv/combinedv3.csv — it's
# the only step that needs to change; every colonization run downstream picks
# up the refined niches automatically without needing to recompute anything
# itself.
#
# Requires microenv_<site>[_h<step>].rds to already exist for every site with
# observations (run_microenv.sh) — reads each site's climate once.
#
# Usage: Rscript scripts/complex_model/characterize_niches.R [height_step]
#   height_step defaults to 0.25 (the production resolution) — pass 0.1 to
#   use the unsuffixed manifests instead, if you have those for every site.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
HEIGHT_STEP <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
manifest_suffix <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""

VARS <- c("temp", "relhum")

niches <- read.csv("data/csv/combinedv3.csv")
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
sites <- sort(unique(niches$Area_or_Site))

# ── Pass 1: per-observation climate, and each site's own vertical climate
# range (pooled afterwards into one shared reference scale for the
# minimum-width floor — see niche_geometry() in get_colonization.R). ─────────
obs_clim_rows  <- list()
site_ranges    <- list()

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

  clim_scalars <- do.call(rbind, lapply(cc$clim_by_height, function(cl) {
    if (is.null(cl)) return(setNames(rep(NA_real_, length(VARS)), VARS))
    vapply(VARS, function(v) mean(cl[[v]], na.rm = TRUE), numeric(1))
  }))
  site_ranges[[s]] <- apply(clim_scalars, 2, function(x) diff(range(x, na.rm = TRUE)))

  obs_site <- niches[niches$Area_or_Site == s, ]
  for (i in seq_len(nrow(obs_site))) {
    h   <- obs_site$Height_m[i]
    sp  <- obs_site$FinalID[i]
    cv  <- clim_scalars[which.min(abs(heights - h)), , drop = TRUE]
    obs_clim_rows[[length(obs_clim_rows) + 1]] <- c(species = sp, cv)
  }
  message(sprintf("  %s: %d observations, %d height tiers", s, nrow(obs_site), length(heights)))
}

if (length(obs_clim_rows) == 0) stop("No usable observations across any site with a microenv file.")

obs_df <- as.data.frame(do.call(rbind, obs_clim_rows), stringsAsFactors = FALSE)
for (v in VARS) obs_df[[v]] <- as.numeric(obs_df[[v]])

# Shared reference scale for the minimum-width floor: the largest per-height
# climate range seen at any single site, so the floor stays meaningful (not
# swamped) regardless of how many sites happen to be included.
min_width_frac <- 0.10
site_range_mat <- do.call(rbind, site_ranges)
ref_range      <- apply(site_range_mat, 2, max, na.rm = TRUE)
min_width      <- min_width_frac * ref_range
message(sprintf("Reference vertical climate range (max across sites): temp=%.2f C, RH=%.2f %%",
                ref_range["temp"], ref_range["relhum"]))

# ── Pass 2: pooled niche geometry per species, across every site it was
# observed at ────────────────────────────────────────────────────────────────
species_ids <- sort(unique(obs_df$species))
niche_cache <- lapply(species_ids, function(sp) {
  clim_vals <- as.matrix(obs_df[obs_df$species == sp, VARS])
  n_obs     <- nrow(clim_vals)
  n_sites   <- length(unique(niches$Area_or_Site[niches$FinalID == sp]))
  message(sprintf("%-30s %3d observations across %d site(s)", sp, n_obs, n_sites))
  niche_geometry(clim_vals, min_width)
})
names(niche_cache) <- species_ids

out_path <- file.path(PROCESSED_DIR, "species_niches.rds")
saveRDS(niche_cache, out_path)
message(sprintf("\nSaved %d species niches to %s", length(niche_cache), out_path))
