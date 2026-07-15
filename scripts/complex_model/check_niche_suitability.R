# check_niche_widths.R
# Diagnostic: shows each species' realized climate niche (temp/RH range) and
# how many of the site's height tiers fall inside it, at a given niche_pad.
# Useful for sanity-checking get_niche()/niche_match() without running a
# full colonization job — read-only, no cluster job needed.
#
# Usage: Rscript scripts/complex_model/check_niche_widths.R <site> [pad1,pad2,...]
#   e.g.: Rscript scripts/complex_model/check_niche_widths.R Maquipucuna 0,0.05
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
site_name <- if (length(args) >= 1) args[1] else "Maquipucuna"
pads      <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else c(0, 0.05)

niches <- read.csv("data/csv/combinedv3.csv")
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
site_obs <- niches[niches$Area_or_Site == site_name, ]

microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s.rds", site_name))
if (!file.exists(microenv_path)) stop("No microenv for ", site_name, " at ", microenv_path)
microenv <- readRDS(microenv_path)
heights  <- microenv_heights(microenv)

message("Building climate cache (reads all ", length(heights), " height files once)...")
cc <- build_clim_cache(microenv)
message("Done.")

for (pad in pads) {
  cat(sprintf("\n=== %s | niche_pad = %.3f ===\n", site_name, pad))
  niche_list <- get_niche(site_obs, heights, cc$clim_by_height, pad = pad)
  for (sp in names(niche_list)) {
    n <- niche_list[[sp]]
    if (is.null(n)) { cat(sp, ": NULL niche (no usable observations)\n"); next }
    n_match <- sum(sapply(cc$clim_by_height, function(cl) {
      if (is.null(cl)) return(FALSE)
      niche_match(list(temp = mean(cl$temp, na.rm = TRUE),
                       relhum = mean(cl$relhum, na.rm = TRUE)), n) >= 1
    }))
    cat(sprintf("%-25s temp[%.2f,%.2f] RH[%.2f,%.2f] -- %d/%d heights match\n",
                sp, n$lo["temp"], n$hi["temp"], n$lo["relhum"], n$hi["relhum"],
                n_match, length(heights)))
  }
}
