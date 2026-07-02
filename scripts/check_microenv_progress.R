# check_microenv_progress.R
# Reports, per site (and optionally per height-tier spacing), how many of the
# expected per-height climate files have been computed on scratch, and
# whether the final microenv_<site>[_h<step>].rds manifest has been saved
# yet. Read-only — safe to run any time.
#
# Usage:
#   export CANOPY_SCRATCH=$(ws_find canopymicroenv)
#   Rscript scripts/check_microenv_progress.R                       # all sites, 0.1m
#   Rscript scripts/check_microenv_progress.R Maquipucuna 0.1,0.25,0.5,1.0  # one site, several height steps
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/paths.R")

args <- commandArgs(trailingOnly = TRUE)
SITES <- if (length(args) >= 1) args[1] else
  c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu")
HEIGHT_STEPS <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else 0.1

scratch <- Sys.getenv("CANOPY_SCRATCH", unset = "")
if (nchar(scratch) == 0 || !dir.exists(scratch)) {
  stop("CANOPY_SCRATCH not set to an existing directory.\n",
       "Run: export CANOPY_SCRATCH=$(ws_find canopymicroenv)")
}

niches <- read.csv("data/csv/combinedv3.csv")
niches <- niches[!is.na(niches$Height_m), ]

cat(sprintf("%-15s %-6s %-14s %-10s %s\n", "Site", "Step", "Heights", "Manifest", "Height dir"))
cat(strrep("-", 80), "\n")

for (s in SITES) {
  hmax <- suppressWarnings(max(niches$Height_m[niches$Area_or_Site == s], na.rm = TRUE))
  for (step in HEIGHT_STEPS) {
    suffix <- if (step != 0.1) sprintf("_h%.2f", step) else ""
    height_dir <- file.path(scratch, sprintf("microenv_%s%s_heights", s, suffix))
    manifest_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", s, suffix))
    manifest_status <- if (file.exists(manifest_path)) "done" else "pending"

    if (!is.finite(hmax)) {
      cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", s, step, "no obs in CSV", manifest_status, height_dir))
      next
    }
    n_expected <- length(seq(0.1, hmax, by = step))

    if (!dir.exists(height_dir)) {
      cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", s, step, sprintf("0/%d", n_expected), manifest_status,
                  paste0(height_dir, " (not created yet)")))
      next
    }
    n_done <- length(list.files(height_dir, pattern = "\\.rds$"))
    cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", s, step, sprintf("%d/%d", n_done, n_expected),
                manifest_status, height_dir))
  }
}
