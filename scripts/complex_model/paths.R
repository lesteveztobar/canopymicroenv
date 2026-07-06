# paths.R
# Project directory constants for canopymicroenv
# All paths derived from BASE_DIR — change only BASE_DIR if project moves
# Per-site subdirectories (era5, dtm, soil, etc.) are built dynamically
# inside the site loop in getmicroenv.R using BASE_DIR as the root
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

BASE_DIR <- "/home/s38leste_hpc/canopymicroenv"

# data
RAW_DIR       <- file.path(BASE_DIR, "data", "raw")
CSV_DIR       <- file.path(BASE_DIR, "data", "csv")
PROCESSED_DIR <- file.path(BASE_DIR, "data", "processed")
PARAMS_DIR    <- file.path(BASE_DIR, "data", "params")

# project
SCRIPTS_DIR   <- file.path(BASE_DIR, "scripts")
OUTPUT_DIR    <- file.path(BASE_DIR, "output")
LOGS_DIR      <- file.path(BASE_DIR, "logs")

# create all directories if they don't exist
# showWarnings = FALSE silently skips dirs that already exist
for (d in c(RAW_DIR, CSV_DIR, PROCESSED_DIR, PARAMS_DIR, SCRIPTS_DIR, OUTPUT_DIR, LOGS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}