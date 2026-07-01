#!/bin/bash
# run_plots.sh — generate every project plot from whatever results exist in
# data/processed/ and geojson_to_csv/. Safe to re-run any time.
# Usage: sh scripts/run_plots.sh   (or bash scripts/run_plots.sh)
# Run from: /home/s38leste_hpc/canopymicroenv/
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/canopymicroenv
Rscript scripts/plot_all.R
