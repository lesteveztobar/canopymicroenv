#!/bin/bash
# microenv_array.sh — submit one job per site to the Marvin cluster
# Usage: sbatch microenv_array.sh [n_months] [height_step]
#   e.g. sbatch microenv_array.sh 12 0.5
# Run from: /home/s38leste_hpc/canopymicroenv/
#
# NOTE: run_microclimate_site.R now uses the canopy-height raster (vhgt.tif)
# as the height ceiling, not the tallest recorded epiphyte observation — the
# ceiling can be much taller than before (e.g. Maquipucuna: ~5.5m -> ~15m).
# At height_step=0.1 that's ~3x more height tiers, and each tier is a real
# microclimf run — consider a coarser step (0.25-0.5m) for practicality.
# Existing unsuffixed microenv_<site>.rds manifests (from the old, truncated
# ceiling) must be deleted or moved first, or the script will see them and
# exit immediately without regenerating anything — see check_microenv_progress.R.
#
# NOTE (2026-07-22): MiradorMindo retired -- its raw GeoJSON export actually
# bundled observations from 3 separate locations. Replaced here by
# LaElenita/MindoMirador/Saloya (see scripts/00_data_conversion/rebuild_combined_csv.py). This
# only takes effect once combinedv3.csv itself has been updated with those
# site names (make_sites() in lib.R groups by Area_or_Site in
# combinedv3.csv, not by this array) -- promote the split rows from
# combined.csv into combinedv3.csv (with species IDs) before submitting.
# Every other pipeline script derives its site list dynamically from the
# observations CSV, so the new sites are picked up automatically. The one
# exception is scripts/legacy/allsites.R, which still hardcodes the old
# 5-site list incl. "MiradorMindo" -- it predates the site split entirely
# and is not part of the active pipeline (see its own header comment).
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

SITES=("Maquipucuna" "Mashpi" "MindoTarabita" "LaElenita" "MindoMirador" "Saloya" "Yanayacu")
N_MONTHS=${1:-12}      # e.g. sbatch microenv_array.sh 12
HEIGHT_STEP=${2:-0.1}  # e.g. sbatch microenv_array.sh 12 0.5

for i in "${!SITES[@]}"; do
  sbatch /home/$USER/canopymicroenv/scripts/01_microclimate/run_microenv.sh "${SITES[$i]}" "$N_MONTHS" "$HEIGHT_STEP"
done
