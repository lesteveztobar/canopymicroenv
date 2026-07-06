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
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

SITES=("Maquipucuna" "Mashpi" "MindoTarabita" "MiradorMindo" "Yanayacu")
N_MONTHS=${1:-12}      # e.g. sbatch microenv_array.sh 12
HEIGHT_STEP=${2:-0.1}  # e.g. sbatch microenv_array.sh 12 0.5

for i in "${!SITES[@]}"; do
  sbatch /home/$USER/canopymicroenv/scripts/run_microenv.sh "${SITES[$i]}" "$N_MONTHS" "$HEIGHT_STEP"
done
