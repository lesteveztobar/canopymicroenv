#!/bin/bash
# height_res_array.sh — regenerate microclimate for ONE site at several
# height-tier spacings, for the height-resolution efficiency diagnostic.
# Each spacing reuses the site's cached weather/DTM/point-model data (those
# don't depend on height spacing) and only redoes the per-height loop, so
# coarser spacings are proportionally cheaper, not free — this still submits
# real, independent cluster jobs.
#
# Usage: sbatch height_res_array.sh <SiteName>
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

SITE=${1:?"Usage: sbatch height_res_array.sh <SiteName>"}
N_MONTHS=12

# 0.1m (55 tiers for Maquipucuna, hObs_max=5.5m) is production resolution and
# should already exist — only test coarser spacings here.
HEIGHT_STEPS=(0.25 0.5 1.0)

for step in "${HEIGHT_STEPS[@]}"; do
  sbatch /home/$USER/canopymicroenv/scripts/run_microenv.sh "$SITE" "$N_MONTHS" "$step"
done
