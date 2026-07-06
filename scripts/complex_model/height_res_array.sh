#!/bin/bash
# height_res_array.sh — regenerate microclimate for ONE site at several
# height-tier spacings, for the height-resolution efficiency diagnostic.
# Each spacing reuses the site's cached weather/DTM/point-model data (those
# don't depend on height spacing) and only redoes the per-height loop, so
# coarser spacings are proportionally cheaper, not free — this still submits
# real, independent cluster jobs.
#
# Overnight setup: 0.1/0.5/1.0m for Maquipucuna, on lm_medium with more
# cpus/mem than run_microenv.sh's own defaults, so the height loop finishes
# well within the 24h window instead of crawling on 4 cores. 0.25m is not
# included here since that's the production run already in progress
# separately (microenv_array.sh) — run_resolution_diagnostics.R will compare
# all four (0.1, 0.25, 0.5, 1.0) once everything's done.
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
HEIGHT_STEPS=(0.1 0.5 1.0)

for step in "${HEIGHT_STEPS[@]}"; do
  sbatch --partition=lm_medium --time=24:00:00 --cpus-per-task=32 --mem=500G \
    /home/$USER/canopymicroenv/scripts/run_microenv.sh "$SITE" "$N_MONTHS" "$step"
done
