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

# Chained sequentially, not submitted independently: these all share one
# ERA5-merged-file cache for $SITE (get_weather(), lib.R)
# with no locking around its extend/rename step, so running them
# concurrently can race on that shared file and silently corrupt it --
# 2026-07-24, exactly this happened to Saloya via the analogous loop in
# run_full_analysis_pipeline.sh (three height-step jobs hit the ERA5
# extend path within the same second; two manifests came back with 0 valid
# height tiers despite the SLURM jobs exiting 0).
prev_step_id=""
for step in "${HEIGHT_STEPS[@]}"; do
  if [ -n "$prev_step_id" ]; then
    dep_arg="--dependency=afterok:${prev_step_id}"
  else
    dep_arg=""
  fi
  id=$(sbatch --parsable $dep_arg --partition=lm_medium --time=24:00:00 --cpus-per-task=32 --mem=500G \
    /home/$USER/canopymicroenv/scripts/01_microclimate/run_microenv.sh "$SITE" "$N_MONTHS" "$step")
  prev_step_id="$id"
done
