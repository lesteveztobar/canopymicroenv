#!/bin/bash
# characterize_niches.sh — precompute each species' realized climate niche
# pooled across every site it was observed at, saving
# data/processed/species_niches.rds for init_colonization() to use.
#
# Rerun this whenever you add observations to data/csv/combinedv3.csv — every
# colonization run downstream picks up the refined niches automatically.
# Requires microenv_<site>[_h<step>].rds for every site with observations
# already generated (run_microenv.sh / microenv_array.sh).
#
# Usage: sbatch scripts/complex_model/characterize_niches.sh [height_step]
#   height_step defaults to 0.25 (the production resolution).
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/canopymicroenv
Rscript scripts/complex_model/characterize_niches.R "${1:-0.25}"
