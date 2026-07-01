#!/bin/bash
# batch_exp.sh — submit one sensitivity-experiment job per site to the Marvin cluster
# Each job sweeps one parameter (5 values, handled inside run_colonization_onesite.R)
# across all timesteps for one site. SITES x EXP = 5 x 6 = 30 jobs.
# Usage: sbatch batch_exp.sh
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

Rscript scripts/make_params.R

SITES=("Maquipucuna" "Mashpi" "MindoTarabita" "MiradorMindo" "Yanayacu")

# EXP and PARAMS are paired by index — PARAMS[k] is the RDS make_params.R
# built for the experiment EXP[k]. Keep these two arrays in sync.
EXP=("pollination_success" "adult_survival_intercept" "germination_probability" "reproduction_cost" "climate_sensitivity_rh" "precipitation_sensitivity")
PARAMS=("p_poll.rds" "beta0A.rds" "p_germ.rds" "cost_repro.rds" "beta_rh.rds" "beta_precip.rds")

PARAMS_DIR="/home/$USER/canopymicroenv/data/params"

for i in "${!SITES[@]}"; do
    for k in "${!EXP[@]}"; do
        sbatch /home/$USER/canopymicroenv/scripts/run_colonization.sh \
            "${SITES[$i]}" \
            "$PARAMS_DIR/${PARAMS[$k]}" \
            "${EXP[$k]}"
    done
done
