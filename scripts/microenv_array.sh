#!/bin/bash
# microenv_array.sh — submit one job per site to the Marvin cluster
# Usage: sbatch microenv_array.sh
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

SITES=("Maquipucuna" "Mashpi" "MindoTarabita" "MiradorMindo" "Yanayacu")
N_MONTHS=${1:-1}  # e.g. sbatch microenv_array.sh 12

for i in "${!SITES[@]}"; do
  sbatch /home/$USER/canopymicroenv/scripts/run_microenv.sh "${SITES[$i]}" "$N_MONTHS"
done
