#!/bin/bash
# hpc_job.sh — run microclimate model for a single site
# Called by hpc_array.sh; do not submit directly.
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
export CANOPY_PYTHON="/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"

SITE=$1      # e.g. "Maquipucuna"
N_MONTHS=${2:-1}  # months of ERA5 to fetch backwards from tme_end (default: 1)

cd /home/$USER/canopymicroenv
$CANOPY_PYTHON -c "import ee; print('ee import OK')" 2>&1
Rscript scripts/run_microclimate_site.R "$SITE" "$N_MONTHS"
