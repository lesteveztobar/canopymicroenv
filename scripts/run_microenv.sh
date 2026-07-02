#!/bin/bash
# run_microenv.sh — run microclimate model for a single site
# Called by microenv_array.sh; do not submit directly.
#
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
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
N_MONTHS=${2:-12}  # months of ERA5 to fetch backwards from tme_end (default: 12, matches run_microclimate_site.R's own default and the rest of the pipeline's annual climate resolution)
HEIGHT_STEP=${3:-0.1}  # m between height tiers (default: 0.1, production resolution)

# Allocate (or find existing) Lustre workspace for per-height temp files.
# ws_allocate is idempotent: calling it again returns the same path.
export CANOPY_SCRATCH=$(ws_allocate canopymicroenv 90)
echo "Scratch workspace: $CANOPY_SCRATCH"

cd /home/$USER/canopymicroenv
$CANOPY_PYTHON -c "import ee; print('ee import OK')" 2>&1
Rscript scripts/run_microclimate_site.R "$SITE" "$N_MONTHS" "$HEIGHT_STEP"
