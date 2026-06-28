#!/bin/bash
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_test_%j.out

module purge
module load GCCcore/13.3.0
module load R
module load Miniforge3/24.1.2-0
export LD_LIBRARY_PATH="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib:/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
export CANOPY_PYTHON="/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"

echo "--- Python test ---"
$CANOPY_PYTHON -c "import ee; import numpy; print('Python OK — ee', ee.__version__, '| numpy', numpy.__version__)" 2>&1

echo "--- R test ---"
cd /home/$USER/canopymicroenv
Rscript scripts/test_env.R