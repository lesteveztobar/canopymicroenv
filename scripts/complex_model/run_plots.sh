#!/bin/bash
# run_plots.sh — generate every project plot from whatever results exist in
# data/processed/ and geojson_to_csv/. Safe to re-run any time.
# Usage: sh scripts/complex_model/run_plots.sh   (or bash scripts/complex_model/run_plots.sh)
# Run from: /home/s38leste_hpc/canopymicroenv/
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"

cd /home/$USER/canopymicroenv
Rscript scripts/complex_model/plot_all.R
