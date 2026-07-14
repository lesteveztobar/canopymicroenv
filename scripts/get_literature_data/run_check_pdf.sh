#!/bin/bash
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_check_pdf_%j.out

# Usage — anything after the script name is forwarded as-is to check_pdf.py:
#   sbatch run_check_pdf.sh                                              # no args: batch over the whole observation_extraction folder
#   sbatch run_check_pdf.sh data/literature/observation_extraction/foo.pdf     # a single paper
#   sbatch run_check_pdf.sh data/literature/observation_extraction/*.pdf --source "Author et al. 2024" --site "Some Site"
#   sbatch run_check_pdf.sh data/literature/observation_extraction/foo.pdf --dump-candidates

module purge
module load Miniforge3/24.1.2-0

cd /home/$USER/canopymicroenv

if [ "$#" -eq 0 ]; then
    set -- data/literature/observation_extraction/*.pdf --outdir data/literature_staging
fi

python3 scripts/get_literature_data/check_pdf.py "$@"
