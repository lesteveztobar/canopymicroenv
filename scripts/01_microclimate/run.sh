#!/bin/bash
# run.sh — single entry point for the microclimate stage.
#
# Dispatches to the underlying scripts, which remain separate files on
# purpose: run_microenv.sh and microenv_array.sh carry their own #SBATCH
# resource directives (partition/time/mem/cpus), which SLURM requires at the
# top of the exact file passed to `sbatch` -- they cannot be folded into a
# shared dispatcher without losing that. Everything else here (the plain R
# entry points) has no SBATCH header and can be dispatched directly.
#
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
HERE="scripts/01_microclimate"

usage() {
  cat <<'EOF'
Usage: scripts/01_microclimate/run.sh <subcommand> [args...]

Direct (no SLURM submission -- assumes the run environment, e.g. an
interactive SLURM allocation or a node with modules/conda already loaded):
  interactive
      Run run_microclimate.R: all sites sequentially, interactive
      (needs a console for the one-time Earth Engine auth step).

  site SITE [N_MONTHS] [HEIGHT_STEP]
      Run run_microclimate_site.R directly for one site (default
      N_MONTHS=12, HEIGHT_STEP=0.1). This is what run_microenv.sh invokes
      under SLURM -- use that (via `submit-site` below) for production runs.

  progress [SITES] [HEIGHT_STEPS]
      Run check_microenv_progress.R (read-only). SITES and HEIGHT_STEPS are
      comma-separated, e.g. `progress Maquipucuna,Mashpi 0.1,0.25,0.5,1.0`.

  fix-dtm
      Run regenerate_missing_dtm.R: re-fetch any site's missing dtm.tif.

SLURM submission (pass-through to the SBATCH-bearing scripts):
  submit-site SITE [N_MONTHS] [HEIGHT_STEP]
      sbatch run_microenv.sh -- single-site production run.

  submit-array [N_MONTHS] [HEIGHT_STEP]
      sbatch microenv_array.sh -- submits one run_microenv.sh job per site.

  help
      Show this message.
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift

case "$CMD" in
  interactive)   Rscript "$HERE/run_microclimate.R" "$@" ;;
  site)          Rscript "$HERE/run_microclimate_site.R" "$@" ;;
  progress)      Rscript "$HERE/check_microenv_progress.R" "$@" ;;
  fix-dtm)       Rscript "$HERE/regenerate_missing_dtm.R" "$@" ;;
  submit-site)   sbatch "$HERE/run_microenv.sh" "$@" ;;
  submit-array)  sbatch "$HERE/microenv_array.sh" "$@" ;;
  help|-h|--help) usage ;;
  *) echo "Unknown subcommand: $CMD" >&2; usage; exit 1 ;;
esac
