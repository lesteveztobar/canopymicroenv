# run.R — single entry point for the simple (no-microclimate) standalone model.
# Dispatches to the underlying scripts, which remain separate files since each
# is a distinct analysis with its own top-level driver logic (baseline run +
# plots, a battery of sensitivity experiments, or a 2D extinction-threshold
# scan) rather than a library function.
#
# Usage: Rscript scripts/simple_model/run.R <baseline|experiments|heatmap>
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
HERE <- "scripts/simple_model"
args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: Rscript scripts/simple_model/run.R <subcommand>\n\n",
    "Subcommands:\n",
    "  baseline     Run simple_colonization.R's default simulation + plots.\n",
    "  experiments  Run simple_experiments.R: OAT sensitivity sweeps,\n",
    "               animations, and the combined overview figure.\n",
    "  heatmap      Run run_extinction_heatmap.R: the p_est x repro_rate\n",
    "               extinction-threshold scan and heatmap figure.\n",
    sep = ""
  )
}

if (length(args) < 1) { usage(); quit(status = 1) }

switch(args[1],
  baseline = {
    # simple_colonization_no_run left unset (default FALSE) so sourcing it
    # runs its own default simulation + plots, as it does when run standalone.
    source(file.path(HERE, "simple_colonization.R"))
  },
  experiments = {
    source(file.path(HERE, "simple_experiments.R"))
  },
  heatmap = {
    source(file.path(HERE, "run_extinction_heatmap.R"))
  },
  { usage(); quit(status = 1) }
)
