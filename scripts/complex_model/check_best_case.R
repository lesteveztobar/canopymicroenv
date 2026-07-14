# check_best_case.R
# Quick first look at a `run_replicated()` colonization result: a text
# summary (did it persist? final abundance per replicate?) followed by the
# full plot suite via plot_default_colonization_run().
#
# Usage: sh scripts/complex_model/run_check_best_case.sh [site] [exp_tag]
#   Defaults to Maquipucuna / best_case_h0.25 (the first post-p_germ-fix
#   run — see report/methods.tex, Persistence validation).
#   Don't invoke this file directly with bare `Rscript` -- it sources
#   plot_functions.R, which loads sf (needs GDAL/PROJ/GEOS/UDUNITS on
#   LD_LIBRARY_PATH); the .sh wrapper sets that up. See run_plots.sh for
#   the same pattern used elsewhere in this pipeline.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/complex_model/plot_functions.R")

args      <- commandArgs(trailingOnly = TRUE)
site_name <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
exp_tag   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "best_case_h0.25"

in_path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
if (!file.exists(in_path)) stop("No results at ", in_path)
result <- readRDS(in_path)

if (is.data.frame(result)) {
  stop("This is a sweep/factorial result, not a run_replicated() result -- ",
       "use plot_colonization_experiment() or plot_factorial_experiment() instead.")
}

cat(sprintf("\n== %s / %s ==\n", site_name, exp_tag))
cat(sprintf("%d replicate(s) requested, %d succeeded\n",
            length(result$runs), sum(!vapply(result$runs, is.null, logical(1)))))

# Per-replicate outcome: final-year abundance and the same "extinct" test
# run_replicated() uses (adult total == 0 for the whole second half of the run).
final_t <- max(result$summary$t)
final <- result$summary[result$summary$t == final_t, ]
final <- final[order(final$rep), c("rep", "totalS", "totalJ", "totalA", "total", "extinct")]
cat(sprintf("\nFinal-year (t=%d) abundance by replicate:\n", final_t))
print(final, row.names = FALSE)

cat(sprintf("\nPersisted: %d/%d replicates (adults present at some point in the second half of the run)\n",
            sum(!final$extinct), nrow(final)))

# Trajectory summary across all years, averaged over replicates -- a quick
# read on trend (growing/declining/flat) without opening a plot.
traj <- aggregate(cbind(totalS, totalJ, totalA, total) ~ t, data = result$summary, FUN = mean)
cat("\nMean trajectory across replicates (selected years):\n")
print(traj[traj$t %in% unique(c(1, round(final_t / 2), final_t)), ], row.names = FALSE)

# Full plot suite: per-replicate abundance PNGs + 3D snapshot/animation for
# replicate 1. Saved under OUTPUT_DIR (see plot_default_colonization_run()
# in plot_functions.R for exact filenames).
cat("\nGenerating plots...\n")
plot_default_colonization_run(site_name, exp_tag = exp_tag)
cat("\nDone -- see", OUTPUT_DIR, "\n")
