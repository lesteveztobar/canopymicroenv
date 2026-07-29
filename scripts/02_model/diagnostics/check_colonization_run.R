# check_colonization_run.R
# Quick first look at a colonization result -- text summary only, no plots
# (for plots, use run_plots.sh / plot_default_colonization_run() /
# plot_factorial_experiment() in plot_functions.R directly). Handles both
# shapes produced by this pipeline:
#   - run_replicated() output (a list): per-replicate final abundance +
#     persistence.
#   - run_experiment()/run_factorial_experiment() output (a data frame, one
#     or more swept parameter columns): persistence rate overall and per
#     swept-parameter level, plus the parameter bounding box among
#     persisting combinations.
#
# Usage: Rscript scripts/02_model/diagnostics/check_colonization_run.R [site] [exp_tag]
#   Defaults to Maquipucuna / best_case_h0.25 if no args given. Only needs
#   paths.R (no sf/GDAL/PROJ/GEOS/UDUNITS), so no module-loading wrapper is
#   required -- run_check_colonization_run.sh still works too, it just loads
#   more than this script actually needs.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")

args      <- commandArgs(trailingOnly = TRUE)
site_name <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
exp_tag   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "best_case_h0.25"

in_path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
if (!file.exists(in_path)) stop("No results at ", in_path)
result <- readRDS(in_path)

if (is.data.frame(result)) {
  cat(sprintf("\n== %s / %s (sweep/factorial) ==\n", site_name, exp_tag))

  fixed_cols <- c("rep", "t", "totalS", "totalJ", "totalA", "total", "extinct")
  swept      <- setdiff(names(result), fixed_cols)
  cat(sprintf("Swept parameter(s): %s\n", paste(swept, collapse = ", ")))

  t_max <- max(result$t)
  final <- result[result$t == t_max, ]
  # One row per unique parameter combination (averaging over reps, if >1).
  form  <- as.formula(paste("cbind(extinct, total) ~", paste(swept, collapse = " + ")))
  combo <- aggregate(form, data = final, FUN = mean)
  combo$persisted <- !as.logical(round(combo$extinct))

  cat(sprintf("%d unique combination(s), final year t=%d\n", nrow(combo), t_max))
  cat(sprintf("\nPersisted: %d/%d (%.1f%%) combinations (adults present in the second half of the run)\n",
              sum(combo$persisted), nrow(combo), 100 * mean(combo$persisted)))

  cat("\nPersistence rate by parameter level (marginal — averaged over the other swept parameters):\n")
  for (v in swept) {
    rate <- tapply(combo$persisted, combo[[v]], mean)
    cat(sprintf("  %s:\n", v))
    print(round(100 * rate, 1))
  }

  persisting <- combo[combo$persisted, ]
  if (nrow(persisting) > 0) {
    cat(sprintf("\nAmong persisting combinations: final total abundance %.0f–%.0f (median %.0f)\n",
                min(persisting$total), max(persisting$total), median(persisting$total)))
    cat("\nParameter bounding box among persisting combinations:\n")
    for (v in swept) {
      cat(sprintf("  %s: %s to %s\n", v,
                  format(min(persisting[[v]])), format(max(persisting[[v]]))))
    }
  } else {
    cat("\nNo combinations persisted.\n")
  }
  quit(save = "no")
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
