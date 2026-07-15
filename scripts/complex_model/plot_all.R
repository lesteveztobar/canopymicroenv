# plot_all.R
# Generates every project plot from whatever results currently exist under
# data/processed/ and geojson_to_csv/. Safe to re-run at any point in the
# pipeline — each function in plot_functions.R skips (with a message) if its
# inputs aren't there yet, instead of erroring.
#
# Usage: sh scripts/complex_model/run_plots.sh   (wraps: Rscript scripts/complex_model/plot_all.R)
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/complex_model/plot_functions.R")

SITES <- c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu")
EXP   <- names(EXP_PARAM_MAP)

cat("== Field site map ==\n")
plot_site_map()

cat("\n== Temperature profiles ==\n")
for (site in SITES) plot_temperature_profile(site)

cat("\n== Niche suitability (per-axis + before/after normalization) ==\n")
for (site in SITES) plot_niche_suitability(site)

cat("\n== Best-fit 3D comparison (simple model) ==\n")
plot_bestfit_3d_comparison()

cat("\n== Colonization sensitivity experiments ==\n")
for (site in SITES) {
  for (tag in EXP) {
    plot_colonization_experiment(site, tag)
  }
}

cat("\n== Reproduction factorial (p_poll x p_germ x p_s1 x n_founders) ==\n")
# reproduction_factorial_v3: levels bracketed realistic.rds -> best_case.rds
# (see report/methods.tex, Full factorial experiment). Only Maquipucuna has
# a v3 result so far -- add more sites here as they land in data/processed/.
plot_factorial_experiment("Maquipucuna", "reproduction_factorial_v3_h0.25")

cat("\n== Default (unswept) colonization runs ==\n")
for (site in SITES) {
  plot_default_colonization_run(site)
}

cat("\n== Persistence validation runs (best-case / realistic) ==\n")
# Only Maquipucuna has post-p_germ-fix runs so far (see report/methods.tex,
# Persistence validation) -- add more sites here as their own
# best_case_h0.25 / realistic_h0.25 results land in data/processed/.
plot_default_colonization_run("Maquipucuna", exp_tag = "best_case_h0.25")
plot_default_colonization_run("Maquipucuna", exp_tag = "realistic_h0.25")

cat("\nAll done. Plots (where inputs existed) are in", OUTPUT_DIR, "\n")
