# plot_bestfit_3d.R
# 3D scatter plot of adult epiphyte abundance at the last timestep
# of the best-fit parameter combination.
#
# Requires best_run to be in the environment (from simple_experiments.R)
# or re-runs it from scratch if not found.
#
# Output: output/bestfit_3d.png
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

library(scatterplot3d)
library(RColorBrewer)

# ── Re-run best fit if not in environment ─────────────────────────────────────
if (!exists("best_run")) {
  source("scripts/simple_colonization.R")
  best_params <- params
  best_params$establishment_prob <- 0.10
  best_params$repro_rate         <- 500
  best_params$survival_A         <- 0.95
  best_params$maxDisp            <- 4
  best_params$n_founders         <- 100
  best_params$carCap             <- 40
  best_params$timesteps          <- 30
  cat("Re-running best fit...\n")
  best_run <- run_simple_colonization(best_params, seed = 42)
}

p       <- best_run$params
T_last  <- p$timesteps

# ── Extract last timestep ─────────────────────────────────────────────────────
abundA_last <- best_run$abundA[, , , T_last]
zone_arr    <- best_run$zone

# Build data frame of occupied voxels only
idx <- which(abundA_last > 0, arr.ind = TRUE)
df  <- data.frame(
  x    = idx[, 1],
  y    = idx[, 2],
  z    = idx[, 3],
  n    = abundA_last[idx],
  zone = zone_arr[idx]
)

# Convert voxel indices to metres
df$x_m <- (df$x - 0.5) * p$resolution
df$y_m <- (df$y - 0.5) * p$resolution
df$z_m <- (df$z - 0.5) * (p$max_height / p$zDim)

# ── Colour by zone (Johansson) ────────────────────────────────────────────────
zone_cols <- c(
  "1" = "#8B4513",   # trunk base   — brown
  "2" = "#A0522D",   # lower trunk  — sienna
  "3" = "#6B8E23",   # upper trunk  — olive
  "4" = "#228B22",   # inner crown  — forest green
  "5" = "#32CD32"    # outer crown  — lime green
)
df$col <- zone_cols[as.character(df$zone)]

# Point size proportional to abundance (capped so plot stays readable)
df$pt_size <- pmin(df$n / max(df$n) * 6 + 1, 7)

# ── Run a "normal" realistic simulation for comparison ────────────────────────
# Use p_est=0.05, repro=150 — above the extinction threshold but not saturating
source("scripts/simple_colonization.R")
base_p <- params
base_p$establishment_prob <- 0.05
base_p$repro_rate         <- 150
base_p$carCap             <- 8
base_p$n_founders         <- 30
base_p$timesteps          <- 30
cat("Running baseline simulation...\n")
base_run <- suppressMessages(capture.output(
  base_out <- run_simple_colonization(base_p, seed = 7),
  type = "output"
))

make_df <- function(run_out, timestep) {
  abA  <- run_out$abundA[, , , timestep]
  zArr <- run_out$zone
  idx  <- which(abA > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)
  data.frame(
    x_m  = (idx[,1] - 0.5) * run_out$params$resolution,
    y_m  = (idx[,2] - 0.5) * run_out$params$resolution,
    z_m  = (idx[,3] - 0.5) * (run_out$params$max_height / run_out$params$zDim),
    n    = abA[idx],
    zone = zArr[idx]
  )
}

df_base <- make_df(base_out, base_p$timesteps)
df_best <- make_df(best_run, T_last)

# Subsample best-fit (very dense); keep all baseline points
set.seed(42)
df_best_plot <- df_best[sample(nrow(df_best), min(nrow(df_best), 5000)), ]
df_base_plot <- df_base   # usually sparse enough to keep all

add_cols <- function(df) {
  df$col <- zone_cols[as.character(df$zone)]
  df
}
df_base_plot <- add_cols(df_base_plot)
df_best_plot <- add_cols(df_best_plot)

# ── Two-panel PNG ─────────────────────────────────────────────────────────────
out_path <- "output/bestfit_3d.png"
png(out_path, width = 3200, height = 1400, res = 180)

layout(matrix(c(1, 2, 3), nrow = 1), widths = c(10, 10, 3))

plot_panel <- function(df_p, title_str, cex_sym = 0.55) {
  par(mar = c(2, 2, 3, 1))
  scatterplot3d(
    x = df_p$x_m, y = df_p$y_m, z = df_p$z_m,
    color       = df_p$col,
    pch         = 16,
    cex.symbols = cex_sym,
    xlab        = "East–West (m)",
    ylab        = "South–North (m)",
    zlab        = "Height (m)",
    main        = title_str,
    angle       = 35,
    scale.y     = 0.6,
    grid        = TRUE,
    box         = FALSE,
    col.axis    = "grey40",
    col.grid    = "grey88",
    col.lab     = "grey20",
    cex.axis    = 0.8,
    cex.lab     = 0.9
  )
}

plot_panel(df_base_plot,
  sprintf("Realistic params  (p_e=0.05, λ=150, carCap=8)\nyear %d — %d adults",
          base_p$timesteps, sum(df_base$n)))

plot_panel(df_best_plot,
  sprintf("Best-fit params  (p_e=0.10, λ=500, carCap=40)\nyear %d — %d adults  [5k sample]",
          T_last, sum(df_best$n)))

# ── Legend panel ──────────────────────────────────────────────────────────────
par(mar = c(2, 0, 3, 1))
plot.new()
legend(
  "center",
  legend    = c("Zone 1 – trunk base", "Zone 2 – lower trunk",
                "Zone 3 – upper trunk", "Zone 4 – inner crown",
                "Zone 5 – outer crown"),
  col       = unname(zone_cols),
  pch       = 16,
  pt.cex    = 1.5,
  bty       = "n",
  cex       = 0.95,
  title     = "Johansson zone",
  title.col = "grey20",
  title.font = 2
)

dev.off()
cat("Saved:", out_path, "\n")
