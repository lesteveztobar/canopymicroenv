# simple_experiments.R
# Sensitivity experiments for the simple 3D colonization model.
# Sources simple_colonization.R (functions only — run guard suppresses auto-exec).
# Each experiment varies one parameter at a time and compares trajectories.
# Includes dispersal animations (side-view GIF via gganimate, 3D via plotly).
#
# Usage: source("scripts/simple_model/simple_experiments.R")
# ─────────────────────────────────────────────────────────────────────────────
options(simple_colonization_no_run = TRUE)
source("scripts/simple_model/simple_colonization.R")
library(ggplot2)
library(patchwork)
library(parallel)

N_CORES <- max(1L, detectCores() - 1L)  # leave one core free
cat(sprintf("Using %d parallel cores\n", N_CORES))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Baseline params ───────────────────────────────────────────────────────────
# repro_rate = effective seeds dispersed per adult per year — already collapses
# flowering rate, fruit number, and pollination, so ~100 is realistic.
# establishment_prob × suitability ≈ P(seed → seedling); with repro_rate=100
# and suit~0.7, setting this to 0.02 gives ~1–3 recruits/adult/yr at good sites.

base_params <- params
base_params$repro_rate         <- 100
base_params$establishment_prob <- 0.02
base_params$survival_A         <- 0.88
base_params$n_founders         <- 80
base_params$timesteps          <- 50

# ── Experiment runner ─────────────────────────────────────────────────────────

run_experiment <- function(param_name, values, base = base_params, n_reps = 3) {
  # Build one job per (value × rep) combination, run all in parallel
  jobs <- expand.grid(val = values, rep = seq_len(n_reps), stringsAsFactors = FALSE)

  rows <- mclapply(seq_len(nrow(jobs)), function(i) {
    val <- jobs$val[i]; rep <- jobs$rep[i]
    p <- base; p[[param_name]] <- val
    capture.output(r <- run_simple_colonization(p, seed = rep * 7 + val * 13))
    T <- p$timesteps
    data.frame(
      param_value = as.character(val),
      rep         = rep,
      t           = 1:T,
      totalA      = r$totalA,
      totalS      = r$totalS,
      total       = r$totalA + r$totalS,
      extinct     = all(r$totalA[(T %/% 2):T] == 0)
    )
  }, mc.cores = N_CORES)

  cat(sprintf("  %s: %d jobs done\n", param_name, nrow(jobs)))
  df <- do.call(rbind, rows)
  df$param_value <- factor(df$param_value, levels = as.character(values))
  df
}

# ── Plot helper: S and A panels side by side ──────────────────────────────────

plot_experiment <- function(df, param_name, title = NULL) {
  title <- title %||% sprintf("Effect of %s", param_name)

  make_panel <- function(y_var, y_lab) {
    mean_df <- aggregate(as.formula(paste(y_var, "~ param_value + t")),
                         data = df, FUN = mean)
    ggplot(df, aes(x = t, y = .data[[y_var]], colour = param_value,
                   group = interaction(param_value, rep))) +
      geom_line(alpha = 0.20, linewidth = 0.4) +
      geom_line(data = mean_df,
                aes(x = t, y = .data[[y_var]], colour = param_value,
                    group = param_value),
                linewidth = 1.3, inherit.aes = FALSE) +
      scale_colour_brewer(palette = "RdYlBu", direction = -1,
                          name = param_name) +
      labs(x = "Year", y = y_lab) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "right")
  }

  p_A <- make_panel("totalA", "Adults")
  p_S <- make_panel("totalS", "Seedlings")

  (p_A + p_S) +
    plot_annotation(
      title   = title,
      caption = "Thick line = mean of 3 replicates; thin lines = individual runs"
    )
}

# ── Animation: side-view dispersal + abundance over time ─────────────────────
# Requires gganimate and gifski. Install once with:
#   install.packages(c("gganimate", "gifski"))

animate_colonization <- function(result, filename = NULL, fps = 4,
                                 what = c("abundance", "dispersal")) {
  what <- match.arg(what)

  if (!requireNamespace("gganimate", quietly = TRUE))
    stop("Install gganimate: install.packages('gganimate')")
  if (!requireNamespace("gifski", quietly = TRUE))
    stop("Install gifski: install.packages('gifski')")

  T       <- result$params$timesteps
  heights <- result$heights
  yDim    <- result$params$yDim

  if (what == "abundance") {
    # Side view: sum adults + seedlings over x → [yDim, zDim] per timestep
    frames <- do.call(rbind, lapply(1:T, function(t) {
      A <- apply(result$abundanceA[,,,t], c(2, 3), sum)  # [yDim, zDim]
      S <- apply(result$abundanceS[,,,t], c(2, 3), sum)
      do.call(rbind, lapply(seq_along(heights), function(z)
        data.frame(y = 1:yDim, height = heights[z], t = t,
                   adults = A[, z], seedlings = S[, z])
      ))
    }))

    p <- ggplot(frames, aes(x = y, y = height)) +
      geom_raster(aes(fill = adults)) +
      scale_fill_gradientn(
        colours  = c("white", "#d1e5f0", "#4393c3", "#08306b"),
        name     = "Adults", na.value = "grey95"
      ) +
      geom_point(data = frames[frames$seedlings > 0, ],
                 aes(size = seedlings), colour = "#4dac26", alpha = 0.5,
                 shape = 16) +
      scale_size_continuous(name = "Seedlings", range = c(0.5, 3)) +
      labs(x = "S → N position (voxels)", y = "Height (m)",
           title = "Colonization dynamics — year {frame_time}") +
      theme_minimal(base_size = 12) +
      gganimate::transition_time(t) +
      gganimate::ease_aes("linear")

  } else {
    # Dispersal: seed rain side-view from saved disp_side [yDim, zDim, T]
    frames <- do.call(rbind, lapply(1:T, function(t)
      do.call(rbind, lapply(seq_along(heights), function(z)
        data.frame(y = 1:yDim, height = heights[z], t = t,
                   seeds = result$disp_side[, z, t])
      ))
    ))

    p <- ggplot(frames, aes(x = y, y = height, fill = seeds)) +
      geom_raster() +
      scale_fill_gradientn(
        colours  = c("white", "#fff7bc", "#fe9929", "#cc4c02"),
        name     = "Seeds", na.value = "grey95"
      ) +
      labs(x = "S → N position (voxels)", y = "Height (m)",
           title = "Seed rain — year {frame_time}") +
      theme_minimal(base_size = 12) +
      gganimate::transition_time(t) +
      gganimate::ease_aes("linear")
  }

  anim <- gganimate::animate(p, fps = fps, nframes = T,
                              width = 600, height = 400, renderer = gganimate::gifski_renderer())
  if (!is.null(filename)) {
    gganimate::anim_save(filename, anim)
    cat(sprintf("Animation saved to %s\n", filename))
  }
  anim
}

# ── 3D interactive plot: adult abundance at a given timestep ──────────────────

plot_3d_abundance <- function(result, t = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Install plotly: install.packages('plotly')")

  T <- result$params$timesteps
  if (is.null(t)) t <- T

  # Build long data frame of occupied adult voxels
  idx <- which(result$abundanceA[,,,t] > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) { message("No adults at t=", t); return(invisible(NULL)) }

  df <- data.frame(
    x       = idx[, 1],
    y       = idx[, 2],
    z       = idx[, 3],
    height  = result$heights[idx[, 3]],
    N       = result$abundanceA[,,,t][idx],
    zone    = result$zone[idx]
  )

  plotly::plot_ly(df, x = ~x, y = ~y, z = ~height,
                  color = ~factor(zone),
                  colors = c("#d73027","#fc8d59","#fee090","#91cf60","#1a9850"),
                  size  = ~N, sizes = c(3, 15),
                  type  = "scatter3d", mode = "markers",
                  marker = list(opacity = 0.75),
                  text  = ~paste0("Zone ", zone, "<br>N=", N,
                                  "<br>h=", round(height, 1), "m")) |>
    plotly::layout(
      title  = sprintf("Adult abundance at t=%d (coloured by Johansson zone)", t),
      scene  = list(
        xaxis = list(title = "x (W→E)"),
        yaxis = list(title = "y (S→N)"),
        zaxis = list(title = "Height (m)")
      ),
      legend = list(title = list(text = "Zone"))
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# BASELINE RUN (for animations)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n====== Baseline run (seed=42) ======\n")
baseline <- run_simple_colonization(base_params, seed = 42)

cat(sprintf("\nBaseline final state: S=%d  A=%d\n",
            baseline$totalS[base_params$timesteps],
            baseline$totalA[base_params$timesteps]))

# Trajectory plot: S and A together
T <- base_params$timesteps
traj_df <- data.frame(
  t     = rep(1:T, 2),
  N     = c(baseline$totalS, baseline$totalA),
  stage = rep(c("Seedling", "Adult"), each = T)
)
p_baseline <- ggplot(traj_df, aes(x = t, y = N, colour = stage)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(Seedling = "#4dac26", Adult = "#08519c")) +
  labs(x = "Year", y = "Abundance", colour = NULL,
       title = "Baseline — seedling and adult trajectory") +
  theme_minimal(base_size = 12)
print(p_baseline)

# 3D interactive at final timestep
print(plot_3d_abundance(baseline))

# Side-view abundance animation
cat("\nRendering abundance animation...\n")
anim_abund <- animate_colonization(baseline, what = "abundance",
                                   filename = "logs/anim_abundance.gif")
print(anim_abund)

# Side-view dispersal animation
cat("\nRendering dispersal animation...\n")
anim_disp <- animate_colonization(baseline, what = "dispersal",
                                  filename = "logs/anim_dispersal.gif")
print(anim_disp)

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENTS
# ─────────────────────────────────────────────────────────────────────────────

cat("\n====== Experiment 1: establishment_prob ======\n")
exp1 <- run_experiment("establishment_prob",
                       values = c(0.005, 0.01, 0.02, 0.05, 0.10))
print(plot_experiment(exp1, "establishment_prob",
                      title = "Exp 1: Effect of establishment probability"))

cat("\n====== Experiment 2: repro_rate ======\n")
exp2 <- run_experiment("repro_rate",
                       values = c(10, 30, 100, 300, 1000))
print(plot_experiment(exp2, "repro_rate",
                      title = "Exp 2: Effect of reproduction rate (seeds/adult/yr)"))

cat("\n====== Experiment 3: survival_A ======\n")
exp3 <- run_experiment("survival_A",
                       values = c(0.70, 0.78, 0.85, 0.90, 0.95))
print(plot_experiment(exp3, "survival_A",
                      title = "Exp 3: Effect of adult annual survival"))

cat("\n====== Experiment 4: maxDisp ======\n")
exp4 <- run_experiment("maxDisp",
                       values = c(2, 4, 6, 10, 15))
print(plot_experiment(exp4, "maxDisp",
                      title = "Exp 4: Effect of dispersal distance (voxels)"))

cat("\n====== Experiment 5: n_founders ======\n")
exp5 <- run_experiment("n_founders",
                       values = c(5, 20, 50, 100, 200))
print(plot_experiment(exp5, "n_founders",
                      title = "Exp 5: Effect of founder population size"))

cat("\n====== Experiment 6: carCap ======\n")
exp6 <- run_experiment("carCap",
                       values = c(2, 5, 10, 20, 50))
print(plot_experiment(exp6, "carCap",
                      title = "Exp 6: Effect of carrying capacity per voxel"))

# ── Combined overview ─────────────────────────────────────────────────────────

p_all <- list(exp1=exp1, exp2=exp2, exp3=exp3, exp4=exp4, exp5=exp5, exp6=exp6)
titles <- c("1: establishment_prob", "2: repro_rate", "3: survival_A",
            "4: maxDisp", "5: n_founders", "6: carCap")
panels <- mapply(function(df, nm, ttl) {
  # adults-only panel for the combined overview (keeps it readable)
  mean_df <- aggregate(totalA ~ param_value + t, data = df, FUN = mean)
  ggplot(df, aes(x = t, y = totalA, colour = param_value,
                 group = interaction(param_value, rep))) +
    geom_line(alpha = 0.15, linewidth = 0.3) +
    geom_line(data = mean_df, aes(x = t, y = totalA, colour = param_value,
                                   group = param_value),
              linewidth = 1.1, inherit.aes = FALSE) +
    scale_colour_brewer(palette = "RdYlBu", direction = -1, name = nm) +
    labs(x = "Year", y = "Adults", title = ttl) +
    theme_minimal(base_size = 9) +
    theme(legend.key.size = unit(0.4, "cm"))
}, p_all, names(p_all), titles, SIMPLIFY = FALSE)

combined <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title   = "Simple colonization model — sensitivity analysis (adults)",
    caption = "Thick = mean of 3 reps; thin = individual runs",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )
print(combined)

# ── Extinction rates ──────────────────────────────────────────────────────────

cat("\n====== Extinction rates ======\n")
for (nm in names(p_all)) {
  df  <- p_all[[nm]]
  ext <- aggregate(extinct ~ param_value, data = df, FUN = mean)
  ext$extinct_pct <- round(ext$extinct * 100)
  cat(sprintf("\n%s:\n", nm))
  print(ext, row.names = FALSE)
}
