# run_experiments_thesis.R
# One-at-a-time sensitivity experiments for the full thesis colonization model.
# Requires microenv_Maquipucuna.rds in data/processed/.
# Runs experiments in parallel across parameter values (one core per value×rep).
#
# Usage: Rscript scripts/run_experiments_thesis.R
#   or:  source("scripts/run_experiments_thesis.R")  (interactive)
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(parallel)
source("scripts/complex_model/paths.R")
source("scripts/complex_model/get_colonization.R")

N_CORES <- max(1L, detectCores() - 1L)
cat(sprintf("Using %d parallel cores\n", N_CORES))

# ── Load site data (shared across all runs) ───────────────────────────────────

microenv_path <- file.path(PROCESSED_DIR, "microenv_Maquipucuna.rds")
if (!file.exists(microenv_path))
  stop("microenv_Maquipucuna.rds not found — run run_microclimate.R first")

cat("Loading microenv...\n")
microenv <- readRDS(microenv_path)
if (is.null(microenv$.weather)) {
  pm <- readRDS(file.path(PROCESSED_DIR, "pointmodel_Maquipucuna.rds"))
  microenv$.weather <- pm[[1]]$weather; rm(pm)
}

niches      <- read.csv("data/csv/combinedv3.csv")
niches      <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                      !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == "Maquipucuna"],
                    na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site        <- list(Site = "Maquipucuna")

forestparams <- list(
  stems_per_ha = 298, mean_hgt = 8.4, sd_hgt = 3.5,
  mean_crown_r = 2.0, sd_crown_r = 0.8, trunk_r = 0.114,
  branch_density = 3.0, epiphyte_footprint_m2 = 0.02
)

# ── Baseline parameters ───────────────────────────────────────────────────────

base_params <- list(
  beta0S  = -0.24, beta0J  =  0.41, beta0A  =  1.73, beta1  =  0.10,
  z_S_min =  0.0,  z_S_max =  1.0,
  z_J_min =  1.0,  z_J_max =  7.0,
  z_A_min =  7.0,  z_A_max = 20.0,
  psi0S        = -3.30, psi0J        = -2.70,
  beta_precip  =  3e-4, beta_rh      =  0.010,
  sigma        =  0.10, delta_z_base =  0.80,
  cost_repro   =  0.50,
  p_poll  = 0.30, p_germ  = 0.001, p_s1 = 0.45,
  canopy_z = mean_canopy, lambda = 1, Ut = 1
)

# ── Logging (suppressed in parallel workers) ──────────────────────────────────

log_msg <- function(msg) message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

# ── Single-run wrapper ────────────────────────────────────────────────────────

run_one <- function(params, tag = "run", timesteps = 20, spinup = 3) {
  tryCatch(
    runcolonization(
      site         = site,
      niches       = niches,
      canopy_grid  = canopy_grid,
      microenv     = microenv,
      timesteps    = timesteps,
      resolution   = 10,
      carCap       = 5,
      maxDisp      = 10,
      spinup       = spinup,
      Visualize    = FALSE,
      parameters   = params,
      forestparams = forestparams
    ),
    error = function(e) { message("ERROR [", tag, "]: ", e$message); NULL }
  )
}

# ── Experiment runner ─────────────────────────────────────────────────────────
# Varies param_name across values; n_reps replicates per value.
# Returns tidy data frame of S, J, A totals over time.

run_experiment <- function(param_name, values, base = base_params,
                           n_reps = 2, timesteps = 20, spinup = 3) {
  jobs <- expand.grid(val = values, rep = seq_len(n_reps),
                      stringsAsFactors = FALSE)
  cat(sprintf("\n====== %s (%d jobs on %d cores) ======\n",
              param_name, nrow(jobs), N_CORES))

  rows <- mclapply(seq_len(nrow(jobs)), function(i) {
    val <- jobs$val[i]; rep <- jobs$rep[i]
    p   <- base; p[[param_name]] <- val
    # Suppress log_msg inside workers
    suppressMessages(
      r <- run_one(p, tag = sprintf("%s=%s rep%d", param_name, val, rep),
                   timesteps = timesteps, spinup = spinup)
    )
    if (is.null(r)) return(NULL)
    T <- timesteps
    data.frame(
      param_value = as.character(val), rep = rep, t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total  = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
  }, mc.cores = N_CORES)

  df <- do.call(rbind, Filter(Negate(is.null), rows))
  df$param_value <- factor(df$param_value, levels = as.character(values))
  cat(sprintf("  Done: %d/%d runs succeeded\n",
              length(Filter(Negate(is.null), rows)), nrow(jobs)))
  df
}

# ── Plot: S, J, A panels ──────────────────────────────────────────────────────

plot_experiment <- function(df, param_name, title = NULL) {
  title <- title %||% sprintf("Effect of %s — Maquipucuna", param_name)

  make_panel <- function(y_var, y_lab, col) {
    mean_df <- aggregate(as.formula(paste(y_var, "~ param_value + t")),
                         data = df, FUN = mean)
    ggplot(df, aes(x = t, y = .data[[y_var]], colour = param_value,
                   group = interaction(param_value, rep))) +
      geom_line(alpha = 0.20, linewidth = 0.4) +
      geom_line(data = mean_df,
                aes(x = t, y = .data[[y_var]], colour = param_value,
                    group = param_value),
                linewidth = 1.2, inherit.aes = FALSE) +
      scale_colour_brewer(palette = "RdYlBu", direction = -1, name = param_name) +
      labs(x = "Year", y = y_lab) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "right")
  }

  (make_panel("totalS", "Seedlings", "#4dac26") +
   make_panel("totalJ", "Juveniles", "#f1a340") +
   make_panel("totalA", "Adults",    "#08519c")) +
    plot_annotation(
      title   = title,
      caption = "Thick = mean of replicates; thin = individual runs"
    )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENTS
# ─────────────────────────────────────────────────────────────────────────────

# ── Exp 1: Pollination success (p_poll) ───────────────────────────────────────
# Epiphytic orchids are strongly pollinator-limited. Zotz & Schmidt (2006):
# "complete pollination would push λ above 1." Range: 0.05–0.60.

exp1 <- run_experiment("p_poll", values = c(0.05, 0.15, 0.30, 0.50, 0.70))
print(plot_experiment(exp1, "p_poll",
  title = "Exp 1: Pollination success — how pollinator-limited is persistence?"))

# ── Exp 2: Adult survival intercept (beta0A) ──────────────────────────────────
# Shifts the logistic survival curve for adults up/down.
# Default -0.24→1.73 gives ~0.85/yr. Range calibrated to 0.65–0.95/yr.

exp2 <- run_experiment("beta0A", values = c(0.50, 1.00, 1.73, 2.50, 3.50))
print(plot_experiment(exp2, "beta0A",
  title = "Exp 2: Adult survival intercept — how sensitive is persistence to mortality?"))

# ── Exp 3: Germination probability (p_germ) ───────────────────────────────────
# Mycorrhizal-gated. McCormick & Jacquemyn (2014): "dispersal not limiting;
# microsite (mycorrhizal) is." Range: 0.0001–0.005.

exp3 <- run_experiment("p_germ", values = c(0.0001, 0.0005, 0.001, 0.003, 0.005))
print(plot_experiment(exp3, "p_germ",
  title = "Exp 3: Germination probability — mycorrhizal microsite availability"))

# ── Exp 4: Cost of reproduction (cost_repro) ─────────────────────────────────
# Zotz (1998): "fruiting individuals showed reduced growth." 1.0 = no cost.

exp4 <- run_experiment("cost_repro", values = c(0.20, 0.40, 0.60, 0.80, 1.00))
print(plot_experiment(exp4, "cost_repro",
  title = "Exp 4: Cost of reproduction — growth penalty after fruiting"))

# ── Exp 5: Climate sensitivity — RH slope (beta_rh) ──────────────────────────
# Scales how much relative humidity drives S→J and J→A transitions.
# Connects the microclimate output to population dynamics.

exp5 <- run_experiment("beta_rh", values = c(0.002, 0.005, 0.010, 0.020, 0.040))
print(plot_experiment(exp5, "beta_rh",
  title = "Exp 5: RH sensitivity of stage transitions — microclimate coupling"))

# ── Exp 6: Precipitation sensitivity (beta_precip) ───────────────────────────
# Zotz & Schmidt (2006): "annual rainfall significantly affected growth."

exp6 <- run_experiment("beta_precip", values = c(1e-4, 2e-4, 3e-4, 6e-4, 1e-3))
print(plot_experiment(exp6, "beta_precip",
  title = "Exp 6: Precipitation sensitivity of stage transitions"))

# ── Combined overview (adults only, all 6) ───────────────────────────────────

all_exp <- list(p_poll=exp1, beta0A=exp2, p_germ=exp3,
                cost_repro=exp4, beta_rh=exp5, beta_precip=exp6)
titles  <- c("1: p_poll", "2: beta0A", "3: p_germ",
             "4: cost_repro", "5: beta_rh", "6: beta_precip")

panels <- mapply(function(df, nm, ttl) {
  mean_df <- aggregate(totalA ~ param_value + t, data = df, FUN = mean)
  ggplot(df, aes(x = t, y = totalA, colour = param_value,
                 group = interaction(param_value, rep))) +
    geom_line(alpha = 0.15, linewidth = 0.3) +
    geom_line(data = mean_df,
              aes(x = t, y = totalA, colour = param_value, group = param_value),
              linewidth = 1.0, inherit.aes = FALSE) +
    scale_colour_brewer(palette = "RdYlBu", direction = -1, name = nm) +
    labs(x = "Year", y = "Adults", title = ttl) +
    theme_minimal(base_size = 9) +
    theme(legend.key.size = unit(0.35, "cm"))
}, all_exp, names(all_exp), titles, SIMPLIFY = FALSE)

print(
  wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title   = "Thesis model sensitivity — Maquipucuna (adults)",
    caption = "All other params at literature defaults (Zotz 2006, McCormick 2014)",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )
)

# ── Extinction rates ──────────────────────────────────────────────────────────

cat("\n====== Extinction rates ======\n")
for (nm in names(all_exp)) {
  ext <- aggregate(extinct ~ param_value, data = all_exp[[nm]], FUN = mean)
  ext$pct <- round(ext$extinct * 100)
  cat(sprintf("\n%s:\n", nm)); print(ext[, c("param_value","pct")], row.names = FALSE)
}

# ── Save results ──────────────────────────────────────────────────────────────

out_path <- file.path(PROCESSED_DIR, "experiments_thesis_Maquipucuna.rds")
saveRDS(all_exp, out_path)
cat(sprintf("\nResults saved to %s\n", out_path))
